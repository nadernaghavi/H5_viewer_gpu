import datashader as ds
import datashader.transfer_functions as tf
import time
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
import pandas as pd
import rapids_singlecell as rsc
import cupy as cp
import scipy.sparse as sp
import cudf
import hvplot.cudf  # enables .hvplot on cudf
import holoviews as hv
import asyncio
from concurrent.futures import ThreadPoolExecutor
import io
from fastapi import FastAPI, UploadFile, File, HTTPException, Request
from fastapi import FastAPI, UploadFile, File, Form, HTTPException
from fastapi.concurrency import run_in_threadpool
from fastapi.middleware.cors import CORSMiddleware
import base64
import traceback
import tempfile
import uvicorn
import nest_asyncio
import threading
from dataclasses import dataclass
import os
import secrets
import gc
import psutil
from math import sqrt
import re
import scipy.sparse as sp


app = FastAPI()
app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://localhost:3000",
        "http://127.0.0.1:3000",
        "http://localhost:5173",
        "http://127.0.0.1:5173",
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

canvas_scatter = ds.Canvas(plot_width=1000, plot_height=1000)
canvas_box = ds.Canvas(plot_width=200, plot_height=1000)

data_obj = None
log_message = ""

class SpatialRNAProcessor:
    def __init__(self, data_file: str, n_top_genes: int = 5):
        self.data_file = data_file
        self.n_top_genes = n_top_genes
        self._pool = ThreadPoolExecutor(max_workers=8)
        self.adata = None
        self.barcodes = None
        self.ncount_spatial = None
        self.variable_genes_index = None
        self.variable_genes = None

    def load_data(self):
        """Load the 10x H5 data into AnnData."""
        #self.adata = sc.read_10x_h5(self.data_file)
        try:
            # 🔥 try 10x first
            self.adata = sc.read_10x_h5(self.data_file)
            self.adata.var_names_make_unique()
            print("Loaded as 10x .h5")
        
        except Exception as e:
            print("10x load failed, trying AnnData...", str(e))
    
            # 🔁 fallback to AnnData
            try:
                self.adata = sc.read_h5ad(self.data_file)
                print("Converting dense matrix to sparse...")
                self.adata.X = sp.csr_matrix(self.adata.X)
                self.adata.X = self.adata.X.astype("float32")
                print("Loaded as AnnData (.h5ad/.h5)")
            except Exception as e2:
                raise RuntimeError(f"Failed to load file as 10x or AnnData: {e2}")
        self.adata.var_names_make_unique()
        print("Loaded as 10x .h5")
        #self.data_type = self.detect_anndata_type_from_first_barcode()
        self.region_num, self.gene_num = self.adata.shape
        self.canvas_size = round(sqrt(self.region_num))
        self.canvas_scatter = ds.Canvas(plot_width=self.canvas_size, plot_height=self.canvas_size)
        print(f"Number of sections in the data: {self.region_num}")
        return self.adata

    def region_num_return(self):
        return self.region_num
        
    def extract_barcodes(self):
        """Extract pixel barcodes."""
        if self.adata is None:
            raise ValueError("Data not loaded. Call load_data() first.")

        self.barcodes = self.adata.obs_names.to_numpy()
        return self.barcodes

    def compute_total_counts(self):
        """Compute total RNA counts per pixel."""
        if self.adata is None:
            raise ValueError("Data not loaded. Call load_data() first.")
        self.ncount_spatial = cp.asarray(self.adata.X.sum(axis=1)).ravel()
        return self.ncount_spatial

    def find_highly_variable_genes(self):
        """Move AnnData to GPU and find top highly variable genes."""
        if self.adata is None:
            raise ValueError("Data not loaded. Call load_data() first.")
        rsc.get.anndata_to_GPU(self.adata)
        # ✅ 1. Normalize total counts per cell
        rsc.pp.normalize_total(
            self.adata,
            target_sum=1e4   # standard (CP10K)
        )
        
        # ✅ 2. Log transform
        rsc.pp.log1p(self.adata)
        rsc.pp.highly_variable_genes(
            self.adata,
            n_top_genes=self.n_top_genes,
            flavor="seurat_v3"
        )

        self.top_gene_idx = self.adata.var.highly_variable.to_numpy()
        self.variable_genes = self.adata.var_names[self.top_gene_idx]
        
        return self.variable_genes

    def min_max_find(self):
        self.min_dict = {}
        self.max_dict = {}
        self.min_dict["gene"] = 0
        self.max_dict["gene"] = 1
        self.min_dict["total"] = 0
        self.max_dict["total"] = 1
        return self.min_dict, self.max_dict 

    def gpu_table_create(self):
        # extract X and Y
        df = pd.DataFrame({"barcode": self.barcodes})

        df[["X", "Y"]] = df["barcode"].str.extract(r"s_\d+um_(\d+)_(\d+)-")
        
        # convert to int
        df["X"] = df["X"].astype(int)
        df["Y"] = df["Y"].astype(int)
        
        # move ncount_spatial to GPU
        y_gpu = cp.asarray(self.ncount_spatial).astype(cp.float32)
        n = len(y_gpu)
    
        # original GPU dataframe
        self.df_gpu = cudf.DataFrame({
            'x': df['X'],
            'y': df['Y'],
            'value': y_gpu,
            'gene_value': y_gpu
        })
    
        # bounds (for original df)
        self.x_max = self.df_gpu['x'].max()
        self.x_min = self.df_gpu['x'].min()
        self.y_max = self.df_gpu['y'].max()
        self.y_min = self.df_gpu['y'].min()
        
        return self.df_gpu

    def scatter_plot(self, df_gpu, var="value", adata=None, gene_name=None, vmin=0, vmax=1):
        vmin = float(vmin)
        stream = cp.cuda.Stream(non_blocking=True)

        with stream:
            if len(df_gpu) == 0:
                raise ValueError("df_gpu is empty")
            n = len(df_gpu)
            temp_canvas_size = round(sqrt(n))
            temp_canvas = ds.Canvas(plot_width=temp_canvas_size, 
                                    plot_height=temp_canvas_size)
            if var == "value":
                df_plot = df_gpu.copy(deep=False)
                value_col = "value"
            else:
                gene_idx = adata.var_names.get_loc(gene_name)
        
                # extract one column on GPU
                col = adata.X[:, gene_idx]
                if hasattr(col, "toarray"):
                    gene_expr = cp.asarray(col.toarray()).ravel()
                else:
                    gene_expr = cp.asarray(col).ravel()
                #if sp.issparse(col):
                #    col = col.tocoo()
                #    gene_expr = cp.zeros(col.shape[0], dtype=cp.float32)
                #    gene_expr[cp.asarray(col.row)] = cp.asarray(col.data, dtype=cp.float32)
                #else:
                #    gene_expr = cp.asarray(col, dtype=cp.float32).ravel()
        
                df_plot = df_gpu.copy(deep=False)
                df_plot["gene_value"] = gene_expr.astype(cp.float32)
                value_col = "gene_value"
        
            v_max_new = float(df_plot[value_col].max())
        
            if not np.isfinite(v_max_new) or v_max_new <= vmin:
                v_max_new = vmin + 1.0
        
            agg = temp_canvas.points(
                df_plot,
                "x",
                "y",
                agg=ds.mean(value_col)
            )

        stream.synchronize()
        img = tf.shade(
            agg,
            cmap=[
                "#000000",
                "#2c105c",
                "#711f81",
                "#b63679",
                "#ee605e",
                "#fdae61",
                "#ffd166"
            ],
            span=(vmin, v_max_new),
            how="log"
        )              
        img = tf.spread(img, px=1)
    
        pil_img = img.to_pil()
    
        with io.BytesIO() as buf:
            pil_img.save(buf, format="PNG")
            return buf.getvalue()


    def gpu_box(self, df_gpu, var="value", adata=None, gene_name=None, vmin=0, vmax=0) -> bytes:
        """
        GPU boxplot stats, then PNG on CPU.
        If there is no valid data after filtering, draw a degenerate boxplot at y=0.
        """
        vmin = float(vmin)
        vmax = float(vmax)
    
        stream = cp.cuda.Stream(non_blocking=True)
    
        with stream:
            if var != "value":
                gene_idx = adata.var_names.get_loc(gene_name)
    
                col = adata.X[:, gene_idx]
                if hasattr(col, "toarray"):
                    gene_expr = cp.asarray(col.toarray()).ravel()
                else:
                    gene_expr = cp.asarray(col).ravel()
                #if sp.issparse(col):
                #    col = col.tocoo()
                #    gene_expr = cp.zeros(col.shape[0], dtype=cp.float32)
                #    gene_expr[cp.asarray(col.row)] = cp.asarray(col.data, dtype=cp.float32)
                #else:
                #    gene_expr = cp.asarray(col, dtype=cp.float32).ravel()
    
                # do not mutate shared df_gpu
                df_local = df_gpu.copy(deep=False)
                df_local["gene_value"] = gene_expr.astype(cp.float32)
    
                x = df_local["gene_value"]
                mask = (x != 0) & (~x.isna())
                x = x[mask].astype("float32")
                title_name = gene_name
            else:
                x = df_gpu["value"].dropna().astype("float32")
                title_name = "all counts"
    
            # remove inf / -inf too
            if len(x) > 0:
                x_cp = x.to_cupy()
                finite_mask = cp.isfinite(x_cp)
                x = x[finite_mask]
    
            # empty case -> collapsed box at zero
            if len(x) == 0:
                min_v = q1 = median = q3 = max_v = 0.0
                lower_whisker = upper_whisker = 0.0
            else:
                q = x.quantile([0.0, 0.25, 0.5, 0.75, 1.0])
                q_cpu = q.to_pandas().values.astype(float)
    
                min_v, q1, median, q3, max_v = q_cpu
    
                # guard against weird quantile output
                if not np.all(np.isfinite([min_v, q1, median, q3, max_v])):
                    min_v = q1 = median = q3 = max_v = 0.0
                    lower_whisker = upper_whisker = 0.0
                else:
                    iqr = q3 - q1
                    lower_whisker = max(min_v, q1 - 1.5 * iqr)
                    upper_whisker = min(max_v, q3 + 1.5 * iqr)
    
        stream.synchronize()
    
        fig, ax = plt.subplots(figsize=(4, 6))
    
        # box
        ax.plot([1, 1], [q1, q3], linewidth=10)
    
        # median
        ax.plot([0.9, 1.1], [median, median], linewidth=3)
    
        # whiskers
        ax.plot([1, 1], [lower_whisker, q1], linewidth=2)
        ax.plot([1, 1], [q3, upper_whisker], linewidth=2)
    
        # whisker caps
        ax.plot([0.95, 1.05], [lower_whisker, lower_whisker], linewidth=2)
        ax.plot([0.95, 1.05], [upper_whisker, upper_whisker], linewidth=2)
    
        ax.set_xlim(0.8, 1.2)
        ax.set_xticks([])
        ax.set_ylabel("")
        ax.set_title(title_name)
    
        # make zero-only case visible
        if min_v == max_v:
            ax.set_ylim(min_v - 1.0, max_v + 1.0)
    
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=120, bbox_inches="tight", pad_inches=0)
        plt.close(fig)
    
        buf.seek(0)
        return buf.getvalue()
        
    async def run_plot_in_thread(self, fn, df_gpu, var, adata, gene_name, vmin, vmax):
        loop = asyncio.get_running_loop()
        return await loop.run_in_executor(self._pool, fn, df_gpu, var, adata, gene_name, vmin,
                                         vmax)

    async def run_all_fig(self, df_gpu, adata, gene_list):

        results = await asyncio.gather(
            self.run_plot_in_thread(self.scatter_plot, df_gpu, "value", None, None, 
                                    self.min_dict["total"], self.max_dict["total"]),
            self.run_plot_in_thread(self.gpu_box, df_gpu, "value", None, None,
                                   self.min_dict["total"], self.max_dict["total"]),
            self.run_plot_in_thread(self.scatter_plot, df_gpu, "gene_value", adata, gene_list[0],
                                   self.min_dict["gene"], self.max_dict["gene"]),
            self.run_plot_in_thread(self.gpu_box, df_gpu, "gene_value", adata, gene_list[0],
                                   self.min_dict["gene"], self.max_dict["gene"]),
            self.run_plot_in_thread(self.scatter_plot, df_gpu, "gene_value", adata, gene_list[1],
                                   self.min_dict["gene"], self.max_dict["gene"]),
            self.run_plot_in_thread(self.gpu_box, df_gpu, "gene_value", adata, gene_list[1],
                                   self.min_dict["gene"], self.max_dict["gene"]),
            self.run_plot_in_thread(self.scatter_plot, df_gpu, "gene_value", adata, gene_list[2],
                                   self.min_dict["gene"], self.max_dict["gene"]),
            self.run_plot_in_thread(self.gpu_box, df_gpu, "gene_value", adata, gene_list[2],
                                   self.min_dict["gene"], self.max_dict["gene"])
        )
    
        return results  # list of 4 PNG byte blobs

    def pixel_to_data(self, px, py, plot_width, plot_height, x_max, y_max):
        
        x = 0 + (px / plot_width) * x_max
        y = y_max - (py / plot_height) * y_max   # image y is top->down
    
        return x, y

    async def run_all_zoom(self, x_start, x_end, y_start, y_end):

        px_min, py_min = self.pixel_to_data(x_start, y_start, self.canvas_size, self.canvas_size,
                                      self.x_max, self.y_max)
        px_max, py_max = self.pixel_to_data(x_end, y_end, self.canvas_size, self.canvas_size,
                                     self.x_max, self.y_max)

        tab_x_min = min(px_min, px_max)
        tab_y_min = min(py_min, py_max)
        tab_x_max = max(px_min, px_max)
        tab_y_max = max(py_min, py_max)
        
        mask_gpu = (
            (self.df_gpu["x"] >= tab_x_min) &
            (self.df_gpu["x"] <= tab_x_max) &
            (self.df_gpu["y"] >= tab_y_min) &
            (self.df_gpu["y"] <= tab_y_max)
        )
        df_gpu_filtered = self.df_gpu[mask_gpu]
        mask_cpu = mask_gpu.to_pandas().to_numpy()
        adata_filtered = self.adata[mask_cpu].copy()
        png_list = await self.run_all_fig(df_gpu_filtered, adata_filtered, self.variable_genes)
        return png_list

    async def run_all(self):
        """Run the full pipeline."""
        self.load_data()
        print("data loaded")
        self.extract_barcodes()
        print("barcodes extracted")
        self.find_highly_variable_genes()
        print("highly variable genes found)")
        self.compute_total_counts()
        print("found total count")
        self.min_max_find()
        print("min max done")
        self.gpu_table_create()
        print("gpu table created")
        self.init_img_list = await self.run_all_fig(self.df_gpu, self.adata, self.variable_genes)
        print("inital images found")
        return self.region_num

    async def init_img_return(self):
        return self.init_img_list

    def var_gene_return(self):
        return self.variable_genes

    

class SingleRNAProcessor:
    def __init__(self, data_file: str, dim_file: str, n_top_genes: int = 5):
        self.data_file = data_file
        self.dim_file = dim_file
        self.n_top_genes = n_top_genes
        self._pool = ThreadPoolExecutor(max_workers=8)
        self.adata = None
        self.barcodes = None
        self.ncount_spatial = None
        self.variable_genes_index = None
        self.variable_genes = None

    def fix_anndata_X_inplace(self):
        if sp.issparse(self.adata.X):
            bad = ~np.isfinite(self.adata.X.data)
            self.adata.X.data[bad] = 0
        else:
            self.adata.X[~np.isfinite(self.adata.X)] = 0

    def load_data(self):
        """Load the 10x H5 data into AnnData."""
        self.adata = sc.read_10x_h5(self.data_file)
        self.fix_anndata_X_inplace()
        self.adata.var_names_make_unique()
        print("Loaded as 10x .h5")
        #self.data_type = self.detect_anndata_type_from_first_barcode()
        self.region_num, self.gene_num = self.adata.shape
        self.canvas_size = round(sqrt(self.region_num))
        self.canvas_scatter = ds.Canvas(plot_width=self.canvas_size, plot_height=self.canvas_size)
        print(f"Number of sections in the data: {self.region_num}")
        return self.adata

    def load_and_match_barcodes_old(self, csv_path, barcodes):
    
        df = pd.read_csv(csv_path)
    
        if "Barcode" not in df.columns:
            raise ValueError("CSV must contain 'Barcode' column")
    
        # make sure everything is string
        df["Barcode"] = df["Barcode"].astype(str)
        barcodes = barcodes.astype(str)
    
        # set index for fast matching
        df = df.set_index("Barcode")
    
        # align to your barcodes (THIS is the key line)
        df_matched = df.reindex(barcodes)
    
        # optional check
        missing = df_matched.isna().any(axis=1).sum()
        if missing > 0:
            print(f"Warning: {missing} barcodes not found")
    
        return df_matched

    def load_and_match_barcodes(self, csv_path, barcodes):
        df = pd.read_csv(csv_path)
        return df
    

    def load_dim_data(self):
        self.dim_data = self.load_and_match_barcodes(self.dim_file, self.barcodes)
    
    def region_num_return(self):
        return self.region_num
        
    def extract_barcodes(self):
        """Extract pixel barcodes."""
        if self.adata is None:
            raise ValueError("Data not loaded. Call load_data() first.")

        self.barcodes = self.adata.obs_names.to_numpy()
        return self.barcodes

    def compute_total_counts(self):
        """Compute total RNA counts per pixel."""
        if self.adata is None:
            raise ValueError("Data not loaded. Call load_data() first.")
        self.ncount_spatial = cp.asarray(self.adata.X.sum(axis=1)).ravel()
        return self.ncount_spatial

    def find_highly_variable_genes_old(self):
        """Move AnnData to GPU and find top highly variable genes."""
        if self.adata is None:
            raise ValueError("Data not loaded. Call load_data() first.")
        print("Starting gpu load")
        rsc.get.anndata_to_GPU(self.adata)
        print("final gpu transfer")
        # ✅ 1. Normalize total counts per cell
        rsc.pp.normalize_total(
            self.adata,
            target_sum=1e4   # standard (CP10K)
        )
        
        # ✅ 2. Log transform
        rsc.pp.log1p(self.adata)
        rsc.pp.highly_variable_genes(
            self.adata,
            n_top_genes=self.n_top_genes,
            flavor="seurat_v3"
        )

        self.top_gene_idx = self.adata.var.highly_variable.to_numpy()
        self.variable_genes = self.adata.var_names[self.top_gene_idx]
        
        return self.variable_genes

    def find_highly_variable_genes(self):
        if self.adata is None:
            raise ValueError("Data not loaded. Call load_data() first.")
    
        rsc.get.anndata_to_GPU(self.adata)
    
        # seurat_v3 expects raw counts, so do HVG first
        rsc.pp.highly_variable_genes(
            self.adata,
            n_top_genes=self.n_top_genes,
            flavor="seurat_v3"
        )
    
        self.top_gene_idx = self.adata.var.highly_variable.to_numpy()
        self.variable_genes = self.adata.var_names[self.top_gene_idx]
    
        # only normalize/log AFTER HVG selection if you want
        rsc.pp.normalize_total(self.adata, target_sum=1e4)
        rsc.pp.log1p(self.adata)
    
        return self.variable_genes

    def min_max_find(self):
        self.min_dict = {}
        self.max_dict = {}
        self.min_dict["gene"] = 0
        self.max_dict["gene"] = 1
        self.min_dict["total"] = 0
        self.max_dict["total"] = 1
        return self.min_dict, self.max_dict 

    def gpu_table_create(self):
        # extract X and Y
        df = pd.DataFrame({"barcode": self.barcodes})

        df[["X", "Y"]] = self.dim_data[["TSNE-1","TSNE-2"]]
        
        # convert to int
        df["X"] = df["X"].astype(int)
        df["Y"] = df["Y"].astype(int)
        
        # move ncount_spatial to GPU
        y_gpu = cp.asarray(self.ncount_spatial).astype(cp.float32)
        n = len(y_gpu)
    
        # original GPU dataframe
        self.df_gpu = cudf.DataFrame({
            'x': df['X'],
            'y': df['Y'],
            'value': y_gpu,
            'gene_value': y_gpu
        })
    
        # bounds (for original df)
        self.x_max = self.df_gpu['x'].max()
        self.x_min = self.df_gpu['x'].min()
        self.y_max = self.df_gpu['y'].max()
        self.y_min = self.df_gpu['y'].min()
        
        return self.df_gpu

    def gpu_table_create_old(self):
        df = pd.DataFrame({"barcode": self.barcodes})
    
        # pull matched coordinates
        df[["X", "Y"]] = self.dim_data[["UMAP-1", "UMAP-2"]]
    
        # force numeric; bad strings become NaN
        df["X"] = pd.to_numeric(df["X"], errors="coerce")
        df["Y"] = pd.to_numeric(df["Y"], errors="coerce")
    
        # keep only rows with finite coordinates
        #good = np.isfinite(df["X"].to_numpy()) & np.isfinite(df["Y"].to_numpy())
    
        #dropped = (~good).sum()
        #if dropped > 0:
        #    print(f"Dropping {dropped} rows with missing/invalid UMAP coords")
    
        #df = df.loc[good].copy()
    
        # keep ncount and adata aligned with filtered rows
        #if isinstance(self.ncount_spatial, cp.ndarray):
        #    self.ncount_spatial = self.ncount_spatial[cp.asarray(good)]
        #else:
        #    self.ncount_spatial = self.ncount_spatial[good]
    
        #self.adata = self.adata[good].copy()
        #self.barcodes = self.barcodes[good]
        #self.dim_data = self.dim_data.loc[self.barcodes]
    
        # do NOT cast to int unless you truly need ints
        df["X"] = df["X"].astype(np.float32)
        df["Y"] = df["Y"].astype(np.float32)
    
        y_gpu = cp.asarray(self.ncount_spatial).astype(cp.float32)
    
        self.df_gpu = cudf.DataFrame({
            "x": df["X"].to_numpy(),
            "y": df["Y"].to_numpy(),
            "value": y_gpu,
            "gene_value": y_gpu
        })
    
        self.x_max = float(self.df_gpu["x"].max())
        self.x_min = float(self.df_gpu["x"].min())
        self.y_max = float(self.df_gpu["y"].max())
        self.y_min = float(self.df_gpu["y"].min())
    
        return self.df_gpu

    def scatter_plot(self, df_gpu, var="value", adata=None, gene_name=None, vmin=0, vmax=1):
        vmin = float(vmin)
        stream = cp.cuda.Stream(non_blocking=True)

        with stream:
            if len(df_gpu) == 0:
                raise ValueError("df_gpu is empty")
            n = len(df_gpu)
            temp_canvas_size = round(sqrt(n))
            temp_canvas = ds.Canvas(plot_width=temp_canvas_size, 
                                    plot_height=temp_canvas_size)
            if var == "value":
                df_plot = df_gpu.copy(deep=False)
                value_col = "value"
            else:
                gene_idx = adata.var_names.get_loc(gene_name)
        
                # extract one column on GPU
                col = adata.X[:, gene_idx]
                if hasattr(col, "toarray"):
                    gene_expr = cp.asarray(col.toarray()).ravel()
                else:
                    gene_expr = cp.asarray(col).ravel()
                #if sp.issparse(col):
                #    col = col.tocoo()
                #    gene_expr = cp.zeros(col.shape[0], dtype=cp.float32)
                #    gene_expr[cp.asarray(col.row)] = cp.asarray(col.data, dtype=cp.float32)
                #else:
                #    gene_expr = cp.asarray(col, dtype=cp.float32).ravel()
        
                df_plot = df_gpu.copy(deep=False)
                df_plot["gene_value"] = gene_expr.astype(cp.float32)
                value_col = "gene_value"
        
            v_max_new = float(df_plot[value_col].max())
        
            if not np.isfinite(v_max_new) or v_max_new <= vmin:
                v_max_new = vmin + 1.0
        
            agg = temp_canvas.points(
                df_plot,
                "x",
                "y",
                agg=ds.mean(value_col)
            )

        stream.synchronize()
        img = tf.shade(
            agg,
            cmap=[
                "#000000",
                "#2c105c",
                "#711f81",
                "#b63679",
                "#ee605e",
                "#fdae61",
                "#ffd166"
            ],
            span=(vmin, v_max_new),
            how="log"
        )              
        img = tf.spread(img, px=1)
    
        pil_img = img.to_pil()
    
        with io.BytesIO() as buf:
            pil_img.save(buf, format="PNG")
            return buf.getvalue()


    def gpu_box(self, df_gpu, var="value", adata=None, gene_name=None, vmin=0, vmax=0) -> bytes:
        """
        GPU boxplot stats, then PNG on CPU.
        If there is no valid data after filtering, draw a degenerate boxplot at y=0.
        """
        vmin = float(vmin)
        vmax = float(vmax)
    
        stream = cp.cuda.Stream(non_blocking=True)
    
        with stream:
            if var != "value":
                gene_idx = adata.var_names.get_loc(gene_name)
    
                col = adata.X[:, gene_idx]
                if hasattr(col, "toarray"):
                    gene_expr = cp.asarray(col.toarray()).ravel()
                else:
                    gene_expr = cp.asarray(col).ravel()
                #if sp.issparse(col):
                #    col = col.tocoo()
                #    gene_expr = cp.zeros(col.shape[0], dtype=cp.float32)
                #    gene_expr[cp.asarray(col.row)] = cp.asarray(col.data, dtype=cp.float32)
                #else:
                #    gene_expr = cp.asarray(col, dtype=cp.float32).ravel()
    
                # do not mutate shared df_gpu
                df_local = df_gpu.copy(deep=False)
                df_local["gene_value"] = gene_expr.astype(cp.float32)
    
                x = df_local["gene_value"]
                mask = (x != 0) & (~x.isna())
                x = x[mask].astype("float32")
                title_name = gene_name
            else:
                x = df_gpu["value"].dropna().astype("float32")
                title_name = "all counts"
    
            # remove inf / -inf too
            if len(x) > 0:
                x_cp = x.to_cupy()
                finite_mask = cp.isfinite(x_cp)
                x = x[finite_mask]
    
            # empty case -> collapsed box at zero
            if len(x) == 0:
                min_v = q1 = median = q3 = max_v = 0.0
                lower_whisker = upper_whisker = 0.0
            else:
                q = x.quantile([0.0, 0.25, 0.5, 0.75, 1.0])
                q_cpu = q.to_pandas().values.astype(float)
    
                min_v, q1, median, q3, max_v = q_cpu
    
                # guard against weird quantile output
                if not np.all(np.isfinite([min_v, q1, median, q3, max_v])):
                    min_v = q1 = median = q3 = max_v = 0.0
                    lower_whisker = upper_whisker = 0.0
                else:
                    iqr = q3 - q1
                    lower_whisker = max(min_v, q1 - 1.5 * iqr)
                    upper_whisker = min(max_v, q3 + 1.5 * iqr)
    
        stream.synchronize()
    
        fig, ax = plt.subplots(figsize=(4, 6))
    
        # box
        ax.plot([1, 1], [q1, q3], linewidth=10)
    
        # median
        ax.plot([0.9, 1.1], [median, median], linewidth=3)
    
        # whiskers
        ax.plot([1, 1], [lower_whisker, q1], linewidth=2)
        ax.plot([1, 1], [q3, upper_whisker], linewidth=2)
    
        # whisker caps
        ax.plot([0.95, 1.05], [lower_whisker, lower_whisker], linewidth=2)
        ax.plot([0.95, 1.05], [upper_whisker, upper_whisker], linewidth=2)
    
        ax.set_xlim(0.8, 1.2)
        ax.set_xticks([])
        ax.set_ylabel("")
        ax.set_title(title_name)
    
        # make zero-only case visible
        if min_v == max_v:
            ax.set_ylim(min_v - 1.0, max_v + 1.0)
    
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=120, bbox_inches="tight", pad_inches=0)
        plt.close(fig)
    
        buf.seek(0)
        return buf.getvalue()
        
    async def run_plot_in_thread(self, fn, df_gpu, var, adata, gene_name, vmin, vmax):
        loop = asyncio.get_running_loop()
        return await loop.run_in_executor(self._pool, fn, df_gpu, var, adata, gene_name, vmin,
                                         vmax)

    async def run_all_fig(self, df_gpu, adata, gene_list):

        results = await asyncio.gather(
            self.run_plot_in_thread(self.scatter_plot, df_gpu, "value", None, None, 
                                    self.min_dict["total"], self.max_dict["total"]),
            self.run_plot_in_thread(self.gpu_box, df_gpu, "value", None, None,
                                   self.min_dict["total"], self.max_dict["total"]),
            self.run_plot_in_thread(self.scatter_plot, df_gpu, "gene_value", adata, gene_list[0],
                                   self.min_dict["gene"], self.max_dict["gene"]),
            self.run_plot_in_thread(self.gpu_box, df_gpu, "gene_value", adata, gene_list[0],
                                   self.min_dict["gene"], self.max_dict["gene"]),
            self.run_plot_in_thread(self.scatter_plot, df_gpu, "gene_value", adata, gene_list[1],
                                   self.min_dict["gene"], self.max_dict["gene"]),
            self.run_plot_in_thread(self.gpu_box, df_gpu, "gene_value", adata, gene_list[1],
                                   self.min_dict["gene"], self.max_dict["gene"]),
            self.run_plot_in_thread(self.scatter_plot, df_gpu, "gene_value", adata, gene_list[2],
                                   self.min_dict["gene"], self.max_dict["gene"]),
            self.run_plot_in_thread(self.gpu_box, df_gpu, "gene_value", adata, gene_list[2],
                                   self.min_dict["gene"], self.max_dict["gene"])
        )
    
        return results  # list of 4 PNG byte blobs

    def pixel_to_data(self, px, py, plot_width, plot_height, x_max, y_max, x_min, y_min):
        
        x = x_min + (px / plot_width) * (x_max-x_min)
        y = y_max - (py / plot_height) * (y_max-y_min)   # image y is top->down
    
        return x, y

    async def run_all_zoom(self, x_start, x_end, y_start, y_end):

        px_min, py_min = self.pixel_to_data(x_start, y_start, self.canvas_size, self.canvas_size,
                                      self.x_max, self.y_max, self.x_min, self.y_min)
        px_max, py_max = self.pixel_to_data(x_end, y_end, self.canvas_size, self.canvas_size,
                                     self.x_max, self.y_max, self.x_min, self.y_min)

        tab_x_min = min(px_min, px_max)
        tab_y_min = min(py_min, py_max)
        tab_x_max = max(px_min, px_max)
        tab_y_max = max(py_min, py_max)
        
        mask_gpu = (
            (self.df_gpu["x"] >= tab_x_min) &
            (self.df_gpu["x"] <= tab_x_max) &
            (self.df_gpu["y"] >= tab_y_min) &
            (self.df_gpu["y"] <= tab_y_max)
        )
        df_gpu_filtered = self.df_gpu[mask_gpu]
        mask_cpu = mask_gpu.to_pandas().to_numpy()
        adata_filtered = self.adata[mask_cpu].copy()
        png_list = await self.run_all_fig(df_gpu_filtered, adata_filtered, self.variable_genes)
        return png_list

    async def run_all(self):
        """Run the full pipeline."""
        self.load_data()
        print("data loaded")
        self.extract_barcodes()
        print("barcodes extracted")
        self.load_dim_data()
        self.find_highly_variable_genes()
        print("highly variable genes found)")
        self.compute_total_counts()
        print("found total count")
        self.min_max_find()
        print("min max done")
        self.gpu_table_create()
        print("gpu table created")
        self.init_img_list = await self.run_all_fig(self.df_gpu, self.adata, self.variable_genes)
        print("inital images found")
        return self.region_num

    async def init_img_return(self):
        return self.init_img_list

    def var_gene_return(self):
        return self.variable_genes
    
@app.post("/run-h5")
async def run_h5(
    file1: UploadFile = File(...),
    file2: UploadFile = File(...),
    user_text: str = Form(...)
):
    global data_obj
    global log_message

    temp_path1 = None
    temp_path2 = None

    try:
        # validate extensions (optional, adjust as needed)
        if not file1.filename.lower().endswith(".h5"):
            raise HTTPException(status_code=400, detail="file1 must be .h5")

        # save file1
        with tempfile.NamedTemporaryFile(delete=False, suffix=".h5") as tmp1:
            content1 = await file1.read()
            tmp1.write(content1)
            temp_path1 = tmp1.name

        if user_text == "RNA":
            # save file2
            with tempfile.NamedTemporaryFile(delete=False, suffix=".csv") as tmp2:
                content2 = await file2.read()
                tmp2.write(content2)
                temp_path2 = tmp2.name

        print("Preprocessing started!")
        print(f"User text: {user_text}")
        log_message += f"Data type is {user_text}"

        # 🔥 your logic — modify constructor if needed
        if user_text == "SPATIAL":
            data_obj = SpatialRNAProcessor(
                temp_path1,
                n_top_genes=5
            )
        else:
            data_obj = SingleRNAProcessor(
                temp_path1,
                temp_path2,
                n_top_genes=5
            )

        region_num = await data_obj.run_all()

        log_temp = f"Total number of regions/cells found in the data: {region_num}\n"
        log_message += log_temp

        print("Preprocessing finished!")

        return {
            "status": "ok",
            "filename1": file1.filename,
            "filename2": file2.filename,
            "text": user_text,
            "count": region_num
        }

    except HTTPException:
        raise

    except Exception as e:
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=str(e))

    finally:
        if temp_path1 and os.path.exists(temp_path1):
            os.remove(temp_path1)
        if temp_path2 and os.path.exists(temp_path2):
            os.remove(temp_path2)

@app.get("/init_render")
async def init_render():
    global data_obj
    global log_message
    try:
        start = time.time()
        init_images = await data_obj.init_img_return()
        
        # return images as base64 strings
        encoded_images = [
            base64.b64encode(png).decode("utf-8")
            for png in init_images
        ]

        gene_list = data_obj.var_gene_return()
        
        images_titles = ["Total Count", "Total Count"] + [g for gene in gene_list for g in (gene, gene)]
        end = time.time()
        log_temp = f"Total init render time: {end - start:.4f} sec\n"
        log_message += log_temp
        print(log_temp)
        return {
            "status": "ok",
            "count": len(encoded_images),
            "images": encoded_images,
            "titles": images_titles
        }
    except HTTPException:
        raise

    except Exception as e:
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/submit-bbox")
async def submit_bbox(request: Request):
    global data_obj
    global log_message
    data = await request.json()   # ← THIS is the key line

    x_start = float(data["x_start"])
    x_end   = float(data["x_end"])
    y_start = float(data["y_start"])
    y_end   = float(data["y_end"])

    print("Received:")
    print(x_start, x_end, y_start, y_end)
    start = time.time()

            # run pipeline on saved file
    png_list = await data_obj.run_all_zoom(x_start, x_end, y_start, y_end)

        # return images as base64 strings
    encoded_images = [
        base64.b64encode(png).decode("utf-8")
        for png in png_list
    ]
    end = time.time()
    log_temp = f"Total zoom render time: {end - start:.4f} sec\n"
    print(log_temp)
    log_message += log_temp
    
    return {
        "status": "ok",
        "count": len(encoded_images),
        "images": encoded_images
    }

@app.get("/func_end")
async def func_end():
        global data_obj
        global log_message
        region_num = data_obj.region_num_return()
        process = psutil.Process(os.getpid())
        mem_mb = process.memory_info().rss / (1024 ** 2)

        log_temp = f"For number of data point: {region_num}, RAM usage is: {mem_mb:.2f} MB \n\n\n"
        print(log_temp)
        log_message += log_temp

        with open("log.txt", "a") as f:   # "a" = append mode
            f.write(log_message)

        del data_obj
        log_message = ""
        gc.collect()

        return {
            "status": "ok",
        }
        



    
