import React, { useRef, useState } from "react";

const API_BASE = "/api";

export default function App() {
  const [selectedMode, setSelectedMode] = useState("");

  const [selectedFile1, setSelectedFile1] = useState(null);
  const [selectedFile2, setSelectedFile2] = useState(null);

  const [figures, setFigures] = useState([]);
  const [titles, setTitles] = useState([]);

  const [loading, setLoading] = useState(false);
  const [bboxLoading, setBboxLoading] = useState(false);
  const [initLoading, setInitLoading] = useState(false);
  const [restoreLoading, setRestoreLoading] = useState(false);
  const [exitLoading, setExitLoading] = useState(false);

  const [error, setError] = useState("");
  const [uploadedFilename1, setUploadedFilename1] = useState("");
  const [uploadedFilename2, setUploadedFilename2] = useState("");

  const [preprocessMessage, setPreprocessMessage] = useState("");
  const [regionCount, setRegionCount] = useState(null);
  const [canShowResults, setCanShowResults] = useState(false);
  const [isZoomed, setIsZoomed] = useState(false);

  const [dragStart, setDragStart] = useState(null);
  const [dragCurrent, setDragCurrent] = useState(null);
  const [isDragging, setIsDragging] = useState(false);

  const fileInputRef1 = useRef(null);
  const fileInputRef2 = useRef(null);
  const firstImageRef = useRef(null);

  const resetSelection = () => {
    setDragStart(null);
    setDragCurrent(null);
    setIsDragging(false);
  };

  const resetToInitialPage = () => {
    setSelectedMode("");
    setSelectedFile1(null);
    setSelectedFile2(null);
    setFigures([]);
    setTitles([]);
    setLoading(false);
    setBboxLoading(false);
    setInitLoading(false);
    setRestoreLoading(false);
    setExitLoading(false);
    setError("");
    setUploadedFilename1("");
    setUploadedFilename2("");
    setPreprocessMessage("");
    setRegionCount(null);
    setCanShowResults(false);
    setIsZoomed(false);
    resetSelection();

    if (fileInputRef1.current) fileInputRef1.current.value = "";
    if (fileInputRef2.current) fileInputRef2.current.value = "";
  };

  const resetUploadStateOnly = () => {
    setFigures([]);
    setTitles([]);
    setUploadedFilename1("");
    setUploadedFilename2("");
    setPreprocessMessage("");
    setRegionCount(null);
    setCanShowResults(false);
    setIsZoomed(false);
    resetSelection();
  };

  const handleModeChange = (event) => {
    const value = event.target.value;

    setSelectedMode(value);
    setSelectedFile1(null);
    setSelectedFile2(null);
    setError("");
    resetUploadStateOnly();

    if (fileInputRef1.current) fileInputRef1.current.value = "";
    if (fileInputRef2.current) fileInputRef2.current.value = "";
  };

  const handleFile1Change = (event) => {
    const file = event.target.files?.[0] || null;

    setError("");
    resetUploadStateOnly();

    if (!file) {
      setSelectedFile1(null);
      return;
    }

    const fileName = file.name.toLowerCase();
    if (!fileName.endsWith(".h5")) {
      setError("The first file must be a valid .h5 file.");
      setSelectedFile1(null);
      if (fileInputRef1.current) fileInputRef1.current.value = "";
      return;
    }

    setSelectedFile1(file);
  };

  const handleFile2Change = (event) => {
    const file = event.target.files?.[0] || null;

    setError("");
    resetUploadStateOnly();

    if (!file) {
      setSelectedFile2(null);
      return;
    }

    const fileName = file.name.toLowerCase();
    if (!fileName.endsWith(".csv")) {
      setError("The second file must be a valid .csv file.");
      setSelectedFile2(null);
      if (fileInputRef2.current) fileInputRef2.current.value = "";
      return;
    }

    setSelectedFile2(file);
  };

  const clamp = (value, min, max) => {
    return Math.max(min, Math.min(max, value));
  };

  const getRelativeCoords = (event) => {
    const img = firstImageRef.current;
    if (!img) return null;

    const rect = img.getBoundingClientRect();

    const x = clamp(event.clientX - rect.left, 0, rect.width);
    const y = clamp(event.clientY - rect.top, 0, rect.height);

    return { x, y, rect };
  };

  const handleSubmit = async (event) => {
    event.preventDefault();

    if (!selectedMode) {
      setError("Please select either Spatial Gene Expression or Single Cell RNA-seq.");
      return;
    }

    if (!selectedFile1) {
      setError("Please choose the required .h5 file first.");
      return;
    }

    if (selectedMode === "RNA" && !selectedFile2) {
      setError("Please choose the UMAP .csv file for RNA mode.");
      return;
    }

    setLoading(true);
    setError("");
    setFigures([]);
    setTitles([]);
    setUploadedFilename1("");
    setUploadedFilename2("");
    setRegionCount(null);
    setCanShowResults(false);
    setIsZoomed(false);
    resetSelection();
    setPreprocessMessage("Preprocessing...");

    try {
      const file1 = selectedFile1;
      const file2 =
        selectedMode === "SPATIAL"
          ? new File(["dummy"], "dummy.h5", { type: "application/x-hdf5" })
          : selectedFile2;

      const formData = new FormData();
      formData.append("file1", file1);
      formData.append("file2", file2);
      formData.append("user_text", selectedMode);

      const startTime = performance.now();

      const response = await fetch(`${API_BASE}/run-h5`, {
        method: "POST",
        headers: {
          "ngrok-skip-browser-warning": "true",
        },
        body: formData,
      });

      const endTime = performance.now();

      const data = await response.json();

      const returnedCount =
        typeof data.count === "number" ? data.count : Number(data.count);

      console.log(`mode: ${selectedMode}`);
      console.log(`regions/cells: ${returnedCount}`);
      console.log(`preprocessing: ${(endTime - startTime).toFixed(2)} ms`);

      if (!response.ok) {
        throw new Error(data.detail || `Backend error: ${response.status}`);
      }

      setUploadedFilename1(data.filename1 || file1.name);
      setUploadedFilename2(
        selectedMode === "RNA" ? data.filename2 || file2.name : "dummy.h5"
      );

      setRegionCount(Number.isFinite(returnedCount) ? returnedCount : null);
      setPreprocessMessage(
        `Preprocessing completed and total number of regions/cells is ${
          Number.isFinite(returnedCount) ? returnedCount : "unknown"
        }`
      );
      setCanShowResults(true);

      setSelectedFile1(null);
      setSelectedFile2(null);

      if (fileInputRef1.current) fileInputRef1.current.value = "";
      if (fileInputRef2.current) fileInputRef2.current.value = "";
    } catch (err) {
      setError(err.message || "Something went wrong while uploading.");
      setPreprocessMessage("");
    } finally {
      setLoading(false);
    }
  };

  const loadInitialRender = async () => {
    setInitLoading(true);
    setError("");
    resetSelection();

    try {
      const startTime = performance.now();

      const response = await fetch(`${API_BASE}/init_render`, {
        method: "GET",
        headers: {
          "ngrok-skip-browser-warning": "true",
        },
      });

      const endTime = performance.now();
      console.log(`initial rendering: ${(endTime - startTime).toFixed(2)} ms`);

      const data = await response.json();

      if (!response.ok) {
        throw new Error(data.detail || `Backend error: ${response.status}`);
      }

      if (!data.images || !Array.isArray(data.images)) {
        throw new Error("Invalid response format. Expected { images: [] }");
      }

      const imageUrls = data.images.map(
        (imgBase64) => `data:image/png;base64,${imgBase64}`
      );

      setFigures(imageUrls);
      setTitles(Array.isArray(data.titles) ? data.titles : []);
      setIsZoomed(false);
    } catch (err) {
      setError(err.message || "Something went wrong while loading results.");
    } finally {
      setInitLoading(false);
    }
  };

  const restoreInitialRender = async () => {
    setRestoreLoading(true);
    setError("");
    resetSelection();

    try {
      const startTime = performance.now();

      const response = await fetch(`${API_BASE}/init_render`, {
        method: "GET",
        headers: {
          "ngrok-skip-browser-warning": "true",
        },
      });

      const endTime = performance.now();
      console.log(`initial rendering: ${(endTime - startTime).toFixed(2)} ms`);

      const data = await response.json();

      if (!response.ok) {
        throw new Error(data.detail || `Backend error: ${response.status}`);
      }

      if (!data.images || !Array.isArray(data.images)) {
        throw new Error("Invalid response format. Expected { images: [] }");
      }

      const imageUrls = data.images.map(
        (imgBase64) => `data:image/png;base64,${imgBase64}`
      );

      setFigures(imageUrls);
      setTitles(Array.isArray(data.titles) ? data.titles : []);
      setIsZoomed(false);
    } catch (err) {
      setError(err.message || "Something went wrong while restoring results.");
    } finally {
      setRestoreLoading(false);
    }
  };

  const submitBBox = async () => {
    const img = firstImageRef.current;
    if (!img || !dragStart || !dragCurrent) return;

    const displayWidth = img.clientWidth;
    const displayHeight = img.clientHeight;
    const naturalWidth = img.naturalWidth;
    const naturalHeight = img.naturalHeight;

    if (!displayWidth || !displayHeight || !naturalWidth || !naturalHeight) {
      setError("Image size is not available.");
      return;
    }

    const x1 = Math.min(dragStart.x, dragCurrent.x);
    const x2 = Math.max(dragStart.x, dragCurrent.x);
    const y1 = Math.min(dragStart.y, dragCurrent.y);
    const y2 = Math.max(dragStart.y, dragCurrent.y);

    if (Math.abs(x2 - x1) < 3 || Math.abs(y2 - y1) < 3) {
      resetSelection();
      return;
    }

    const scaleX = naturalWidth / displayWidth;
    const scaleY = naturalHeight / displayHeight;

    const payload = {
      x_start: x1 * scaleX,
      x_end: x2 * scaleX,
      y_start: y1 * scaleY,
      y_end: y2 * scaleY,
    };

    setBboxLoading(true);
    setError("");

    try {
      const startTime = performance.now();

      const response = await fetch(`${API_BASE}/submit-bbox`, {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          "ngrok-skip-browser-warning": "true",
        },
        body: JSON.stringify(payload),
      });

      const endTime = performance.now();
      console.log(`zoom rendering: ${(endTime - startTime).toFixed(2)} ms`);

      const data = await response.json();

      if (!response.ok) {
        throw new Error(data.detail || `Backend error: ${response.status}`);
      }

      if (!data.images || !Array.isArray(data.images)) {
        throw new Error("Invalid response format from bbox endpoint.");
      }

      const imageUrls = data.images.map(
        (imgBase64) => `data:image/png;base64,${imgBase64}`
      );

      setFigures(imageUrls);
      setIsZoomed(true);
      resetSelection();
    } catch (err) {
      setError(err.message || "Something went wrong while submitting bbox.");
    } finally {
      setBboxLoading(false);
    }
  };

  const handleExit = async () => {
    setExitLoading(true);
    setError("");

    try {
      const response = await fetch(`${API_BASE}/func_end`, {
        method: "GET",
        headers: {
          "ngrok-skip-browser-warning": "true",
        },
      });

      const data = await response.json();

      if (!response.ok) {
        throw new Error(data.detail || `Backend error: ${response.status}`);
      }

      if (data.status !== "ok") {
        throw new Error("Backend did not return success.");
      }

      resetToInitialPage();
    } catch (err) {
      setError(err.message || "Something went wrong while exiting.");
    } finally {
      setExitLoading(false);
    }
  };

  const handleMouseDown = (event) => {
    if (
      bboxLoading ||
      loading ||
      initLoading ||
      restoreLoading ||
      figures.length === 0 ||
      isZoomed
    ) {
      return;
    }

    const coords = getRelativeCoords(event);
    if (!coords) return;

    setDragStart({ x: coords.x, y: coords.y });
    setDragCurrent({ x: coords.x, y: coords.y });
    setIsDragging(true);
  };

  const handleMouseMove = (event) => {
    if (!isDragging || isZoomed) return;

    const coords = getRelativeCoords(event);
    if (!coords) return;

    setDragCurrent({ x: coords.x, y: coords.y });
  };

  const handleMouseUp = async () => {
    if (!isDragging || isZoomed) return;
    setIsDragging(false);
    await submitBBox();
  };

  const handleMouseLeave = async () => {
    if (!isDragging || isZoomed) return;
    setIsDragging(false);
    await submitBBox();
  };

  const selectionBox = (() => {
    if (!dragStart || !dragCurrent || isZoomed) return null;

    const left = Math.min(dragStart.x, dragCurrent.x);
    const top = Math.min(dragStart.y, dragCurrent.y);
    const width = Math.abs(dragCurrent.x - dragStart.x);
    const height = Math.abs(dragCurrent.y - dragStart.y);

    return { left, top, width, height };
  })();

  return (
    <div style={styles.page}>
      <div style={styles.container}>
        <h1 style={styles.title}>Gene Expression Figure Viewer</h1>

        {figures.length === 0 && (
          <form onSubmit={handleSubmit} style={styles.formColumn}>
            <div style={styles.modeSection}>
              <p style={styles.modeTitle}>Choose input type</p>

              <label style={styles.radioLabel}>
                <input
                  type="radio"
                  name="mode"
                  value="SPATIAL"
                  checked={selectedMode === "SPATIAL"}
                  onChange={handleModeChange}
                  disabled={loading || initLoading || exitLoading}
                />
                <span>Spatial Gene Expression</span>
              </label>

              <label style={styles.radioLabel}>
                <input
                  type="radio"
                  name="mode"
                  value="RNA"
                  checked={selectedMode === "RNA"}
                  onChange={handleModeChange}
                  disabled={loading || initLoading || exitLoading}
                />
                <span>Single Cell RNA-seq</span>
              </label>
            </div>

            {selectedMode && (
              <div style={styles.form}>
                <div style={styles.fileGroup}>
                  <label style={styles.fileLabel}>Upload H5 file</label>
                  <input
                    ref={fileInputRef1}
                    type="file"
                    accept=".h5,application/x-hdf5"
                    onChange={handleFile1Change}
                    style={styles.input}
                    disabled={loading || initLoading || exitLoading}
                  />
                </div>

                {selectedMode === "RNA" && (
                  <div style={styles.fileGroup}>
                    <label style={styles.fileLabel}>Upload TSNE CSV file</label>
                    <input
                      ref={fileInputRef2}
                      type="file"
                      accept=".csv,text/csv"
                      onChange={handleFile2Change}
                      style={styles.input}
                      disabled={loading || initLoading || exitLoading}
                    />
                  </div>
                )}

                <button
                  type="submit"
                  disabled={
                    loading ||
                    !selectedMode ||
                    !selectedFile1 ||
                    (selectedMode === "RNA" && !selectedFile2)
                  }
                  style={styles.button}
                >
                  {loading ? "Uploading..." : "Upload and Analyze"}
                </button>
              </div>
            )}
          </form>
        )}

        {selectedMode && figures.length === 0 && (
          <p style={styles.fileText}>
            Selected mode: <strong>{selectedMode}</strong>
          </p>
        )}

        {selectedFile1 && figures.length === 0 && (
          <p style={styles.fileText}>
            Selected file 1: <strong>{selectedFile1.name}</strong>
          </p>
        )}

        {selectedFile2 && figures.length === 0 && (
          <p style={styles.fileText}>
            Selected file 2: <strong>{selectedFile2.name}</strong>
          </p>
        )}

        {uploadedFilename1 && (
          <p style={styles.fileText}>
            Uploaded file 1: <strong>{uploadedFilename1}</strong>
          </p>
        )}

        {uploadedFilename2 && (
          <p style={styles.fileText}>
            Uploaded file 2: <strong>{uploadedFilename2}</strong>
          </p>
        )}

        {preprocessMessage && (
          <p style={styles.statusText}>{preprocessMessage}</p>
        )}

        {loading && <p style={styles.loadingText}>Preprocessing...</p>}

        {canShowResults && figures.length === 0 && (
          <div style={styles.centerRow}>
            <button
              type="button"
              onClick={loadInitialRender}
              disabled={initLoading}
              style={styles.button}
            >
              {initLoading ? "Loading results..." : "Show Results"}
            </button>
          </div>
        )}

        {figures.length > 0 && (
          <>
            <p style={styles.helpText}>
              Drag on the first image to select a region and update the figures.
            </p>

            {bboxLoading && (
              <p style={styles.loadingText}>Submitting selected region...</p>
            )}

            {isZoomed && (
              <div style={styles.centerRow}>
                <button
                  type="button"
                  onClick={restoreInitialRender}
                  disabled={restoreLoading}
                  style={styles.secondaryButton}
                >
                  {restoreLoading ? "Restoring..." : "Back to Initial Render"}
                </button>
              </div>
            )}

            <div style={styles.grid}>
              {figures.map((figure, index) => {
                const isOddNumberFigure = (index + 1) % 2 === 1;

                return (
                  <div
                    key={index}
                    style={{
                      ...styles.card,
                      ...(isOddNumberFigure
                        ? styles.oddCard
                        : styles.evenCard),
                    }}
                  >
                    {titles[index] && (
                      <div style={styles.figureTitle}>{titles[index]}</div>
                    )}

                    {index === 0 ? (
                      <div
                        style={styles.imageWrapper}
                        onMouseDown={handleMouseDown}
                        onMouseMove={handleMouseMove}
                        onMouseUp={handleMouseUp}
                        onMouseLeave={handleMouseLeave}
                      >
                        <img
                          ref={firstImageRef}
                          src={figure}
                          alt={titles[index] || `Figure ${index + 1}`}
                          style={{
                            ...styles.image,
                            cursor:
                              bboxLoading || isZoomed ? "default" : "crosshair",
                          }}
                          draggable={false}
                        />

                        {selectionBox && (
                          <div
                            style={{
                              ...styles.selectionBox,
                              left: `${selectionBox.left}px`,
                              top: `${selectionBox.top}px`,
                              width: `${selectionBox.width}px`,
                              height: `${selectionBox.height}px`,
                            }}
                          />
                        )}
                      </div>
                    ) : (
                      <img
                        src={figure}
                        alt={titles[index] || `Figure ${index + 1}`}
                        style={styles.image}
                      />
                    )}
                  </div>
                );
              })}
            </div>

            <div style={styles.exitSection}>
              <button
                type="button"
                onClick={handleExit}
                disabled={exitLoading}
                style={styles.exitButton}
              >
                {exitLoading ? "Exiting..." : "Exit"}
              </button>
            </div>
          </>
        )}

        {error && <p style={styles.error}>{error}</p>}
      </div>
    </div>
  );
}

const styles = {
  page: {
    minHeight: "100vh",
    backgroundColor: "#f5f7fb",
    padding: "30px",
    fontFamily: "Arial, sans-serif",
  },
  container: {
    maxWidth: "1400px",
    margin: "0 auto",
    backgroundColor: "#ffffff",
    borderRadius: "12px",
    padding: "24px",
    boxShadow: "0 4px 20px rgba(0,0,0,0.08)",
  },
  title: {
    marginBottom: "20px",
    textAlign: "center",
  },
  formColumn: {
    display: "flex",
    flexDirection: "column",
    gap: "18px",
    alignItems: "center",
    marginBottom: "16px",
  },
  modeSection: {
    display: "flex",
    flexDirection: "column",
    gap: "10px",
    alignItems: "flex-start",
    padding: "16px",
    border: "1px solid #ddd",
    borderRadius: "10px",
    minWidth: "320px",
    backgroundColor: "#fafafa",
  },
  modeTitle: {
    margin: 0,
    fontWeight: "bold",
  },
  radioLabel: {
    display: "flex",
    alignItems: "center",
    gap: "8px",
    cursor: "pointer",
  },
  form: {
    display: "flex",
    gap: "12px",
    alignItems: "center",
    justifyContent: "center",
    flexWrap: "wrap",
    marginBottom: "16px",
  },
  fileGroup: {
    display: "flex",
    flexDirection: "column",
    gap: "6px",
    alignItems: "flex-start",
  },
  fileLabel: {
    fontWeight: "bold",
    fontSize: "14px",
  },
  input: {
    padding: "8px",
  },
  button: {
    padding: "10px 16px",
    border: "none",
    borderRadius: "8px",
    cursor: "pointer",
    backgroundColor: "#1f6feb",
    color: "white",
    fontWeight: "bold",
  },
  secondaryButton: {
    padding: "10px 16px",
    border: "none",
    borderRadius: "8px",
    cursor: "pointer",
    backgroundColor: "#2da44e",
    color: "white",
    fontWeight: "bold",
  },
  exitButton: {
    padding: "12px 22px",
    border: "none",
    borderRadius: "8px",
    cursor: "pointer",
    backgroundColor: "#d1242f",
    color: "white",
    fontWeight: "bold",
  },
  centerRow: {
    display: "flex",
    justifyContent: "center",
    marginTop: "14px",
    marginBottom: "14px",
  },
  fileText: {
    textAlign: "center",
    marginBottom: "12px",
  },
  statusText: {
    textAlign: "center",
    marginBottom: "12px",
    color: "#1a1a1a",
    fontWeight: "bold",
  },
  helpText: {
    textAlign: "center",
    marginBottom: "12px",
    color: "#333",
    fontWeight: "bold",
  },
  loadingText: {
    textAlign: "center",
    marginBottom: "12px",
    color: "#1f6feb",
    fontWeight: "bold",
  },
  error: {
    color: "red",
    textAlign: "center",
    marginTop: "16px",
    marginBottom: "8px",
  },
  grid: {
    display: "grid",
    gridTemplateColumns: "repeat(12, 1fr)",
    gap: "16px",
    marginTop: "24px",
  },
  card: {
    backgroundColor: "#fafafa",
    border: "1px solid #ddd",
    borderRadius: "10px",
    padding: "12px",
    textAlign: "center",
  },
  oddCard: {
    gridColumn: "span 6",
  },
  evenCard: {
    gridColumn: "span 6",
  },
  figureTitle: {
    marginBottom: "8px",
    fontWeight: "bold",
  },
  imageWrapper: {
    position: "relative",
    display: "inline-block",
    width: "100%",
    userSelect: "none",
  },
  image: {
    width: "100%",
    height: "260px",
    objectFit: "contain",
    borderRadius: "6px",
    backgroundColor: "#fff",
    display: "block",
  },
  selectionBox: {
    position: "absolute",
    border: "2px solid #ff0000",
    backgroundColor: "rgba(255, 0, 0, 0.15)",
    pointerEvents: "none",
  },
  exitSection: {
    display: "flex",
    justifyContent: "center",
    marginTop: "28px",
  },
};