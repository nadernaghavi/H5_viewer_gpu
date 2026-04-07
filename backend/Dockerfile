FROM rapidsai/rapidsai-core:23.06-cuda11.8-runtime-ubuntu22.04-py3.10

WORKDIR /app

SHELL ["/bin/bash", "-c"]

COPY requirements.txt /app/requirements.txt

RUN source /opt/conda/etc/profile.d/conda.sh && \
    conda activate rapids && \
    pip install --no-cache-dir -r /app/requirements.txt && \
    pip install --no-cache-dir rapids-singlecell

COPY . /app

EXPOSE 8000

# 🔥 THIS LINE FIXES EVERYTHING
ENTRYPOINT ["/bin/bash", "-lc"]

CMD ["source /opt/conda/etc/profile.d/conda.sh && conda activate rapids && exec uvicorn app:app --host 0.0.0.0 --port 8000"]