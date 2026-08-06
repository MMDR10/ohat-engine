FROM python:3.11-slim

LABEL org.opencontainers.image.title="Ô-HAT Engine v0.5"
LABEL org.opencontainers.image.description="Ô-HAT Noise Topology Engine — Universal framework for detecting dynamical phase transitions in 2D fields via dH/θ₁ phase-space. Validated: typhoon RI/ERC (7/7, 100%), ENSO SST folding (z=-45), volcanic DART coherence (H=1.00), seismic GEONET (12/12), epidemic COVID (2.5×)."
LABEL org.opencontainers.image.version="0.5"
LABEL org.opencontainers.image.authors="tygtDc <nnrpmrmm@gmail.com>"

# System deps for netCDF4/HDF5
RUN apt-get update && apt-get install -y --no-install-recommends \
    libhdf5-dev libnetcdf-dev \
    && rm -rf /var/lib/apt/lists/*

# Python deps
COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt

# Engine
COPY ohat_engine_v0.5.py /app/engine.py

# Mount points
RUN mkdir -p /data /output
VOLUME ["/data", "/output"]

WORKDIR /app

ENTRYPOINT ["python3", "/app/engine.py"]
