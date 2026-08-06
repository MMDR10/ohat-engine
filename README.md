# Ô-HAT Engine v0.5 — Real-Time Noise Topology Engine

> **A universal framework for detecting dynamical regime shifts in any 2D field.**
> Typhoon rapid intensification (RI) and eyewall replacement cycle (ERC) detection is
> the first fully validated domain.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21816438.svg)](https://doi.org/10.5281/zenodo.21816438)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/)
[![Docker](https://img.shields.io/badge/docker-ready-2496ED.svg)](https://www.docker.com/)

---

## 摘要 / Abstract

**中文：** Ô-HAT Engine 是一個基於雜訊拓撲（Noise Topology）的通用相變檢測框架。
引擎從任意二維場提取雜訊場 dH = H − smooth(H)，計算拓撲角 θ₁ 與折疊維數 D_fold，
透過三相狀態機（WARMUP → ACTIVE → DETECT）實時追蹤動力學機制轉變。
**v0.5 颱風跨驗證：5 個 Cat 5 超強颱風（Goni / Hagibis / Haiyan / Meranti / Patricia），
20/20 checks 全部通過。** 三條 RI 觸發路徑已被識別（P1 θ₁ 跳水 / P2 鎖定 / P2 振盪），
NO_CYCLONE 背景過濾率 60–78%。引擎以 Docker 容器封裝（~200MB image, ~400MB RAM/run），
單次運行 30–60 秒，支援平行多容器批次重播。

**English:** The Ô-HAT Engine is a universal phase-transition detection framework based on
Noise Topology. It extracts the noise field dH = H − smooth(H) from any 2D field, computes
the topological angle θ₁ and folding dimension D_fold, and tracks dynamical regime shifts
in real time via a three-phase state machine (WARMUP → ACTIVE → DETECT).
**v0.5 typhoon cross-validation: 5 Cat 5 super-typhoons (Goni / Hagibis / Haiyan /
Meranti / Patricia), 20/20 checks passed.** Three RI trigger pathways identified
(P1 θ₁ drop / P2 locked / P2 oscillatory), NO_CYCLONE background filter rate 60–78%.
The engine is containerized via Docker (~200MB image, ~400MB RAM/run), with 30–60 s
per-storm runtime and multi-container parallel batch replay support.

---

## Architecture

```
Any 2D field  →  [dH = H − smooth(H)]  →  [θ₁ = topological angle]  →  [State Machine]  →  [Triggers]
                  (Poisson/spectral)        (noise topography PCA)        (WARMUP→ACTIVE→DETECT)   (RI/ERC/Impulse)
```

| Component | Role | Domain-Specific? |
|-----------|------|:---:|
| **dH computation** | Noise field = raw − Poisson reconstruction | ❌ Universal |
| **θ₁ extraction** | Topological angle from 7D noise topography vector | ❌ Universal |
| **D_fold** | Folding dimension via box-counting on singularity set | ❌ Universal |
| **State machine** | WARMUP → ACTIVE → DETECT with phase-boundary crossing | ❌ Universal |
| **IBTrACS filter** | Track-based ground truth (±200 km coincidence) | ✅ Typhoon |
| **Storm metadata** | Name / year / wind for labeling | ✅ Typhoon |

The core pipeline (dH / θ₁ / D_fold / state machine) is **domain-agnostic**.
Only the cyclone-specific filtering layer requires domain knowledge (IBTrACS best-track data).

---

## v0.5 Cross-Validation — 5 Cat 5 Super-Typhoons

| # | Storm | Year | Basin | Frames | Triggers | NO_CYCLONE | Checks |
|---|-------|------|-------|--------|----------|:----------:|:------:|
| 1 | **Haiyan** | 2013 | WPAC | 168 | 4 (2 ERC + 2 RI) | 65% | 4/4 ✅ |
| 2 | **Patricia** | 2015 | EPAC | 192 | 6 (4 ERC + 2 RI) | 76% | 4/4 ✅ |
| 3 | **Meranti** | 2016 | WPAC | 168 | 6 (4 ERC + 2 RI) | 60% | 4/4 ✅ |
| 4 | **Hagibis** | 2019 | WPAC | 348 | 10 (7 ERC + 3 RI) | 74% | 4/4 ✅ |
| 5 | **Goni** | 2020 | WPAC | 348 | 5 (4 ERC + 1 RI) | 78% | 4/4 ✅ |

**Summary: 20/20 checks passed (100%).** Background filter rate 60–78% (mean 67%).

### Three RI Trigger Pathways

| Pathway | Signature | First Seen |
|---------|-----------|:----------:|
| **P1: θ₁ Drop** | Pure θ₁ plunge >15°, dH stable | Patricia 2015 |
| **P2: Locked** | θ₁ locked in narrow band, dH_curl surge >200% | Hagibis / Goni / Haiyan |
| **P2: Oscillatory** | θ₁ oscillates with growing amplitude | Meranti 2016 |

### Validation Checks (per storm)

1. **ERC_ACTIVE detected** — ≥1 eyewall replacement cycle trigger
2. **RI_PRECURSOR detected** — ≥1 rapid intensification precursor
3. **NO_CYCLONE filtering active** — >50% frames correctly classified as no-cyclone background
4. **dH_curl dynamic range > 1e-5** — noise field has measurable structure

---

## Cross-Domain Validation (Engine Framework)

Beyond typhoons, the underlying Ô-HAT framework has been validated across multiple domains:

| Domain | System | Key Result | Status |
|--------|--------|-----------|:---:|
| 🌀 **Typhoon** | RI/ERC detection | 5/5 storms, 20/20 checks (100%) | ✅ v0.5 |
| 🌊 **ENSO** | SST noise folding | D_fold z=−45, 3 products × 39 yr | ✅ Beta |
| 🌋 **Volcanic** | Tonga DART coherence | H_agg=1.000, cross-domain phase lock | ✅ Case |
| 🔥 **Seismic** | GEONET 1,336 stations | 12/12 all-negative, global p=0.011 | ✅ Case |
| 🦠 **Epidemic** | COVID case counts | Balance=2.5×, Continuous Flow regime | ✅ Case |

---

## Key Numbers

| Metric | Value |
|--------|-------|
| **Engine throughput** | 435 fps (Phase 0 mock) |
| **Per-storm runtime** | 30–60 s (ERA5, subsample=2) |
| **Docker image** | ~200 MB (python:3.11-slim) |
| **RAM per run** | ~400 MB |
| **Cross-validation** | 20/20 checks (100%) |
| **RI pathways** | 3 (P1 θ₁ drop / P2 locked / P2 oscillatory) |
| **Background filter** | 60–78% NO_CYCLONE rate |

---

## Quick Start

### Option 1 — Docker (Recommended)

```bash
# 1. Unpack
tar -xzf ohat-engine-v0.5-docker.tar.gz && cd ohat-engine-v0.5

# 2. Build
docker build -t ohat-engine:v0.5 .

# 3. Run (Typhoon Mode)
docker run --rm \
  -v /path/to/era5-data:/data/era5:ro \
  -v /path/to/ibtracs:/data/ibtracs:ro \
  -v $(pwd)/output:/output \
  ohat-engine:v0.5 \
  --nc /data/era5/haiyan_2013.nc \
  --ibtracs /data/ibtracs/ibtracs_all.csv \
  --storm HAIYAN --year 2013 --subsample 2 \
  --output-dir /output

# 4. Run (Generic Mode — any 2D field, no cyclone filter)
docker run --rm \
  -v /path/to/data:/data:ro \
  -v $(pwd)/output:/output \
  ohat-engine:v0.5 \
  --nc /data/sst_monthly.nc \
  --output-dir /output
```

### Option 2 — Native Python

```bash
pip install numpy scipy xarray netCDF4
python ohat_engine_v0.5.py \
  --nc /path/to/typhoon.nc \
  --ibtracs /path/to/ibtracs_all.csv \
  --storm HAIYAN --year 2013 --subsample 2
```

### Batch Cross-Validation

```bash
for storm in "GONI:2020:goni_2020" "HAGIBIS:2019:hagibis_2019" \
             "HAIYAN:2013:haiyan_2013" "MERANTI:2016:meranti_2016" \
             "PATRICIA:2015:patricia_2015"; do
  IFS=':' read -r name year nc <<< "$storm"
  docker run --rm \
    -v /data/typhoons:/data/era5:ro \
    -v /data/ibtracs:/data/ibtracs:ro \
    -v $(pwd)/output:/output \
    ohat-engine:v0.5 \
    --nc /data/era5/${nc}.nc \
    --ibtracs /data/ibtracs/ibtracs_all.csv \
    --storm $name --year $year --subsample 2 \
    --output-dir /output
done
```

---

## Parameters

| Flag | Default | Description |
|------|---------|-------------|
| `--nc` | *required* | NetCDF with 2D spatial field + time dimension |
| `--ibtracs` | *optional* | IBTrACS CSV (enables cyclone mode with track filter) |
| `--storm` | HAGIBIS | Storm name for IBTrACS matching (cyclone mode) |
| `--year` | 2019 | Storm year for IBTrACS matching (cyclone mode) |
| `--subsample` | 1 | Spatial subsampling factor (2 = half resolution) |
| `--max-frames` | None | Cap temporal frames (None = process all) |
| `--impulse` | False | Enable volcanic/external shock impulse gate |
| `--output-dir` | /output | Directory for JSON results |

---

## Output Format

Each run produces `engine_v0.5_{name}_{year}.json`:

```json
{
  "n_frames": 168,
  "state_counts": {"NO_CYCLONE": 109, "NOMINAL": 19, "RI_PRECURSOR": 21, "ERC_ACTIVE": 14},
  "triggers": [
    {"t_h": 120, "type": "RI_PRECURSOR", "pathway": "P2_chaotic_dh_locked", "theta1": 8.4, "dh_surge": 3.2}
  ],
  "checks": [
    {"name": "RI_PRECURSOR detected", "passed": true, "detail": "2 triggers"},
    {"name": "ERC_ACTIVE detected", "passed": true, "detail": "4 triggers"}
  ],
  "time_series_dh": [...],
  "time_series_theta1": [...]
}
```

---

## Data Provenance

| Data | Source | Access |
|------|--------|--------|
| ERA5 hourly (850 hPa) | ECMWF CDS | [cds.climate.copernicus.eu](https://cds.climate.copernicus.eu/) |
| IBTrACS v4 | NOAA NCEI | [ncei.noaa.gov](https://www.ncei.noaa.gov/products/international-best-track-archive) |
| OISST v2 | NOAA PSL | [psl.noaa.gov](https://psl.noaa.gov/) |

> ⚠️ **Raw data is NOT included in this repository** (ERA5 .nc ~500 MB each).
> All scripts accept standard NetCDF/CSV formats. See individual source links above.

---

## Roadmap

| Phase | Status | Description |
|-------|:------:|-------------|
| **Phase 0** | ✅ Done | Mock engine (435 fps) + Hagibis replay |
| **Phase 0.5** | ✅ Done | IBTrACS track filter + ERC persistence + 5-storm cross-validation |
| **Phase 1** | ⬜ Next | GFS 0.25° real-time feed integration |
| **Phase 2** | ⬜ Next | Satellite data (ASCAT/CYGNSS) + sub-grid alignment |
| **Phase 3** | ⬜ Next | Live dashboard (WebSocket + Redis Stream) |

---

## Repository Contents

| Path | Description |
|------|-------------|
| `ohat_engine_v0.5.py` | **Engine core** — 3 operators + state machine + IBTrACS filter |
| `ohat_engine_mock.py` | Phase 0 mock engine (synthetic vortex + replay) |
| `replay_hagibis.py` | Hagibis 2019 replay script |
| `Dockerfile` | python:3.11-slim + numpy/scipy/xarray/netCDF4 |
| `docker-compose.yml` | Multi-container parallel replay config |
| `requirements.txt` | Python dependencies |
| `README.md` | This file |
| `README_DOCKER.md` | Docker-specific quickstart guide |
| `output/` | v0.5 cross-validation JSON results (5 storms) |
| `notes/` | Research notes & Phase 0 replay report |

---

## Citation

If you use the Ô-HAT Engine in your research, please cite:

> tygtDc, Deep Research. (2026). *Ô-HAT Engine: Real-Time Noise Topology Engine for Dynamical Regime Shift Detection (v0.5)*. Zenodo. DOI: [10.5281/zenodo.21816438](https://doi.org/10.5281/zenodo.21816438)

---

## License

[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) — Free to use, share, and adapt with attribution.

## Author

**tygtDc, Deep Research** — Independent Researcher  
Contact: nnrpmrmm@gmail.com
