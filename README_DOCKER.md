# Ô-HAT Engine v0.5 — Docker Quickstart

🌀 **Noise Topology Engine** — Detect phase transitions in any 2D field via dH/θ₁ phase-space.

Not a typhoon tool. A **universal framework** for detecting dynamical regime shifts.
Typhoon is the first validated domain; the engine is domain-agnostic.

## Architecture

```
Any 2D field → [dH = H − smooth(H)] → [θ₁ = topological angle] → [State Machine] → [Triggers]
```

| Component | What It Does | Domain-Specific? |
|-----------|-------------|:---:|
| dH computation | Compute noise field = raw − Poisson reconstruction | ❌ Universal |
| θ₁ extraction | Topological angle from noise topography vector | ❌ Universal |
| State machine | WARMUP → ACTIVE → DETECT states | ❌ Universal |
| Trigger detection | Phase-space boundary crossing events | ❌ Universal |
| IBTrACS filter | Track-based ground truth for cyclones | ✅ Typhoon |
| Storm metadata | Name/year/wind for labeling | ✅ Typhoon |

## Validated Domains

| Domain | System | Key Result | Status |
|--------|--------|-----------|:---:|
| 🌀 **Typhoon** | RI/ERC detection | 7/7 storms, 28/28 checks (100%) | ✅ v0.5 |
| 🌊 **ENSO** | SST noise folding | D_fold z=−45, 3 products, 39yr | ✅ Beta |
| 🌋 **Volcanic** | Tonga DART coherence | H_agg=1.000, cross-domain lock | ✅ Case |
| 🔥 **Seismic** | GEONET 1336 stations | 12/12 all-negative, global p=0.011 | ✅ Case |
| 🦠 **Epidemic** | COVID case counts | Balance=2.5×, Continuous Flow | ✅ Case |

## 1. Build

```bash
tar -xzf ohat-engine-v0.5-docker.tar.gz && cd ohat-engine-v0.5
docker build -t ohat-engine:v0.5 .
```

## 2. Run — Typhoon Mode

```bash
docker run --rm \
  -v /path/to/era5-data:/data/era5:ro \
  -v /path/to/ibtracs:/data/ibtracs:ro \
  -v $(pwd)/output:/output \
  ohat-engine:v0.5 \
  --nc /data/era5/haiyan_2013.nc \
  --ibtracs /data/ibtracs/ibtracs_all.csv \
  --storm HAIYAN --year 2013 --subsample 2 \
  --output-dir /output
```

## 3. Run — Generic Mode (any 2D field)

The engine works on any NetCDF with 2D spatial + temporal dims:

```bash
docker run --rm \
  -v /path/to/data:/data:ro \
  -v $(pwd)/output:/output \
  ohat-engine:v0.5 \
  --nc /data/sst_monthly.nc \
  --output-dir /output
```

When no `--ibtracs` is given, the engine runs in pure noise-topology mode:
dH computation, θ₁ tracking, state machine transitions — no cyclone-specific filtering.

## 4. Batch — 7 Typhoon Cross-Validation

```bash
for storm in \
  "GONI:2020:goni_2020" \
  "HAGIBIS:2019:hagibis_2019" \
  "HAIYAN:2013:haiyan_2013" \
  "MERANTI:2016:meranti_2016" \
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

## Parameters

| Flag | Default | Description |
|------|---------|-------------|
| `--nc` | *required* | NetCDF with 2D spatial field |
| `--ibtracs` | *optional* | IBTrACS CSV (enables cyclone mode) |
| `--storm` | HAGIBIS | Storm name (cyclone mode only) |
| `--year` | 2019 | Storm year (cyclone mode only) |
| `--subsample` | 1 | Spatial subsampling factor |
| `--max-frames` | None | Cap temporal frames |
| `--impulse` | False | Enable volcanic/external shock gate |
| `--output-dir` | /output | Results directory |

## Output

Each run produces `engine_v0.5_{name}_{year}.json`:

```json
{
  "n_frames": 168,
  "state_counts": {"NO_CYCLONE": 109, "NOMINAL": 19, "RI_PRECURSOR": 21, ...},
  "triggers": [
    {"t_h": 120, "type": "RI_PRECURSOR", "pathway": "P2_chaotic_dh_locked", ...}
  ],
  "checks": {"RI_PRECURSOR": true, "ERC_ACTIVE": true, ...},
  "time_series_dh": [...],
  "time_series_theta1": [...]
}
```

## Resource

| Metric | Value |
|--------|-------|
| Image size | ~200MB |
| RAM per run | ~400MB |
| Runtime (subsample=2) | ~30–60s per storm |
| Dependencies | numpy, scipy, xarray, netCDF4 (auto-installed) |
