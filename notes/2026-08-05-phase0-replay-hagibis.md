# DR_REPORT_META
> **歸檔時間：** 2026-08-05 17:13:17 UTC
> **Agent ID：** tygtDc (DR)
> **報告類型：** research
> **目標 Collection：** dr_research

---

## 🔬 DR Reflection（研究後自我評估）

| 維度 | 自評 | 備註 |
|:-----|:----:|:-----|
| **信源質量** | ⬜ ⭐⭐⭐⭐⭐ | |
| **交叉驗證** | ⬜ ⭐⭐⭐⭐⭐ | |
| **分析深度** | ⬜ ⭐⭐⭐⭐⭐ | |
| **不確定性標註** | ⬜ ⭐⭐⭐⭐⭐ | |
| **整體信心** | ⬜ ⭐⭐⭐⭐⭐ | |

### 💡 做得好的
- 

### 🔧 可改進的
- 

### 🧠 關鍵洞察
- 

---
# Phase 0 Engine Replay — Hagibis 2019 驗證報告

**Date:** 2026-08-05  
**Agent:** tygtDc (DR)  
**Status:** ✅ Phase 0 完成，發現結構性 false positive 問題，已記錄 fix 方向

---

## 目標

用真實 ERA5 Hagibis 2019 數據（850hPa, 61 日, 696 frames）replay Ô-HAT mock engine，驗證：
1. Pipeline 技術可行性（3 算子 streaming + state machine）
2. Trigger 正確性（對比 Phase 2 已知結果）

---

## 技術路線

### 數據處理

Hagibis ERA5 netCDF 只有 `vo`（相對渦度）同 `d`（輻散），冇 u10/v10。需要 spectral Poisson solver 重建風場：

1. FFT along x → tridiagonal Helmholtz solver along y → 得 streamfunction ψ 同 velocity potential χ
2. u = -∂ψ/∂y + ∂χ/∂x, v = ∂ψ/∂x + ∂χ/∂y

### Engine 適配

| 算子 | 原設計 | Phase 0 適配 |
|------|--------|-------------|
| dH_curl | 從 u/v 計 vorticity | **直接用 ERA5 raw `vo`**（避開 Poisson 平滑） |
| θ₁ | 從 u/v 計 PCA | 用 reconstructed u/v（方向/朝向不受 magnitude 影響） |
| D_fold | 從 u/v 計 vorticity magnitude | 用 ERA5 raw `vo` magnitude |

### Cyclone Presence Filter

第一版 replay 發現 ERC_ACTIVE 幾乎長開（157/174 frames）— θ₁ range 15°–45° 係背景大氣正常變異。加入：

- `CYCLONE_VORT_MIN = 1e-4`：grid max \|vo\| > 1e-4 先當有 TC
- `DH_ABS_MIN = 1e-10`：dH_curl 必須過 noise floor

---

## 結果

### Performance

```
174 frames (subsample=4, original 696)
Grid: 101×201 (850hPa)
Poisson reconstruction: ~95s
Pipeline runtime: ~20s
Total: 115s
Throughput: ~1.5 frames/s (reconstruction dominates)
```

### Trigger Summary

| Trigger | Count | 評價 |
|---------|-------|------|
| ERC_ACTIVE | 156/174 | ❌ 太多，θ₁ range 15°–45° 在 24h window 內幾乎 always true |
| RI_PRECURSOR | 6 | ❌ 實際只有 1 次真正 RI（Oct 6-7） |
| IMPULSE_RESPONSE | 2 | ❌ Hagibis 非火山事件 |

### Operator 動態範圍

| Operator | Range | 評價 |
|----------|-------|------|
| dH_curl | −4.35e-05 ~ +4.35e-05 | ✅ ERA5 raw vo 產生物理合理值 |
| θ₁ | 4.7° ~ 38.9° | ✅ 有動態範圍，但非 TC 天氣都有高變異 |
| D_fold | 1.47 ~ 1.48 | ⚠️ 幾乎無變化（grid 太粗？） |

---

## 發現嘅結構性問題

### P1: ERC threshold 對連續監測太寬

θ₁ 在 24h window 內 range 15°–45° 係大氣正常變異，非 ERC 專有信號。61 日 continuous run 幾乎 always trigger。

**Fix:** 加 TC-track spatial coincidence（IBTrACS best track ±200km）+ 連續 ≥12h persistence。

### P2: % surge 對近零 dH 太敏感

dH 從 10⁻¹⁰ 升到 10⁻⁹ = 1000% surge，但絕對值仍然好細。RI/IMPULSE trigger 需要 absolute magnitude threshold。

**Fix:** 已完成（`RI_DH_ABS_SURGE = 1e-10`），但需要根據 Phase 2 實測值 calibrate。

### P3: Poisson reconstruction 平滑效應

Spectral Poisson solver 將 vo/d 重建為 u/v 時會丟失高波數分量，導致 θ₁ 敏感度下降。直接用 ERA5 u10/v10（如有）可以避開呢個問題。

**Fix:** 正式部署用 GFS u10/v10（唔使 reconstruction）。

### P4: State machine 需要 spatial context

目前 state machine 只睇 operator time series，唔知 vortex 喺邊。需要 spatial coincidence filter。

**Fix:** 對接 IBTrACS / automated vortex tracker。

---

## Phase 0 結論

**Engine pipeline 技術可行。** 3 算子 streaming + state machine 架構正確。但 continuous monitoring 需要以下 fix 先可以上線：

| Priority | Fix | Effort |
|----------|-----|--------|
| 🔴 P1 | IBTrACS track coincidence filter | 1 session |
| 🔴 P2 | Absolute dH threshold calibration | 0.5 session |
| 🟡 P3 | ERA5 u10/v10 直接數據源 | 已有（dolphin-watch） |
| 🟡 P4 | ERC persistence（≥12h sustained） | 0.5 session |

**Phase 0.5 路線：** 加 IBTrACS filter + persistence → replay Dolphin/Hagibis/Tonga → 確認 trigger 正確 → Docker 封裝 → Phase 1 GFS live feed。

---

## 產出

- Engine: `projects/ohat_engine/ohat_engine_mock.py`（v0.2, 含 cyclone filter）
- Replay: `projects/ohat_engine/replay_hagibis.py`
- Output: `projects/ohat_engine/output/replay_hagibis_2019.json`
