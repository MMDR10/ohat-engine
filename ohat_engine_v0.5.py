#!/usr/bin/env python3
"""
🧪 Ô-HAT Real-Time Engine — v0.5 (Phase 0.5 Fixes)
====================================================
Phase 0.5 upgrades over v0.2 mock engine:

  P1 ✅ IBTrACS track coincidence filter (±200km, ±3h)
  P2 ✅ Absolute dH thresholds calibrated from Phase 2 S-Curve
  P4 ✅ ERC persistence ≥12h sustained before trigger
  P3 ⚠️  uv-reconstruction retained (no direct u10/v10 for most typhoons)

Key calibration values (Phase 2 S-Curve):
  - Mature TC: dH_curl < −1.5e-5 s⁻¹
  - GFS saturation floor: ~−2.4e-5 s⁻¹
  - Hagibis dH_curl range: −4.35e-5 to +4.35e-5 (ERA5 vo-based)

Usage:
  python ohat_engine_v0.5.py --case hagibis
  python ohat_engine_v0.5.py --case dolphin
  python ohat_engine_v0.5.py --nc data/typhoons/hagibis_2019.nc --ibtracs projects/cyclone/data/ibtracs_all.csv
"""
import sys, json, os, argparse, time, csv
import numpy as np
from scipy import stats
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple
from datetime import datetime, timedelta

OUT_DIR = Path(__file__).parent / 'output'
OUT_DIR.mkdir(exist_ok=True)

# ═══════════════════════════════════════════════════════════════
# 1. CONFIG — Phase 2 Calibrated + Phase 0.5 Fixes
# ═══════════════════════════════════════════════════════════════

@dataclass
class EngineConfig:
    # ── Spatial ──
    core_radius_km: float = 30.0
    shell_radius_km: float = 100.0
    dx_km: float = 0.25  # ERA5 0.25°
    
    # ── θ₁ window ──
    theta_window_h: float = 6.0
    theta_fallback: float = 90.0
    
    # ── D_fold window ──
    dfold_window_frames: int = 6
    dfold_threshold_pct: float = 80
    
    # ── Cyclone presence pre-filter (from v0.2) ──
    CYCLONE_VORT_MIN: float = 1e-4   # max |vo| in domain > 1e-4 = TC present
    DH_ABS_MIN: float = 1e-10        # |dH_curl| noise floor
    
    # ── P2 Fix: Calibrated absolute dH thresholds ──
    # Phase 2 S-Curve: mature TC dH < −1.5e-5, saturation ~−2.4e-5
    # For RI detection, dH_curl must enter mature zone (not just % surge)
    DH_MATURE_THRESHOLD: float = -1.5e-5   # dH below this = organized TC
    DH_RI_ABS_DROP: float = 0.5e-5          # absolute dH drop for RI (>0.5e-5)
    DH_IMPULSE_ABS_DROP: float = 0.3e-5     # absolute dH drop for impulse
    
    # ── RI Detection: Two Pathways ──
    # Pathway 1 (Dolphin-style): θ₁ drops BEFORE RI → precursor
    RI_THETA_DROP_DEG: float = 15.0
    RI_THETA_WINDOW_H: float = 15.0
    # Pathway 2 (Chaotic RI, A+ Fix): dH plummets + θ₁ volatile → rapid organization
    RI_DH_ABS_PLUMMET: float = 5.0e-5       # dH_curl drops by >5e-5 in window
    RI_THETA_VOLATILITY_DEG: float = 8.0    # θ₁ std > 8° = oscillatory chaos (Hagibis)
    RI_THETA_LOCKED_MIN_DEG: float = 35.0   # θ₁ sustained > 35° = locked disorder (Goni)
    RI_CHAOTIC_WINDOW_H: float = 24.0       # window for chaotic RI detection
    
    # ── ERC Oscillation (Phase 2 validated: Hagibis 15°–45°) ──
    ERC_THETA_OSC_MIN: float = 15.0
    ERC_THETA_OSC_MAX: float = 45.0
    ERC_WINDOW_H: float = 24.0
    ERC_PERSISTENCE_H: float = 12.0         # P4 Fix: ≥12h sustained
    
    # ── Impulse Response (Phase 2 validated: Tonga θ₁↓10.5° in 6h) ──
    # IMPULSE mode is ONLY for volcanic/external shock cases — gated by flag
    IMPULSE_ENABLED: bool = False
    IMPULSE_THETA_DROP: float = 10.0
    IMPULSE_WINDOW_H: float = 6.0
    IMPULSE_DH_SURGE: float = 100.0
    
    # ── P1 Fix: IBTrACS spatial coincidence ──
    IBTRACS_MAX_DIST_KM: float = 200.0      # max distance from IBTrACS position
    IBTRACS_MAX_DT_H: float = 3.0           # max time offset from IBTrACS time


cfg = EngineConfig()


# ═══════════════════════════════════════════════════════════════
# 2. IBTrACS TRACK LOADER (P1 Fix)
# ═══════════════════════════════════════════════════════════════

def load_ibtracs_track(csv_path: str, storm_name: str, year: int) -> Optional[List[Dict]]:
    """
    Load IBTrACS best-track for a specific storm.
    Returns list of dicts: {time: datetime, lat: float, lon: float, wind: float, pres: float}
    """
    track = []
    with open(csv_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = (row.get('NAME') or '').strip().upper()
            iso_time = (row.get('ISO_TIME') or '').strip()
            if not iso_time or not name:
                continue
            if name != storm_name.upper():
                continue
            try:
                t = datetime.fromisoformat(iso_time)
            except ValueError:
                continue
            if t.year != year:
                continue
            
            lat_str = row.get('LAT', '') or row.get('USA_LAT', '')
            lon_str = row.get('LON', '') or row.get('USA_LON', '')
            wind_str = row.get('WMO_WIND', '') or row.get('USA_WIND', '')
            pres_str = row.get('WMO_PRES', '') or row.get('USA_PRES', '')
            
            try:
                lat = float(lat_str)
                lon = float(lon_str)
            except (ValueError, TypeError):
                continue
            
            wind = float(wind_str) if wind_str and wind_str.strip() else None
            pres = float(pres_str) if pres_str and pres_str.strip() else None
            
            track.append({
                'time': t,
                'lat': lat,
                'lon': lon,
                'wind': wind,
                'pres': pres
            })
    
    track.sort(key=lambda x: x['time'])
    return track if track else None


def find_ibtracs_match(t_hours: float, t0: datetime, cy: float, cx: float,
                        lats: np.ndarray, lons: np.ndarray,
                        track: List[Dict]) -> Optional[Dict]:
    """
    Check if current frame has a matching IBTrACS position within spatial+temporal window.
    Returns the matched IBTrACS point or None.
    
    Parameters:
      t_hours: hours since ERA5 start
      t0: ERA5 start datetime
      cy, cx: vortex center in PIXEL coordinates
      lats, lons: grid coordinate arrays (1D)
      track: IBTrACS track list
    """
    frame_time = t0 + timedelta(hours=float(t_hours))
    
    # Convert pixel center to lat/lon
    cy_int = int(np.clip(cy, 0, len(lats) - 1))
    cx_int = int(np.clip(cx, 0, len(lons) - 1))
    center_lat = float(lats[cy_int])
    center_lon = float(lons[cx_int])
    
    best_dist = float('inf')
    best_point = None
    
    for pt in track:
        dt_h = abs((frame_time - pt['time']).total_seconds() / 3600)
        if dt_h > cfg.IBTRACS_MAX_DT_H:
            continue
        
        # Haversine distance
        dlat = np.radians(pt['lat'] - center_lat)
        dlon = np.radians(pt['lon'] - center_lon)
        a = np.sin(dlat/2)**2 + np.cos(np.radians(center_lat)) * np.cos(np.radians(pt['lat'])) * np.sin(dlon/2)**2
        dist_km = 6371.0 * 2 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))
        
        if dist_km < best_dist:
            best_dist = dist_km
            best_point = pt
    
    if best_dist <= cfg.IBTRACS_MAX_DIST_KM and best_point is not None:
        return {
            'dist_km': float(best_dist),
            'dt_h': float(abs((frame_time - best_point['time']).total_seconds() / 3600)),
            'ibtracs_lat': best_point['lat'],
            'ibtracs_lon': best_point['lon'],
            'ibtracs_wind': best_point['wind'],
            'ibtracs_pres': best_point['pres'],
        }
    # Track loaded but no match → return empty dict (distinguishes from "no track")
    return {}


# ═══════════════════════════════════════════════════════════════
# 3. CORE OPERATORS (unchanged from v0.2)
# ═══════════════════════════════════════════════════════════════

def compute_dh_curl_vo(vo_frame: np.ndarray, core_mask: np.ndarray, 
                        shell_mask: np.ndarray) -> float:
    """dH_curl from raw ERA5 vorticity (avoids Poisson smoothing)."""
    core_vo = float(np.mean(vo_frame[core_mask])) if core_mask.sum() > 0 else 0.0
    shell_vo = float(np.mean(vo_frame[shell_mask])) if shell_mask.sum() > 0 else 0.0
    return shell_vo - core_vo


def compute_theta1(u: np.ndarray, v: np.ndarray, core_mask: np.ndarray) -> float:
    """θ₁ = arctan(√(λ_min/λ_max)) from PCA of (u,v) in core region."""
    u_c = u[core_mask]; v_c = v[core_mask]
    if len(u_c) < 4:
        return cfg.theta_fallback
    
    X = np.column_stack([u_c, v_c])
    X -= X.mean(axis=0)
    cov = (X.T @ X) / (len(X) - 1)
    
    try:
        ev = np.linalg.eigvalsh(cov)
    except np.linalg.LinAlgError:
        return cfg.theta_fallback
    
    if ev[1] < 1e-30:
        return 0.0
    ratio = max(ev[0], 0) / ev[1]
    return float(np.degrees(np.arctan(np.sqrt(ratio))))


def compute_dfold(field: np.ndarray, threshold_pct: float = 80) -> float:
    """D_fold: box-counting dimension of singularity set."""
    ny, nx = field.shape
    thresh = np.percentile(np.abs(field), threshold_pct)
    mask = np.abs(field) >= thresh
    
    y_idx, x_idx = np.where(mask)
    if len(y_idx) < 10:
        return float('nan')
    
    max_box = max(nx, ny) // 2
    sizes = [2**k for k in range(0, int(np.log2(max_box)) + 1)]
    sizes = [s for s in sizes if s <= max_box]
    
    counts = []
    for s in sizes:
        bins_x = np.arange(0, nx + s, s)
        bins_y = np.arange(0, ny + s, s)
        h, _, _ = np.histogram2d(x_idx, y_idx, bins=[bins_x, bins_y])
        counts.append(np.sum(h > 0))
    
    valid = [i for i, c in enumerate(counts) if c > 0]
    if len(valid) < 3:
        return float('nan')
    
    log_s = np.log([1.0 / sizes[i] for i in valid])
    log_n = np.log([counts[i] for i in valid])
    slope, _, _, _, _ = stats.linregress(log_s, log_n)
    return float(slope)


def make_masks(ny: int, nx: int, cx: float, cy: float) -> Tuple[np.ndarray, np.ndarray]:
    """Core and shell masks centered at (cx, cy) in pixel coordinates."""
    y, x = np.meshgrid(np.arange(ny), np.arange(nx), indexing='ij')
    r = np.sqrt((x - cx)**2 + (y - cy)**2)
    core_px = max(3, cfg.core_radius_km / (cfg.dx_km * 111.0))
    shell_px = max(6, cfg.shell_radius_km / (cfg.dx_km * 111.0))
    core_mask = r <= core_px
    shell_mask = (r > core_px) & (r <= shell_px)
    return core_mask, shell_mask


def find_center(u: np.ndarray, v: np.ndarray) -> Tuple[float, float]:
    """Vortex center: centroid of lowest 20% wind speed region."""
    speed = np.sqrt(u**2 + v**2)
    thresh = np.percentile(speed, 20)
    mask = speed <= thresh
    if mask.sum() == 0:
        return float(u.shape[0] // 2), float(u.shape[1] // 2)
    y_idx, x_idx = np.where(mask)
    return float(np.mean(y_idx)), float(np.mean(x_idx))


def find_center_vo(vo: np.ndarray) -> Tuple[float, float]:
    """Vortex center from raw ERA5 vorticity: max |vo| = cyclone core.
    More accurate than reconstructed-wind center for IBTrACS matching."""
    vo_abs = np.abs(vo)
    # Blur slightly to avoid noise spikes
    from scipy import ndimage
    vo_smooth = ndimage.gaussian_filter(vo_abs, sigma=2.0)
    cy, cx = np.unravel_index(np.argmax(vo_smooth), vo_smooth.shape)
    return float(cy), float(cx)


# ═══════════════════════════════════════════════════════════════
# 4. IMPROVED STATE MACHINE (P2 + P4 Fixes)
# ═══════════════════════════════════════════════════════════════

class StateMachine:
    def __init__(self):
        self.history: List[Dict] = []
        self.triggers: List[Dict] = []
        self.erc_start_frame: Optional[int] = None  # P4: track ERC onset
        self._last_erc_trigger_t: Optional[float] = None
        self._last_ri_trigger_t: Optional[float] = None  # A+: RI cooldown
        self.RI_COOLDOWN_H: float = 48.0  # don't re-trigger RI within 48h
    
    def evaluate(self, frame_idx: int, t_h: float,
                 theta1: float, dh_curl: float, dfold: float,
                 core_vorticity: float = 0.0,
                 ibtracs_match: Optional[Dict] = None,
                 ibtracs_loaded: bool = False) -> Tuple[str, Dict]:
        """
        Evaluate state. Returns (state_str, metadata_dict).
        
        P1 IBTrACS logic:
          - ibtracs_loaded=False → fall back to vorticity-only TC detection
          - ibtracs_loaded=True, ibtracs_match={}  → IBTrACS loaded, no spatial match → NO_CYCLONE
          - ibtracs_loaded=True, ibtracs_match populated → confirmed TC → evaluate triggers
        """
        meta = {
            'frame': frame_idx, 't_h': t_h,
            'theta1': theta1, 'dh_curl': dh_curl, 'dfold': dfold,
            'ibtracs': ibtracs_match
        }
        
        self.history.append({
            'frame': frame_idx, 't_h': t_h,
            'theta1': theta1, 'dh_curl': dh_curl, 'dfold': dfold,
            'core_vorticity': core_vorticity,
            'ibtracs_match': ibtracs_match
        })
        
        # Warmup
        if len(self.history) < 6 or (t_h - self.history[0]['t_h']) < 6.0:
            return 'WARMUP', meta
        
        # P1: IBTrACS spatial coincidence check
        if not ibtracs_loaded:
            # No IBTrACS data → fall back to vorticity filter
            if abs(core_vorticity) < cfg.CYCLONE_VORT_MIN:
                return 'NO_CYCLONE', meta
            if abs(dh_curl) < cfg.DH_ABS_MIN:
                return 'NO_CYCLONE', meta
        elif not ibtracs_match or ibtracs_match.get('dist_km', 999) > cfg.IBTRACS_MAX_DIST_KM:
            # IBTrACS loaded but no spatial match → no TC at this location/time
            return 'NO_CYCLONE', meta
        # else: confirmed TC by IBTrACS → proceed to trigger evaluation
        
        states = []
        
        # Check RI precursor (A+: dual pathway + cooldown)
        ri = self._check_ri_precursor()
        if ri:
            states.append('RI_PRECURSOR')
            t_now = self.history[-1]['t_h']
            if self._last_ri_trigger_t is None or (t_now - self._last_ri_trigger_t) > self.RI_COOLDOWN_H:
                self.triggers.append(ri)
                self._last_ri_trigger_t = t_now
        
        # Check ERC oscillation (P4: persistence ≥12h)
        erc = self._check_erc()
        if erc:
            states.append('ERC_ACTIVE')
            if self._last_erc_trigger_t is None or (t_h - self._last_erc_trigger_t) > cfg.ERC_PERSISTENCE_H:
                self.triggers.append(erc)
                self._last_erc_trigger_t = t_h
        
        # Check impulse response (B Fix: gated — only for volcanic/external shock)
        if cfg.IMPULSE_ENABLED:
            imp = self._check_impulse()
            if imp:
                states.append('IMPULSE')
                self.triggers.append(imp)
        
        if not states:
            states.append('NOMINAL')
        
        meta['state'] = '|'.join(states)
        return meta['state'], meta
    
    def _window(self, h: float) -> List[Dict]:
        if not self.history:
            return []
        t_now = self.history[-1]['t_h']
        return [r for r in self.history if r['t_h'] >= t_now - h]
    
    def _check_ri_precursor(self) -> Optional[Dict]:
        """
        Two RI pathways:
          Pathway 1 (Dolphin): θ₁ drops >15° + dH enters mature zone  — precursor
          Pathway 2 (Chaotic, A+): dH plummets >1e-4 + θ₁ volatile >8°  — rapid organization
        """
        # ── Pathway 1: θ₁-drop precursor (Dolphin-style) ──
        w1 = self._window(cfg.RI_THETA_WINDOW_H)
        if len(w1) >= 2:
            th_start = w1[0]['theta1']
            th_end = w1[-1]['theta1']
            th_drop = th_start - th_end
            dh_now = w1[-1]['dh_curl']
            dh_start = w1[0]['dh_curl']
            dh_drop_abs = dh_start - dh_now
            
            if (th_drop > cfg.RI_THETA_DROP_DEG and 
                dh_now < cfg.DH_MATURE_THRESHOLD and
                dh_drop_abs > cfg.DH_RI_ABS_DROP):
                return {
                    'type': 'RI_PRECURSOR',
                    'pathway': 'P1_theta_drop',
                    't_h': self.history[-1]['t_h'],
                    'theta_drop_deg': float(th_drop),
                    'dh_now': float(dh_now),
                    'dh_drop_abs': float(dh_drop_abs),
                    'confidence': 'empirical_baseline_v0.5'
                }
        
        # ── Pathway 2: Chaotic RI (A+ Fix) ──
        # dH plummets (core rapidly organizing) + θ₁ chaos (oscillatory OR locked)
        w2 = self._window(cfg.RI_CHAOTIC_WINDOW_H)
        if len(w2) >= 3:
            dhs = [r['dh_curl'] for r in w2]
            ths = [r['theta1'] for r in w2]
            
            dh_start2 = dhs[0]
            dh_end2 = dhs[-1]
            dh_plummet = dh_start2 - dh_end2  # positive = dH becoming more negative
            
            th_std = float(np.std(ths))
            th_min = float(np.min(ths))
            
            # θ₁ chaos: either oscillatory (std > 8°) or locked disorder (all > 35°)
            th_chaotic = (th_std > cfg.RI_THETA_VOLATILITY_DEG or 
                         th_min > cfg.RI_THETA_LOCKED_MIN_DEG)
            
            if (dh_plummet > cfg.RI_DH_ABS_PLUMMET and
                dh_end2 < cfg.DH_MATURE_THRESHOLD and
                th_chaotic):
                chaos_type = 'oscillatory' if th_std > cfg.RI_THETA_VOLATILITY_DEG else 'locked'
                return {
                    'type': 'RI_PRECURSOR',
                    'pathway': f'P2_chaotic_dh_{chaos_type}',
                    't_h': self.history[-1]['t_h'],
                    'dh_plummet': float(dh_plummet),
                    'dh_end': float(dh_end2),
                    'theta_std_deg': float(th_std),
                    'theta_min_deg': float(th_min),
                    'confidence': 'empirical_baseline_v0.5'
                }
        
        return None
    
    def _check_erc(self) -> Optional[Dict]:
        """P4 Fix: ERC requires ≥12h sustained θ₁ oscillation + organized TC structure."""
        w = self._window(cfg.ERC_WINDOW_H)
        if len(w) < 4:
            self.erc_start_frame = None
            return None
        
        ths = [r['theta1'] for r in w]
        th_range = max(ths) - min(ths)
        
        if not (cfg.ERC_THETA_OSC_MIN < th_range < cfg.ERC_THETA_OSC_MAX):
            self.erc_start_frame = None
            return None
        
        # P2 Fix: ERC requires organized TC structure (dH must be negative)
        dh_recent = np.median([r['dh_curl'] for r in w[-4:]])
        if dh_recent > 0:
            return None  # not organized enough for ERC
        
        # Check persistence: θ₁ must have been oscillating for ≥12h
        w_long = self._window(cfg.ERC_PERSISTENCE_H + cfg.ERC_WINDOW_H)
        if len(w_long) < 4:
            return None
        
        # Check multiple sub-windows for sustained oscillation
        sustained = True
        window_step_h = 6.0
        n_checks = max(1, int(cfg.ERC_PERSISTENCE_H / window_step_h))
        for i in range(n_checks):
            offset_h = i * window_step_h
            sub = [r for r in w_long 
                   if r['t_h'] >= (self.history[-1]['t_h'] - cfg.ERC_WINDOW_H - offset_h)
                   and r['t_h'] <= (self.history[-1]['t_h'] - offset_h)]
            if len(sub) >= 4:
                sub_ths = [r['theta1'] for r in sub]
                sub_range = max(sub_ths) - min(sub_ths)
                if sub_range < cfg.ERC_THETA_OSC_MIN * 0.7:  # 70% of min threshold
                    sustained = False
                    break
        
        if not sustained:
            return None
        
        return {
            'type': 'ERC_ACTIVE',
            't_h': self.history[-1]['t_h'],
            'theta_range_deg': float(th_range),
            'persistence_h': cfg.ERC_PERSISTENCE_H,
            'confidence': 'empirical_baseline_v0.5'
        }
    
    def _check_impulse(self) -> Optional[Dict]:
        """Tonga-type impulse: rapid θ₁ drop + dH_curl surge."""
        w = self._window(cfg.IMPULSE_WINDOW_H)
        if len(w) < 2:
            return None
        
        th_drop = w[0]['theta1'] - w[-1]['theta1']
        dh_now = w[-1]['dh_curl']
        dh_start = w[0]['dh_curl']
        dh_drop_abs = dh_start - dh_now
        
        if th_drop > cfg.IMPULSE_THETA_DROP and dh_drop_abs > cfg.DH_IMPULSE_ABS_DROP:
            return {
                'type': 'IMPULSE_RESPONSE',
                't_h': self.history[-1]['t_h'],
                'theta_drop_deg': float(th_drop),
                'dh_drop_abs': float(dh_drop_abs),
                'confidence': 'empirical_baseline_v0.5'
            }
        return None


# ═══════════════════════════════════════════════════════════════
# 5. POISSON SOLVER (for ERA5 vo/d → u/v reconstruction)
# ═══════════════════════════════════════════════════════════════

def solve_poisson_spectral(field: np.ndarray) -> np.ndarray:
    """Solve ∇²ψ = field via FFT (periodic x) + tridiagonal (Dirichlet y)."""
    ny, nx = field.shape
    f_hat = np.fft.rfft(field, axis=1)
    kx = 2 * np.pi * np.fft.rfftfreq(nx)
    psi_hat = np.zeros_like(f_hat, dtype=complex)
    
    for j in range(f_hat.shape[1]):
        kj2 = kx[j]**2
        if kj2 < 1e-30:
            psi_hat[:, j] = _solve_poisson_1d(f_hat[:, j].real)
        else:
            psi_hat[:, j] = _solve_helmholtz_1d(f_hat[:, j], kj2)
    
    return np.fft.irfft(psi_hat, n=nx, axis=1)


def _solve_poisson_1d(f: np.ndarray) -> np.ndarray:
    ny = len(f)
    A = np.zeros((ny-2, ny-2))
    for i in range(ny-2):
        A[i, i] = -2.0
        if i > 0: A[i, i-1] = 1.0
        if i < ny-3: A[i, i+1] = 1.0
    rhs = f[1:ny-1]
    try:
        psi_inner = np.linalg.solve(A, rhs)
    except np.linalg.LinAlgError:
        psi_inner = np.zeros(ny-2)
    psi = np.zeros(ny)
    psi[1:ny-1] = psi_inner
    return psi


def _solve_helmholtz_1d(f_hat: np.ndarray, k2: float) -> np.ndarray:
    ny = len(f_hat)
    A = np.zeros((ny-2, ny-2), dtype=complex)
    for i in range(ny-2):
        A[i, i] = -2.0 - k2
        if i > 0: A[i, i-1] = 1.0
        if i < ny-3: A[i, i+1] = 1.0
    rhs = f_hat[1:ny-1]
    try:
        psi_inner = np.linalg.solve(A, rhs)
    except np.linalg.LinAlgError:
        psi_inner = np.zeros(ny-2, dtype=complex)
    psi = np.zeros(ny, dtype=complex)
    psi[1:ny-1] = psi_inner
    return psi


def reconstruct_winds(vo: np.ndarray, d: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Reconstruct u, v from vorticity & divergence via spectral Poisson."""
    psi = solve_poisson_spectral(vo)
    chi = solve_poisson_spectral(d)
    dpsi_dy, dpsi_dx = np.gradient(psi)
    dchi_dx, dchi_dy = np.gradient(chi)
    u = -dpsi_dy + dchi_dx
    v = dpsi_dx + dchi_dy
    return u, v


# ═══════════════════════════════════════════════════════════════
# 6. MAIN PIPELINE
# ═══════════════════════════════════════════════════════════════

def run_pipeline(nc_path: str, ibtracs_path: Optional[str] = None,
                 storm_name: str = "HAGIBIS", storm_year: int = 2019,
                 subsample: int = 1, max_frames: Optional[int] = None) -> Dict:
    """
    Full Phase 0.5 pipeline:
      1. Load ERA5 NetCDF (vo, d)
      2. Reconstruct winds (u, v)
      3. Load IBTrACS track (P1)
      4. Frame-by-frame: 3 operators + state machine with P1/P2/P4 fixes
    """
    import xarray as xr
    
    t0_wall = time.time()
    
    # ── 6a. Load ERA5 ──
    ds = xr.open_dataset(nc_path)
    vo_all = ds.vo.isel(pressure_level=0).values  # 850hPa
    d_all = ds.d.isel(pressure_level=0).values
    times = ds.valid_time.values
    lats = ds.latitude.values
    lons = ds.longitude.values
    
    t0_dt = times[0].astype('datetime64[h]').astype(float)
    t_hours = np.array([(t.astype('datetime64[h]').astype(float) - t0_dt) for t in times])
    t0_datetime = times[0].astype('M8[s]').astype(datetime) if hasattr(times[0], 'astype') else datetime(2019,9,1)
    
    n_full = vo_all.shape[0]
    ds.close()
    
    # ── 6b. Subsample ──
    if subsample > 1:
        vo_all = vo_all[::subsample]
        d_all = d_all[::subsample]
        t_hours = t_hours[::subsample]
    if max_frames:
        vo_all = vo_all[:max_frames]
        d_all = d_all[:max_frames]
        t_hours = t_hours[:max_frames]
    
    n_frames, ny, nx = vo_all.shape
    
    print(f"\n{'═'*60}")
    print(f"  🧪 Ô-HAT Engine v0.5 — {storm_name} {storm_year}")
    print(f"  Frames: {n_frames}  Grid: {ny}×{nx}  Duration: {t_hours[-1]:.0f}h")
    print(f"  Fixes: P1(IBTrACS ±{cfg.IBTRACS_MAX_DIST_KM:.0f}km) "
          f"P2(dH<{cfg.DH_MATURE_THRESHOLD:.1e}) "
          f"P4(ERC≥{cfg.ERC_PERSISTENCE_H:.0f}h)")
    print(f"{'═'*60}")
    
    # ── 6c. Load IBTrACS (P1) ──
    track = None
    if ibtracs_path and os.path.exists(ibtracs_path):
        track = load_ibtracs_track(ibtracs_path, storm_name, storm_year)
        if track:
            print(f"  📍 IBTrACS track loaded: {len(track)} points, "
                  f"{track[0]['time']} → {track[-1]['time']}")
            print(f"     Wind range: {min(p['wind'] for p in track if p['wind']):.0f} – "
                  f"{max(p['wind'] for p in track if p['wind']):.0f} kts")
        else:
            print(f"  ⚠️  IBTrACS: no track found for {storm_name} {storm_year}")
    else:
        print(f"  ⚠️  IBTrACS: no file provided — using vorticity-only TC detection")
    
    # ── 6d. Reconstruct winds ──
    print(f"  🔄 Reconstructing wind field (spectral Poisson)...")
    t_rec = time.time()
    u_all = np.zeros_like(vo_all)
    v_all = np.zeros_like(vo_all)
    
    for i in range(n_frames):
        if i % max(n_frames // 5, 1) == 0:
            print(f"    ... frame {i}/{n_frames}")
        u_all[i], v_all[i] = reconstruct_winds(vo_all[i], d_all[i])
    
    print(f"    ✅ Reconstruction done in {time.time()-t_rec:.1f}s")
    
    # ── 6e. Streaming pipeline ──
    sm = StateMachine()
    results = []
    ibtracs_hits = 0
    ibtracs_misses = 0
    
    t_pipe = time.time()
    
    for i in range(n_frames):
        u = u_all[i]; v = v_all[i]; vo = vo_all[i]
        t_h = t_hours[i]
        
        # Find center: use raw vo for accuracy with IBTrACS
        if track:
            cy, cx = find_center_vo(vo)
        else:
            cy, cx = find_center(u, v)
        core_mask, shell_mask = make_masks(ny, nx, cx, cy)
        
        # dH_curl: from ERA5 raw vo (P3 workaround)
        dh_curl = compute_dh_curl_vo(vo, core_mask, shell_mask)
        
        # θ₁: from reconstructed winds
        theta1 = compute_theta1(u, v, core_mask)
        
        # Cyclone vorticity for fallback detection
        cyclone_vort = float(np.max(np.abs(vo)))
        
        # D_fold (every 6 frames)
        if i % cfg.dfold_window_frames == 0 or i < 2:
            zeta_mag = np.abs(
                (np.roll(v, -1, axis=1) - np.roll(v, 1, axis=1)) / 2.0 -
                (np.roll(u, -1, axis=0) - np.roll(u, 1, axis=0)) / 2.0
            )
            dfold = compute_dfold(zeta_mag, cfg.dfold_threshold_pct)
        else:
            dfold = float('nan')
        
        # IBTrACS match (P1)
        ib_match = None
        if track:
            ib_match = find_ibtracs_match(t_h, t0_datetime, cy, cx, lats, lons, track)
            if ib_match and 'dist_km' in ib_match:
                ibtracs_hits += 1
            else:
                ibtracs_misses += 1
        
        # State machine
        state, meta = sm.evaluate(i, t_h, theta1, dh_curl, dfold, cyclone_vort, 
                                  ib_match if ib_match else {}, 
                                  ibtracs_loaded=(track is not None))
        
        results.append({
            'frame': i, 't_h': float(t_h),
            'theta1': theta1, 'dh_curl': dh_curl, 'dfold': dfold,
            'state': state,
            'ibtracs': ib_match,
            'center_y': float(cy), 'center_x': float(cx)
        })
        
        if i % max(n_frames // 20, 1) == 0:
            flag = '⚠️' if any(s in state for s in ['RI_PRECURSOR','ERC_ACTIVE','IMPULSE']) else ' '
            if 'NO_CYCLONE' in state: flag = '·'
            if 'IBTRACS_MISS' in state: flag = '○'
            ib_str = f"IB={ib_match['dist_km']:.0f}km" if (ib_match and 'dist_km' in ib_match) else 'IB=---'
            df_str = f"{dfold:.2f}" if not (isinstance(dfold, float) and np.isnan(dfold)) else "  ---"
            print(f"  [{flag}] t={t_h:6.0f}h  θ₁={theta1:5.1f}°  dH={dh_curl:+.2e}  "
                  f"D={df_str}  {ib_str}  {state}")
    
    pipe_elapsed = time.time() - t_pipe
    
    # ── 6f. Summary ──
    triggers = sm.triggers
    theta_arr = np.array([r['theta1'] for r in results])
    dh_arr = np.array([r['dh_curl'] for r in results])
    
    trigger_types = {}
    for t in triggers:
        tt = t['type']
        trigger_types[tt] = trigger_types.get(tt, 0) + 1
    
    # Count frames in each state
    state_counts = {}
    for r in results:
        for s in r['state'].split('|'):
            state_counts[s] = state_counts.get(s, 0) + 1
    
    print(f"\n{'─'*60}")
    print(f"  📊 PIPELINE COMPLETE")
    print(f"  ⏱️  Reconstruction: {t_rec - t0_wall:.1f}s  "
          f"Pipeline: {pipe_elapsed:.1f}s  "
          f"Total: {time.time()-t0_wall:.1f}s")
    print(f"  📡 Frames: {n_frames}  Throughput: {n_frames/pipe_elapsed:.1f} fps")
    if track:
        print(f"  📍 IBTrACS: {ibtracs_hits} hits / {ibtracs_misses} misses "
              f"({100*ibtracs_hits/max(ibtracs_hits+ibtracs_misses,1):.0f}% hit rate)")
    print(f"  📏 θ₁: {theta_arr.min():.1f}° – {theta_arr.max():.1f}° "
          f"(μ={theta_arr.mean():.1f}°, σ={theta_arr.std():.1f}°)")
    print(f"  📏 dH_curl: {dh_arr.min():.2e} – {dh_arr.max():.2e}")
    print(f"  🎯 State distribution: {state_counts}")
    print(f"  ⚡ Trigger events: {len(triggers)} {trigger_types}")
    
    # Validate against Phase 2
    print(f"\n{'─'*60}")
    print(f"  🔍 PHASE 2 VALIDATION")
    print(f"{'─'*60}")
    
    checks = []
    
    # 1. ERC check: should fire during Hagibis ERC period (Oct 6-8 for Hagibis 2019)
    erc_triggers = [t for t in triggers if t['type'] == 'ERC_ACTIVE']
    if storm_name.upper() == 'HAGIBIS':
        erc_ok = 1 <= len(erc_triggers) <= 10  # should fire but not 156 times!
        checks.append(('ERC_ACTIVE count reasonable (1-10)', erc_ok, 
                      f'{len(erc_triggers)} triggers (v0.2 had 156)'))
    else:
        checks.append(('ERC_ACTIVE detected', len(erc_triggers) > 0,
                      f'{len(erc_triggers)} triggers'))
    
    # 2. RI check
    ri_triggers = [t for t in triggers if t['type'] == 'RI_PRECURSOR']
    checks.append(('RI_PRECURSOR detected', len(ri_triggers) > 0,
                  f'{len(ri_triggers)} triggers'))
    
    # 3. NO_CYCLONE frames dominate outside storm period
    nc_frames = state_counts.get('NO_CYCLONE', 0)
    checks.append(('NO_CYCLONE filtering active', nc_frames > n_frames * 0.3,
                  f'{nc_frames}/{n_frames} frames ({100*nc_frames/n_frames:.0f}%)'))
    
    # 4. dH_curl range
    dh_range = float(dh_arr.max() - dh_arr.min())
    checks.append(('dH_curl dynamic range > 1e-5', dh_range > 1e-5,
                  f'range={dh_range:.2e}'))
    
    passed = sum(1 for _, ok, _ in checks if ok)
    for name, ok, detail in checks:
        print(f"  {'✅' if ok else '❌'} {name}: {detail}")
    print(f"  Result: {passed}/{len(checks)} checks")
    
    # ── 6g. Save ──
    out = {
        'engine_version': '0.5',
        'fixes': ['P1_IBTrACS_spatial', 'P2_absolute_dH', 'P4_ERC_persistence'],
        'case': f'{storm_name.lower()}_{storm_year}',
        'n_frames': n_frames,
        'duration_h': float(t_hours[-1]),
        'checks': [{'name': n, 'passed': o, 'detail': d} for n, o, d in checks],
        'n_passed': passed,
        'n_total': len(checks),
        'triggers': triggers,
        'state_counts': state_counts,
        'ibtracs_stats': {'hits': ibtracs_hits, 'misses': ibtracs_misses} if track else None,
        'summary': {
            'theta_min': float(theta_arr.min()),
            'theta_max': float(theta_arr.max()),
            'theta_mean': float(theta_arr.mean()),
            'theta_std': float(theta_arr.std()),
            'dh_min': float(dh_arr.min()),
            'dh_max': float(dh_arr.max()),
            'dh_mean': float(dh_arr.mean()),
            'dh_std': float(dh_arr.std()),
            'triggers': trigger_types,
        },
        'config': {
            'ibtracs_max_dist_km': cfg.IBTRACS_MAX_DIST_KM,
            'ibtracs_max_dt_h': cfg.IBTRACS_MAX_DT_H,
            'dh_mature_threshold': cfg.DH_MATURE_THRESHOLD,
            'dh_ri_abs_drop': cfg.DH_RI_ABS_DROP,
            'erc_persistence_h': cfg.ERC_PERSISTENCE_H,
        },
        'results': results[-20:]  # last 20 frames for inspection
    }
    
    out_path = OUT_DIR / f'engine_v0.5_{storm_name.lower()}_{storm_year}.json'
    json.dump(out, open(out_path, 'w'), indent=2, default=str)
    print(f"\n  💾 Saved: {out_path}")
    
    return out


# ═══════════════════════════════════════════════════════════════
# 7. CLI
# ═══════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description='Ô-HAT Engine v0.5 — Phase 0.5 Fixes')
    parser.add_argument('--nc', type=str, 
                       default='data/typhoons/hagibis_2019.nc')
    parser.add_argument('--ibtracs', type=str,
                       default='projects/cyclone/data/ibtracs_all.csv')
    parser.add_argument('--storm', type=str, default='HAGIBIS')
    parser.add_argument('--year', type=int, default=2019)
    parser.add_argument('--subsample', type=int, default=1)
    parser.add_argument('--max-frames', type=int, default=None)
    parser.add_argument('--impulse', action='store_true', 
                       help='Enable IMPULSE mode (volcanic/external shock detection)')
    parser.add_argument('--output-dir', type=str, default=None,
                       help='Output directory for engine results (default: script dir/output/)')
    args = parser.parse_args()
    
    if args.output_dir:
        global OUT_DIR
        OUT_DIR = Path(args.output_dir)
        OUT_DIR.mkdir(exist_ok=True)
    
    if args.impulse:
        cfg.IMPULSE_ENABLED = True
        print("  🌋 IMPULSE mode ENABLED (volcanic/external shock detection)")
    
    run_pipeline(
        nc_path=args.nc,
        ibtracs_path=args.ibtracs,
        storm_name=args.storm,
        storm_year=args.year,
        subsample=args.subsample,
        max_frames=args.max_frames
    )


if __name__ == '__main__':
    main()
