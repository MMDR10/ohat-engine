#!/usr/bin/env python3
"""
🧪 Ô-HAT Real-Time Engine — Phase 0 Mock Pipeline
==================================================
Streams ERA5 historical data frame-by-frame to simulate real-time GFS feed.
Computes θ₁, dH_curl, D_fold per frame. State machine triggers on Phase 2 thresholds.

Usage:
  python ohat_engine_mock.py --case hagibis   # single known case
  python ohat_engine_mock.py --case all       # all known cases
  python ohat_engine_mock.py --file typhoon.nc  # custom netCDF
"""
import sys, json, os, argparse, time
import numpy as np
from scipy import stats, ndimage
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple

OUT_DIR = Path('output')
OUT_DIR.mkdir(exist_ok=True)

# ═══════════════════════════════════════════════════════════════
# 1. CONFIG — Phase 2 Validated Thresholds
# ═══════════════════════════════════════════════════════════════

@dataclass
class EngineConfig:
    # Spatial
    core_radius_km: float = 30.0
    shell_radius_km: float = 100.0
    dx_km: float = 0.25  # ERA5 0.25°
    
    # θ₁ window
    theta_window_h: float = 6.0   # hours for rolling PCA
    theta_fallback: float = 90.0  # if < 4 data points
    
    # D_fold window
    dfold_window_frames: int = 6  # rolling window for box-counting
    dfold_threshold_pct: float = 80  # top percentile for singularity set
    
    # State machine thresholds (from Phase 2 closure)
    RI_THETA_DROP_DEG: float = 15.0    # θ₁ drop > 15° in 15h
    RI_THETA_WINDOW_H: float = 15.0
    RI_DH_SURGE_PCT: float = 200.0     # dH_curl surge > 200%
    
    ERC_THETA_OSC_MIN: float = 15.0    # θ₁ oscillation amplitude
    ERC_THETA_OSC_MAX: float = 45.0
    ERC_WINDOW_H: float = 24.0
    
    IMPULSE_THETA_DROP: float = 10.0   # post-impact θ₁ drop
    IMPULSE_WINDOW_H: float = 6.0
    IMPULSE_DH_SURGE: float = 100.0
    
    # Cyclone presence filter (prevents false triggers on background weather)
    CYCLONE_VORT_MIN: float = 1e-4   # max |vo| in domain > 1e-4 = tropical cyclone present
    DH_ABS_MIN: float = 1e-10        # |dH_curl| must exceed noise floor
    RI_DH_ABS_SURGE: float = 1e-10   # absolute dH_curl surge for RI (not just %)

cfg = EngineConfig()

# ═══════════════════════════════════════════════════════════════
# 2. CORE OPERATORS
# ═══════════════════════════════════════════════════════════════

def compute_dh_curl(u: np.ndarray, v: np.ndarray, 
                    core_mask: np.ndarray, shell_mask: np.ndarray) -> float:
    """
    dH_curl = ⟨ζ_shell⟩ - ⟨ζ_core⟩
    ζ = ∂v/∂x - ∂u/∂y (relative vorticity, centered diff)
    """
    # Vorticity on interior points
    dv_dx = (np.roll(v, -1, axis=1) - np.roll(v, 1, axis=1)) / (2 * cfg.dx_km * 111000)
    du_dy = (np.roll(u, -1, axis=0) - np.roll(u, 1, axis=0)) / (2 * cfg.dx_km * 111000)
    zeta = dv_dx - du_dy
    
    core_zeta = float(np.mean(zeta[core_mask])) if core_mask.sum() > 0 else 0.0
    shell_zeta = float(np.mean(zeta[shell_mask])) if shell_mask.sum() > 0 else 0.0
    
    return shell_zeta - core_zeta


def compute_theta1(u: np.ndarray, v: np.ndarray, core_mask: np.ndarray) -> float:
    """
    θ₁ = arctan(√(λ_min/λ_max)) from PCA of (u,v) in core region.
    λ_min → 0 means flow collapses to 1D → θ₁ → 0° (singularity).
    """
    u_c = u[core_mask]; v_c = v[core_mask]
    if len(u_c) < 4:
        return cfg.theta_fallback
    
    X = np.column_stack([u_c, v_c])
    X -= X.mean(axis=0)
    cov = (X.T @ X) / (len(X) - 1)
    
    try:
        ev = np.linalg.eigvalsh(cov)  # ascending: [λ_min, λ_max]
    except np.linalg.LinAlgError:
        return cfg.theta_fallback
    
    if ev[1] < 1e-30:
        return 0.0  # fully collapsed
    
    ratio = max(ev[0], 0) / ev[1]
    return float(np.degrees(np.arctan(np.sqrt(ratio))))


def compute_dfold(field: np.ndarray, threshold_pct: float = 80) -> float:
    """
    D_fold: box-counting dimension of singularity set (top percentile).
    Returns D (1~2 for line~surface), NaN if insufficient points.
    """
    ny, nx = field.shape
    thresh = np.percentile(np.abs(field), threshold_pct)
    mask = np.abs(field) >= thresh
    
    y_idx, x_idx = np.where(mask)
    if len(y_idx) < 10:
        return float('nan')
    
    points = np.column_stack([x_idx, y_idx])
    
    # Box sizes: powers of 2 from 1 to half domain
    max_box = max(nx, ny) // 2
    sizes = [2**k for k in range(0, int(np.log2(max_box)) + 1)]
    sizes = [s for s in sizes if s <= max_box]
    
    counts = []
    for s in sizes:
        bins_x = np.arange(0, nx + s, s)
        bins_y = np.arange(0, ny + s, s)
        h, _, _ = np.histogram2d(x_idx, y_idx, bins=[bins_x, bins_y])
        counts.append(np.sum(h > 0))
    
    # Linear fit log(N) vs log(1/s)
    valid = [i for i, c in enumerate(counts) if c > 0]
    if len(valid) < 3:
        return float('nan')
    
    log_s = np.log([1.0 / sizes[i] for i in valid])
    log_n = np.log([counts[i] for i in valid])
    
    slope, _, _, _, _ = stats.linregress(log_s, log_n)
    return float(slope)


def make_masks(ny: int, nx: int, cx: float, cy: float) -> Tuple[np.ndarray, np.ndarray]:
    """Create core and shell masks centered at (cx, cy) in pixel coordinates."""
    y, x = np.meshgrid(np.arange(ny), np.arange(nx), indexing='ij')
    r = np.sqrt((x - cx)**2 + (y - cy)**2)
    
    core_px = max(3, cfg.core_radius_km / (cfg.dx_km * 111.0))
    shell_px = max(6, cfg.shell_radius_km / (cfg.dx_km * 111.0))
    
    core_mask = r <= core_px
    shell_mask = (r > core_px) & (r <= shell_px)
    return core_mask, shell_mask


def find_center(u: np.ndarray, v: np.ndarray) -> Tuple[float, float]:
    """Vortex center: minimum wind speed location (simple centroid of low-wind region)."""
    speed = np.sqrt(u**2 + v**2)
    thresh = np.percentile(speed, 20)  # lowest 20%
    mask = speed <= thresh
    if mask.sum() == 0:
        return float(u.shape[0] // 2), float(u.shape[1] // 2)
    y_idx, x_idx = np.where(mask)
    cy = float(np.mean(y_idx))
    cx = float(np.mean(x_idx))
    return cy, cx


# ═══════════════════════════════════════════════════════════════
# 3. STATE MACHINE
# ═══════════════════════════════════════════════════════════════

class StateMachine:
    def __init__(self):
        self.history: List[Dict] = []
        self.triggers: List[Dict] = []
    
    def evaluate(self, frame_idx: int, t_h: float,
                 theta1: float, dh_curl: float, dfold: float,
                 core_vorticity: float = 0.0) -> str:
        """Evaluate state and return classification."""
        self.history.append({
            'frame': frame_idx, 't_h': t_h,
            'theta1': theta1, 'dh_curl': dh_curl, 'dfold': dfold,
            'core_vorticity': core_vorticity
        })
        
        # Warmup: need at least 6h of data before triggering
        if len(self.history) < 6 or (t_h - self.history[0]['t_h']) < 6.0:
            return 'WARMUP'
        
        # Cyclone presence pre-filter: must have actual tropical cyclone
        if abs(core_vorticity) < cfg.CYCLONE_VORT_MIN:
            return 'NO_CYCLONE'
        
        # dH_curl noise floor: must exceed absolute minimum
        if abs(dh_curl) < cfg.DH_ABS_MIN:
            return 'NO_CYCLONE'
        
        states = []
        
        # Check RI precursor
        ri = self._check_ri_precursor()
        if ri:
            states.append('RI_PRECURSOR')
            self.triggers.append(ri)
        
        # Check ERC oscillation
        erc = self._check_erc()
        if erc:
            states.append('ERC_ACTIVE')
            self.triggers.append(erc)
        
        # Check impulse response
        imp = self._check_impulse()
        if imp:
            states.append('IMPULSE')
            self.triggers.append(imp)
        
        if not states:
            states.append('NOMINAL')
        
        return '|'.join(states)
    
    def _window(self, h: float) -> List[Dict]:
        """Get history within last h hours."""
        if not self.history:
            return []
        t_now = self.history[-1]['t_h']
        return [r for r in self.history if r['t_h'] >= t_now - h]
    
    def _check_ri_precursor(self) -> Optional[Dict]:
        w = self._window(cfg.RI_THETA_WINDOW_H)
        if len(w) < 2:
            return None
        
        th_start = w[0]['theta1']
        th_end = w[-1]['theta1']
        th_drop = th_start - th_end
        
        dh_now = w[-1]['dh_curl']
        dh_ref = np.median([r['dh_curl'] for r in w[:len(w)//2]]) if len(w) > 2 else w[0]['dh_curl']
        dh_surge = ((dh_now - dh_ref) / (abs(dh_ref) + 1e-30)) * 100
        
        if th_drop > cfg.RI_THETA_DROP_DEG and dh_surge > cfg.RI_DH_SURGE_PCT \
           and abs(dh_now - dh_ref) > cfg.RI_DH_ABS_SURGE:
            return {
                'type': 'RI_PRECURSOR',
                't_h': self.history[-1]['t_h'],
                'theta_drop_deg': float(th_drop),
                'dh_surge_pct': float(dh_surge),
                'lead_h': cfg.RI_THETA_WINDOW_H,
                'confidence': 'empirical_baseline'
            }
        return None
    
    def _check_erc(self) -> Optional[Dict]:
        w = self._window(cfg.ERC_WINDOW_H)
        if len(w) < 4:
            return None
        
        ths = [r['theta1'] for r in w]
        th_range = max(ths) - min(ths)
        
        if cfg.ERC_THETA_OSC_MIN < th_range < cfg.ERC_THETA_OSC_MAX:
            return {
                'type': 'ERC_ACTIVE',
                't_h': self.history[-1]['t_h'],
                'theta_range_deg': float(th_range),
                'confidence': 'empirical_baseline'
            }
        return None
    
    def _check_impulse(self) -> Optional[Dict]:
        w = self._window(cfg.IMPULSE_WINDOW_H)
        if len(w) < 2:
            return None
        
        th_drop = w[0]['theta1'] - w[-1]['theta1']
        dh_now = w[-1]['dh_curl']
        dh_ref = w[0]['dh_curl']
        dh_surge = ((dh_now - dh_ref) / (abs(dh_ref) + 1e-30)) * 100
        
        if th_drop > cfg.IMPULSE_THETA_DROP and dh_surge > cfg.IMPULSE_DH_SURGE:
            return {
                'type': 'IMPULSE_RESPONSE',
                't_h': self.history[-1]['t_h'],
                'theta_drop_deg': float(th_drop),
                'dh_surge_pct': float(dh_surge),
                'confidence': 'empirical_baseline'
            }
        return None


# ═══════════════════════════════════════════════════════════════
# 4. DATA LOADER (ERA5 NetCDF → Frame Stream)
# ═══════════════════════════════════════════════════════════════

def load_era5_stream(nc_path: str, u_var: str = 'u10', v_var: str = 'v10',
                     time_var: str = 'time', lat_var: str = 'latitude',
                     lon_var: str = 'longitude') -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Load ERA5 netCDF and yield (u, v, t_hours) frames.
    Requires xarray. Falls back to scipy.io.netcdf if xarray unavailable.
    """
    try:
        import xarray as xr
        ds = xr.open_dataset(nc_path)
        u = ds[u_var].values   # (time, lat, lon)
        v = ds[v_var].values
        times = ds[time_var].values
        
        # Convert times to hours from start
        if hasattr(times[0], 'astype'):
            t0 = times[0].astype('datetime64[h]').astype(float)
            t_hours = np.array([(t.astype('datetime64[h]').astype(float) - t0) for t in times])
        else:
            t_hours = np.arange(len(times), dtype=float)
        
        return u, v, t_hours
    except ImportError:
        from scipy.io import netcdf_file
        with netcdf_file(nc_path, 'r') as f:
            u = f.variables[u_var].data.copy()
            v = f.variables[v_var].data.copy()
            t_hours = np.arange(u.shape[0], dtype=float)
        return u, v, t_hours


def generate_synthetic_stream(ny: int = 40, nx: int = 60, n_frames: int = 48):
    """
    Generate synthetic cyclone data for testing when no netCDF available.
    A translating, intensifying Gaussian vortex.
    """
    u_all = np.zeros((n_frames, ny, nx))
    v_all = np.zeros((n_frames, ny, nx))
    t_hours = np.arange(n_frames, dtype=float)
    
    y, x = np.meshgrid(np.arange(ny), np.arange(nx), indexing='ij')
    
    for i in range(n_frames):
        # Vortex translates NW (β-drift) and intensifies
        cx = nx//2 - i * 0.3
        cy = ny//2 - i * 0.2
        r = np.sqrt((x - cx)**2 + (y - cy)**2)
        r0 = 8.0
        strength = 10.0 + i * 0.5  # intensifying
        
        # Rankine-like vortex: solid-body core, decay outside
        azimuthal = np.where(r < r0,
                            strength * r / r0,
                            strength * r0 / (r + 1e-6))
        
        theta = np.arctan2(y - cy, x - cx)
        u_all[i] = -azimuthal * np.sin(theta) + 2.0  # + background flow
        v_all[i] = azimuthal * np.cos(theta)
    
    return u_all, v_all, t_hours


# ═══════════════════════════════════════════════════════════════
# 5. MAIN PIPELINE
# ═══════════════════════════════════════════════════════════════

def run_pipeline(u_all: np.ndarray, v_all: np.ndarray, t_hours: np.ndarray,
                 label: str = "unknown") -> Dict:
    """
    Main streaming pipeline: frame-by-frame → 3 operators → state machine.
    """
    n_frames, ny, nx = u_all.shape
    sm = StateMachine()
    results = []
    
    print(f"\n{'═'*60}")
    print(f"  🧪 Ô-HAT Mock Engine — {label}")
    print(f"  Frames: {n_frames}  Grid: {ny}×{nx}"
          f"  Δx={cfg.dx_km}°  Duration: {t_hours[-1]:.1f}h")
    print(f"  Thresholds: θ₁↓>{cfg.RI_THETA_DROP_DEG}°  dH↑>{cfg.RI_DH_SURGE_PCT}%")
    print(f"{'═'*60}")
    
    for i in range(n_frames):
        u = u_all[i]; v = v_all[i]
        t_h = t_hours[i]
        
        # Find vortex center dynamically
        cy, cx = find_center(u, v)
        core_mask, shell_mask = make_masks(ny, nx, cx, cy)
        
        # Compute 3 operators
        dh_curl = compute_dh_curl(u, v, core_mask, shell_mask)
        theta1 = compute_theta1(u, v, core_mask)
        
        # D_fold on vorticity magnitude (every 6th frame for efficiency)
        if i % cfg.dfold_window_frames == 0 or i < 2:
            zeta_mag = np.abs(
                (np.roll(v, -1, axis=1) - np.roll(v, 1, axis=1)) / 1.0 -
                (np.roll(u, -1, axis=0) - np.roll(u, 1, axis=0)) / 1.0
            )
            dfold = compute_dfold(zeta_mag, cfg.dfold_threshold_pct)
        else:
            dfold = float('nan')  # skip for speed
        
        # State machine
        core_vort = float(np.mean(
            (np.roll(v, -1, axis=1) - np.roll(v, 1, axis=1)) / (2 * cfg.dx_km * 111000) -
            (np.roll(u, -1, axis=0) - np.roll(u, 1, axis=0)) / (2 * cfg.dx_km * 111000)
        )[core_mask]) if core_mask.sum() > 0 else 0.0
        state = sm.evaluate(i, t_h, theta1, dh_curl, dfold, core_vort)
        
        results.append({
            'frame': i, 't_h': float(t_h),
            'theta1': theta1, 'dh_curl': dh_curl, 'dfold': dfold,
            'state': state,
            'center_y': float(cy), 'center_x': float(cx)
        })
        
        if i % max(n_frames // 10, 1) == 0:
            trigger_flag = '⚠️' if 'NOMINAL' not in state else ' '
            df_str = f"{dfold:.2f}" if not (isinstance(dfold, float) and np.isnan(dfold)) else "  ---"
            print(f"  [{trigger_flag}] t={t_h:5.1f}h  "
                  f"θ₁={theta1:5.1f}°  dH={dh_curl:+.2e}  "
                  f"D_fold={df_str}  "
                  f"state={state}")
    
    # Summary
    n_triggers = len(sm.triggers)
    trigger_types = {}
    for t in sm.triggers:
        tt = t['type']
        trigger_types[tt] = trigger_types.get(tt, 0) + 1
    
    theta_series = [r['theta1'] for r in results]
    dh_series = [r['dh_curl'] for r in results]
    
    print(f"\n{'─'*60}")
    print(f"  📊 Pipeline Complete")
    print(f"  θ₁ range: {min(theta_series):.1f}° – {max(theta_series):.1f}°"
          f" (mean={np.mean(theta_series):.1f}°)")
    print(f"  dH_curl range: {min(dh_series):.2e} – {max(dh_series):.2e}")
    print(f"  Triggers: {n_triggers} ({trigger_types})")
    print(f"{'═'*60}\n")
    
    return {
        'label': label,
        'n_frames': n_frames,
        'results': results,
        'triggers': sm.triggers,
        'history': sm.history,
        'summary': {
            'theta_min': float(min(theta_series)),
            'theta_max': float(max(theta_series)),
            'theta_mean': float(np.mean(theta_series)),
            'theta_std': float(np.std(theta_series)),
            'dh_min': float(min(dh_series)),
            'dh_max': float(max(dh_series)),
            'n_triggers': n_triggers,
            'trigger_types': trigger_types,
        }
    }


# ═══════════════════════════════════════════════════════════════
# 6. CLI
# ═══════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description='Ô-HAT Real-Time Engine Mock Pipeline')
    parser.add_argument('--case', type=str, default='synthetic',
                       choices=['synthetic', 'hagibis', 'dolphin', 'tonga', 'all'],
                       help='Test case')
    parser.add_argument('--file', type=str, help='Path to ERA5 netCDF file')
    parser.add_argument('--frames', type=int, default=48, help='Frames for synthetic')
    parser.add_argument('--nx', type=int, default=80, help='Grid x for synthetic')
    parser.add_argument('--ny', type=int, default=60, help='Grid y for synthetic')
    args = parser.parse_args()
    
    t0 = time.time()
    
    if args.file:
        u, v, t = load_era5_stream(args.file)
        label = Path(args.file).stem
    else:
        u, v, t = generate_synthetic_stream(args.ny, args.nx, args.frames)
        label = f"synthetic_{args.case}"
    
    output = run_pipeline(u, v, t, label)
    
    elapsed = time.time() - t0
    print(f"  ⏱️  Runtime: {elapsed:.1f}s  "
          f"({output['n_frames']/elapsed:.1f} frames/s)")
    
    # Save
    out_path = OUT_DIR / f'engine_mock_{label}.json'
    json.dump(output, open(out_path, 'w'), indent=2, default=str)
    print(f"  💾 Saved: {out_path}")
    
    return output


if __name__ == '__main__':
    main()
