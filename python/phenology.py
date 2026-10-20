"""Core functions: harmonic regression of NDVI time series and phenology metrics.

A first-order harmonic model is fitted to each field's monthly NDVI series:

    NDVI(t) = b0 + b1 cos(2 pi t) + b2 sin(2 pi t)

with t in fractional years. From the coefficients:

* amplitude            A   = sqrt(b1^2 + b2^2)
* baseline                 = b0
* normalised amplitude     = A / b0   (seasonal variability)
* phase                phi = atan2(b2, b1), wrapped to [0, 2 pi)
* peak day of year         = phi / (2 pi) * 365
* SOS / EOS                = first / last day the fitted curve exceeds an
                             NDVI threshold (0.4; 0.5 used as a fallback, as
                             in the original analysis), evaluated on a
                             1000-step grid over one year.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

GRID_SIZE = 1000


def fit_harmonic(ndvi: np.ndarray, t: np.ndarray) -> tuple[np.ndarray, float]:
    """Ordinary least squares fit; returns coefficients (b0, b1, b2) and R^2."""
    X = np.column_stack([np.ones_like(t), np.cos(2 * np.pi * t), np.sin(2 * np.pi * t)])
    ok = np.isfinite(ndvi)
    beta, *_ = np.linalg.lstsq(X[ok], ndvi[ok], rcond=None)
    resid = ndvi[ok] - X[ok] @ beta
    r2 = 1 - resid @ resid / np.sum((ndvi[ok] - ndvi[ok].mean()) ** 2)
    return beta, float(r2)


def season_bounds(beta: np.ndarray, threshold: float = 0.4, fallback: float = 0.5) -> tuple[float, float]:
    """Start and end of season (day of year) from the fitted annual curve."""
    grid = np.linspace(0, 1, GRID_SIZE)
    curve = beta[0] + beta[1] * np.cos(2 * np.pi * grid) + beta[2] * np.sin(2 * np.pi * grid)
    above = np.flatnonzero(curve > threshold)
    if above.size == 0:
        above = np.flatnonzero(curve > fallback)
    if above.size == 0:
        return np.nan, np.nan
    return grid[above[0]] * 365, grid[above[-1]] * 365


def harmonic_metrics(ndvi: np.ndarray, t: np.ndarray, threshold: float = 0.4) -> dict:
    beta, r2 = fit_harmonic(ndvi, t)
    amplitude = float(np.hypot(beta[1], beta[2]))
    phase = float(np.arctan2(beta[2], beta[1]) % (2 * np.pi))
    sos, eos = season_bounds(beta, threshold)
    return {
        "mean_ndvi": float(np.nanmean(ndvi)),
        "ndvi_peak": float(np.nanmax(ndvi)),
        "harmonic_peak": float(beta[0] + amplitude),
        "baseline": float(beta[0]),
        "amplitude": amplitude,
        "normalized_amplitude": amplitude / beta[0],
        "phase_rad": phase,
        "harmonic_r2": r2,
        "peak_doy": phase / (2 * np.pi) * 365,
        "sos_doy": sos,
        "eos_doy": eos,
    }


def field_metrics(long_df: pd.DataFrame, field_col: str = "field_id", threshold: float = 0.4) -> pd.DataFrame:
    """Apply :func:`harmonic_metrics` to every field of a long table with
    columns ``field_col, year, month, ndvi``."""
    rows = []
    for fid, g in long_df.groupby(field_col):
        t = (g["year"] + (g["month"] - 1) / 12).to_numpy(float)
        rows.append({field_col: fid, **harmonic_metrics(g["ndvi"].to_numpy(float), t, threshold)})
    return pd.DataFrame(rows)


# --- Categorical helpers -------------------------------------------------------

def compass_class(deg: float) -> str | None:
    """Four-sector compass class (N, E, S, W) for an angle in degrees."""
    if pd.isna(deg):
        return None
    deg = deg % 360
    if deg < 45 or deg >= 315:
        return "North"
    if deg < 135:
        return "East"
    if deg < 225:
        return "South"
    return "West"


def usda_texture(sand: float, silt: float, clay: float) -> str | None:
    """USDA soil texture class from particle-size fractions (%)."""
    if any(pd.isna(v) for v in (sand, silt, clay)):
        return None
    if silt + 1.5 * clay < 15:
        return "Sand"
    if silt + 2 * clay < 30:
        return "Loamy Sand"
    if clay >= 40:
        if silt >= 40:
            return "Silty Clay"
        return "Sandy Clay" if sand > 45 else "Clay"
    if clay >= 35 and sand > 45:
        return "Sandy Clay"
    if clay >= 27:
        if sand <= 20:
            return "Silty Clay Loam"
        return "Clay Loam" if sand <= 45 else "Sandy Clay Loam"
    if clay >= 20 and silt < 28 and sand > 45:
        return "Sandy Clay Loam"
    if silt >= 80 and clay < 12:
        return "Silt"
    if silt >= 50:
        return "Silt Loam"
    if clay >= 7 and silt >= 28 and sand <= 52:
        return "Loam"
    return "Sandy Loam"


def wind_from_uv(u: np.ndarray, v: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Wind speed (m/s) and meteorological direction (degrees the wind blows
    FROM) from ERA5 u/v components, following ECMWF conventions."""
    speed = np.hypot(u, v)
    direction = np.mod(180 + np.degrees(np.arctan2(u, v)), 360)
    return speed, direction
