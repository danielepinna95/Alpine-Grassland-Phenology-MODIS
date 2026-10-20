"""Step 1: reshape the GEE export and build per-field climate/terrain/soil attributes.

Input : the wide CSV exported by gee/extract_field_data.js (one row per field,
        monthly columns <VAR>_<YEAR>_<MONTH>).
Output: <outdir>/timeseries_long.csv   one row per field x month
        <outdir>/field_attributes.csv  one row per field

Per-field climate summaries over the whole period (2001-2023 in the paper):
  mean_temperature_k     mean of monthly 2 m air temperature (K)
  net_radiation_j_m2     mean of monthly surface net solar radiation (J m-2)
  precip_total_m         total precipitation over the period (m)
  mean_annual_precip_mm  precip_total_m / n_years * 1000
  snowfall_total_m_we    total snowfall over the period (m water equivalent)
  wind_speed_m_s         mean of monthly wind speed from u/v
  wind_direction_deg     circular mean of monthly wind direction (deg, meteorological)

Usage: python python/01_prepare.py field_timeseries_wide.csv data/
"""

from __future__ import annotations

import argparse
import sys
import re
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from phenology import compass_class, usda_texture, wind_from_uv

SERIES = ["NDVI", "AirTemp", "Precipitation", "NetRad", "Snow", "UWind", "VWind"]
COL_RE = re.compile(r"^(%s)_(\d{4})_(\d{1,2})$" % "|".join(SERIES))
NDVI_SCALE = 1e-4  # MOD13Q1 NDVI is stored x 10,000


def to_long(wide: pd.DataFrame) -> pd.DataFrame:
    series_cols = [c for c in wide.columns if COL_RE.match(c)]
    long = wide[["field_id"] + series_cols].melt(id_vars="field_id", var_name="col", value_name="value")
    parts = long["col"].str.extract(COL_RE)
    long["var"], long["year"], long["month"] = parts[0], parts[1].astype(int), parts[2].astype(int)
    long["value"] = pd.to_numeric(long["value"], errors="coerce")
    out = long.pivot_table(index=["field_id", "year", "month"], columns="var", values="value").reset_index()
    out.columns.name = None
    out = out.rename(columns={"NDVI": "ndvi", "AirTemp": "air_temp_k", "Precipitation": "precip_m",
                              "NetRad": "net_rad_j_m2", "Snow": "snowfall_m_we",
                              "UWind": "u_wind", "VWind": "v_wind"})
    out["ndvi"] = out["ndvi"] * NDVI_SCALE
    out["wind_speed"], out["wind_direction"] = wind_from_uv(out["u_wind"].to_numpy(), out["v_wind"].to_numpy())
    return out.sort_values(["field_id", "year", "month"])


def circular_mean_deg(angles: pd.Series) -> float:
    """Mean of angles in degrees (an arithmetic mean fails across 0/360)."""
    a = np.radians(angles.dropna().to_numpy())
    return float(np.degrees(np.arctan2(np.sin(a).mean(), np.cos(a).mean())) % 360) if a.size else np.nan


def climate_summary(long: pd.DataFrame) -> pd.DataFrame:
    n_years = long["year"].nunique()
    g = long.groupby("field_id")
    s = pd.DataFrame({
        "mean_temperature_k": g["air_temp_k"].mean(),
        "net_radiation_j_m2": g["net_rad_j_m2"].mean(),
        "precip_total_m": g["precip_m"].sum(),
        "snowfall_total_m_we": g["snowfall_m_we"].sum(),
        "wind_speed_m_s": g["wind_speed"].mean(),
        "wind_direction_deg": g["wind_direction"].agg(circular_mean_deg),
    })
    s["mean_annual_precip_mm"] = s["precip_total_m"] / n_years * 1000
    s["wind_direction_class"] = s["wind_direction_deg"].map(compass_class)
    return s.reset_index()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("wide_csv", type=Path)
    ap.add_argument("outdir", type=Path)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    wide = pd.read_csv(args.wide_csv, low_memory=False)
    long = to_long(wide)
    long.to_csv(args.outdir / "timeseries_long.csv", index=False)

    static = wide[[c for c in wide.columns if not COL_RE.match(c) and not c.startswith(("system:", ".geo"))]].copy()
    static["aspect_class"] = static["meanAspect"].map(compass_class)
    if {"soil_clay", "soil_sand"} <= set(static.columns):
        static["soil_silt"] = 100 - static["soil_clay"] - static["soil_sand"]
        static["soil_texture"] = [usda_texture(s, si, c) for s, si, c in
                                  zip(static["soil_sand"], static["soil_silt"], static["soil_clay"])]

    attrs = static.merge(climate_summary(long), on="field_id", how="left")
    attrs.to_csv(args.outdir / "field_attributes.csv", index=False)
    print(f"{long['field_id'].nunique()} fields, {len(long)} monthly rows -> {args.outdir}")


if __name__ == "__main__":
    main()
