"""Step 3: Kendall's tau correlations and Kruskal-Wallis tests.

* Kendall's tau (scipy) between every environmental driver and every NDVI
  metric, with 95% bootstrap confidence intervals (1,000 resamples).
* Kruskal-Wallis H test with epsilon-squared effect size for categorical
  drivers (aspect class, prevailing wind direction, soil texture).

Input : data/field_metrics.csv
Output: outputs/kendall_tau.csv, outputs/kruskal_wallis.csv

Usage: python python/03_correlation_tests.py data/field_metrics.csv outputs/
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kendalltau, kruskal

METRICS = ["mean_ndvi", "ndvi_peak", "normalized_amplitude", "sos_doy", "peak_doy", "eos_doy"]
DRIVERS = ["meanElevation", "meanSlope", "mean_temperature_k", "net_radiation_j_m2", "mean_annual_precip_mm",
           "snowfall_total_m_we", "wind_speed_m_s", "soil_soc", "soil_ph", "soil_wc33", "soil_bd",
           "soil_clay", "soil_silt", "soil_sand"]
CATEGORICAL = ["aspect_class", "wind_direction_class", "soil_texture"]


def bootstrap_tau_ci(x: np.ndarray, y: np.ndarray, n_boot: int, rng: np.random.Generator) -> tuple[float, float]:
    n = len(x)
    taus = [kendalltau(x[i], y[i])[0] for i in (rng.integers(0, n, n) for _ in range(n_boot))]
    return tuple(np.nanpercentile(taus, [2.5, 97.5]))


def kendall_table(df: pd.DataFrame, n_boot: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    for d in DRIVERS:
        for m in METRICS:
            sub = df[[d, m]].dropna()
            tau, p = kendalltau(sub[d], sub[m])
            lo, hi = bootstrap_tau_ci(sub[d].to_numpy(), sub[m].to_numpy(), n_boot, rng)
            rows.append({"driver": d, "metric": m, "n": len(sub), "tau": tau, "p_value": p, "ci_low": lo, "ci_high": hi})
    return pd.DataFrame(rows)


def kruskal_table(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for g in CATEGORICAL:
        for m in METRICS:
            sub = df[[g, m]].dropna()
            groups = [v[m].to_numpy() for _, v in sub.groupby(g) if len(v) > 1]
            if len(groups) < 2:
                continue
            h, p = kruskal(*groups)
            n, k = len(sub), len(groups)
            rows.append({"factor": g, "metric": m, "groups": k, "n": n, "H": h, "p_value": p,
                         "epsilon_sq": (h - k + 1) / (n - k)})
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("metrics_csv", type=Path)
    ap.add_argument("outdir", type=Path)
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.metrics_csv)
    kt = kendall_table(df, args.n_boot, args.seed)
    kw = kruskal_table(df)
    kt.to_csv(args.outdir / "kendall_tau.csv", index=False)
    kw.to_csv(args.outdir / "kruskal_wallis.csv", index=False)

    pd.set_option("display.width", 160)
    print(kt.pivot(index="driver", columns="metric", values="tau").round(2).loc[DRIVERS, METRICS])
    print()
    print(kw.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
