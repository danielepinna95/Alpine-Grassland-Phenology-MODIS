"""Step 2: harmonic regression per field and phenology metrics.

Input : timeseries_long.csv (field_id, year, month, ndvi, ...) from step 1
        field_attributes.csv (optional) to join terrain, climate and soil
Output: field_metrics.csv, one row per field

Usage: python python/02_phenology_metrics.py data/timeseries_long.csv data/field_metrics_new.csv \
           --attributes data/field_attributes.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from phenology import field_metrics


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("timeseries", type=Path)
    ap.add_argument("output", type=Path)
    ap.add_argument("--attributes", type=Path, help="per-field attributes to join (from step 1)")
    ap.add_argument("--threshold", type=float, default=0.4, help="NDVI threshold for SOS/EOS (default 0.4)")
    args = ap.parse_args()

    long = pd.read_csv(args.timeseries)
    metrics = field_metrics(long, threshold=args.threshold)
    if args.attributes:
        metrics = pd.read_csv(args.attributes).merge(metrics, on="field_id", how="inner")
    metrics.to_csv(args.output, index=False)
    print(f"Wrote metrics for {len(metrics)} fields to {args.output}")


if __name__ == "__main__":
    main()
