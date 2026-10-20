# Alpine grassland phenology from 23 years of MODIS NDVI

Code and data for:

> Pinna, D., Sozzi, M., Pornaro, C., Macolino, S., Pezzuolo, A., & Marinello, F. (2025).
> **Modelling vegetation phenology in alpine grasslands of North-East Italy integrating microclimatic and meteorological drivers: a 23-year time-series analysis.**
> *Computers and Electronics in Agriculture*, 232, 110154. https://doi.org/10.1016/j.compag.2025.110154

The study covers 889 permanent pastures and meadows (> 6 ha, 828–2,580 m a.s.l.) in the Province of Trento, Italy. For each field it extracts the seasonal NDVI cycle from the 2001–2023 MODIS series and relates it to topography, climate and soil:

1. **Harmonic regression** of monthly NDVI → mean NDVI, peak NDVI, seasonal variability (normalised amplitude), day of peak, start (SOS) and end (EOS) of the growing season.
2. **Kendall's τ** (with bootstrap CIs) and **Kruskal-Wallis** tests between these metrics and their drivers.
3. **Generalised Additive Models** that predict each metric from 13 terrain, climate and soil drivers (test R² 0.70–0.91; SOS/EOS to within ≈ 11 days).

## Pipeline

```
gee/extract_field_data.js        Google Earth Engine: MODIS NDVI, ERA5-Land, NASADEM, OpenLandMap soil per field
python/01_prepare.py             wide GEE export -> monthly long table + per-field climate/terrain/soil attributes
python/02_phenology_metrics.py   harmonic regression per field -> phenology metrics
python/03_correlation_tests.py   Kendall's tau + bootstrap CIs, Kruskal-Wallis + epsilon²
R/04_gam_models.R                one GAM per metric (mgcv, cubic regression splines, REML)
python/phenology.py              shared functions (harmonic fit, SOS/EOS, USDA texture, wind from u/v)
```

## Quick start

```bash
pip install -r requirements.txt

# Statistics and models on the paper's field table
python python/03_correlation_tests.py data/field_metrics.csv outputs
Rscript R/04_gam_models.R               # needs mgcv, caret, dplyr, ggplot2

# Full pipeline on the 20-field example export
python python/01_prepare.py data/example_gee_export.csv outputs
python python/02_phenology_metrics.py outputs/timeseries_long.csv outputs/field_metrics_example.csv \
       --attributes outputs/field_attributes.csv
```

To run it on your own area, upload your field polygons to Earth Engine (one `field_id` per feature). Set `fields` in `gee/extract_field_data.js`, export the table and start from step 01.

## Reproducibility

| Check | Result |
|---|---|
| Harmonic metrics (step 02) vs. the paper's field table | identical for all 889 fields (max. difference < 1e-7) |
| Kendall's τ, elevation vs mean NDVI / NDVI peak / peak day / EOS | −0.76 / −0.65 / 0.47 / −0.72, as in the paper |
| GAM test R² (mean NDVI, NDVI peak, SOS, peak day, EOS) | 0.90, 0.69, 0.85, 0.72, 0.77 (paper: 0.91, 0.70, 0.88, 0.72, 0.78; the gap comes from the random split) |

## Notes on the cleaned code

- `data/field_metrics.csv` is the table used in the paper. Its terrain and climate columns come from the original processing, which used arithmetic means.
- The cleaned pipeline averages aspect and wind direction as **circular means**. An arithmetic mean of angles breaks across 0°/360°: a field facing 350° and 10° averages to 180°. As a result, aspect and wind direction from the new pipeline can differ from the published table, and north-facing fields can now appear.
- Soil texture is classified with the full USDA texture triangle.
- `gee/extract_field_data.js` is a consolidated rewrite of the original GEE scripts. It exports everything in one table and uses mean reducers throughout. Test it on a few polygons before a large export.

## Data

See [`data/README.md`](data/README.md). All inputs are public: MODIS MOD13Q1, ERA5-Land, NASADEM, OpenLandMap, and the field boundaries of the Province of Trento ([SIAT](https://siat.provincia.tn.it)).

## Citation and license

Please cite the paper above (see [`CITATION.cff`](CITATION.cff)). Code: MIT. Data: CC BY 4.0.

*Funding: Agritech National Research Center, European Union NextGenerationEU (PNRR – Mission 4, Component 2, Investment 1.4 – D.D. 1032 17/06/2022, CN00000022).*
