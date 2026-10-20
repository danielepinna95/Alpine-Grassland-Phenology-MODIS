# Data

## `field_metrics.csv` — the paper's dataset (889 fields)

One row per field. Permanent pastures and meadows > 6 ha in the Province of Trento, 2001–2023.

| Column | Description |
|---|---|
| `field_id` | Field identifier (`ID_New` in the original analysis) |
| `source_id` | Original polygon identifier |
| `sub_area` | Sub-region of the province: GIUD, NE, NO, SE |
| `area_ha` | Field area (ha) |
| `meanElevation` | Mean elevation, NASADEM (m a.s.l.) |
| `meanSlope` | Mean slope, NASADEM (degrees, `ee.Terrain.slope`) |
| `meanAspect` | Mean aspect, NASADEM (degrees; arithmetic mean, see main README) |
| `aspect_class` | Aspect in four sectors (North / East / South / West) |
| `mean_ndvi` | Mean monthly MODIS NDVI over 2001–2023 |
| `ndvi_peak` | Maximum monthly NDVI over 2001–2023 |
| `harmonic_peak` | Peak of the fitted harmonic curve (baseline + amplitude) |
| `baseline` | Harmonic intercept b0 |
| `amplitude` | Harmonic amplitude √(b1² + b2²) |
| `normalized_amplitude` | amplitude / baseline (seasonal variability) |
| `phase_rad` | Phase of the harmonic (radians, 0–2π) |
| `harmonic_r2` | R² of the harmonic fit |
| `peak_doy` | Day of year of the NDVI peak |
| `sos_doy`, `eos_doy` | Start / end of season: first / last day the fitted curve exceeds NDVI 0.4 |
| `mean_temperature_k` | Mean of monthly ERA5-Land 2 m temperature (K) |
| `net_radiation_j_m2` | Mean of monthly ERA5-Land surface net solar radiation (J m⁻²) |
| `precip_total_m` | Total ERA5-Land precipitation 2001–2023 (m) |
| `mean_annual_precip_mm` | Mean annual precipitation (mm) = precip_total_m / 23 × 1000 |
| `snowfall_total_m_we` | Total ERA5-Land snowfall 2001–2023 (m water equivalent) |
| `wind_speed_m_s` | Mean of monthly ERA5-Land 10 m wind speed from u/v (m s⁻¹) |
| `wind_direction_deg` | Mean prevailing wind direction (degrees, meteorological convention) |
| `wind_direction_class` | Wind direction in four sectors |
| `soil_soc` | Soil organic carbon, OpenLandMap topsoil (band b0) |
| `soil_ph` | Soil pH in H₂O, OpenLandMap topsoil |
| `soil_wc33` | Volumetric water content at 33 kPa (field capacity), OpenLandMap topsoil (%) |
| `soil_bd` | Fine-earth bulk density, OpenLandMap topsoil (OpenLandMap units) |
| `soil_clay`, `soil_sand`, `soil_silt` | Particle-size fractions, OpenLandMap topsoil (%); silt = 100 − clay − sand |
| `soil_texture` | USDA texture class |

Soil variables are reported in the units produced by the GEE export. See the OpenLandMap dataset pages for scaling factors.

## `example_gee_export.csv` — demo input (20 fields)

A wide table in the format written by `gee/extract_field_data.js`. It holds monthly `NDVI_<year>_<month>` (MODIS ×10,000), `AirTemp_`, `Precipitation_`, `NetRad_`, `Snow_`, `UWind_` and `VWind_` columns for 20 real fields from the study, plus terrain and soil attributes. Use it to try steps 01–02.

## License

CC BY 4.0. Please cite Pinna et al. (2025), *Computers and Electronics in Agriculture* 232, 110154.
