// =============================================================================
// Google Earth Engine: per-field NDVI, climate, topography and soil (2001-2023)
// =============================================================================
// Paste into the GEE Code Editor. Set `fields` to your FeatureCollection of
// grassland polygons (each with a unique `field_id` property) and run.
// The export is one CSV row per field with monthly columns named
// <VAR>_<YEAR>_<MONTH>, read by python/01_prepare.py.
//
// Datasets
//   NDVI        MODIS/061/MOD13Q1                        (250 m, 16-day -> monthly mean, scaled x 10,000)
//   Climate     ECMWF/ERA5_LAND/MONTHLY_AGGR             (~11 km, monthly)
//   Topography  NASA/NASADEM_HGT/001                     (30 m)
//   Soil        OpenLandMap, band b0 = 0 cm (topsoil)    (250 m)
//
// Climate variables are averaged over each polygon (ee.Reducer.mean).

var fields = ee.FeatureCollection('users/YOUR_USER/YOUR_FIELDS');   // <-- edit
var START_YEAR = 2001;
var END_YEAR   = 2023;

// --- Topography ------------------------------------------------------------------
// Aspect is an angle, so it is averaged as a circular mean (via sin/cos);
// an arithmetic mean would turn north-facing fields (350 deg, 10 deg) into 180 deg.
var dem     = ee.Image('NASA/NASADEM_HGT/001').select('elevation');
var aspRad  = ee.Terrain.aspect(dem).multiply(Math.PI / 180);
var terrain = dem.addBands(ee.Terrain.slope(dem))
                 .addBands(aspRad.sin().rename('aspect_sin'))
                 .addBands(aspRad.cos().rename('aspect_cos'));

var addTerrain = function (f) {
  var s = terrain.reduceRegion({reducer: ee.Reducer.mean(), geometry: f.geometry(), scale: 30, maxPixels: 1e9});
  var aspect = ee.Number(s.get('aspect_sin')).atan2(ee.Number(s.get('aspect_cos')))
                 .multiply(180 / Math.PI).add(360).mod(360);
  return f.set({meanElevation: s.get('elevation'), meanSlope: s.get('slope'), meanAspect: aspect});
};

// --- Soil (OpenLandMap, topsoil) ---------------------------------------------------
var soil = ee.Image.cat([
  ee.Image('OpenLandMap/SOL/SOL_ORGANIC-CARBON_USDA-6A1C_M/v02').select('b0').rename('soil_soc'),
  ee.Image('OpenLandMap/SOL/SOL_PH-H2O_USDA-4C1A2A_M/v02').select('b0').rename('soil_ph'),
  ee.Image('OpenLandMap/SOL/SOL_WATERCONTENT-33KPA_USDA-4B1C_M/v01').select('b0').rename('soil_wc33'),
  ee.Image('OpenLandMap/SOL/SOL_BULKDENS-FINEEARTH_USDA-4A1H_M/v02').select('b0').rename('soil_bd'),
  ee.Image('OpenLandMap/SOL/SOL_CLAY-WFRACTION_USDA-3A1A1A_M/v02').select('b0').rename('soil_clay'),
  ee.Image('OpenLandMap/SOL/SOL_SAND-WFRACTION_USDA-3A1A1A_M/v02').select('b0').rename('soil_sand')
]);

var addSoil = function (f) {
  var s = soil.reduceRegion({reducer: ee.Reducer.mean(), geometry: f.geometry(), scale: 250, maxPixels: 1e9});
  return f.set(s);
};

// --- Monthly image stacks -------------------------------------------------------------
var months = ee.List.sequence(0, (END_YEAR - START_YEAR + 1) * 12 - 1).map(function (i) {
  return ee.Date.fromYMD(START_YEAR, 1, 1).advance(i, 'month');
});

var monthlyStack = function (collection, bands, newNames, scale) {
  var images = months.map(function (d) {
    d = ee.Date(d);
    var suffix = ee.String('_').cat(d.format('YYYY')).cat('_').cat(ee.Number(d.get('month')).format());
    var img = ee.ImageCollection(collection).filterDate(d, d.advance(1, 'month')).select(bands).mean();
    return img.rename(ee.List(newNames).map(function (n) { return ee.String(n).cat(suffix); }));
  });
  return ee.ImageCollection.fromImages(images).toBands()
    .rename(ee.ImageCollection.fromImages(images).toBands().bandNames()
      .map(function (b) { return ee.String(b).replace('^[0-9]+_', ''); }));
};

var ndviStack = monthlyStack('MODIS/061/MOD13Q1', ['NDVI'], ['NDVI'], 250);
var era5Stack = monthlyStack('ECMWF/ERA5_LAND/MONTHLY_AGGR',
  ['temperature_2m', 'total_precipitation_sum', 'surface_net_solar_radiation_sum', 'snowfall_sum',
   'u_component_of_wind_10m', 'v_component_of_wind_10m'],
  ['AirTemp', 'Precipitation', 'NetRad', 'Snow', 'UWind', 'VWind'], 1000);

var addSeries = function (f) {
  var n = ndviStack.reduceRegion({reducer: ee.Reducer.mean(), geometry: f.geometry(), scale: 250, maxPixels: 1e9});
  var c = era5Stack.reduceRegion({reducer: ee.Reducer.mean(), geometry: f.geometry(), scale: 1000, maxPixels: 1e9});
  return f.set(n).set(c);
};

var out = fields.map(addTerrain).map(addSoil).map(addSeries);

Export.table.toDrive({
  collection: out.map(function (f) { return f.setGeometry(null); }),
  description: 'field_timeseries_wide',
  fileFormat: 'CSV'
});
