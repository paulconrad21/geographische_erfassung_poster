# Datenhinweise

Dieser Ordner enthaelt kleine Eingangsdaten fuer das Projekt. Grosse Rasterdaten, insbesondere Sentinel-2-Mosaike als `.tif`, werden wegen der Dateigroesse **nicht** auf GitHub gespeichert.

## Erwartete Dateien

Das aktuelle Skript `scripts/skript_v3.R` erwartet folgende Datenpfade:

```text
data/training_polygons.shp
data/S2_NDS_2025_Q_20m_mosaic.tif
```

Die Trainingspolygone liegen als Shapefile mit Sidecar-Dateien vor. Das Raster `S2_NDS_2025_Q_20m_mosaic.tif` muss lokal erzeugt bzw. aus Google Earth Engine exportiert werden.

## Sentinel-2-Raster aus Google Earth Engine

Das Raster sollte die im Skript verwendeten Layernamen enthalten:

```text
NDVI_Q1, NDVI_Q2, NDVI_Q3, NDVI_Q4
NDWI_Q1, NDWI_Q2, NDWI_Q3, NDWI_Q4
NDMI_Q1, NDMI_Q2, NDMI_Q3, NDMI_Q4
BSI_Q1,  BSI_Q2,  BSI_Q3,  BSI_Q4
```

Empfohlener Export:

- Gebiet: Niedersachsen
- Sensor: Sentinel-2 Surface Reflectance
- Zeitraum: 2025
- Aufloesung: 20 m
- Komposit: wolkenreduziertes Quartalsmosaik, z.B. Median
- Exportformat: GeoTIFF
- lokaler Dateiname: `data/S2_NDS_2025_Q_20m_mosaic.tif`

## Beispiel-Logik fuer Google Earth Engine

Der folgende Code ist als Vorlage gedacht. Die Variable `nds` muss in Google Earth Engine durch eine eigene Niedersachsen-Grenze ersetzt werden, z.B. ein importiertes Shapefile/Asset.

```javascript
// Beispiel: eigene Niedersachsen-Geometrie/Asset einsetzen
var nds = /* color: #d63000 */ ee.Geometry.Polygon([]);

var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(nds)
  .filterDate('2025-01-01', '2025-12-31')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 40));

function maskS2clouds(image) {
  var scl = image.select('SCL');
  var mask = scl.neq(3)   // cloud shadow
    .and(scl.neq(8))      // cloud medium probability
    .and(scl.neq(9))      // cloud high probability
    .and(scl.neq(10))     // cirrus
    .and(scl.neq(11));    // snow/ice
  return image.updateMask(mask);
}

function addIndices(image) {
  var ndvi = image.normalizedDifference(['B8', 'B4']).rename('NDVI');
  var ndwi = image.normalizedDifference(['B3', 'B8']).rename('NDWI');
  var ndmi = image.normalizedDifference(['B8', 'B11']).rename('NDMI');
  var bsi = image.expression(
    '((swir + red) - (nir + blue)) / ((swir + red) + (nir + blue))',
    {
      swir: image.select('B11'),
      red: image.select('B4'),
      nir: image.select('B8'),
      blue: image.select('B2')
    }
  ).rename('BSI');
  return image.addBands([ndvi, ndwi, ndmi, bsi]);
}

function quarterlyComposite(start, end, suffix) {
  var composite = s2
    .filterDate(start, end)
    .map(maskS2clouds)
    .map(addIndices)
    .select(['NDVI', 'NDWI', 'NDMI', 'BSI'])
    .median()
    .clip(nds);

  return composite.rename([
    'NDVI_' + suffix,
    'NDWI_' + suffix,
    'NDMI_' + suffix,
    'BSI_' + suffix
  ]);
}

var q1 = quarterlyComposite('2025-01-01', '2025-04-01', 'Q1');
var q2 = quarterlyComposite('2025-04-01', '2025-07-01', 'Q2');
var q3 = quarterlyComposite('2025-07-01', '2025-10-01', 'Q3');
var q4 = quarterlyComposite('2025-10-01', '2026-01-01', 'Q4');

var stack = q1.addBands(q2).addBands(q3).addBands(q4);

Export.image.toDrive({
  image: stack,
  description: 'S2_NDS_2025_Q_20m_mosaic',
  folder: 'GEE_exports',
  fileNamePrefix: 'S2_NDS_2025_Q_20m_mosaic',
  region: nds,
  scale: 20,
  crs: 'EPSG:25832',
  maxPixels: 1e13
});
```

## Hinweis zu GitHub

Die `.gitignore` schliesst grosse Rasterdateien wie `.tif` aus. Das ist Absicht. Wer das Projekt reproduzieren moechte, muss das Raster lokal mit dem oben beschriebenen Namen in `data/` ablegen.

