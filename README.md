# Niedermoor-Klassifikation Niedersachsen

Dieses Uni-Projekt untersucht, mit welchen Sentinel-2-basierten Variablen sich Niedermoore in Niedersachsen von anderen Landnutzungsklassen unterscheiden lassen. Ziel ist eine binaere Klassifikation: **Niedermoor** vs. **kein Niedermoor**.

Als Trainingsdaten dienen kuratierte Trainingspolygone mit beiden Klassen. Das aktuelle Hauptskript ist `scripts/skript_v3.R`. Es extrahiert Rasterwerte aus den Trainingspolygonen, berechnet Zusatzfeatures, trainiert ein Random-Forest-Modell mit raeumlicher Cross Validation, waehlt einen F1-optimierten Schwellenwert und erzeugt Karten fuer Moor-Wahrscheinlichkeit, binaere Klassifikation und Area of Applicability.

## Projektstruktur

```text
.
+-- data/                  # kleine Eingangsdaten und Hinweise zu grossen Rasterdaten
+-- docs/                  # Methodik-Grafiken und Erklaerungsdokumente
+-- output/                # automatisch erzeugte Ergebnisse, nicht fuer Git gedacht
+-- scripts/               # R-Skripte der Analyse
+-- .gitignore
+-- GITHUB_ANLEITUNG.md
+-- README.md
```

## Aktueller Stand

Die Version 3 der Pipeline ist der aktuelle stabile Zwischenstand. Sie nutzt die vorhandenen Sentinel-2-Quartalsvariablen, prueft die Trainingsdaten, bewertet das Modell mit Spatial Cross Validation und nutzt den F1-Score als Hauptkriterium fuer den Klassifikationsschwellenwert.

Neue Variablen, Indizes oder Zeitraeume koennen spaeter ergaenzt und gegen diese V3-Basis verglichen werden.

## Daten

Grosse Rasterdaten werden nicht direkt in GitHub gespeichert. Hinweise zum benoetigten Raster und zum Export ueber Google Earth Engine stehen in `data/README_data.md`.

## Reproduzierbarkeit

Die wichtigsten Ergebnisse und Einstellungen werden beim Ausfuehren des Skripts automatisch in `output/v3/` gespeichert, unter anderem Metriken, Feature Importance, Kontrollplots, Karten und `run_info.txt` mit Paket- und Sitzungsinformationen.
