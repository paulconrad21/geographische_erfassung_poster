# ============================================================
# Niedermoor-Klassifikation mit Sentinel-2
# Version 3: saubere Checks, Spatial CV, F1-Threshold und AOA
# ============================================================

# Ziel dieser Version:
# - gleiche aktuellen Variablen wie v2 verwenden
# - Trainingsdaten und Klassen explizit pruefen
# - Spatial-CV-Vorhersagen speichern
# - finalen Schwellenwert anhand des F1-Scores aus der Spatial CV waehlen
# - Metriken, Threshold-Kurve, Feature Importance und AOA nachvollziehbar speichern

setwd(this.path::here())

# ------------------------------------------------------------
# 0. Pakete laden
# ------------------------------------------------------------

library(terra)
library(sf)
library(caret)
library(blockCV)
library(CAST)
library(pROC)
library(ggplot2)
library(viridis)
library(this.path)

# ------------------------------------------------------------
# 1. Pfade und Einstellungen
# ------------------------------------------------------------

raster_path <- "data/S2_NDS_2025_Q_20m_mosaic.tif"
shape_path  <- "data/training_polygons.shp"
output_dir  <- "output/v3"
plot_dir    <- file.path(output_dir, "plots")

# Namen der Klassenspalte und der beiden Klassen so, wie sie im Shapefile stehen.
# Danach werden sie intern in kurze, modellfreundliche Labels umcodiert.
class_col   <- "Class"
moor_label  <- "Niedermoor"
other_label <- "Kein Niedermoor"

positive_class <- "moor"
negative_class <- "other"

# block_size ist die Kantenlaenge der raeumlichen CV-Bloecke in Metern,
# weil das Raster in einem metrischen Koordinatensystem liegt. Groessere
# Bloecke machen die Validierung strenger, weil Trainings- und Testflaechen
# raeumlich staerker getrennt werden.
block_size <- 80000
k_folds    <- 5
seed_value <- 42

# Threshold wird nicht manuell festgelegt, sondern auf den Spatial-CV-
# Vorhersagen ueber den F1-Score optimiert. Die Karte kann spaeter trotzdem
# mit einem konservativeren Threshold interpretiert werden, falls zu viele
# Flaechen als Niedermoor markiert werden.
threshold_grid <- seq(0.01, 0.99, by = 0.01)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

set.seed(seed_value)

# ------------------------------------------------------------
# 2. Hilfsfunktionen
# ------------------------------------------------------------

# Kleine Schutzfunktion fuer Kennzahlen wie Precision oder Recall:
# Falls im Nenner 0 steht, wird NA statt Inf/Fehler zurueckgegeben.
safe_div <- function(numerator, denominator) {
  ifelse(denominator == 0, NA_real_, numerator / denominator)
}

# Gibt im R-Log klare Abschnittsueberschriften aus, damit man spaeter
# schneller sieht, an welcher Stelle der Pipeline etwas passiert ist.
message_step <- function(text) {
  cat("\n")
  cat("============================================================\n")
  cat(text, "\n")
  cat("============================================================\n")
}

# Speichert einen terra-Rasterplot als PNG und zeigt ihn direkt zusaetzlich
# in R/RStudio an. So entstehen automatisch Poster-/Kontrollgrafiken.
save_raster_plot <- function(x, filename, main, col = viridis::viridis(100)) {
  png(file.path(plot_dir, filename), width = 1800, height = 1200, res = 180)
  terra::plot(x, col = col, main = main)
  dev.off()

  terra::plot(x, col = col, main = main)
}

# Speichert "normale" base-R-Plots. substitute() und parent.frame() sorgen
# dafuer, dass der Plot-Ausdruck exakt in der Umgebung ausgewertet wird,
# in der die Funktion aufgerufen wurde.
save_base_plot <- function(filename, expr) {
  plot_expr <- substitute(expr)
  plot_env <- parent.frame()

  png(file.path(plot_dir, filename), width = 1800, height = 1200, res = 180)
  eval(plot_expr, envir = plot_env)
  dev.off()

  eval(plot_expr, envir = plot_env)
}

# Einheitliche Speicherfunktion fuer ggplot-Grafiken.
# print() zeigt den Plot im Plotfenster, ggsave() schreibt ihn in den Ordner.
save_ggplot <- function(plot_object, filename, width = 8, height = 5) {
  print(plot_object)
  ggplot2::ggsave(
    filename = file.path(plot_dir, filename),
    plot = plot_object,
    width = width,
    height = height,
    dpi = 300
  )
}

# Berechnet die wichtigsten Klassifikationsmetriken fuer einen gegebenen
# Wahrscheinlichkeits-Schwellenwert. Aus "prob_moor" wird also erst hier
# eine harte Klasse: ab threshold = Niedermoor, darunter = other.
calc_metrics <- function(obs, prob_moor, threshold) {
  pred <- ifelse(prob_moor >= threshold, positive_class, negative_class)
  pred <- factor(pred, levels = c(positive_class, negative_class))
  obs  <- factor(obs,  levels = c(positive_class, negative_class))

  # Confusion-Matrix-Zaehler:
  # TP = richtig als Niedermoor erkannt, FP = faelschlich als Niedermoor,
  # TN = richtig als other erkannt, FN = verpasstes Niedermoor.
  tp <- sum(pred == positive_class & obs == positive_class, na.rm = TRUE)
  fp <- sum(pred == positive_class & obs == negative_class, na.rm = TRUE)
  tn <- sum(pred == negative_class & obs == negative_class, na.rm = TRUE)
  fn <- sum(pred == negative_class & obs == positive_class, na.rm = TRUE)

  precision <- safe_div(tp, tp + fp)
  recall    <- safe_div(tp, tp + fn)
  specificity <- safe_div(tn, tn + fp)
  f1 <- safe_div(2 * precision * recall, precision + recall)

  data.frame(
    threshold = threshold,
    TP = tp,
    FP = fp,
    TN = tn,
    FN = fn,
    Accuracy = safe_div(tp + tn, tp + fp + tn + fn),
    Balanced_Accuracy = mean(c(recall, specificity), na.rm = TRUE),
    Precision = precision,
    Recall = recall,
    Sensitivity = recall,
    Specificity = specificity,
    F1 = f1
  )
}

# Testet viele moegliche Thresholds und waehlt den mit dem besten F1-Score.
# F1 ist hier sinnvoll, weil er Precision und Recall gemeinsam bewertet:
# also einerseits wenige False Positives, andererseits wenige verpasste Moore.
choose_f1_threshold <- function(obs, prob_moor, thresholds) {
  metric_list <- lapply(thresholds, function(th) calc_metrics(obs, prob_moor, th))
  threshold_metrics <- do.call(rbind, metric_list)

  # Primaeres Kriterium: maximaler F1.
  # Bei Gleichstand: hoeherer Recall fuer Niedermoor, danach Threshold naeher an 0.5.
  threshold_metrics$distance_to_05 <- abs(threshold_metrics$threshold - 0.5)
  threshold_metrics <- threshold_metrics[order(
    -threshold_metrics$F1,
    -threshold_metrics$Recall,
    threshold_metrics$distance_to_05
  ), ]

  list(
    best = threshold_metrics[1, setdiff(names(threshold_metrics), "distance_to_05")],
    all  = threshold_metrics[, setdiff(names(threshold_metrics), "distance_to_05")]
  )
}

# Feature Engineering fuer die Trainings-Tabelle:
# Aus den Quartals-Indizes werden Mittelwerte, Amplituden und saisonale
# Differenzen gebildet. Diese Variablen beschreiben nicht nur den Zustand
# in einem Quartal, sondern auch die jahreszeitliche Dynamik.
add_features_df <- function(df) {
  df$NDVI_mean <- rowMeans(df[, c("NDVI_Q1", "NDVI_Q2", "NDVI_Q3", "NDVI_Q4")], na.rm = TRUE)
  df$NDWI_mean <- rowMeans(df[, c("NDWI_Q1", "NDWI_Q2", "NDWI_Q3", "NDWI_Q4")], na.rm = TRUE)
  df$NDMI_mean <- rowMeans(df[, c("NDMI_Q1", "NDMI_Q2", "NDMI_Q3", "NDMI_Q4")], na.rm = TRUE)
  df$BSI_mean  <- rowMeans(df[, c("BSI_Q1",  "BSI_Q2",  "BSI_Q3",  "BSI_Q4")],  na.rm = TRUE)

  df$NDVI_amp <- apply(df[, c("NDVI_Q1", "NDVI_Q2", "NDVI_Q3", "NDVI_Q4")], 1, max, na.rm = TRUE) -
    apply(df[, c("NDVI_Q1", "NDVI_Q2", "NDVI_Q3", "NDVI_Q4")], 1, min, na.rm = TRUE)

  df$NDMI_amp <- apply(df[, c("NDMI_Q1", "NDMI_Q2", "NDMI_Q3", "NDMI_Q4")], 1, max, na.rm = TRUE) -
    apply(df[, c("NDMI_Q1", "NDMI_Q2", "NDMI_Q3", "NDMI_Q4")], 1, min, na.rm = TRUE)

  df$NDVI_Q3_minus_Q1 <- df$NDVI_Q3 - df$NDVI_Q1
  df$NDMI_Q3_minus_Q1 <- df$NDMI_Q3 - df$NDMI_Q1

  df
}

# Dasselbe Feature Engineering wie oben, aber fuer den Rasterstack.
# Wichtig: Training und Karte muessen exakt dieselben Predictor bekommen,
# sonst lernt das Modell etwas anderes, als spaeter kartiert wird.
add_features_raster <- function(x) {
  x$NDVI_mean <- mean(x[[c("NDVI_Q1", "NDVI_Q2", "NDVI_Q3", "NDVI_Q4")]], na.rm = TRUE)
  x$NDWI_mean <- mean(x[[c("NDWI_Q1", "NDWI_Q2", "NDWI_Q3", "NDWI_Q4")]], na.rm = TRUE)
  x$NDMI_mean <- mean(x[[c("NDMI_Q1", "NDMI_Q2", "NDMI_Q3", "NDMI_Q4")]], na.rm = TRUE)
  x$BSI_mean  <- mean(x[[c("BSI_Q1",  "BSI_Q2",  "BSI_Q3",  "BSI_Q4")]],  na.rm = TRUE)

  x$NDVI_amp <- max(x[[c("NDVI_Q1", "NDVI_Q2", "NDVI_Q3", "NDVI_Q4")]], na.rm = TRUE) -
    min(x[[c("NDVI_Q1", "NDVI_Q2", "NDVI_Q3", "NDVI_Q4")]], na.rm = TRUE)

  x$NDMI_amp <- max(x[[c("NDMI_Q1", "NDMI_Q2", "NDMI_Q3", "NDMI_Q4")]], na.rm = TRUE) -
    min(x[[c("NDMI_Q1", "NDMI_Q2", "NDMI_Q3", "NDMI_Q4")]], na.rm = TRUE)

  x$NDVI_Q3_minus_Q1 <- x$NDVI_Q3 - x$NDVI_Q1
  x$NDMI_Q3_minus_Q1 <- x$NDMI_Q3 - x$NDMI_Q1

  x
}

# ------------------------------------------------------------
# 3. Raster und Trainingspolygone laden
# ------------------------------------------------------------

message_step("1) Raster und Trainingspolygone laden")

# Raster = flaechendeckende Predictor-Layer fuer Niedersachsen.
# train_sf = Polygone mit bekannter Klasse, aus denen Trainingswerte
# fuer das Modell extrahiert werden.
features <- terra::rast(raster_path)
train_sf <- sf::st_read(shape_path)

print(features)
print(names(features))
print(sf::st_crs(train_sf))
print(terra::crs(features))

cat("\nRasterdimensionen:\n")
print(dim(features))
cat("\nRasterlayer:\n")
print(names(features))

save_base_plot(
  "01_raster_overview_first_layers.png",
  {
    terra::plot(
      features[[1:min(4, terra::nlyr(features))]],
      col = viridis::viridis(100),
      main = names(features)[1:min(4, terra::nlyr(features))]
    )
  }
)

# Geometrien werden repariert und in dasselbe Koordinatensystem wie das
# Raster gebracht. Ohne gleiche CRS wuerden Extraktion und Plot raeumlich
# nicht korrekt zusammenpassen.
train_sf <- sf::st_make_valid(train_sf)
train_sf <- sf::st_transform(train_sf, terra::crs(features))

# ------------------------------------------------------------
# 4. Daten- und Klassenchecks
# ------------------------------------------------------------

message_step("2) Daten- und Klassenchecks")

# Diese Checks lassen das Skript frueh und klar abbrechen, falls die
# Trainingsdaten anders aussehen als erwartet. Das ist besser, als spaeter
# mit falsch codierten Klassen ein scheinbar gutes Modell zu trainieren.
if (!class_col %in% names(train_sf)) {
  stop("Die Klassenspalte fehlt im Trainingsdatensatz: ", class_col)
}

observed_classes <- sort(unique(as.character(train_sf[[class_col]])))
expected_classes <- sort(c(moor_label, other_label))

unexpected_classes <- setdiff(observed_classes, expected_classes)
missing_classes <- setdiff(expected_classes, observed_classes)

if (length(unexpected_classes) > 0) {
  stop("Unerwartete Klassen gefunden: ", paste(unexpected_classes, collapse = ", "))
}

if (length(missing_classes) > 0) {
  stop("Erwartete Klassen fehlen: ", paste(missing_classes, collapse = ", "))
}

train_sf$class_bin <- ifelse(
  as.character(train_sf[[class_col]]) == moor_label,
  positive_class,
  negative_class
)
# Reihenfolge der Faktor-Level ist wichtig: caret/pROC sollen "moor" als
# positive Klasse behandeln, weil dafuer Recall, Precision und F1 berechnet
# und optimiert werden.
train_sf$class_bin <- factor(train_sf$class_bin, levels = c(positive_class, negative_class))

class_counts <- as.data.frame(table(train_sf$class_bin))
names(class_counts) <- c("class", "n")
print(class_counts)
write.csv(class_counts, file.path(output_dir, "training_class_counts.csv"), row.names = FALSE)

save_base_plot(
  "02_training_polygons_overview.png",
  {
    terra::plot(features[[1]], col = gray.colors(100), main = "Trainingspolygone auf Raster")
    plot(
      sf::st_geometry(train_sf),
      add = TRUE,
      border = ifelse(train_sf$class_bin == positive_class, "#2f7d5a", "#c36b4c"),
      lwd = 2
    )
    legend(
      "topright",
      legend = c("Niedermoor", "Kein Niedermoor"),
      col = c("#2f7d5a", "#c36b4c"),
      lwd = 3,
      bg = "white"
    )
  }
)

class_plot <- ggplot(class_counts, aes(x = class, y = n, fill = class)) +
  geom_col(width = 0.65) +
  scale_fill_manual(values = c(moor = "#2f7d5a", other = "#c36b4c")) +
  theme_minimal() +
  labs(
    title = "Trainingsdaten: Klassenverteilung",
    x = "Klasse",
    y = "Anzahl Polygone"
  ) +
  theme(legend.position = "none")

save_ggplot(class_plot, "03_training_class_counts.png", width = 6, height = 4)

# ------------------------------------------------------------
# 5. Rasterwerte pro Polygon extrahieren
# ------------------------------------------------------------

message_step("3) Rasterwerte pro Polygon extrahieren")

# Diese Layer muessen im Raster vorhanden sein, weil sowohl Training als
# auch spaetere Rastervorhersage darauf aufbauen.
required_bands <- c(
  "NDVI_Q1", "NDVI_Q2", "NDVI_Q3", "NDVI_Q4",
  "NDWI_Q1", "NDWI_Q2", "NDWI_Q3", "NDWI_Q4",
  "NDMI_Q1", "NDMI_Q2", "NDMI_Q3", "NDMI_Q4",
  "BSI_Q1",  "BSI_Q2",  "BSI_Q3",  "BSI_Q4"
)

missing_bands_raster <- setdiff(required_bands, names(features))
if (length(missing_bands_raster) > 0) {
  stop("Diese benoetigten Rasterlayer fehlen: ", paste(missing_bands_raster, collapse = ", "))
}

train_vect <- terra::vect(train_sf)

# Pro Polygon wird der Median jedes Rasterlayers extrahiert.
# Der Median ist robuster als der Mittelwert, wenn einzelne Pixel im Polygon
# Ausreisser enthalten, z.B. durch Randpixel, Wolkenreste oder Mischpixel.
train_df_raw <- terra::extract(
  features,
  train_vect,
  fun = median,
  na.rm = TRUE
)

train_df_raw$class_bin <- train_sf$class_bin

# Polygone mit fehlenden Rasterwerten werden entfernt. Danach muessen
# train_df und train_sf synchron gekuerzt werden, damit Tabellenzeilen,
# Geometrien und CV-Folds weiterhin dieselben Trainingsfaelle meinen.
valid_rows <- complete.cases(train_df_raw)
if (sum(!valid_rows) > 0) {
  warning(sum(!valid_rows), " Trainingspolygone wurden wegen fehlender Rasterwerte entfernt.")
}

train_df <- train_df_raw[valid_rows, ]
train_sf <- train_sf[valid_rows, ]

train_df$ID <- NULL
train_df$class_bin <- droplevels(train_df$class_bin)
train_sf$class_bin <- droplevels(train_sf$class_bin)

cat("\nTrainingsfaelle nach Rasterextraktion:\n")
print(table(train_df$class_bin))
print(paste("n =", nrow(train_df)))

write.csv(train_df, file.path(output_dir, "training_extracted_medians.csv"), row.names = FALSE)

# ------------------------------------------------------------
# 6. Zusatzfeatures berechnen
# ------------------------------------------------------------

message_step("4) Zusatzfeatures berechnen und Feature-Verteilungen plotten")

# Zweiter Sicherheitscheck: Die Extraktion muss alle benoetigten Baender
# geliefert haben, bevor daraus Zusatzfeatures berechnet werden.
missing_bands_training <- setdiff(required_bands, names(train_df))
if (length(missing_bands_training) > 0) {
  stop("Diese benoetigten Trainingslayer fehlen: ", paste(missing_bands_training, collapse = ", "))
}

train_df <- add_features_df(train_df)

# Auch nach dem Feature Engineering koennen theoretisch NA-Werte entstehen.
# Diese werden entfernt, damit caret spaeter keine stillen Probleme bekommt.
valid_rows_2 <- complete.cases(train_df)
if (sum(!valid_rows_2) > 0) {
  warning(sum(!valid_rows_2), " Trainingspolygone wurden nach Feature Engineering entfernt.")
}

train_df <- train_df[valid_rows_2, ]
train_sf <- train_sf[valid_rows_2, ]
train_df$class_bin <- droplevels(train_df$class_bin)
train_sf$class_bin <- droplevels(train_sf$class_bin)

predictor_names <- setdiff(names(train_df), "class_bin")
write.csv(train_df, file.path(output_dir, "training_features_final.csv"), row.names = FALSE)

cat("\nTrainingsfaelle nach Feature Engineering:\n")
print(table(train_df$class_bin))
print(paste("n =", nrow(train_df)))
cat("\nAnzahl Predictor:\n")
print(length(predictor_names))

plot_features <- intersect(
  c("NDVI_Q3", "NDWI_Q3", "NDMI_Q3", "BSI_Q3", "NDVI_mean", "NDWI_mean", "NDMI_mean", "BSI_mean"),
  predictor_names
)

# Kontrollplot: Wenn sich die Boxplots der beiden Klassen stark ueberlappen,
# ist ein einzelner Predictor allein wahrscheinlich nicht trennscharf.
# Der Random Forest kann aber Kombinationen mehrerer Predictor nutzen.
feature_long <- do.call(
  rbind,
  lapply(plot_features, function(feature_name) {
    data.frame(
      feature = feature_name,
      value = train_df[[feature_name]],
      class_bin = train_df$class_bin
    )
  })
)

feature_distribution_plot <- ggplot(feature_long, aes(x = class_bin, y = value, fill = class_bin)) +
  geom_boxplot(alpha = 0.78, outlier.alpha = 0.45) +
  facet_wrap(~ feature, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c(moor = "#2f7d5a", other = "#c36b4c")) +
  theme_minimal() +
  labs(
    title = "Kontrollplot: Predictor-Verteilungen je Klasse",
    x = "Klasse",
    y = "Wert"
  ) +
  theme(legend.position = "none")

save_ggplot(feature_distribution_plot, "04_feature_distributions_by_class.png", width = 11, height = 6)

# ------------------------------------------------------------
# 7. Predictor-Korrelationen dokumentieren
# ------------------------------------------------------------

message_step("5) Predictor-Korrelationen berechnen")

# Stark korrelierte Predictor tragen sehr aehnliche Informationen.
# Das ist fuer Random Forest nicht automatisch falsch, aber wichtig fuer
# die Interpretation der Feature Importance: Bedeutung kann sich zwischen
# korrelierten Variablen aufteilen oder scheinbar "wandern".
cor_mat <- cor(train_df[, predictor_names], use = "pairwise.complete.obs")
cor_idx <- which(abs(cor_mat) >= 0.90 & upper.tri(cor_mat), arr.ind = TRUE)

high_correlations <- data.frame(
  feature_1 = rownames(cor_mat)[cor_idx[, "row"]],
  feature_2 = colnames(cor_mat)[cor_idx[, "col"]],
  correlation = cor_mat[cor_idx]
)

high_correlations <- high_correlations[order(-abs(high_correlations$correlation)), ]
write.csv(high_correlations, file.path(output_dir, "high_predictor_correlations.csv"), row.names = FALSE)

print("Staerkste Predictor-Korrelationen:")
print(head(high_correlations, 15))

cor_df <- as.data.frame(as.table(cor_mat))
names(cor_df) <- c("feature_1", "feature_2", "correlation")

correlation_plot <- ggplot(cor_df, aes(x = feature_1, y = feature_2, fill = correlation)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#3b6ea8",
    mid = "white",
    high = "#c04e3f",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  theme_minimal() +
  labs(
    title = "Predictor-Korrelationen",
    x = NULL,
    y = NULL,
    fill = "r"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    axis.text.y = element_text(size = 7)
  )

save_ggplot(correlation_plot, "05_predictor_correlation_matrix.png", width = 9, height = 8)

# ------------------------------------------------------------
# 8. Spatial Cross Validation erzeugen
# ------------------------------------------------------------

message_step("6) Spatial Cross Validation erzeugen")

# Spatial CV teilt die Trainingsdaten raeumlich auf. Dadurch wird getestet,
# ob das Modell auch in anderen Gegenden funktioniert und nicht nur sehr
# aehnliche Nachbarflaechen wiedererkennt. Das ist bei Fernerkundungsdaten
# methodisch deutlich aussagekraeftiger als zufaellige CV.
spatial_folds <- blockCV::cv_spatial(
  x = train_sf,
  column = "class_bin",
  k = k_folds,
  size = block_size,
  selection = "random",
  iteration = 100,
  progress = TRUE
)

folds <- spatial_folds$folds_list

# caret erwartet die Trainings- und Test-Indizes in zwei Listen:
# index = Trainingszeilen je Fold, indexOut = Testzeilen je Fold.
caret_index <- list()
caret_indexOut <- list()

for (fold_id in seq_along(folds)) {
  caret_index[[fold_id]]    <- folds[[fold_id]][[1]]
  caret_indexOut[[fold_id]] <- folds[[fold_id]][[2]]
}

fold_class_counts <- do.call(rbind, lapply(seq_along(caret_indexOut), function(i) {
  tab <- table(train_df$class_bin[caret_indexOut[[i]]])
  data.frame(
    fold = i,
    moor = as.integer(tab[positive_class]),
    other = as.integer(tab[negative_class])
  )
}))

print(fold_class_counts)
write.csv(fold_class_counts, file.path(output_dir, "spatial_cv_fold_class_counts.csv"), row.names = FALSE)

# Dieser Plot prueft, ob in jedem Test-Fold beide Klassen vorkommen.
# Fehlt eine Klasse, koennen Metriken wie AUC, Recall oder F1 unbrauchbar
# oder gar nicht berechenbar werden.
fold_plot_df <- rbind(
  data.frame(
    fold = fold_class_counts$fold,
    class = "moor",
    n = fold_class_counts$moor
  ),
  data.frame(
    fold = fold_class_counts$fold,
    class = "other",
    n = fold_class_counts$other
  )
)

fold_plot <- ggplot(fold_plot_df, aes(x = factor(fold), y = n, fill = class)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c(moor = "#2f7d5a", other = "#c36b4c")) +
  theme_minimal() +
  labs(
    title = "Spatial-CV: Klassen pro Test-Fold",
    x = "Fold",
    y = "Anzahl Trainingspolygone",
    fill = "Klasse"
  )

save_ggplot(fold_plot, "06_spatial_cv_fold_class_counts.png", width = 7, height = 4.5)

if (any(fold_class_counts$moor == 0 | fold_class_counts$other == 0)) {
  stop("Mindestens eine Test-Falte enthaelt nicht beide Klassen. Bitte k_folds oder block_size anpassen.")
}

# ------------------------------------------------------------
# 9. Random Forest mit Spatial CV trainieren
# ------------------------------------------------------------

message_step("7) Random Forest trainieren")

p <- length(predictor_names)

# In v3 wird nur mtry automatisch verglichen. mtry steuert, wie viele
# Predictor pro Split zufaellig ausprobiert werden. Die Anzahl der Baeume
# ist hier fest auf 500 gesetzt; das ist fuer Random Forest meist stabil,
# ohne das Training unnoetig gross zu machen.
tune_grid <- expand.grid(
  mtry = unique(c(2, floor(sqrt(p)), floor(p / 3), floor(p / 2)))
)

# Die vorher erzeugten Spatial-CV-Indizes werden hier an caret uebergeben.
# twoClassSummary berechnet u.a. ROC/AUC; classProbs = TRUE ist dafuer
# noetig, weil AUC Wahrscheinlichkeiten statt harter Klassen bewertet.
ctrl <- caret::trainControl(
  method = "cv",
  index = caret_index,
  indexOut = caret_indexOut,
  classProbs = TRUE,
  summaryFunction = caret::twoClassSummary,
  savePredictions = "final",
  allowParallel = TRUE
)

set.seed(seed_value)

final_rf <- caret::train(
  x = train_df[, predictor_names],
  y = train_df$class_bin,
  method = "rf",
  metric = "ROC",
  tuneGrid = tune_grid,
  trControl = ctrl,
  importance = TRUE,
  ntree = 500
)

print(final_rf)
print(final_rf$bestTune)

saveRDS(final_rf, file.path(output_dir, "final_rf_caret_model.rds"))

# ------------------------------------------------------------
# 10. Spatial-CV-Auswertung und F1-Threshold
# ------------------------------------------------------------

message_step("8) Spatial-CV auswerten und F1-Threshold bestimmen")

cv_pred <- final_rf$pred

# Falls mehrere mtry-Werte getestet wurden, enthaelt final_rf$pred alle
# CV-Vorhersagen. Fuer die finale Auswertung werden nur die Vorhersagen
# des besten mtry verwendet.
if ("mtry" %in% names(cv_pred)) {
  cv_pred_best <- cv_pred[cv_pred$mtry == final_rf$bestTune$mtry, ]
} else {
  cv_pred_best <- cv_pred
}

# AUC bewertet die Rangfolge der Wahrscheinlichkeiten ueber alle moeglichen
# Thresholds hinweg. Ein hoher AUC-Wert bedeutet: echte Niedermoore bekommen
# meist hoehere Wahrscheinlichkeiten als andere Flaechen.
roc_obj <- pROC::roc(
  response = cv_pred_best$obs,
  predictor = cv_pred_best[[positive_class]],
  levels = c(negative_class, positive_class),
  quiet = TRUE
)

cv_auc <- as.numeric(pROC::auc(roc_obj))

# Der F1-optimale Threshold wird nur aus den Spatial-CV-Vorhersagen bestimmt.
# Dadurch wird der Schwellenwert auf unabhaengiger wirkenden Testdaten
# gewaehlt und nicht direkt auf dem final trainierten Modell "zurechtgelegt".
threshold_result <- choose_f1_threshold(
  obs = cv_pred_best$obs,
  prob_moor = cv_pred_best[[positive_class]],
  thresholds = threshold_grid
)

best_threshold_metrics <- threshold_result$best
threshold_metrics <- threshold_result$all
best_threshold <- best_threshold_metrics$threshold

# Jetzt werden die Wahrscheinlichkeiten mit dem gewaehlten Threshold in
# harte Klassen umgewandelt, damit Confusion Matrix, Precision, Recall usw.
# berechnet werden koennen.
cv_pred_best$pred_threshold <- ifelse(
  cv_pred_best[[positive_class]] >= best_threshold,
  positive_class,
  negative_class
)
cv_pred_best$pred_threshold <- factor(
  cv_pred_best$pred_threshold,
  levels = c(positive_class, negative_class)
)

conf_mat <- caret::confusionMatrix(
  data = cv_pred_best$pred_threshold,
  reference = cv_pred_best$obs,
  positive = positive_class
)

print(paste("Spatial-CV AUC:", round(cv_auc, 3)))
print(paste("F1-optimaler Threshold:", round(best_threshold, 2)))
print(best_threshold_metrics)
print(conf_mat)

write.csv(cv_pred_best, file.path(output_dir, "spatial_cv_predictions.csv"), row.names = FALSE)
write.csv(threshold_metrics, file.path(output_dir, "threshold_metrics_f1.csv"), row.names = FALSE)
write.csv(best_threshold_metrics, file.path(output_dir, "best_threshold_metrics.csv"), row.names = FALSE)
write.csv(as.data.frame(conf_mat$table), file.path(output_dir, "confusion_matrix_best_f1_threshold.csv"), row.names = FALSE)

summary_metrics <- data.frame(
  feature_set = "S2_quarterly_NDVI_NDWI_NDMI_BSI_plus_summary_features",
  n_training = nrow(train_df),
  n_predictors = length(predictor_names),
  k_folds = k_folds,
  block_size_m = block_size,
  seed = seed_value,
  best_mtry = final_rf$bestTune$mtry,
  auc = cv_auc,
  best_threshold = best_threshold,
  f1 = best_threshold_metrics$F1,
  precision = best_threshold_metrics$Precision,
  recall = best_threshold_metrics$Recall,
  specificity = best_threshold_metrics$Specificity,
  balanced_accuracy = best_threshold_metrics$Balanced_Accuracy,
  accuracy = best_threshold_metrics$Accuracy
)

write.csv(summary_metrics, file.path(output_dir, "model_summary_metrics.csv"), row.names = FALSE)

metric_plot_df <- data.frame(
  metric = c("F1", "Precision", "Recall", "Specificity", "Balanced Accuracy", "Accuracy"),
  value = c(
    best_threshold_metrics$F1,
    best_threshold_metrics$Precision,
    best_threshold_metrics$Recall,
    best_threshold_metrics$Specificity,
    best_threshold_metrics$Balanced_Accuracy,
    best_threshold_metrics$Accuracy
  )
)

metric_plot <- ggplot(metric_plot_df, aes(x = reorder(metric, value), y = value, fill = metric)) +
  geom_col(width = 0.68) +
  coord_flip() +
  scale_fill_viridis_d(option = "C", end = 0.85) +
  ylim(0, 1) +
  theme_minimal() +
  labs(
    title = "Modellmetriken beim F1-optimalen Threshold",
    subtitle = paste0("Threshold = ", round(best_threshold, 2), ", AUC = ", round(cv_auc, 3)),
    x = NULL,
    y = "Wert"
  ) +
  theme(legend.position = "none")

save_ggplot(metric_plot, "07_model_metrics_best_threshold.png", width = 8, height = 4.8)

# ------------------------------------------------------------
# 11. Plots fuer Auswertung
# ------------------------------------------------------------

# ROC-Kurve: zeigt den Zielkonflikt zwischen Trefferquote fuer Niedermoor
# und False-Positive-Rate ueber viele Thresholds. Die gestrichelte Linie
# waere ein Zufallsmodell; je weiter die Kurve darueber liegt, desto besser.
roc_df <- data.frame(
  sensitivity = roc_obj$sensitivities,
  specificity = roc_obj$specificities,
  threshold = roc_obj$thresholds
)

roc_plot <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey55") +
  geom_line(linewidth = 1.1, color = "#4b6f9e") +
  theme_minimal() +
  labs(
    title = paste0("Spatial-CV ROC, AUC = ", round(cv_auc, 3)),
    x = "False Positive Rate",
    y = "True Positive Rate"
  )

# Threshold-Kurve: zeigt direkt, bei welchem Wahrscheinlichkeitswert der
# F1-Score maximal wird. Das ist der im Skript verwendete Karten-Threshold.
threshold_plot <- ggplot(threshold_metrics, aes(x = threshold, y = F1)) +
  geom_line(linewidth = 1.1, color = "#7a6aa8") +
  geom_vline(xintercept = best_threshold, linetype = "dashed", color = "#c36b4c") +
  geom_point(
    data = best_threshold_metrics,
    aes(x = threshold, y = F1),
    color = "#c36b4c",
    size = 3
  ) +
  theme_minimal() +
  labs(
    title = paste0("F1 nach Threshold, bester Threshold = ", round(best_threshold, 2)),
    x = "Threshold fuer Niedermoor-Wahrscheinlichkeit",
    y = "F1"
  )

save_ggplot(roc_plot, "08_spatial_cv_roc.png", width = 7, height = 5)
save_ggplot(threshold_plot, "09_threshold_f1_curve.png", width = 7, height = 5)

ggsave(file.path(output_dir, "spatial_cv_roc.png"), roc_plot, width = 7, height = 5, dpi = 300)
ggsave(file.path(output_dir, "threshold_f1_curve.png"), threshold_plot, width = 7, height = 5, dpi = 300)

# ------------------------------------------------------------
# 12. Feature Importance
# ------------------------------------------------------------

message_step("9) Feature Importance berechnen und plotten")

# Feature Importance zeigt, welche Predictor dem Random Forest beim Splitten
# besonders geholfen haben. Sie ist aber keine Kausalitaetsanalyse; bei
# korrelierten Features sollte man sie immer zusammen mit der Korrelations-
# matrix und fachlichem Wissen interpretieren.
importance_obj <- caret::varImp(final_rf, scale = FALSE)
importance_raw <- importance_obj$importance

if ("Overall" %in% names(importance_raw)) {
  importance_values <- importance_raw$Overall
} else {
  importance_values <- rowMeans(importance_raw, na.rm = TRUE)
}

importance_df <- data.frame(
  feature = rownames(importance_raw),
  importance = importance_values
)

importance_df <- importance_df[order(-importance_df$importance), ]

importance_plot <- ggplot(importance_df, aes(x = reorder(feature, importance), y = importance, fill = importance)) +
  geom_col() +
  coord_flip() +
  scale_fill_viridis_c(option = "C", end = 0.9) +
  theme_minimal() +
  labs(
    title = "Feature Importance",
    x = "Feature",
    y = "Importance"
  )

print("Top 15 Feature Importance:")
print(head(importance_df, 15))

save_ggplot(importance_plot, "10_feature_importance.png", width = 8, height = 10)
ggsave(file.path(output_dir, "feature_importance.png"), importance_plot, width = 8, height = 10, dpi = 300)
write.csv(importance_df, file.path(output_dir, "feature_importance.csv"), row.names = FALSE)

# ------------------------------------------------------------
# 13. Finale Raster-Prediction mit F1-Threshold
# ------------------------------------------------------------

message_step("10) Finale Moor-Wahrscheinlichkeit fuer Niedersachsen berechnen")

# Vor der Rastervorhersage werden exakt dieselben Zusatzfeatures auf dem
# Rasterstack berechnet wie vorher in der Trainingstabelle.
features <- add_features_raster(features)

# Harte Kontrolle: Das Modell darf nur auf Predictor angewendet werden,
# die im Raster mit identischem Namen vorhanden sind. Zusatzlayer sind okay,
# fehlende Modell-Predictor waeren ein echter Fehler.
missing_predictors <- setdiff(predictor_names, names(features))
extra_predictors   <- setdiff(names(features), predictor_names)

print("Fehlende Predictor im Raster:")
print(missing_predictors)

print("Zusaetzliche Rasterlayer, die nicht verwendet werden:")
print(extra_predictors)

if (length(missing_predictors) > 0) {
  stop("Diese Modell-Predictor fehlen im Raster: ", paste(missing_predictors, collapse = ", "))
}

features_model <- features[[predictor_names]]

# terra::predict wendet das trainierte Modell Pixel fuer Pixel auf ganz
# Niedersachsen an. type = "prob" erzeugt eine Wahrscheinlichkeitskarte;
# index = positive_class waehlt daraus die Spalte fuer "moor".
moor_prob <- terra::predict(
  object = features_model,
  model = final_rf,
  type = "prob",
  index = positive_class,
  na.rm = TRUE
)

names(moor_prob) <- "moor_probability"

terra::writeRaster(
  moor_prob,
  file.path(output_dir, "moor_probability_rf.tif"),
  overwrite = TRUE
)

save_raster_plot(
  moor_prob,
  "11_moor_probability_rf.png",
  paste0("Moor-Wahrscheinlichkeit, RF, Threshold = ", round(best_threshold, 2)),
  col = viridis::viridis(100, option = "C")
)

# Aus der kontinuierlichen Wahrscheinlichkeitskarte wird mit dem F1-Threshold
# eine binaere Karte. Werte TRUE/1 bedeuten: als Niedermoor klassifiziert.
moor_binary <- moor_prob >= best_threshold
names(moor_binary) <- "moor_binary"

terra::writeRaster(
  moor_binary,
  file.path(output_dir, "moor_binary_rf_f1_threshold.tif"),
  overwrite = TRUE
)

save_raster_plot(
  moor_binary,
  "12_moor_binary_rf_f1_threshold.png",
  paste0("Binaere Moor-Klassifikation, Threshold = ", round(best_threshold, 2)),
  col = c("#f4e8c1", "#2f7d5a")
)

# ------------------------------------------------------------
# 14. Area of Applicability mit CAST, speicherschonend
# ------------------------------------------------------------

message_step("11) Area of Applicability berechnen, speicherschonend")

# AOA ist sehr speicherintensiv. Deshalb wird sie hier auf einer
# groberen Rasteraufloesung berechnet. Die normale Moor-Wahrscheinlichkeitskarte
# bleibt davon unberuehrt.
# Interpretation: AOA sagt nicht, ob eine Vorhersage richtig ist. Sie sagt,
# ob ein Pixel im Merkmalsraum aehnlich zu den Trainingsdaten ist. Ausserhalb
# der AOA sollte man Modellwerte deshalb deutlich vorsichtiger lesen.

aoa_aggregation_factor <- 5   # 5 bedeutet: 20 m -> ca. 100 m
aoa_use_lpd <- FALSE          # LPD braucht zusaetzlich viel Speicher

cat("\nAOA-Einstellungen:\n")
print(paste("Aggregation factor:", aoa_aggregation_factor))
print(paste("LPD berechnen:", aoa_use_lpd))

features_model_aoa <- terra::aggregate(
  features_model,
  fact = aoa_aggregation_factor,
  fun = mean,
  na.rm = TRUE
)

cat("\nAOA-Raster nach Aggregation:\n")
print(features_model_aoa)

# CAST::aoa berechnet u.a. den Dissimilarity Index (DI). Hohe DI-Werte
# bedeuten: die Pixelkombination der Predictor liegt weiter weg von dem,
# was das Modell in den Trainingsdaten gesehen hat.
AOA <- CAST::aoa(
  newdata = features_model_aoa,
  model = final_rf,
  LPD = aoa_use_lpd,
  verbose = FALSE
)

terra::writeRaster(
  AOA$DI,
  file.path(output_dir, "AOA_DI_100m.tif"),
  overwrite = TRUE
)

terra::writeRaster(
  AOA$AOA,
  file.path(output_dir, "AOA_binary_100m.tif"),
  overwrite = TRUE
)

save_raster_plot(
  AOA$DI,
  "13_AOA_dissimilarity_index_100m.png",
  "Area of Applicability: Dissimilarity Index, ca. 100 m",
  col = viridis::viridis(100, option = "magma")
)

save_raster_plot(
  AOA$AOA,
  "14_AOA_binary_100m.png",
  "Area of Applicability, ca. 100 m",
  col = c("#d9d9d9", "#2f7d5a")
)

# Nur schreiben, falls LPD berechnet wurde.
if (!is.null(AOA$LPD)) {
  terra::writeRaster(
    AOA$LPD,
    file.path(output_dir, "AOA_LPD_100m.tif"),
    overwrite = TRUE
  )

  save_raster_plot(
    AOA$LPD,
    "15_AOA_local_point_density_100m.png",
    "Area of Applicability: Local Point Density, ca. 100 m",
    col = viridis::viridis(100, option = "plasma")
  )
}

aoa_freq <- as.data.frame(terra::freq(AOA$AOA))
print("AOA-Haeufigkeiten:")
print(aoa_freq)

write.csv(
  aoa_freq,
  file.path(output_dir, "AOA_frequency_100m.csv"),
  row.names = FALSE
)


# AOA wieder auf die Aufloesung der Moor-Wahrscheinlichkeitskarte bringen,
# damit sie als Maske verwendet werden kann.
# Dabei wird "near" verwendet, weil AOA eine kategoriale Maske ist
# (innerhalb/ausserhalb) und nicht weich interpoliert werden soll.
AOA_fullres <- terra::resample(
  AOA$AOA,
  moor_prob,
  method = "near"
)

# Maskierte Wahrscheinlichkeitskarte: zeigt nur noch Pixel innerhalb der AOA.
# Das ist oft die wissenschaftlich vorsichtigere Karte fuer Interpretation
# und Poster, weil stark extrapolierte Bereiche ausgeblendet werden.
moor_prob_AOA <- terra::mask(
  moor_prob,
  AOA_fullres,
  maskvalue = 0
)

terra::writeRaster(
  AOA_fullres,
  file.path(output_dir, "AOA_binary_resampled_to_fullres.tif"),
  overwrite = TRUE
)

terra::writeRaster(
  moor_prob_AOA,
  file.path(output_dir, "moor_probability_within_AOA.tif"),
  overwrite = TRUE
)

save_raster_plot(
  moor_prob_AOA,
  "16_moor_probability_within_AOA.png",
  "Moor-Wahrscheinlichkeit innerhalb der AOA",
  col = viridis::viridis(100, option = "C")
)

# ------------------------------------------------------------
# 15. Reproduktions-Metadaten speichern
# ------------------------------------------------------------

# run_info dokumentiert die wichtigsten Einstellungen und Paketversionen.
# Das hilft enorm, wenn ihr spaeter fuer Poster, Abgabe oder GitHub genau
# nachvollziehen wollt, mit welcher Umgebung die Ergebnisse erzeugt wurden.
sink(file.path(output_dir, "run_info.txt"))
cat("Niedermoor-Klassifikation v3\n")
cat("Datum:", as.character(Sys.time()), "\n")
cat("Raster:", raster_path, "\n")
cat("Training:", shape_path, "\n")
cat("Output:", output_dir, "\n")
cat("Plots:", plot_dir, "\n")
cat("Feature set: S2 quarterly NDVI, NDWI, NDMI, BSI plus summary features\n")
cat("Spatial CV folds:", k_folds, "\n")
cat("Spatial CV block size:", block_size, "m\n")
cat("Seed:", seed_value, "\n")
cat("F1-optimaler Threshold:", best_threshold, "\n")
cat("Spatial-CV AUC:", cv_auc, "\n")
cat("AOA aggregation factor:", aoa_aggregation_factor, "\n")
cat("AOA LPD:", aoa_use_lpd, "\n")
cat("\nSession info:\n")
print(sessionInfo())
sink()

print("Fertig: v3-Modell, F1-Threshold, Metriken, Karten, AOA und farbige Kontrollplots wurden gespeichert.")
