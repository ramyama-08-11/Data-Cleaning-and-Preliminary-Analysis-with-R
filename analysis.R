# Week 1: Data cleaning and preliminary analysis
# Run from the project root with: Rscript analysis.R

options(stringsAsFactors = FALSE, scipen = 999)

input_file <- file.path("data", "penguins.csv")
output_dir <- "outputs"
if (!file.exists(input_file)) stop("Dataset not found. Run this script from the project root.")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

penguins_raw <- read.csv(input_file, na.strings = c("NA", ""), check.names = FALSE)

cat("=== SOURCE STRUCTURE ===\n")
cat("Rows:", nrow(penguins_raw), " Columns:", ncol(penguins_raw), "\n")
str(penguins_raw)
cat("\n=== SOURCE SUMMARY ===\n")
print(summary(penguins_raw))

numeric_traits <- c("bill_length_mm", "bill_depth_mm", "flipper_length_mm", "body_mass_g")
category_fields <- c("species", "island", "sex")

missing_counts <- data.frame(
  variable = names(penguins_raw),
  missing_n = vapply(penguins_raw, function(x) sum(is.na(x)), integer(1)),
  missing_percent = round(vapply(penguins_raw, function(x) mean(is.na(x)) * 100, numeric(1)), 2)
)
write.csv(missing_counts, file.path(output_dir, "missingness.csv"), row.names = FALSE)
cat("\n=== MISSINGNESS ===\n")
print(missing_counts)

# IQR screening is performed on observed values. Flags are retained; no valid
# biological observations are deleted solely because they fall beyond a fence.
outlier_summary <- do.call(rbind, lapply(numeric_traits, function(field) {
  values <- penguins_raw[[field]]
  q1 <- quantile(values, 0.25, na.rm = TRUE, names = FALSE)
  q3 <- quantile(values, 0.75, na.rm = TRUE, names = FALSE)
  spread <- q3 - q1
  lower <- q1 - 1.5 * spread
  upper <- q3 + 1.5 * spread
  flags <- !is.na(values) & (values < lower | values > upper)
  data.frame(variable = field, lower_fence = lower, upper_fence = upper,
             flagged_n = sum(flags), flagged_rows = paste(which(flags), collapse = ";"))
}))
write.csv(outlier_summary, file.path(output_dir, "outlier_screen.csv"), row.names = FALSE)
cat("\n=== IQR OUTLIER SCREEN (REVIEW, NOT AUTO-DELETE) ===\n")
print(outlier_summary[, c("variable", "lower_fence", "upper_fence", "flagged_n")])

penguins_clean <- penguins_raw
for (field in numeric_traits) {
  missing_flag <- is.na(penguins_clean[[field]])
  median_value <- median(penguins_clean[[field]], na.rm = TRUE)
  penguins_clean[[paste0(field, "_was_missing")]] <- missing_flag
  penguins_clean[[field]][missing_flag] <- median_value
}
for (field in category_fields) {
  penguins_clean[[field]] <- as.character(penguins_clean[[field]])
  penguins_clean[[field]][is.na(penguins_clean[[field]])] <- "Unknown"
  penguins_clean[[field]] <- factor(penguins_clean[[field]])
}
penguins_clean$year <- as.integer(penguins_clean$year)
write.csv(penguins_clean, file.path(output_dir, "penguins_clean.csv"), row.names = FALSE)

# Min-max scaled fields are additional features; original units remain intact.
penguins_scaled <- penguins_clean
for (field in numeric_traits) {
  bounds <- range(penguins_clean[[field]])
  penguins_scaled[[paste0(field, "_minmax")]] <-
    (penguins_clean[[field]] - bounds[1]) / (bounds[2] - bounds[1])
}
write.csv(penguins_scaled, file.path(output_dir, "penguins_scaled.csv"), row.names = FALSE)

# Explicit full one-hot encoding retains a column for every level, including
# the Unknown category. Drop a reference column later if a model requires it.
encoded <- do.call(cbind, lapply(category_fields, function(field) {
  factor_value <- penguins_clean[[field]]
  indicators <- sapply(levels(factor_value), function(level) as.integer(factor_value == level))
  colnames(indicators) <- paste(field, levels(factor_value), sep = "_")
  indicators
}))
encoded <- as.data.frame(encoded, check.names = FALSE)
model_data <- cbind(penguins_scaled[, c(numeric_traits, "year")], as.data.frame(encoded))
write.csv(model_data, file.path(output_dir, "model_matrix.csv"), row.names = FALSE)
cat("\n=== ENCODING ===\n")
cat("One-hot predictor columns:", ncol(encoded), "\n")
cat("Encoded fields:", paste(colnames(encoded), collapse = ", "), "\n")

# Descriptive statistics, including SD and observed sample size.
descriptive <- do.call(rbind, lapply(numeric_traits, function(field) {
  x <- penguins_raw[[field]]
  data.frame(variable = field, n = sum(!is.na(x)), missing_n = sum(is.na(x)),
             mean = mean(x, na.rm = TRUE), median = median(x, na.rm = TRUE),
             sd = sd(x, na.rm = TRUE), min = min(x, na.rm = TRUE),
             q1 = quantile(x, 0.25, na.rm = TRUE, names = FALSE),
             q3 = quantile(x, 0.75, na.rm = TRUE, names = FALSE),
             max = max(x, na.rm = TRUE))
}))
write.csv(descriptive, file.path(output_dir, "descriptive_statistics.csv"), row.names = FALSE)
cat("\n=== DESCRIPTIVE STATISTICS (RAW OBSERVED VALUES) ===\n")
print(descriptive, row.names = FALSE, digits = 4)

species_counts <- as.data.frame(table(penguins_raw$species), stringsAsFactors = FALSE)
names(species_counts) <- c("species", "n")
write.csv(species_counts, file.path(output_dir, "species_counts.csv"), row.names = FALSE)

species_means <- aggregate(penguins_raw[numeric_traits],
                           by = list(species = penguins_raw$species),
                           FUN = function(x) mean(x, na.rm = TRUE))
write.csv(species_means, file.path(output_dir, "species_means.csv"), row.names = FALSE)
cat("\n=== SPECIES COUNTS ===\n")
print(species_counts, row.names = FALSE)
cat("\n=== SPECIES MEANS ===\n")
print(species_means, row.names = FALSE, digits = 4)

correlation <- cor(penguins_raw[numeric_traits], use = "pairwise.complete.obs", method = "pearson")
write.csv(round(correlation, 3), file.path(output_dir, "correlations.csv"))
cat("\n=== PEARSON CORRELATIONS (PAIRWISE COMPLETE) ===\n")
print(round(correlation, 3))

complete_n <- sum(complete.cases(penguins_raw))
cat("\nComplete cases:", complete_n, "of", nrow(penguins_raw), "\n")
cat("Rows retained after imputation:", nrow(penguins_clean), "\n")
cat("Remaining missing cells in cleaned data:", sum(is.na(penguins_clean)), "\n")

# Plot 1: distribution and species comparison.
png(file.path(output_dir, "figure_1_distributions.png"), width = 1500, height = 900, res = 150)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), cex = 0.9)
plot_colors <- c("#287271", "#D98C35", "#6A8EAE")
for (i in seq_along(numeric_traits)) {
  field <- numeric_traits[i]
  boxplot(penguins_raw[[field]] ~ penguins_raw$species, col = plot_colors,
          xlab = "Species", ylab = field, main = gsub("_", " ", field))
}
dev.off()

# Plot 2: morphology relationship with a linear trend guide.
png(file.path(output_dir, "figure_2_bill_flipper.png"), width = 1200, height = 850, res = 150)
par(mar = c(5, 5, 3, 1))
species_levels <- c("Adelie", "Chinstrap", "Gentoo")
point_colors <- plot_colors[match(penguins_raw$species, species_levels)]
plot(penguins_raw$bill_length_mm, penguins_raw$flipper_length_mm,
     pch = 19, col = adjustcolor(point_colors, alpha.f = 0.72),
     xlab = "Bill length (mm)", ylab = "Flipper length (mm)",
     main = "Bill length and flipper length")
legend("topleft", legend = species_levels, col = plot_colors, pch = 19, bty = "n")
abline(lm(flipper_length_mm ~ bill_length_mm, data = penguins_raw), col = "#263238", lwd = 2)
dev.off()

# Plot 3: missingness by variable.
png(file.path(output_dir, "figure_3_missingness.png"), width = 1200, height = 700, res = 150)
par(mar = c(5, 10, 3, 1))
barplot(missing_counts$missing_n, names.arg = missing_counts$variable, horiz = TRUE,
        las = 1, col = "#D98C35", xlab = "Missing records", main = "Missing values by field")
dev.off()

cat("\nFigures written to outputs/figure_*.png\n")
