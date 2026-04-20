# ============================================================
# Script 04: Spearman Rank Correlation Analysis
# Candrea et al. - Gut Microbiota Comparative Analysis
# Biomedicines 2025
# ============================================================
# Description:
# Spearman rank correlation analysis between microbial taxa
# relative abundance and clinical/functional parameters.
# Multiple testing correction using Benjamini-Hochberg FDR.
# Correlation strength interpreted per Schober et al. (2018):
#   |r| < 0.10 negligible, 0.10-0.39 weak, 0.40-0.69 moderate,
#   0.70-0.89 strong, >= 0.90 very strong.
#
# Reference: Schober P, et al. Anesth Analg. 2018;126(5):1763-1768.
#            DOI: 10.1213/ANE.0000000000002864
#
# Input:
#   - Adriana Rusu_Baza de date pacienti.xlsx
#
# Output:
#   - spearman_results.csv
#   - Correlogram.png (Supplementary Figure S1)
#   - Supplementary_Figure_S2.png (scatter plots A/B/C)
# ============================================================

library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(corrplot)

# ============================================================
# LOAD AND PREPARE DATA
# ============================================================

read_sheet <- function(sheet_name) {
  df <- read_excel("Adriana Rusu_Baza de date pacienti.xlsx",
                   sheet = sheet_name, col_names = FALSE)
  headers    <- as.character(df[2, ])
  headers[1] <- "Subject"
  headers[2] <- "Group"
  data       <- df[3:nrow(df), ]
  colnames(data) <- headers
  return(data)
}

hc  <- read_sheet("Healthy control")
ibs <- read_sheet("Irritable bowl syndrome")
ai  <- read_sheet("Autoimmunity")
ag  <- read_sheet("Anxiety")

combined        <- bind_rows(hc, ibs, ai, ag)
combined$Group  <- factor(combined$Group, levels = c("HC", "IBS", "AG", "AI"))

# Variables of interest
microbiome_vars <- c("Actinobacteria**", "Bifidobacterium spp.**",
                     "Bacteroides spp.**", "Firmicutes**",
                     "Bacteroidetes**", "Prevotella spp.**")

clinical_vars   <- c("Lactate production", "Butyrate production (%)",
                     "Acetate / propionate", "BMI (kg/m2)",
                     "Vitamin D (ng/mL)", "Total cholesterol (mg/dL)",
                     "LDL cholesterol (mg/dL)", "HDL cholesterol (mg/dL)")

all_vars <- c(microbiome_vars, clinical_vars)
combined[all_vars] <- lapply(combined[all_vars], as.numeric)

# Rename for clarity
rename_map <- c(
  "Actinobacteria**"          = "Actinobacteria",
  "Bifidobacterium spp.**"    = "Bifidobacterium",
  "Bacteroides spp.**"        = "Bacteroides",
  "Firmicutes**"              = "Firmicutes",
  "Bacteroidetes**"           = "Bacteroidetes",
  "Prevotella spp.**"         = "Prevotella",
  "Lactate production"        = "Lactate",
  "Butyrate production (%)"   = "Butyrate",
  "Acetate / propionate"      = "Acetate.Propionate",
  "BMI (kg/m2)"               = "BMI",
  "Vitamin D (ng/mL)"         = "Vitamin.D",
  "Total cholesterol (mg/dL)" = "Total.Cholesterol",
  "LDL cholesterol (mg/dL)"   = "LDL",
  "HDL cholesterol (mg/dL)"   = "HDL"
)

for (old_name in names(rename_map)) {
  if (old_name %in% colnames(combined)) {
    colnames(combined)[colnames(combined) == old_name] <- rename_map[old_name]
  }
}

micro_clean <- c("Actinobacteria", "Bifidobacterium", "Bacteroides",
                 "Firmicutes", "Bacteroidetes", "Prevotella")
clin_clean  <- c("Lactate", "Butyrate", "Acetate.Propionate", "BMI",
                 "Vitamin.D", "Total.Cholesterol", "LDL", "HDL")

cat("Data prepared. n =", nrow(combined), "\n")

# ============================================================
# SPEARMAN CORRELATIONS + BH-FDR CORRECTION
# ============================================================

results <- data.frame()

for (micro in micro_clean) {
  for (clin in clin_clean) {
    x    <- combined[[micro]]
    y    <- combined[[clin]]
    mask <- !is.na(x) & !is.na(y)
    if (sum(mask) >= 5) {
      test <- cor.test(x[mask], y[mask], method = "spearman", exact = FALSE)
      results <- rbind(results, data.frame(
        Microbiome = micro,
        Clinical   = clin,
        n          = sum(mask),
        Spearman_r = round(test$estimate, 3),
        p_value    = round(test$p.value, 4)
      ))
    }
  }
}

# BH-FDR correction
results$p_adj_FDR <- round(p.adjust(results$p_value, method = "BH"), 4)

# Interpretation per Schober et al. 2018
results$Strength <- cut(abs(results$Spearman_r),
                        breaks = c(0, 0.10, 0.40, 0.70, 0.90, 1.00),
                        labels = c("Negligible", "Weak", "Moderate", "Strong", "Very strong"),
                        include.lowest = TRUE)

write.csv(results, "spearman_results.csv", row.names = FALSE)

cat("\nSignificant correlations (FDR < 0.05):\n")
print(results[results$p_adj_FDR < 0.05,
              c("Microbiome", "Clinical", "n", "Spearman_r", "p_value", "p_adj_FDR", "Strength")])

# ============================================================
# CORRELOGRAM (Supplementary Figure S1)
# ============================================================

cor_data   <- combined[, c(micro_clean, clin_clean)]
cor_matrix <- cor(cor_data, method = "spearman", use = "pairwise.complete.obs")

p_matrix <- matrix(NA, nrow = ncol(cor_data), ncol = ncol(cor_data))
rownames(p_matrix) <- colnames(cor_data)
colnames(p_matrix) <- colnames(cor_data)

for (i in 1:ncol(cor_data)) {
  for (j in 1:ncol(cor_data)) {
    test <- cor.test(cor_data[[i]], cor_data[[j]], method = "spearman", exact = FALSE)
    p_matrix[i, j] <- test$p.value
  }
}

# Clean names for plot
colnames(cor_matrix) <- gsub("Acetate.Propionate", "Acetate/Propionate", colnames(cor_matrix))
rownames(cor_matrix) <- gsub("Acetate.Propionate", "Acetate/Propionate", rownames(cor_matrix))
colnames(cor_matrix) <- gsub("Total.Cholesterol",  "Total Cholesterol",  colnames(cor_matrix))
rownames(cor_matrix) <- gsub("Total.Cholesterol",  "Total Cholesterol",  rownames(cor_matrix))
colnames(cor_matrix) <- gsub("Vitamin.D",          "Vitamin D",          colnames(cor_matrix))
rownames(cor_matrix) <- gsub("Vitamin.D",          "Vitamin D",          rownames(cor_matrix))
colnames(p_matrix)   <- colnames(cor_matrix)
rownames(p_matrix)   <- rownames(cor_matrix)

png("Correlogram.png", width = 2800, height = 2800, res = 300)
corrplot(cor_matrix,
         method       = "color",
         type         = "full",
         order        = "hclust",
         p.mat        = p_matrix,
         sig.level    = 0.05,
         insig        = "blank",
         addCoef.col  = "black",
         number.cex   = 0.8,
         tl.col       = "black",
         tl.srt       = 45,
         tl.cex       = 1.0,
         col          = colorRampPalette(c("#E63946", "white", "#457B9D"))(200),
         mar          = c(0, 0, 1, 0))
dev.off()
cat("Supplementary Figure S1 saved: Correlogram.png\n")

# ============================================================
# SCATTER PLOTS (Supplementary Figure S2)
# ============================================================

group_colors <- c("HC" = "#E76F51", "IBS" = "#2A9D8F", "AG" = "#9B5DE5", "AI" = "#00B4D8")

p_a <- ggplot(combined, aes(x = Actinobacteria, y = Lactate, color = Group)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "gray40", linetype = "dashed") +
  scale_color_manual(values = group_colors) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top", size = 3) +
  theme_bw() +
  theme(legend.title = element_blank(), panel.grid.minor = element_blank()) +
  labs(x = "Actinobacteria (%)", y = "Lactate production (%)")

p_b <- ggplot(combined, aes(x = Bifidobacterium, y = Lactate, color = Group)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "gray40", linetype = "dashed") +
  scale_color_manual(values = group_colors) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top", size = 3) +
  theme_bw() +
  theme(legend.title = element_blank(), panel.grid.minor = element_blank()) +
  labs(x = "Bifidobacterium spp. (%)", y = "Lactate production (%)")

p_c <- ggplot(combined, aes(x = Bacteroides, y = Acetate.Propionate, color = Group)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "gray40", linetype = "dashed") +
  scale_color_manual(values = group_colors) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top", size = 3) +
  theme_bw() +
  theme(legend.title = element_blank(), panel.grid.minor = element_blank()) +
  labs(x = "Bacteroides spp. (%)", y = "Acetate/Propionate ratio")

combined_plot <- ggarrange(p_a, p_b, p_c,
                           ncol = 3, nrow = 1,
                           labels = c("A", "B", "C"),
                           common.legend = TRUE,
                           legend = "right")

ggsave("Supplementary_Figure_S2.png", combined_plot, width = 14, height = 5, dpi = 300)
cat("Supplementary Figure S2 saved: Supplementary_Figure_S2.png\n")
cat("\nSpearman analysis complete.\n")


ggsave("Supplementary_Figure_S2.png", combined_plot, width = 14, height = 5, dpi = 300)
cat("Supplementary Figure S2 saved: Supplementary_Figure_S2.png\n")
cat("\nSpearman analysis complete.\n")
