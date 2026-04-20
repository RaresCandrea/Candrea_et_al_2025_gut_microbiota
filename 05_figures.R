# ============================================================
# Script 05: Figure Generation
# Candrea et al. - Gut Microbiota Comparative Analysis
# Biomedicines 2025
# ============================================================
# Description:
# Generates all manuscript figures (Figures 1-3 and 6).
# Note: Figure 4 (alpha diversity) and Figure 5 (PCoA) were
# generated using MicrobiomeAnalyst v2.0 (Lu et al., 2023).
#
# Input:
#   - De analizat.xlsx (sheets: Phylum, Genus)
#   - Adriana Rusu_Baza de date pacienti.xlsx
#
# Output:
#   - Figure_1A.png/pdf  (Phylum-level stacked bar plot)
#   - Figure_1B.png/pdf  (Actinobacteria boxplot)
#   - Figure_2.png/pdf   (Genus stacked bar plot top 10)
#   - Figure_3A.png/pdf  (Bacteroides spp. boxplot)
#   - Figure_3B.png/pdf  (Bifidobacterium spp. boxplot)
#   - Figure_6.png/pdf   (Spearman scatter plots A/B/C)
# ============================================================

library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# Group colors and order consistent across all figures
group_colors <- c("HC" = "#E76F51", "IBS" = "#2A9D8F", "AG" = "#9B5DE5", "AI" = "#00B4D8")
group_levels <- c("HC", "AG", "AI", "IBS")

# ============================================================
# LOAD AND CLEAN DATA
# ============================================================

clean_taxa_data <- function(sheet_name) {
  df        <- read_excel("De analizat.xlsx", sheet = sheet_name)
  taxa_cols <- colnames(df)[!colnames(df) %in% c("Subject", "Group")]
  df[taxa_cols] <- lapply(df[taxa_cols], as.numeric)
  clean_names   <- trimws(gsub(" \\(%\\)|\\(%\\)", "", taxa_cols))
  colnames(df)[colnames(df) %in% taxa_cols] <- clean_names
  df$Group <- factor(df$Group, levels = group_levels)
  return(df)
}

phylum_data <- clean_taxa_data("Phylum")
genus_data  <- clean_taxa_data("Genus")

phylum_taxa <- colnames(phylum_data)[!colnames(phylum_data) %in% c("Subject", "Group")]
genus_taxa  <- colnames(genus_data)[!colnames(genus_data)  %in% c("Subject", "Group")]

cat("Data loaded. Phylum taxa:", length(phylum_taxa),
    "| Genus taxa:", length(genus_taxa), "\n")

# ============================================================
# FIGURE 1A — Phylum stacked bar plot
# ============================================================

# Normalize per sample to 100%
phylum_norm <- phylum_data
phylum_norm[, phylum_taxa] <- phylum_data[, phylum_taxa] /
  rowSums(phylum_data[, phylum_taxa], na.rm = TRUE) * 100

# Main phyla — consistent with Table 2
main_phyla <- c("Firmicutes", "Bacteroidetes", "Proteobacteria",
                "Actinobacteria", "Verrucomicrobia")

mean_phylum <- phylum_norm %>%
  group_by(Group) %>%
  summarise(across(all_of(main_phyla), mean, na.rm = TRUE), .groups = "drop")

mean_phylum$Other <- 100 - rowSums(mean_phylum[, main_phyla])

mean_phylum_long <- mean_phylum %>%
  pivot_longer(cols = c(all_of(main_phyla), "Other"),
               names_to = "Taxon", values_to = "Abundance")

taxa_order_phylum <- c("Other", "Verrucomicrobia", "Actinobacteria",
                       "Proteobacteria", "Bacteroidetes", "Firmicutes")
mean_phylum_long$Taxon <- factor(mean_phylum_long$Taxon, levels = taxa_order_phylum)

n_phylum <- phylum_data %>% count(Group)
group_labels_phylum <- setNames(
  paste0(n_phylum$Group, "\n(n=", n_phylum$n, ")"), n_phylum$Group)

phylum_colors <- c(
  "Other"           = "#CCCCCC",
  "Verrucomicrobia" = "#F9C74F",
  "Actinobacteria"  = "#9B5DE5",
  "Proteobacteria"  = "#F77F00",
  "Bacteroidetes"   = "#457B9D",
  "Firmicutes"      = "#E63946"
)

p1a <- ggplot(mean_phylum_long, aes(x = Group, y = Abundance, fill = Taxon)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  scale_fill_manual(values = phylum_colors) +
  scale_x_discrete(labels = group_labels_phylum) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 101)) +
  theme_bw() +
  theme(axis.text.x    = element_text(size = 10),
        axis.text.y    = element_text(size = 10),
        legend.text    = element_text(size = 9, face = "italic"),
        legend.key.size = unit(0.4, "cm"),
        legend.title   = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Study group", y = "Relative abundance (%)")

ggsave("Figure_1A.pdf", p1a, width = 7, height = 5, dpi = 300)
ggsave("Figure_1A.png", p1a, width = 7, height = 5, dpi = 300)
cat("Figure 1A saved.\n")

# ============================================================
# FIGURE 1B — Actinobacteria boxplot
# ============================================================

actin_data       <- phylum_data[, c("Group", "Actinobacteria")]
colnames(actin_data)[2] <- "Abundance"

p1b <- ggplot(actin_data, aes(x = Group, y = Abundance, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = group_colors) +
  theme_bw() +
  theme(legend.position = "none", panel.grid.minor = element_blank()) +
  labs(x = NULL, y = "Relative Abundance (%)")

ggsave("Figure_1B.pdf", p1b, width = 5, height = 5, dpi = 300)
ggsave("Figure_1B.png", p1b, width = 5, height = 5, dpi = 300)
cat("Figure 1B saved.\n")

# ============================================================
# FIGURE 2 — Genus stacked bar plot (top 10 + Other)
# ============================================================

genus_norm <- genus_data
genus_norm[, genus_taxa] <- genus_data[, genus_taxa] /
  rowSums(genus_data[, genus_taxa], na.rm = TRUE) * 100

global_means_genus <- colMeans(genus_norm[, genus_taxa], na.rm = TRUE)
top10 <- names(sort(global_means_genus, decreasing = TRUE))[1:10]

mean_genus <- genus_norm %>%
  group_by(Group) %>%
  summarise(across(all_of(genus_taxa), mean, na.rm = TRUE), .groups = "drop")

mean_genus$Other <- 100 - rowSums(mean_genus[, top10])

mean_genus_long <- mean_genus %>%
  pivot_longer(cols = c(all_of(top10), "Other"),
               names_to = "Taxon", values_to = "Abundance")

mean_genus_long$Taxon <- factor(mean_genus_long$Taxon,
                                levels = c("Other", rev(top10)))
mean_genus_long$Group <- factor(mean_genus_long$Group, levels = group_levels)

n_genus <- genus_data %>% count(Group)
group_labels_genus <- setNames(
  paste0(n_genus$Group, "\n(n=", n_genus$n, ")"), n_genus$Group)

taxon_colors_genus <- c(
  "Other" = "#CCCCCC",
  setNames(c("#E63946","#457B9D","#2A9D8F","#E9C46A","#F4A261",
             "#264653","#9B5DE5","#00B4D8","#90BE6D","#F77F00"), top10))

p2 <- ggplot(mean_genus_long, aes(x = Group, y = Abundance, fill = Taxon)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  scale_fill_manual(values = taxon_colors_genus) +
  scale_x_discrete(labels = group_labels_genus) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 101)) +
  theme_bw() +
  theme(axis.text.x    = element_text(size = 10),
        legend.text    = element_text(size = 8, face = "italic"),
        legend.key.size = unit(0.4, "cm"),
        panel.grid.minor = element_blank()) +
  labs(x = "Study group", y = "Relative abundance (%)", fill = "Genus")

ggsave("Figure_2.pdf", p2, width = 9, height = 6, dpi = 300)
ggsave("Figure_2.png", p2, width = 9, height = 6, dpi = 300)
cat("Figure 2 saved.\n")

# ============================================================
# FIGURE 3A — Bacteroides spp. boxplot
# ============================================================

bact_data       <- genus_data[, c("Group", "Bacteroides spp.")]
colnames(bact_data)[2] <- "Abundance"

p3a <- ggplot(bact_data, aes(x = Group, y = Abundance, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = group_colors) +
  theme_bw() +
  theme(legend.position = "none", panel.grid.minor = element_blank()) +
  labs(x = NULL, y = "Relative Abundance (%)")

ggsave("Figure_3A.pdf", p3a, width = 5, height = 5, dpi = 300)
ggsave("Figure_3A.png", p3a, width = 5, height = 5, dpi = 300)
cat("Figure 3A saved.\n")

# ============================================================
# FIGURE 3B — Bifidobacterium spp. boxplot
# ============================================================

bifi_data       <- genus_data[, c("Group", "Bifidobacterium spp.")]
colnames(bifi_data)[2] <- "Abundance"

p3b <- ggplot(bifi_data, aes(x = Group, y = Abundance, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = group_colors) +
  theme_bw() +
  theme(legend.position = "none", panel.grid.minor = element_blank()) +
  labs(x = NULL, y = "Relative Abundance (%)")

ggsave("Figure_3B.pdf", p3b, width = 5, height = 5, dpi = 300)
ggsave("Figure_3B.png", p3b, width = 5, height = 5, dpi = 300)
cat("Figure 3B saved.\n")

# ============================================================
# FIGURE 6 — Spearman scatter plots A/B/C
# ============================================================

read_sheet_clinical <- function(sheet_name) {
  df         <- read_excel("Adriana Rusu_Baza de date pacienti.xlsx",
                            sheet = sheet_name, col_names = FALSE)
  headers    <- as.character(df[2, ])
  headers[1] <- "Subject"
  headers[2] <- "Group"
  data       <- df[3:nrow(df), ]
  colnames(data) <- headers
  return(data)
}

hc  <- read_sheet_clinical("Healthy control")
ibs <- read_sheet_clinical("Irritable bowl syndrome")
ai  <- read_sheet_clinical("Autoimmunity")
ag  <- read_sheet_clinical("Anxiety")

combined       <- bind_rows(hc, ibs, ai, ag)
combined$Group <- factor(combined$Group, levels = group_levels)

vars_needed <- c("Actinobacteria**", "Bifidobacterium spp.**",
                 "Bacteroides spp.**", "Lactate production",
                 "Acetate / propionate")
combined[vars_needed] <- lapply(combined[vars_needed], as.numeric)

colnames(combined)[colnames(combined) == "Actinobacteria**"]       <- "Actinobacteria"
colnames(combined)[colnames(combined) == "Bifidobacterium spp.**"] <- "Bifidobacterium"
colnames(combined)[colnames(combined) == "Bacteroides spp.**"]     <- "Bacteroides"
colnames(combined)[colnames(combined) == "Lactate production"]     <- "Lactate"
colnames(combined)[colnames(combined) == "Acetate / propionate"]   <- "Acetate.Propionate"

p6a <- ggplot(combined, aes(x = Actinobacteria, y = Lactate, color = Group)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "gray40", linetype = "dashed") +
  scale_color_manual(values = group_colors) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top", size = 3) +
  theme_bw() +
  theme(legend.title = element_blank(), panel.grid.minor = element_blank()) +
  labs(x = "Actinobacteria (%)", y = "Lactate production (%)")

p6b <- ggplot(combined, aes(x = Bifidobacterium, y = Lactate, color = Group)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "gray40", linetype = "dashed") +
  scale_color_manual(values = group_colors) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top", size = 3) +
  theme_bw() +
  theme(legend.title = element_blank(), panel.grid.minor = element_blank()) +
  labs(x = "Bifidobacterium spp. (%)", y = "Lactate production (%)")

p6c <- ggplot(combined, aes(x = Bacteroides, y = Acetate.Propionate, color = Group)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "gray40", linetype = "dashed") +
  scale_color_manual(values = group_colors) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top", size = 3) +
  theme_bw() +
  theme(legend.title = element_blank(), panel.grid.minor = element_blank()) +
  labs(x = "Bacteroides spp. (%)", y = "Acetate/Propionate ratio")

figure_6 <- ggarrange(p6a, p6b, p6c,
                      ncol = 3, nrow = 1,
                      labels = c("A", "B", "C"),
                      common.legend = TRUE,
                      legend = "right")

ggsave("Figure_6.pdf", figure_6, width = 14, height = 5, dpi = 300)
ggsave("Figure_6.png", figure_6, width = 14, height = 5, dpi = 300)
cat("Figure 6 saved.\n")

cat("\nAll figures saved successfully.\n")
cat("Note: Figure 4 (alpha diversity) and Figure 5 (PCoA) were\n")
cat("generated using MicrobiomeAnalyst v2.0 (Lu et al., 2023).\n")
