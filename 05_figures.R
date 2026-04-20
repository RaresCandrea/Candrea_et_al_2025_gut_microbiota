# ============================================================
# Script 05: Figure Generation
# Candrea et al. - Gut Microbiota Comparative Analysis
# Biomedicines 2025
# ============================================================
# Description:
# Generates all manuscript figures (Figures 1-6).
#
# Input:
#   - De analizat.xlsx
#   - Adriana Rusu_Baza de date pacienti.xlsx
#
# Output:
#   - Figure_1A.png/pdf  (Phylum dot plot)
#   - Figure_1B.png/pdf  (Actinobacteria boxplot)
#   - Figure_2.png/pdf   (Genus stacked bar plot)
#   - Figure_3A.png/pdf  (Bacteroides boxplot)
#   - Figure_3B.png/pdf  (Bifidobacterium boxplot)
#   - Figure_6.png/pdf   (Spearman scatter plots A/B/C)
# ============================================================

library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(RColorBrewer)

# Group colors consistent across all figures
group_colors <- c("HC" = "#E76F51", "IBS" = "#2A9D8F", "AG" = "#9B5DE5", "AI" = "#00B4D8")
group_levels <- c("HC", "AG", "AI", "IBS")

# ============================================================
# LOAD DATA
# ============================================================

phylum_data  <- read_excel("De analizat.xlsx", sheet = "Phylum")
genus_data   <- read_excel("De analizat.xlsx", sheet = "Genus")

# Clean column names and convert to numeric
clean_data <- function(df) {
  taxa_cols <- colnames(df)[!colnames(df) %in% c("Subject", "Group")]
  df[taxa_cols] <- lapply(df[taxa_cols], as.numeric)
  clean_names <- gsub(" \\(%\\)|\\(%\\)", "", taxa_cols)
  colnames(df)[colnames(df) %in% taxa_cols] <- clean_names
  df$Group <- factor(df$Group, levels = group_levels)
  return(df)
}

phylum_data <- clean_data(phylum_data)
genus_data  <- clean_data(genus_data)

phylum_taxa <- colnames(phylum_data)[!colnames(phylum_data) %in% c("Subject", "Group")]
genus_taxa  <- colnames(genus_data)[!colnames(genus_data)  %in% c("Subject", "Group")]

# ============================================================
# FIGURE 1A 
# ============================================================

library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)

# Incarca datele
df <- read_excel("De analizat.xlsx", sheet = "Phylum")

# Curata coloanele
phylum_taxa <- colnames(df)[!colnames(df) %in% c("Subject", "Group")]
df[phylum_taxa] <- lapply(df[phylum_taxa], as.numeric)

# Curata numele coloanelor
clean_names <- gsub(" \\(%\\)|\\(%\\)", "", phylum_taxa)
clean_names <- trimws(clean_names)
colnames(df)[colnames(df) %in% phylum_taxa] <- clean_names
phylum_taxa <- clean_names

# Normalizeaza la 100% per sample
df[phylum_taxa] <- df[phylum_taxa] / rowSums(df[phylum_taxa], na.rm = TRUE) * 100

# Phyla principale - exact cele din tabel
main_phyla <- c("Firmicutes", "Bacteroidetes", "Proteobacteria", 
                "Actinobacteria", "Verrucomicrobia")

# Calculeaza media per grup
df$Group <- factor(df$Group, levels = c("HC", "AG", "AI", "IBS"))

mean_df <- df %>%
  group_by(Group) %>%
  summarise(across(all_of(main_phyla), mean, na.rm = TRUE))

mean_df$Other <- 100 - rowSums(mean_df[, main_phyla])

mean_long <- mean_df %>%
  pivot_longer(cols = c(all_of(main_phyla), "Other"),
               names_to = "Taxon", values_to = "Abundance")

# Ordine taxa in legenda
taxa_order <- c("Other", "Verrucomicrobia", "Actinobacteria",
                "Proteobacteria", "Bacteroidetes", "Firmicutes")
mean_long$Taxon <- factor(mean_long$Taxon, levels = taxa_order)

# Label grupuri cu n
n_labels <- df %>% count(Group)
group_labels <- setNames(
  paste0(n_labels$Group, "\n(n=", n_labels$n, ")"),
  n_labels$Group
)

# Culori consistente cu Figure 2
phylum_colors <- c(
  "Other"           = "#CCCCCC",
  "Verrucomicrobia" = "#F9C74F",
  "Actinobacteria"  = "#9B5DE5",
  "Proteobacteria"  = "#F77F00",
  "Bacteroidetes"   = "#457B9D",
  "Firmicutes"      = "#E63946"
)

# Figura
p_phylum <- ggplot(mean_long, aes(x = Group, y = Abundance, fill = Taxon)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  scale_fill_manual(values = phylum_colors) +
  scale_x_discrete(labels = group_labels) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 101)) +
  theme_bw() +
  theme(
    axis.text.x  = element_text(size = 10),
    axis.text.y  = element_text(size = 10),
    legend.text  = element_text(size = 9, face = "italic"),
    legend.key.size = unit(0.4, "cm"),
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(x = "Study group", y = "Relative abundance (%)")

ggsave("Figure_1A_phylum_stacked.pdf", p_phylum, width = 7, height = 5, dpi = 300)
ggsave("Figure_1A_phylum_stacked.png", p_phylum, width = 7, height = 5, dpi = 300)
cat("Figure 1A saved!\n")

# ============================================================
# FIGURE 1B — Actinobacteria boxplot
# ============================================================

actin_data <- phylum_data[, c("Group", "Actinobacteria")]
colnames(actin_data)[2] <- "Abundance"
actin_data$Group <- factor(actin_data$Group, levels = group_levels)

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

# Normalize to 100% of genus-level abundance per sample
genus_norm <- genus_data
genus_norm[, genus_taxa] <- genus_data[, genus_taxa] /
  rowSums(genus_data[, genus_taxa], na.rm = TRUE) * 100

# Top 10 genera by global mean
global_means <- colMeans(genus_norm[, genus_taxa], na.rm = TRUE)
top10 <- names(sort(global_means, decreasing = TRUE))[1:10]

mean_df <- genus_norm %>%
  group_by(Group) %>%
  summarise(across(all_of(genus_taxa), mean, na.rm = TRUE))

mean_df$Other <- 100 - rowSums(mean_df[, top10])

all_taxa   <- c(top10, "Other")
mean_long  <- mean_df %>%
  pivot_longer(cols = all_of(all_taxa), names_to = "Taxon", values_to = "Abundance")

top10_ordered <- c("Other", rev(top10))
mean_long$Taxon <- factor(mean_long$Taxon, levels = top10_ordered)
mean_long$Group <- factor(mean_long$Group, levels = group_levels)

# Group labels with n
n_labels    <- genus_data %>% count(Group)
group_labels <- setNames(paste0(n_labels$Group, "\n(n=", n_labels$n, ")"), n_labels$Group)

taxon_colors <- c("Other" = "#CCCCCC",
                  setNames(c("#E63946","#457B9D","#2A9D8F","#E9C46A","#F4A261",
                             "#264653","#9B5DE5","#00B4D8","#90BE6D","#F77F00"), top10))

p2 <- ggplot(mean_long, aes(x = Group, y = Abundance, fill = Taxon)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  scale_fill_manual(values = taxon_colors) +
  scale_x_discrete(labels = group_labels) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 101)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        legend.text = element_text(size = 8, face = "italic"),
        legend.key.size = unit(0.4, "cm"),
        panel.grid.minor = element_blank()) +
  labs(x = "Study group", y = "Relative abundance (%)", fill = "Genus")

ggsave("Figure_2.pdf", p2, width = 9, height = 6, dpi = 300)
ggsave("Figure_2.png", p2, width = 9, height = 6, dpi = 300)
cat("Figure 2 saved.\n")

# ============================================================
# FIGURE 3A — Bacteroides boxplot
# ============================================================

bact_data       <- genus_data[, c("Group", "Bacteroides spp.")]
colnames(bact_data)[2] <- "Abundance"
bact_data$Group <- factor(bact_data$Group, levels = group_levels)

p3a <- ggplot(bact_data, aes(x = Group, y = Abundance, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(limits = c(0, 60)) +
  theme_bw() +
  theme(legend.position = "none", panel.grid.minor = element_blank()) +
  labs(x = NULL, y = "Relative Abundance (%)")

ggsave("Figure_3A.pdf", p3a, width = 5, height = 5, dpi = 300)
ggsave("Figure_3A.png", p3a, width = 5, height = 5, dpi = 300)
cat("Figure 3A saved.\n")

# ============================================================
# FIGURE 3B — Bifidobacterium boxplot
# ============================================================

bifi_data       <- genus_data[, c("Group", "Bifidobacterium spp.")]
colnames(bifi_data)[2] <- "Abundance"
bifi_data$Group <- factor(bifi_data$Group, levels = group_levels)

p3b <- ggplot(bifi_data, aes(x = Group, y = Abundance, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(limits = c(0, 25)) +
  theme_bw() +
  theme(legend.position = "none", panel.grid.minor = element_blank()) +
  labs(x = NULL, y = "Relative Abundance (%)")

ggsave("Figure_3B.pdf", p3b, width = 5, height = 5, dpi = 300)
ggsave("Figure_3B.png", p3b, width = 5, height = 5, dpi = 300)
cat("Figure 3B saved.\n")

# ============================================================
# FIGURE 6 — Spearman scatter plots (run after 04_Spearman)
# ============================================================
# Note: Requires 'combined' object from 04_Spearman_correlations.R
# Source that script first, or load data here independently.

cat("\nAll figures saved successfully.\n")
cat("Note: Figure 4 (alpha diversity) and Figure 5 (PCoA) were\n")
cat("generated using MicrobiomeAnalyst v2.0 (Lu et al., 2023).\n")
