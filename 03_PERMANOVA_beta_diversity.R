# ============================================================
# Script 03: PERMANOVA Beta Diversity Analysis
# Candrea et al. - Gut Microbiota Comparative Analysis
# Biomedicines 2025
# ============================================================
# Description:
# Performs PERMANOVA (permutational multivariate analysis of variance)
# on genus-level relative abundance data using three distance metrics:
# Bray-Curtis, Jaccard, and Jensen-Shannon.
# Pairwise comparisons with Benjamini-Hochberg FDR correction.
# PCoA variance explained calculated via eigenvalue decomposition.
#
# Input:
#   - linda_otu_genus.csv  (from 01_data_preparation.R)
#   - linda_meta_phylum.csv (metadata)
#
# Output:
#   - permanova_results.csv
#   - pcoa_variance_explained.csv
#
# Software: vegan package
# ============================================================

library(vegan)
library(philentropy)

# ============================================================
# LOAD DATA
# ============================================================

otu  <- read.csv("linda_otu_genus.csv", row.names = 1, check.names = FALSE)
meta <- read.csv("linda_meta_phylum.csv", row.names = 1)

# Transpose: samples on rows
otu_t <- t(otu)

# Normalize to proportions (row sums to 1)
otu_prop <- otu_t / rowSums(otu_t)

groups <- meta$Group
cat("Samples:", nrow(otu_prop), "| Groups:", table(groups), "\n")

# ============================================================
# DISTANCE MATRICES
# ============================================================

cat("Calculating distance matrices...\n")

# Bray-Curtis
bc_dist  <- vegdist(otu_prop, method = "bray")

# Jaccard
jac_dist <- vegdist(otu_prop, method = "jaccard")

# Jensen-Shannon
js_dist  <- distance(otu_prop, method = "jensen-shannon")
js_dist_matrix <- as.dist(js_dist)

# ============================================================
# GLOBAL PERMANOVA (999 permutations)
# ============================================================

cat("\nRunning global PERMANOVA (999 permutations)...\n")

set.seed(42)
perm_bc  <- adonis2(bc_dist  ~ groups, permutations = 999)
perm_jac <- adonis2(jac_dist ~ groups, permutations = 999)
perm_js  <- adonis2(js_dist_matrix ~ groups, permutations = 999)

cat("\nBray-Curtis PERMANOVA:\n");  print(perm_bc)
cat("\nJaccard PERMANOVA:\n");      print(perm_jac)
cat("\nJensen-Shannon PERMANOVA:\n"); print(perm_js)

# ============================================================
# PCoA VARIANCE EXPLAINED (eigenvalue decomposition)
# ============================================================

pcoa_variance <- function(dist_matrix, name) {
  
  pcoa_res <- cmdscale(dist_matrix, k = 2, eig = TRUE)
  pos_eig  <- pcoa_res$eig[pcoa_res$eig > 0]
  pct1 <- pos_eig[1] / sum(pos_eig) * 100
  pct2 <- pos_eig[2] / sum(pos_eig) * 100
  
  cat(sprintf("%s: PCoA1 = %.1f%%, PCoA2 = %.1f%%\n", name, pct1, pct2))
  return(c(PCoA1 = round(pct1, 1), PCoA2 = round(pct2, 1)))
}

cat("\nPCoA Variance Explained:\n")
var_bc  <- pcoa_variance(bc_dist,         "Bray-Curtis")
var_jac <- pcoa_variance(jac_dist,        "Jaccard")
var_js  <- pcoa_variance(js_dist_matrix,  "Jensen-Shannon")

# Save variance explained
var_df <- data.frame(
  Distance = c("Bray-Curtis", "Jaccard", "Jensen-Shannon"),
  PCoA1_pct = c(var_bc[1], var_jac[1], var_js[1]),
  PCoA2_pct = c(var_bc[2], var_jac[2], var_js[2])
)
write.csv(var_df, "pcoa_variance_explained.csv", row.names = FALSE)

# ============================================================
# PAIRWISE PERMANOVA WITH BH-FDR CORRECTION
# ============================================================

cat("\nRunning pairwise PERMANOVA (Bray-Curtis)...\n")

groups_unique <- unique(groups)
pairs <- combn(groups_unique, 2, simplify = FALSE)

pairwise_results <- data.frame()

for (pair in pairs) {
  g1 <- pair[1]; g2 <- pair[2]
  mask <- groups %in% c(g1, g2)
  sub_dist <- as.dist(as.matrix(bc_dist)[mask, mask])
  sub_groups <- groups[mask]
  
  set.seed(42)
  res <- adonis2(sub_dist ~ sub_groups, permutations = 999)
  
  pairwise_results <- rbind(pairwise_results, data.frame(
    Comparison = paste0(g1, "_vs_", g2),
    Pseudo_F   = round(res$F[1], 3),
    R2         = round(res$R2[1], 3),
    p_value    = res$`Pr(>F)`[1]
  ))
}

# FDR correction
pairwise_results$p_adj_FDR <- round(p.adjust(pairwise_results$p_value, method = "BH"), 4)
pairwise_results <- pairwise_results[order(pairwise_results$p_value), ]

cat("\nPairwise PERMANOVA results (Bray-Curtis, BH-FDR corrected):\n")
print(pairwise_results)

# Save all results
permanova_summary <- data.frame(
  Distance  = c("Bray-Curtis", "Jaccard", "Jensen-Shannon"),
  Pseudo_F  = c(round(perm_bc$F[1], 3), round(perm_jac$F[1], 3), round(perm_js$F[1], 3)),
  R2        = c(round(perm_bc$R2[1], 3), round(perm_jac$R2[1], 3), round(perm_js$R2[1], 3)),
  p_value   = c(perm_bc$`Pr(>F)`[1], perm_jac$`Pr(>F)`[1], perm_js$`Pr(>F)`[1])
)

write.csv(permanova_summary,  "permanova_global_results.csv",   row.names = FALSE)
write.csv(pairwise_results,   "permanova_pairwise_results.csv", row.names = FALSE)

cat("\nPERMANOVA analysis complete.\n")
cat("Results saved to: permanova_global_results.csv, permanova_pairwise_results.csv\n")
