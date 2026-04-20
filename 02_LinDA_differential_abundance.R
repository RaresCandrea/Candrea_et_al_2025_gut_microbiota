# Script 02: LinDA Alternative Differential Abundance Analysis
# Candrea et al. - Gut Microbiota Comparative Analysis
# Biomedicines 2025
#
# Run 01_data_preparation.R first to generate input files.
#
# ============================================================
# LinDA Analysis - Phylum, Genus, Species Level
# Candrea et al. - Gut Microbiota Comparative Analysis
# ============================================================

# Install package if not already installed
if (!requireNamespace("MicrobiomeStat", quietly = TRUE)) {
  install.packages("MicrobiomeStat")
}
library(MicrobiomeStat)

# ============================================================
# INSTRUCTIONS:
# 1. Go to posit.cloud and create a free account
# 2. Create a new RStudio project
# 3. Upload these files to the project:
#    - linda_analysis.R
#    - linda_otu_phylum.csv
#    - linda_otu_genus.csv
#    - linda_otu_species.csv
#    - linda_meta_phylum.csv (same for all levels)
# 4. Open linda_analysis.R and click "Source" or run line by line
# ============================================================

run_linda_analysis <- function(otu_file, meta_file, level_name) {
  
  cat("\n", rep("=", 60), "\n")
  cat("LinDA Analysis —", level_name, "Level\n")
  cat(rep("=", 60), "\n")
  
  # Read data
  otu <- read.csv(otu_file, row.names = 1, check.names = FALSE)
  meta <- read.csv(meta_file, row.names = 1)
  
  # Ensure samples match
  common_samples <- intersect(colnames(otu), rownames(meta))
  otu <- otu[, common_samples]
  meta <- meta[common_samples, , drop = FALSE]
  
  cat("Samples:", ncol(otu), "| Taxa:", nrow(otu), "\n")
  cat("Groups:", table(meta$Group), "\n\n")
  
  # Set HC as reference group
  meta$Group <- relevel(factor(meta$Group), ref = "HC")
  
  # Run LinDA
  linda.obj <- linda(
    feature.dat = as.matrix(otu),
    meta.dat = meta,
    formula = "~Group",
    feature.dat.type = "proportion",
    prev.filter = 0.1,
    p.adj.method = "BH",
    alpha = 0.05,
    verbose = FALSE
  )
  
  # Extract results for each group vs HC
  results_list <- list()
  
  for (group_name in names(linda.obj$output)) {
    res <- linda.obj$output[[group_name]]
    res$Taxon <- rownames(res)
    res$Comparison <- group_name
    res$Level <- level_name
    results_list[[group_name]] <- res
  }
  
  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL
  
  # Show significant results
  cat("Nominal p < 0.05:\n")
  sig_nom <- results_df[results_df$pvalue < 0.05, 
                         c("Taxon", "Comparison", "log2FoldChange", "pvalue", "padj")]
  if (nrow(sig_nom) > 0) {
    print(sig_nom[order(sig_nom$pvalue), ], digits = 4)
  } else {
    cat("None\n")
  }
  
  cat("\nFDR-adjusted p < 0.05:\n")
  sig_fdr <- results_df[results_df$padj < 0.05,
                          c("Taxon", "Comparison", "log2FoldChange", "pvalue", "padj")]
  if (nrow(sig_fdr) > 0) {
    print(sig_fdr[order(sig_fdr$padj), ], digits = 4)
  } else {
    cat("None\n")
  }
  
  # Save full results
  output_file <- paste0("linda_results_", tolower(level_name), ".csv")
  write.csv(results_df[, c("Level", "Taxon", "Comparison", 
                             "log2FoldChange", "lfcSE", "pvalue", "padj")],
            output_file, row.names = FALSE)
  cat("\nFull results saved to:", output_file, "\n")
  
  return(results_df)
}

# ============================================================
# RUN ANALYSIS FOR ALL THREE LEVELS
# ============================================================

# Phylum Level
results_phylum <- run_linda_analysis(
  otu_file = "linda_otu_phylum.csv",
  meta_file = "linda_meta_phylum.csv",
  level_name = "Phylum"
)

# Genus Level
results_genus <- run_linda_analysis(
  otu_file = "linda_otu_genus.csv",
  meta_file = "linda_meta_genus.csv",
  level_name = "Genus"
)

# Species Level
results_species <- run_linda_analysis(
  otu_file = "linda_otu_species.csv",
  meta_file = "linda_meta_species.csv",
  level_name = "Species"
)

# ============================================================
# COMBINE ALL RESULTS
# ============================================================

all_results <- rbind(results_phylum, results_genus, results_species)
write.csv(all_results[, c("Level", "Taxon", "Comparison", 
                           "log2FoldChange", "lfcSE", "pvalue", "padj")],
          "linda_results_ALL_LEVELS.csv", row.names = FALSE)

cat("\n\nAll results saved to: linda_results_ALL_LEVELS.csv\n")
cat("Analysis complete.\n")
