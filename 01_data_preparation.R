# ============================================================
# Script 01: Data Preparation
# Candrea et al. - Gut Microbiota Comparative Analysis
# Biomedicines 2025
# ============================================================
# Description:
# Prepares input CSV files from the raw Excel dataset for use
# in LinDA, PERMANOVA, Spearman correlation, and figure scripts.
#
# Input:
#   - De analizat.xlsx (sheets: Phylum, Genus, Species)
#
# Output:
#   - linda_otu_phylum.csv
#   - linda_otu_genus.csv
#   - linda_otu_species.csv
#   - linda_meta_phylum.csv (metadata, same for all levels)
#   - linda_meta_genus.csv
#   - linda_meta_species.csv
# ============================================================

library(readxl)

prepare_level <- function(sheet_name, level_name) {
  
  cat("Preparing", level_name, "level data...\n")
  
  df <- read_excel("De analizat.xlsx", sheet = sheet_name)
  
  # Identify taxa columns
  taxa_cols <- colnames(df)[!colnames(df) %in% c("Subject", "Group")]
  
  # Clean taxa names
  taxa_names <- gsub(" \\(%\\)|\\(%\\)", "", taxa_cols)
  taxa_names <- trimws(taxa_names)
  taxa_names <- gsub(" ", "_", taxa_names)
  taxa_names <- gsub("[^A-Za-z0-9_]", "", taxa_names)
  
  subjects <- df$Subject
  groups   <- df$Group
  
  # OTU table: taxa on rows, samples on columns
  otu <- as.data.frame(t(df[, taxa_cols]))
  colnames(otu) <- subjects
  rownames(otu) <- taxa_names
  
  # Metadata
  meta <- data.frame(Group = groups, row.names = subjects)
  
  # Save
  otu_file  <- paste0("linda_otu_",  tolower(level_name), ".csv")
  meta_file <- paste0("linda_meta_", tolower(level_name), ".csv")
  
  write.csv(otu,  otu_file,  row.names = TRUE)
  write.csv(meta, meta_file, row.names = TRUE)
  
  cat("  Saved:", otu_file, "(", nrow(otu), "taxa x", ncol(otu), "samples )\n")
  cat("  Saved:", meta_file, "\n")
  cat("  Groups:", table(meta$Group), "\n\n")
}

# Prepare all three levels
prepare_level("Phylum",  "Phylum")
prepare_level("Genus",   "Genus")
prepare_level("Species", "Species")

cat("Data preparation complete.\n")
cat("Files ready for LinDA, PERMANOVA, and Spearman analyses.\n")
