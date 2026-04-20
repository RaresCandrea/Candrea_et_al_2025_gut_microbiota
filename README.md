# Gut Microbiota Comparative Analysis — Candrea et al.

## Comparative Analysis of Gut Microbiota Patterns in Irritable Bowel Syndrome, Anxiety, and Autoimmune Disorders

**Journal:** Biomedicines  
**Authors:** Adelin-Rareș Candrea, Laura Ioana Gavrilaș, Andrei Mocan, Adriana Rusu, Gianina Crișan

---

## Repository Contents

This repository contains all R scripts used for the statistical analyses and figure generation reported in the manuscript.

| Script | Description |
|--------|-------------|
| `01_data_preparation.R` | Prepare input files from raw Excel data for all analyses |
| `02_LinDA_differential_abundance.R` | LinDA alternative differential abundance analysis (phylum, genus, species) |
| `03_PERMANOVA_beta_diversity.R` | PERMANOVA beta diversity analysis (Bray-Curtis, Jaccard, Jensen-Shannon) |
| `04_Spearman_correlations.R` | Spearman rank correlation analysis with FDR correction |
| `05_figures.R` | Generation of all manuscript figures (Figures 1-6) |

---

## Data

The anonymized dataset (relative abundance values at phylum, genus, and species level; short-chain fatty acid measurements; clinical and anthropometric parameters for all 59 participants) is available as **Supplementary Table S3** in the manuscript.

Input files required to run the scripts:
- `De analizat.xlsx` — microbiota relative abundance data (Phylum, Genus, Species sheets)
- `Adriana Rusu_Baza de date pacienti.xlsx` — full clinical dataset

**Note:** Raw 16S amplicon sequencing files (FASTQ) were not transferred to the research team, as gut microbiota analyses were performed by a certified external laboratory (Bioclinica, Romania, in partnership with Ganzimmun Diagnostics, Germany). The research team had access exclusively to processed relative abundance outputs.

---

## Software Requirements

- R version 4.5.3
- MicrobiomeStat v1.18 (LinDA analysis)
- vegan (PERMANOVA analysis)
- ggplot2, ggpubr (figures)
- readxl, dplyr, tidyr (data preparation)
- corrplot (correlation heatmap)
- pheatmap (heatmap figures)

Install all required packages:

```r
install.packages(c("MicrobiomeStat", "vegan", "ggplot2", "ggpubr",
                   "readxl", "dplyr", "tidyr", "corrplot", "pheatmap",
                   "FSA", "RColorBrewer"))
```

---

## How to Reproduce

1. Clone this repository
2. Download the anonymized dataset from Supplementary Table S3
3. Place the data files in the same directory as the scripts
4. Run scripts in order: 01 → 02 → 03 → 04 → 05

---

## Citation

If you use these scripts, please cite:

> Candrea A.R., Gavrilaș L.I., Mocan A., Rusu A., Crișan G. Comparative Analysis of Gut Microbiota Patterns in Irritable Bowel Syndrome, Anxiety, and Autoimmune Disorders. *Biomedicines* 2025.
