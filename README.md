# metacontam-figures

Analysis code for figures in the Metacontam manuscript.

> Jo J, Lee H, Baek JW, Lee S, Singh V, Shoaie S, Mardinoglu A, Choi J, Lee S. 
> **Metacontam: A Negative Control-Free Decontamination Method for Metagenomic Analysis.** 
> *bioRxiv* 2026.04.26.720876; doi: https://doi.org/10.64898/2026.04.26.720876

For the Metacontam tool itself: https://github.com/Lifemining-lab/Metacontam

## Repository Structure

```
metacontam-figures/
├── Rcode/                  # R scripts (1_Env.R → 12_Mouse_experiment.R)
├── Python_code/            # Tissue microbiome classification notebook
└── inputs/                 # Input data files
```

## Usage

### 1. Set the input directory

All scripts read input data from `BASE_INPUT`. Set this to the location of the `inputs/` folder before running:

```r
Sys.setenv(BASE_INPUT = "/path/to/metacontam-figures/inputs")
```

### 2. Run R scripts in order

Scripts must be executed sequentially — `1_Env.R` loads shared objects used by the others.

```r
source("Rcode/1_Env.R")
source("Rcode/2_Benchmark.R")
source("Rcode/3_ANI.R")
# ... continue through 12_Mouse_experiment.R
```

### 3. Run the Python notebook

For tissue microbiome classification:

```bash
jupyter notebook Python_code/Tissue_microbiomoe_classification.ipynb
```

## Script-to-Figure Mapping

| Script | Description | Figures |
|--------|-------------|---------|
| `1_Env.R` | Environment setup, input paths, count matrices | — |
| `2_Benchmark.R` | Performance comparison (Decontam, Squeegee, Metacontam) and stepwise (Blacklist → Network → ANI) | Fig 2d, 3d, 3e, 4a; Supp 5, 8, 10a |
| `3_ANI.R` | ANI distribution: contaminant vs non-contaminant (boxplots + Wilcoxon) | Supp Fig 2 |
| `4_ANI2.R` | Per-species ANI comparison (*S. constellatus* vs *M. radiodurans*) | Fig 2b, 2c |
| `5_Network.R` | Co-occurrence network with Louvain communities | Fig 2a |
| `6_MI_prevalence_abundance.R` | Abundance & prevalence distribution of detected contaminants vs ground truth | Fig 2e, 2f |
| `7_diversity.R` | Alpha diversity (Shannon, Simpson) before/after decontamination | Supp Fig 3 |
| `8_Skin_negative.R` | PCoA + PERMANOVA across extraction kits | Fig 3a; Supp Fig 4 |
| `9_Skin_replicate.R` | Pairwise Jaccard distance between standard/microbiome kit replicates | Fig 3b, 3c |
| `10_Nasal.R` | Nasal contaminant detection (*Moraxella* case study) | Fig 4b, 4c, 4d |
| `11_Groundtruth.R` | Whitelist ratio across prevalence thresholds | Supp Fig 9 |
| `12_Mouse_experiment.R` | BALF/Fecal/Kit contamination tracing in mouse experiment | Supp Fig 10b, 10c, 10d |
| `Tissue_microbiomoe_classification.ipynb` | Tissue microbiome classification | Fig 5b, 5c; Supp Fig 6a, 6b, 6c, 6d |

## Requirements

**R (≥ 4.0)**

```r
install.packages(c(
  "dplyr", "tidyr", "tibble", "purrr", "readr", "stringr",
  "ggplot2", "ggsci", "ggpubr", "patchwork", "cowplot", "scales",
  "igraph", "ggraph", "tidygraph",
  "vegan", "ape", "data.table", "openxlsx"
))

# Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("NetCoMi")
```

**Python (≥ 3.8)**

```bash
pip install jupyter pandas numpy scikit-learn matplotlib seaborn
```

## Citation

If you use this code, please cite the Metacontam preprint:

> Jo J, Lee H, Baek JW, Lee S, Singh V, Shoaie S, Mardinoglu A, Choi J, Lee S. 
> *bioRxiv* 2026.04.26.720876; doi: https://doi.org/10.64898/2026.04.26.720876

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

## Contact

For questions or issues, please open a [GitHub Issue](https://github.com/Lifemining-lab/metacontam-figures/issues).
 
 
