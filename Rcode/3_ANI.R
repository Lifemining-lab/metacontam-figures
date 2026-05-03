# =========================
# Boxplot of ANI by Ground Truth — with ALL simulations
# =========================

options(stringsAsFactors = FALSE)

# ---- Packages ----
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  library(ggsci)
})

# =========================
# Paths (staged inputs root)
# =========================
BASE_INPUT <- Sys.getenv("BASE_INPUT", unset = "/media/junwoojo/18T/Submission/Rcode/inputs")

# =========================
# Load GT_species from TSV
#   File: {BASE_INPUT}/config/GT_species.tsv
# =========================
GT_path <- file.path(BASE_INPUT, "config", "GT_species.tsv")
if (!file.exists(GT_path)) stop("Missing GT_species.tsv: ", GT_path)

GT_df <- read.table(GT_path, sep = "\t", header = TRUE,
                    check.names = FALSE, stringsAsFactors = FALSE)
GT_df$group <- as.character(GT_df$group)
GT_df$taxid <- as.character(GT_df$taxid)
GT_species <- split(GT_df$taxid, GT_df$group)

# =========================
# Palette
# =========================
nejm_cols <- pal_nejm("default")(8)
cols_groundtruth <- c(
  "Non-contaminant" = nejm_cols[5],
  "Contaminant"     = nejm_cols[6]
)

# =========================
# Optional quick peek
# =========================
peek_cols <- function(path, n = 5) {
  stopifnot(file.exists(path))
  df <- read.table(path, sep = "\t", header = TRUE,
                   check.names = FALSE, stringsAsFactors = FALSE)
  cat("File:", path, "\n")
  cat("Columns:", paste(colnames(df), collapse = ", "), "\n\n")
  print(utils::head(df, n))
  invisible(df)
}

# =========================
# Core plotting function
# =========================
ani_boxplot <- function(x,
                        highlight_taxid,
                        title_text,
                        ani_col = NULL,
                        id_col  = NULL,
                        y_limits = c(0.94, 1.02),
                        median_y = NULL,
                        tag = NULL,
                        jitter_size = 1.5,
                        jitter_alpha = 1.0) {
  
  # 1) Read input (path or data.frame)
  df <- if (is.character(x)) {
    stopifnot(file.exists(x))
    read.table(x, sep = "\t", header = TRUE,
               check.names = FALSE, stringsAsFactors = FALSE)
  } else if (is.data.frame(x)) {
    x
  } else {
    stop("x must be a file path or a data.frame")
  }
  cols <- colnames(df)
  
  # 2) Auto-detect ID column
  if (is.null(id_col)) {
    id_candidates <- c(
      "Taxid","taxid","TaxID","tax_id","NCBI_taxid","NCBI TaxID",
      "taxon_id","taxon","Taxon", cols[1]
    )
    id_col <- id_candidates[id_candidates %in% cols][1]
  }
  if (is.na(id_col) || !(id_col %in% cols)) {
    stop(sprintf("ID column '%s' not found. Available: %s", id_col, paste(cols, collapse = ", ")))
  }
  
  # 3) Auto-detect ANI column
  if (is.null(ani_col)) {
    ani_candidates <- c(
      "Mean.Pairwise.ANI","Mean-Pairwise-ANI","Mean_Pairwise_ANI","mean_pairwise_ani",
      "ANI_mean","ANI.Mean","ANI","meanANI","pairwise_ani_mean",
      "MeanANI","Mean.ANI","ANI.mean"
    )
    ani_col <- ani_candidates[ani_candidates %in% cols][1]
    
    # Fallback: choose a numeric column with values mostly in [0.85, 1.05]
    if (is.na(ani_col)) {
      num_cols <- cols[sapply(df, is.numeric)]
      if (length(num_cols) > 0) {
        score <- sapply(num_cols, function(cn) {
          v <- suppressWarnings(as.numeric(df[[cn]]))
          sum(is.finite(v) & v >= 0.85 & v <= 1.05, na.rm = TRUE)
        })
        if (any(score > 0)) ani_col <- names(sort(score, decreasing = TRUE))[1]
      }
    }
  }
  if (is.na(ani_col) || !(ani_col %in% cols)) {
    stop(sprintf("ANI column not found. Available: %s", paste(cols, collapse = ", ")))
  }
  
  # 4) Preprocess
  df <- df %>%
    mutate(
      Taxid_chr   = as.character(.data[[id_col]]),
      ANI_value   = suppressWarnings(as.numeric(.data[[ani_col]])),
      groundtruth = ifelse(Taxid_chr %in% highlight_taxid, "Contaminant", "Non-contaminant"),
      groundtruth = factor(groundtruth, levels = c("Non-contaminant","Contaminant"))
    )
  
  # 5) Wilcoxon test (Contaminant < Non-contaminant)
  ng <- nlevels(droplevels(df$groundtruth))
  pval_df <- NULL
  pval_df_num <- NULL
  
  if (ng == 2) {
    wil <- suppressWarnings(
      wilcox.test(ANI_value ~ groundtruth, data = df, alternative = "less", exact = FALSE)
    )
    pval <- wil$p.value
    p_signif <- if (pval <= 0.001) "***" else if (pval <= 0.01) "**" else if (pval <= 0.05) "*" else "ns"
    
    cat(if (!is.null(tag)) paste0("[", tag, "] ") else "",
        "Wilcoxon p-value:", format(pval, digits = 3), "\tsignif:", p_signif, "\n", sep = "")
    
    # Star position
    y_star <- min(max(df$ANI_value, na.rm = TRUE) + 0.003, y_limits[2] - 0.001)
    # Numeric p-value position (above the star)
    y_num  <- min(y_star + 0.003, y_limits[2] - 0.0005)
    
    pval_df <- data.frame(
      group1     = "Contaminant",
      group2     = "Non-contaminant",
      y.position = y_star,
      p.adj      = pval,
      p.signif   = p_signif,
      stringsAsFactors = FALSE
    )
    
    pval_df_num <- transform(
      pval_df,
      y.position = y_num,
      p.label = paste0("p=", format(p.adj, digits = 2, scientific = TRUE))
    )
  } else {
    cat(if (!is.null(tag)) paste0("[", tag, "] ") else "",
        "Wilcoxon p-value: NA (only one group)\n", sep = "")
  }
  
  # 6) Median labels
  if (is.null(median_y)) median_y <- y_limits[2] - 0.005
  
  med_df <- df %>%
    group_by(groundtruth) %>%
    summarise(med_val = median(ANI_value, na.rm = TRUE), .groups = "drop") %>%
    mutate(x_pos = as.numeric(groundtruth), y_pos = median_y)
  
  # 7) Plot
  p <- ggboxplot(
    df, x = "groundtruth", y = "ANI_value",
    color = "groundtruth",
    add = "jitter", add.params = list(size = jitter_size, alpha = jitter_alpha)
  ) +
    scale_color_manual(values = cols_groundtruth, name = NULL) +
    { if (!is.null(pval_df)) ggpubr::stat_pvalue_manual(
      pval_df, label = "p.signif", tip.length = 0.04, size = 5
    ) } +
    { if (!is.null(pval_df_num)) ggpubr::stat_pvalue_manual(
      pval_df_num, label = "p.label", tip.length = 0, size = 4, bracket.size = 0
    ) } +
    geom_text(
      data = med_df,
      aes(x = x_pos, y = y_pos, label = round(med_val, 3)),
      inherit.aes = FALSE, size = 5
    ) +
    scale_y_continuous(limits = y_limits, breaks = seq(0.95, 1.00, 0.05)) +
    ggtitle(title_text) +
    theme_pubr() +
    theme(
      legend.position    = "none",
      plot.title         = element_text(size = 15, face = "bold", hjust = 0.5),
      axis.title.x       = element_blank(),
      axis.title.y       = element_blank(),
      axis.text.x        = element_text(size = 15, angle = 45, hjust = 1),
      axis.text.y        = element_text(size = 15),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.border       = element_rect(color = "black", fill = NA, linewidth = 1)
    )
  
  return(p)
}

# ============================================================
# Dataset paths (staged outputs)
# ============================================================
MI_output              <- file.path(BASE_INPUT, "Metacontam_output_folder", "metacontam_MI_out")
Skin_microbiome_output <- file.path(BASE_INPUT, "Metacontam_output_folder", "metacontam_microbiome_out")
Skin_standard_output   <- file.path(BASE_INPUT, "Metacontam_output_folder", "metacontam_standard_out")
Nasal_output           <- file.path(BASE_INPUT, "Metacontam_output_folder", "metacontam_nasal_out")

# Simulation outputs (staged)
path_Simul_1_50_output  <- file.path(BASE_INPUT, "Metacontam_output_folder", "Simulation", "metacontam_1_50_out")
path_Simul_1_30_output  <- file.path(BASE_INPUT, "Metacontam_output_folder", "Simulation", "metacontam_1_30_out")
path_Simul_1_10_output  <- file.path(BASE_INPUT, "Metacontam_output_folder", "Simulation", "metacontam_1_10_out")

path_Simul_0.5_50_output <- file.path(BASE_INPUT, "Metacontam_output_folder", "Simulation", "metacontam_0.5_50_out")
path_Simul_0.5_30_output <- file.path(BASE_INPUT, "Metacontam_output_folder", "Simulation", "metacontam_0.5_30_out")
path_Simul_0.5_10_output <- file.path(BASE_INPUT, "Metacontam_output_folder", "Simulation", "metacontam_0.5_10_out")

path_Simul_0.1_50_output <- file.path(BASE_INPUT, "Metacontam_output_folder", "Simulation", "metacontam_0.1_50_out")
path_Simul_0.1_30_output <- file.path(BASE_INPUT, "Metacontam_output_folder", "Simulation", "metacontam_0.1_30_out")
path_Simul_0.1_10_output <- file.path(BASE_INPUT, "Metacontam_output_folder", "Simulation", "metacontam_0.1_10_out")

# Final_prediction paths (ANI table)
path_MI_predict               <- file.path(MI_output,              "Final_prediction.txt")
path_Skin_microbiome_predict  <- file.path(Skin_microbiome_output, "Final_prediction.txt")
path_Skin_standard_predict    <- file.path(Skin_standard_output,   "Final_prediction.txt")
path_Nasal_predict            <- file.path(Nasal_output,           "Final_prediction.txt")

path_Simul_1_50_predict   <- file.path(path_Simul_1_50_output,  "Final_prediction.txt")
path_Simul_1_30_predict   <- file.path(path_Simul_1_30_output,  "Final_prediction.txt")
path_Simul_1_10_predict   <- file.path(path_Simul_1_10_output,  "Final_prediction.txt")

path_Simul_0.5_50_predict <- file.path(path_Simul_0.5_50_output, "Final_prediction.txt")
path_Simul_0.5_30_predict <- file.path(path_Simul_0.5_30_output, "Final_prediction.txt")
path_Simul_0.5_10_predict <- file.path(path_Simul_0.5_10_output, "Final_prediction.txt")

path_Simul_0.1_50_predict <- file.path(path_Simul_0.1_50_output, "Final_prediction.txt")
path_Simul_0.1_30_predict <- file.path(path_Simul_0.1_30_output, "Final_prediction.txt")
path_Simul_0.1_10_predict <- file.path(path_Simul_0.1_10_output, "Final_prediction.txt")

# Check required files (fail fast)
need_files <- c(
  path_Skin_microbiome_predict, path_Skin_standard_predict, path_Nasal_predict,
  path_Simul_1_50_predict, path_Simul_1_30_predict, path_Simul_1_10_predict,
  path_Simul_0.5_50_predict, path_Simul_0.5_30_predict, path_Simul_0.5_10_predict,
  path_Simul_0.1_50_predict, path_Simul_0.1_30_predict, path_Simul_0.1_10_predict
)
missing_files <- need_files[!file.exists(need_files)]
if (length(missing_files) > 0) {
  stop("Missing prediction file(s):\n- ", paste(missing_files, collapse = "\n- "))
}

# ============================================================
# Build plots
# ============================================================
PRED_MI    <- path_MI_predict
PRED_SKMIC <- path_Skin_microbiome_predict
PRED_SKSTD <- path_Skin_standard_predict
PRED_NASAL <- path_Nasal_predict

plots_real <- list(
  # ani_boxplot(PRED_MI, GT_species$MI_permissive, "MI",
  #             ani_col = "Mean-Pairwise-ANI", id_col = "Taxid", tag = "MI"),
  ani_boxplot(PRED_SKMIC, GT_species$Skin_microbiome, "Skin (microbiome kit)",
              ani_col = "Mean-Pairwise-ANI", id_col = "Taxid", tag = "Skin (microbiome)"),
  ani_boxplot(PRED_SKSTD, GT_species$skin_standard, "Skin (standard kit)",
              ani_col = "Mean-Pairwise-ANI", id_col = "Taxid", tag = "Skin (standard)"),
  ani_boxplot(PRED_NASAL, GT_species$nasal, "Nasal",
              ani_col = "Mean-Pairwise-ANI", id_col = "Taxid", tag = "Nasal")
)

plots_sim <- list(
  ani_boxplot(path_Simul_1_50_predict,  GT_species$simulation, "1% Spike-in (n=50)",
              ani_col = "Mean-Pairwise-ANI", id_col = "Taxid", tag = "1% (50)"),
  ani_boxplot(path_Simul_1_30_predict,  GT_species$simulation, "1% Spike-in (n=30)",
              ani_col = "Mean-Pairwise-ANI", id_col = "Taxid", tag = "1% (30)"),
  ani_boxplot(path_Simul_1_10_predict,  GT_species$simulation, "1% Spike-in (n=10)",
              ani_col = "Mean-Pairwise-ANI", id_col = "Taxid", tag = "1% (10)"),
  
  ani_boxplot(path_Simul_0.5_50_predict, GT_species$simulation, "0.5% Spike-in (n=50)",
              ani_col = "Mean-Pairwise-ANI", id_col = "Taxid", tag = "0.5% (50)"),
  ani_boxplot(path_Simul_0.5_30_predict, GT_species$simulation, "0.5% Spike-in (n=30)",
              ani_col = "Mean-Pairwise-ANI", id_col = "Taxid", tag = "0.5% (30)"),
  ani_boxplot(path_Simul_0.5_10_predict, GT_species$simulation, "0.5% Spike-in (n=10)",
              ani_col = "Mean-Pairwise-ANI", id_col = "Taxid", tag = "0.5% (10)"),
  
  ani_boxplot(path_Simul_0.1_50_predict, GT_species$simulation, "0.1% Spike-in (n=50)",
              ani_col = "Mean-Pairwise-ANI", id_col = "Taxid", tag = "0.1% (50)"),
  ani_boxplot(path_Simul_0.1_30_predict, GT_species$simulation, "0.1% Spike-in (n=30)",
              ani_col = "Mean-Pairwise-ANI", id_col = "Taxid", tag = "0.1% (30)"),
  ani_boxplot(path_Simul_0.1_10_predict, GT_species$simulation, "0.1% Spike-in (n=10)",
              ani_col = "Mean-Pairwise-ANI", id_col = "Taxid", tag = "0.1% (10)")
)

# Combine all panels
p_all <- wrap_plots(c(plots_real, plots_sim), ncol = 3)
p_all

# Clean axes (do not change plot geometry; only remove labels)
p_clean <- p_all &
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.title.y = element_blank()
    # axis.ticks.x = element_blank()  # Uncomment if you also want to hide x ticks
  )

p_clean

# ---- Save (optional; commented) ----
# ggsave("/media/junwoojo/18T/Metacontam_Figure/all_ani_with_simulations.pdf",
#        plot = p_clean, width = 14, height = 18, units = "in", dpi = 300)

Supp_Fig2 <- p_clean
