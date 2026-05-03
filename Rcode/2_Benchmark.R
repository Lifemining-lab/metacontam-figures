# =========================================================
# Figure code (clean) - based on GT_species.tsv
#   - highlight_list -> GT_species
#   - ggsave() calls remain commented out
# =========================================================

options(stringsAsFactors = FALSE)

# =========================
# Libraries
# =========================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggsci)
  library(patchwork)
  library(readr)
  library(tibble)
})

# =========================
# Paths / Load GT_species
# =========================
BASE_INPUT <- Sys.getenv("BASE_INPUT", unset = "/media/junwoojo/18T/Submission/Rcode/inputs")

GT_path <- file.path(BASE_INPUT, "config", "GT_species.tsv")
if (!file.exists(GT_path)) stop("Missing GT_species.tsv: ", GT_path)

GT_df <- read.table(
  GT_path,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
GT_df$group <- as.character(GT_df$group)
GT_df$taxid <- as.character(GT_df$taxid)

GT_species <- split(GT_df$taxid, GT_df$group)

# =========================
# Required object checks
# =========================
need_objs <- c(
  # Metacontam predictions
  "Skin_standard_predict", "Skin_microbiome_predict", "MI_predict", "Nasal_predict", "Oral_predict",
  "Simul_1_50_predict", "Simul_0.5_50_predict", "Simul_0.1_50_predict",
  "Simul_0.5_30_predict", "Simul_0.5_10_predict",
  
  # Blacklist predictions
  "skin_standard_black_predict", "skin_microbiome_black_predict", "MI_black_predict", "nasal_black_predict",
  
  # Network community predictions
  "MI_community", "skin_standard_community", "skin_microbiome_community", "nasal_community",
  
  # Squeegee predictions
  "squeegee_skin_microbiome_predict", "squeegee_MI_predict", "squeegee_nasal_predict", "squeegee_nasal_oral_predict",
  
  # Decontam MI
  "decontam_MI_predict"
)

missing <- need_objs[!vapply(need_objs, exists, logical(1), inherits = TRUE)]
if (length(missing) > 0) {
  stop(
    "These objects are missing in the environment.\n",
    "Run your main setup script first (the one that creates predictions/communities).\n\n",
    "Missing:\n- ", paste(missing, collapse = "\n- ")
  )
}

# Squeegee skin-standard typo alias support
if (!exists("squegee_skin_standard_predict", inherits = TRUE) &&
    exists("squeegee_skin_standard_predict", inherits = TRUE)) {
  squegee_skin_standard_predict <- get("squeegee_skin_standard_predict", inherits = TRUE)
}
if (exists("squegee_skin_standard_predict", inherits = TRUE) &&
    !exists("squeegee_skin_standard_predict", inherits = TRUE)) {
  squeegee_skin_standard_predict <- get("squegee_skin_standard_predict", inherits = TRUE)
}
if (!exists("squegee_skin_standard_predict", inherits = TRUE) &&
    !exists("squeegee_skin_standard_predict", inherits = TRUE)) {
  stop("Either 'squegee_skin_standard_predict' or 'squeegee_skin_standard_predict' must exist in the environment.")
}

# =========================
# Helpers
# =========================
compute_f1_from_taxid <- function(true_taxid, predicted_taxid) {
  true_set <- unique(as.character(true_taxid))
  pred_set <- unique(as.character(predicted_taxid))
  
  TP <- length(intersect(true_set, pred_set))
  FP <- length(setdiff(pred_set, true_set))
  FN <- length(setdiff(true_set, pred_set))
  
  precision <- if ((TP + FP) > 0) TP / (TP + FP) else 0
  recall    <- if ((TP + FN) > 0) TP / (TP + FN) else 0
  f1        <- if ((precision + recall) > 0) 2 * precision * recall / (precision + recall) else 0
  
  list(
    precision = round(precision, 3),
    recall    = round(recall, 3),
    f1        = round(f1, 3),
    TP = TP, FP = FP, FN = FN,
    precision_num = TP,
    precision_den = TP + FP,
    recall_num    = TP,
    recall_den    = TP + FN,
    f1_num        = 2 * TP,
    f1_den        = 2 * TP + FP + FN
  )
}

read_taxid_file <- function(path) {
  stopifnot(file.exists(path))
  x <- readr::read_lines(path)
  x <- x[!grepl("^\\s*$", x)]
  x <- x[!grepl("^\\s*#", x)]
  x <- sub("\\t.*$", "", x)
  x <- trimws(x)
  unique(x[nzchar(x) & tolower(x) != "taxid"])
}

safe_f1 <- function(truth, pred_name) {
  if (exists(pred_name, inherits = TRUE)) {
    compute_f1_from_taxid(truth, get(pred_name, inherits = TRUE))$f1
  } else {
    NA_real_
  }
}

# =========================
# Colors / Theme
# =========================
nejm_cols <- pal_nejm("default")(8)
FAMILY_COLORS <- c("Blacklist" = nejm_cols[2], "Network" = nejm_cols[3], "ANI" = nejm_cols[4])
METHOD_COLORS <- c("Decontam" = nejm_cols[2], "Squeegee" = nejm_cols[3], "Metacontam" = nejm_cols[4])

theme_f1 <- theme_minimal(base_size = 15) +
  theme(
    axis.title.x       = element_blank(),
    axis.title.y       = element_text(size = 15, face = "bold"),
    axis.text.x        = element_text(size = 16),
    axis.text.y        = element_text(size = 16),
    legend.title       = element_blank(),
    legend.text        = element_text(size = 15),
    legend.position    = c(0.5, 0.97),
    legend.direction   = "horizontal",
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.9),
    plot.title         = element_text(size = 10, hjust = 0.5)
  )

# =========================================================
# 1) Performance calculation
# Standard labels:
# "Skin (standard)", "Skin (microbiome)",
# "MI (permissive)", "MI (strict)", "Nasal", "Oral"
# =========================================================

# --- Blacklist (F1 only) ---
Blacklist_performance <- list(
  "Skin (standard)"   = compute_f1_from_taxid(GT_species$skin_standard,   skin_standard_black_predict)$f1,
  "Skin (microbiome)" = compute_f1_from_taxid(GT_species$Skin_microbiome, skin_microbiome_black_predict)$f1,
  "MI (permissive)"   = compute_f1_from_taxid(GT_species$MI_permissive,   MI_black_predict)$f1,
  "MI (strict)"       = compute_f1_from_taxid(GT_species$MI_strict,       MI_black_predict)$f1,
  "Nasal"             = compute_f1_from_taxid(GT_species$nasal,           nasal_black_predict)$f1
)

# --- Metacontam (ANI / Original) ---
Skin_standard_performance   <- compute_f1_from_taxid(GT_species$skin_standard,   Skin_standard_predict)
Skin_microbiome_performance <- compute_f1_from_taxid(GT_species$Skin_microbiome, Skin_microbiome_predict)
MI_permissive_performance   <- compute_f1_from_taxid(GT_species$MI_permissive,   MI_predict)
MI_strict_performance       <- compute_f1_from_taxid(GT_species$MI_strict,       MI_predict)
Nasal_performance           <- compute_f1_from_taxid(GT_species$nasal,           Nasal_predict)
Oral_performance            <- compute_f1_from_taxid(GT_species$oral,            Oral_predict)

# --- Squeegee ---
squeegee_Skin_standard_performance   <- compute_f1_from_taxid(GT_species$skin_standard,   squegee_skin_standard_predict)
squeegee_Skin_microbiome_performance <- compute_f1_from_taxid(GT_species$Skin_microbiome, squeegee_skin_microbiome_predict)
squeegee_MI_permissive_performance   <- compute_f1_from_taxid(GT_species$MI_permissive,   squeegee_MI_predict)
squeegee_MI_strict_performance       <- compute_f1_from_taxid(GT_species$MI_strict,       squeegee_MI_predict)
squeegee_Nasal_performance           <- compute_f1_from_taxid(GT_species$nasal,           squeegee_nasal_predict)
squeegee_Nasal_fromNO_performance    <- compute_f1_from_taxid(GT_species$nasal,           squeegee_nasal_oral_predict)
squeegee_Oral_fromNO_performance     <- compute_f1_from_taxid(GT_species$oral,            squeegee_nasal_oral_predict)

# --- Decontam ---
decontam_MI_permissive_performance <- compute_f1_from_taxid(GT_species$MI_permissive, decontam_MI_predict)
decontam_MI_strict_performance     <- compute_f1_from_taxid(GT_species$MI_strict,     decontam_MI_predict)

decon_dir <- file.path(BASE_INPUT, "Decontam_output_folder")

if (!exists("decontam_skin_standard_predict", inherits = TRUE)) {
  decontam_skin_standard_predict <- read_taxid_file(file.path(decon_dir, "skin_standard_contaminant_taxids.txt"))
}
if (!exists("decontam_skin_microbiome_predict", inherits = TRUE)) {
  decontam_skin_microbiome_predict <- read_taxid_file(file.path(decon_dir, "skin_microbiome_contaminant_taxids.txt"))
}
if (!exists("decontam_nasal_predict", inherits = TRUE)) {
  decontam_nasal_predict <- read_taxid_file(file.path(decon_dir, "nasal_contaminant_taxids.txt"))
}
if (!exists("decontam_oral_predict", inherits = TRUE)) {
  decontam_oral_predict <- read_taxid_file(file.path(decon_dir, "oral_contaminant_taxids.txt"))
}

decontam_Skin_standard_performance   <- compute_f1_from_taxid(GT_species$skin_standard,   decontam_skin_standard_predict)
decontam_Skin_microbiome_performance <- compute_f1_from_taxid(GT_species$Skin_microbiome, decontam_skin_microbiome_predict)
decontam_Nasal_performance           <- compute_f1_from_taxid(GT_species$nasal,           decontam_nasal_predict)
decontam_Oral_performance            <- compute_f1_from_taxid(GT_species$oral,            decontam_oral_predict)

# =========================================================
# 2) Dataset-wise benchmark plots (Precision / Recall / F1)
# =========================================================
plot_benchmark_metrics <- function(results_list, title_text = "") {
  metrics_df <- do.call(rbind, lapply(names(results_list), function(tool) {
    res <- results_list[[tool]]
    data.frame(
      Tool   = tool,
      Metric = c("Precision", "Recall", "F1"),
      Value  = c(res$precision, res$recall, res$f1)
    )
  })) %>% as.data.frame()
  
  metrics_df$Tool   <- factor(metrics_df$Tool, levels = unique(names(results_list)))
  metrics_df$Metric <- factor(metrics_df$Metric, levels = c("Precision", "Recall", "F1"))
  
  ggplot(metrics_df, aes(x = Tool, y = Value, fill = Metric)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.65) +
    geom_text(
      aes(label = round(Value, 3)),
      position = position_dodge(width = 0.7),
      vjust = -0.5,
      size = 10
    ) +
    scale_fill_viridis_d() +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(title = title_text, x = NULL, y = "") +
    theme_f1
}

# --- Individual dataset plots ---
p_bench_MI_perm <- plot_benchmark_metrics(
  list(
    Decontam   = decontam_MI_permissive_performance,
    Squeegee   = squeegee_MI_permissive_performance,
    Metacontam = MI_permissive_performance
  ),
  title_text = "MI (permissive): Precision / Recall / F1"
)

p_bench_MI_strict <- plot_benchmark_metrics(
  list(
    Decontam   = decontam_MI_strict_performance,
    Squeegee   = squeegee_MI_strict_performance,
    Metacontam = MI_strict_performance
  ),
  title_text = "MI (strict): Precision / Recall / F1"
)

p_bench_skin_microbiome <- plot_benchmark_metrics(
  list(
    Decontam   = decontam_Skin_microbiome_performance,
    Squeegee   = squeegee_Skin_microbiome_performance,
    Metacontam = Skin_microbiome_performance
  ),
  title_text = "Skin (microbiome): Precision / Recall / F1"
)

p_bench_skin_standard <- plot_benchmark_metrics(
  list(
    Decontam   = decontam_Skin_standard_performance,
    Squeegee   = squeegee_Skin_standard_performance,
    Metacontam = Skin_standard_performance
  ),
  title_text = "Skin (standard): Precision / Recall / F1"
)

p_bench_nasal_sq_nasal <- plot_benchmark_metrics(
  list(
    Decontam                = decontam_Nasal_performance,
    "Squeegee (nasal-only)" = squeegee_Nasal_performance,
    Metacontam              = Nasal_performance
  ),
  title_text = "Nasal: Precision / Recall / F1 (Squeegee: nasal-only)"
)

p_bench_nasal_sq_no <- plot_benchmark_metrics(
  list(
    Decontam                = decontam_Nasal_performance,
    "Squeegee (nasal+oral)" = squeegee_Nasal_fromNO_performance,
    Metacontam              = Nasal_performance
  ),
  title_text = "Nasal: Precision / Recall / F1 (Squeegee: nasal+oral)"
)

p_bench_oral <- plot_benchmark_metrics(
  list(
    Decontam                = decontam_Oral_performance,
    "Squeegee (nasal+oral)" = squeegee_Oral_fromNO_performance,
    Metacontam              = Oral_performance
  ),
  title_text = "Oral: Precision / Recall / F1 (Squeegee trained on Nasal+Oral)"
)

# --- Simulation ---
Simulation_1_50_performance   <- compute_f1_from_taxid(GT_species$simulation, Simul_1_50_predict)
Simulation_0.5_50_performance <- compute_f1_from_taxid(GT_species$simulation, Simul_0.5_50_predict)
Simulation_0.1_50_performance <- compute_f1_from_taxid(GT_species$simulation, Simul_0.1_50_predict)

p_bench_sim_spike <- plot_benchmark_metrics(
  list(
    "1%"   = Simulation_1_50_performance,
    "0.5%" = Simulation_0.5_50_performance,
    "0.1%" = Simulation_0.1_50_performance
  ),
  title_text = "Simulation (n = 50): Spike %"
) + scale_y_continuous(limits = c(0, 1.1), expand = c(0, 0))

Simulation_0.5_30_performance <- compute_f1_from_taxid(GT_species$simulation, Simul_0.5_30_predict)
Simulation_0.5_10_performance <- compute_f1_from_taxid(GT_species$simulation, Simul_0.5_10_predict)

p_bench_sim_samples <- plot_benchmark_metrics(
  list(
    "50 samples" = Simulation_0.5_50_performance,
    "30 samples" = Simulation_0.5_30_performance,
    "10 samples" = Simulation_0.5_10_performance
  ),
  title_text = "Simulation (0.5%): Sample Count"
) + scale_y_continuous(limits = c(0, 1.1), expand = c(0, 0))

# =========================================================
# 2.5) MI-only stepwise benchmark (Blacklist / Network / ANI)
# =========================================================
mi_stepwise_results <- list(
  "Blacklist" = compute_f1_from_taxid(GT_species$MI_permissive, MI_black_predict),
  "Network"   = compute_f1_from_taxid(GT_species$MI_permissive, MI_community),
  "ANI"       = compute_f1_from_taxid(GT_species$MI_permissive, MI_predict)
)

p_bench_MI_stepwise <- plot_benchmark_metrics(
  mi_stepwise_results,
  title_text = "MI (permissive): Blacklist vs Network vs ANI (Precision / Recall / F1)"
)

# =========================================================
# 3) Stepwise plots by dataset
# =========================================================
stepwise_results_by_dataset <- list(
  "MI" = list(
    "Blacklist" = compute_f1_from_taxid(GT_species$MI_permissive, MI_black_predict),
    "Network"   = compute_f1_from_taxid(GT_species$MI_permissive, MI_community),
    "ANI"       = compute_f1_from_taxid(GT_species$MI_permissive, MI_predict)
  ),
  "Skin (microbiome)" = list(
    "Blacklist" = compute_f1_from_taxid(GT_species$Skin_microbiome, skin_microbiome_black_predict),
    "Network"   = compute_f1_from_taxid(GT_species$Skin_microbiome, skin_microbiome_community),
    "ANI"       = compute_f1_from_taxid(GT_species$Skin_microbiome, Skin_microbiome_predict)
  ),
  "Skin (standard)" = list(
    "Blacklist" = compute_f1_from_taxid(GT_species$skin_standard, skin_standard_black_predict),
    "Network"   = compute_f1_from_taxid(GT_species$skin_standard, skin_standard_community),
    "ANI"       = compute_f1_from_taxid(GT_species$skin_standard, Skin_standard_predict)
  ),
  "Nasal" = list(
    "Blacklist" = compute_f1_from_taxid(GT_species$nasal, nasal_black_predict),
    "Network"   = compute_f1_from_taxid(GT_species$nasal, nasal_community),
    "ANI"       = compute_f1_from_taxid(GT_species$nasal, Nasal_predict)
  )
)

p_f1_stepwise_MI <- plot_benchmark_metrics(
  stepwise_results_by_dataset[["MI"]],
  title_text = "MI: Blacklist / Network / ANI"
)

p_f1_stepwise_skin_microbiome <- plot_benchmark_metrics(
  stepwise_results_by_dataset[["Skin (microbiome)"]],
  title_text = "Skin (microbiome): Blacklist / Network / ANI"
)

p_f1_stepwise_skin_standard <- plot_benchmark_metrics(
  stepwise_results_by_dataset[["Skin (standard)"]],
  title_text = "Skin (standard): Blacklist / Network / ANI"
)

p_f1_stepwise_nasal <- plot_benchmark_metrics(
  stepwise_results_by_dataset[["Nasal"]],
  title_text = "Nasal: Blacklist / Network / ANI"
)

# =========================================================
# 4) Tool-level F1 (Decontam / Squeegee / Metacontam)
# =========================================================
Original_performance <- list(
  "Skin (standard)"   = Skin_standard_performance$f1,
  "Skin (microbiome)" = Skin_microbiome_performance$f1,
  "MI (permissive)"   = MI_permissive_performance$f1,
  "MI (strict)"       = MI_strict_performance$f1,
  "Nasal"             = Nasal_performance$f1
)

Squeegee_f1 <- list(
  "Skin (standard)"   = squeegee_Skin_standard_performance$f1,
  "Skin (microbiome)" = squeegee_Skin_microbiome_performance$f1,
  "MI (permissive)"   = squeegee_MI_permissive_performance$f1,
  "MI (strict)"       = squeegee_MI_strict_performance$f1,
  "Nasal"             = squeegee_Nasal_performance$f1
)

Decontam_f1 <- list(
  "MI (strict)"     = decontam_MI_strict_performance$f1,
  "MI (permissive)" = decontam_MI_permissive_performance$f1
)

df_f1_methods <- bind_rows(
  tibble::tibble(Dataset = names(Original_performance), F1 = unlist(Original_performance), Method = "Metacontam"),
  tibble::tibble(Dataset = names(Squeegee_f1),         F1 = unlist(Squeegee_f1),         Method = "Squeegee"),
  tibble::tibble(Dataset = names(Decontam_f1),         F1 = unlist(Decontam_f1),         Method = "Decontam")
) %>%
  mutate(
    Dataset = factor(Dataset, levels = c("MI (permissive)", "MI (strict)", "Skin (microbiome)", "Skin (standard)", "Nasal")),
    Method  = factor(Method,  levels = c("Decontam", "Squeegee", "Metacontam"))
  )

p_f1_methods <- ggplot(df_f1_methods, aes(x = Dataset, y = F1, fill = Method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(
    aes(label = round(F1, 3)),
    position = position_dodge(width = 0.8),
    vjust = -0.5,
    size = 5
  ) +
  scale_fill_manual(values = METHOD_COLORS) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "F1 by Tool across Datasets", y = "F1 Score") +
  theme_f1

# =========================================================
# 5) Step-wise F1 summary
# Left panel: MI-only stepwise metric plot
# Right panel: F1-only comparison across datasets
# =========================================================
Network_performance <- list(
  "Skin (standard)"   = compute_f1_from_taxid(GT_species$skin_standard,   skin_standard_community)$f1,
  "Skin (microbiome)" = compute_f1_from_taxid(GT_species$Skin_microbiome, skin_microbiome_community)$f1,
  "MI (permissive)"   = compute_f1_from_taxid(GT_species$MI_permissive,   MI_community)$f1,
  "MI (strict)"       = compute_f1_from_taxid(GT_species$MI_strict,       MI_community)$f1,
  "Nasal"             = compute_f1_from_taxid(GT_species$nasal,           nasal_community)$f1
)

df_f1_stepwise_summary <- dplyr::bind_rows(
  tibble::tibble(MethodFamily = "Blacklist", Dataset = names(Blacklist_performance), F1 = unlist(Blacklist_performance)),
  tibble::tibble(MethodFamily = "Network",   Dataset = names(Network_performance),   F1 = unlist(Network_performance)),
  tibble::tibble(MethodFamily = "ANI",       Dataset = names(Original_performance),  F1 = unlist(Original_performance))
)

df_f1_stepwise <- df_f1_stepwise_summary %>%
  dplyr::filter(!(Dataset %in% c("MI (strict)", "Oral"))) %>%
  dplyr::mutate(
    Dataset = dplyr::if_else(Dataset == "MI (permissive)", "MI", Dataset),
    MethodFamily = factor(MethodFamily, levels = c("Blacklist", "Network", "ANI"))
  ) %>%
  tidyr::drop_na(F1)

ds_levels <- c("MI", "Skin (microbiome)", "Skin (standard)", "Nasal")
df_f1_stepwise$Dataset <- factor(
  df_f1_stepwise$Dataset,
  levels = intersect(ds_levels, unique(df_f1_stepwise$Dataset))
)

p_f1_stepwise_summary <- ggplot2::ggplot(
  df_f1_stepwise,
  ggplot2::aes(x = Dataset, y = F1, fill = MethodFamily)
) +
  ggplot2::geom_col(
    position = ggplot2::position_dodge(width = 0.8),
    width = 0.7
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = round(F1, 3)),
    position = ggplot2::position_dodge(width = 0.8),
    vjust = -0.5,
    size = 5
  ) +
  ggplot2::scale_fill_manual(values = FAMILY_COLORS) +
  ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  ggplot2::labs(title = "Step-wise F1 across Methods", y = "F1 Score") +
  theme_f1

# Final stepwise panel requested by the user
p_f1_stepwise <- p_f1_stepwise_MI | p_f1_stepwise_summary

# =========================================================
# 6) Panel assembly
# =========================================================
panel_benchmarks <-
  (p_bench_MI_perm | p_bench_MI_strict) /
  (p_bench_skin_microbiome | p_bench_skin_standard) /
  (p_bench_nasal_sq_nasal | p_bench_nasal_sq_no) /
  (p_bench_oral | p_bench_sim_spike) /
  (p_bench_sim_samples | patchwork::plot_spacer())

panel_summary <- p_f1_stepwise / p_f1_methods + plot_layout(heights = c(1, 1))

# =========================================================
# 7) Console output of fractions (Precision / Recall / F1)
# =========================================================
print_metrics_fractions <- function(name, res) {
  cat(sprintf(
    "[%s]\n  Precision: %d/%d, Recall: %d/%d, F1: %d/%d\n",
    name,
    res$precision_num, res$precision_den,
    res$recall_num,    res$recall_den,
    res$f1_num,        res$f1_den
  ))
}

print_block <- function(title, named_list) {
  cat("\n==============================\n")
  cat(sprintf("%s\n", title))
  cat("==============================\n")
  for (nm in names(named_list)) {
    print_metrics_fractions(nm, named_list[[nm]])
  }
}

print_block("Original (Metacontam)", list(
  "Skin (standard)"   = Skin_standard_performance,
  "Skin (microbiome)" = Skin_microbiome_performance,
  "MI (permissive)"   = MI_permissive_performance,
  "MI (strict)"       = MI_strict_performance,
  "Nasal"             = Nasal_performance,
  "Oral"              = Oral_performance
))

print_block("Squeegee (dataset-specific)", list(
  "Skin (standard)"    = squeegee_Skin_standard_performance,
  "Skin (microbiome)"  = squeegee_Skin_microbiome_performance,
  "MI (permissive)"    = squeegee_MI_permissive_performance,
  "MI (strict)"        = squeegee_MI_strict_performance,
  "Nasal (nasal-only)" = squeegee_Nasal_performance
))

print_block("Squeegee (nasal+oral trained)", list(
  "Nasal (from NO)" = squeegee_Nasal_fromNO_performance,
  "Oral  (from NO)" = squeegee_Oral_fromNO_performance
))

print_block("Decontam", list(
  "Skin (standard)"   = decontam_Skin_standard_performance,
  "Skin (microbiome)" = decontam_Skin_microbiome_performance,
  "MI (permissive)"   = decontam_MI_permissive_performance,
  "MI (strict)"       = decontam_MI_strict_performance,
  "Nasal"             = decontam_Nasal_performance,
  "Oral"              = decontam_Oral_performance
))

# =========================================================
# 8) Add Nasal Squeegee bars to the Oral panel
# =========================================================
p_bench_oral_nasal_merge <- plot_benchmark_metrics(
  list(
    Decontam                        = decontam_Oral_performance,
    "Squeegee (nasal+oral) [Oral]"  = squeegee_Oral_fromNO_performance,
    "Squeegee (nasal+oral) [Nasal]" = squeegee_Nasal_fromNO_performance,
    Metacontam                      = Oral_performance
  ),
  title_text = "Oral: Precision / Recall / F1 (Squeegee trained on Nasal+Oral) + Nasal Squeegee bars"
)

# =========================================================
# Preview
# =========================================================
# panel_benchmarks
# panel_summary
# p_f1_stepwise
# p_f1_stepwise_summary
# p_f1_methods
# p_bench_MI_stepwise
# p_bench_oral_nasal_merge

Figure2d   <- p_bench_MI_perm
Figure3d   <- p_bench_skin_standard
Figure3e   <- p_bench_skin_microbiome
Figure4a   <- p_bench_nasal_sq_nasal
Supp_Fig5  <- p_bench_oral_nasal_merge
Supp_Fig8 <- p_f1_stepwise
Supp_Fig10   <- p_bench_sim_samples
Supp_Fig10   <- p_bench_sim_spike
