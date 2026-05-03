# =========================
# S. constellatus vs M. radiodurans ANI comparison
#   - Left plot labels: mean
#   - Right plot labels: median
#   - Bracket/star moved slightly lower
#   - Value labels placed higher to avoid overlap
#   - Uses staged paths (BASE_INPUT)
#   - Uses GT_species (loaded from {BASE_INPUT}/config/GT_species.tsv)
# =========================

options(stringsAsFactors = FALSE)

# ===== Global style =====
BORDER_LW   <- 0.6
BOX_LW      <- 0.9
BRACKET_LW  <- 0.2
DOT_SIZE    <- 2.5
DOT_ALPHA   <- 1.0
STAR_SIZE   <- 5

# ===== FIXED label/bracket geometry (same on BOTH panels) =====
# Order: bracket -> star -> value
Y_LIMITS          <- c(0.94, 1.025)

# bracket / star lower
BRACKET_Y         <- Y_LIMITS[2] - 0.012
BRACKET_TIP_DOWN  <- 0.002
STAR_OFFSET       <- 0.0015

# value labels higher and independent
VALUE_LABEL_Y     <- Y_LIMITS[2] - 0.003

suppressPackageStartupMessages({
  library(stringr)
  library(dplyr)
  library(ggpubr)
  library(ggplot2)
  library(ggsci)
  library(patchwork)
  library(cowplot)
})

# ---------- 0) Settings ----------
ALT <- "less"  # Wilcoxon direction
BASE_INPUT <- Sys.getenv("BASE_INPUT", unset = "/media/junwoojo/18T/Submission/Rcode/inputs")

# ---------- Resolve input table path ----------
mi_pair_path <- local({
  candidates <- c(
    file.path(BASE_INPUT, "Metacontam_output_folder", "metacontam_MI_out", "merged_IS_compare_Table.tsv"),
    "/media/junwoojo/18T/Metacontam_dataset/Metacontam_output_folder/metacontam_MI_out/merged_IS_compare_Table.tsv"
  )
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) {
    stop(
      "Missing merged_IS_compare_Table.tsv. Tried:\n- ",
      paste(candidates, collapse = "\n- ")
    )
  }
  hit
})

# ---------- 1) Helpers ----------
read_tsv_smart <- function(path) {
  stopifnot(file.exists(path))
  read.table(path, sep = "\t", header = TRUE,
             check.names = FALSE, stringsAsFactors = FALSE)
}

get_col <- function(df, candidates, required = FALSE) {
  hit <- candidates[candidates %in% names(df)][1]
  if (is.na(hit) && required) {
    stop(sprintf("Required column not found. Tried: %s", paste(candidates, collapse = ", ")))
  }
  hit
}

extract_ani_by_taxid <- function(df, target_taxid, min_count = 50) {
  col_scaf  <- get_col(df, c("scaffold", "contig", "sequence", "name"))
  col_count <- get_col(df, c("compared_bases_count", "aligned_bases", "Compared.Bases.Count"))
  col_ani   <- get_col(df, c("conANI", "Mean-Pairwise-ANI", "Mean.Pairwise.ANI", "ANI",
                             "ANI_mean", "MeanANI", "ANI.mean"), required = TRUE)
  
  df[[col_ani]] <- suppressWarnings(as.numeric(df[[col_ani]]))
  
  if (!is.na(col_count)) {
    df <- df %>% filter(.data[[col_count]] > min_count)
  }
  
  # (A) Extract taxid from scaffold like "taxid|XXXX"
  if (!is.na(col_scaf)) {
    tx <- str_extract(df[[col_scaf]], "(?<=taxid\\|)\\d+")
    out <- df %>%
      mutate(.target = tx) %>%
      filter(.target == as.character(target_taxid)) %>%
      pull(.data[[col_ani]]) %>%
      stats::na.omit()
    if (length(out) > 0) return(out)
  }
  
  # (B) If there is an explicit Taxid column
  col_tid <- get_col(df, c("Taxid", "taxid", "TaxID", "tax_id", "NCBI_taxid", "taxon_id"))
  if (!is.na(col_tid)) {
    out <- df %>%
      filter(as.character(.data[[col_tid]]) == as.character(target_taxid)) %>%
      pull(.data[[col_ani]]) %>%
      stats::na.omit()
    return(out)
  }
  
  stop("No scaffold/Taxid column found to identify taxid.")
}

# ----- manual bracket drawer -----
add_fixed_bracket <- function(p, x1, x2, label_text,
                              y = BRACKET_Y,
                              tip_down = BRACKET_TIP_DOWN,
                              star_offset = STAR_OFFSET) {
  p +
    geom_segment(aes(x = x1, xend = x2, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = BRACKET_LW,
                 lineend = "butt", linejoin = "mitre") +
    geom_segment(aes(x = x1, xend = x1, y = y - tip_down, yend = y),
                 inherit.aes = FALSE, linewidth = BRACKET_LW,
                 lineend = "butt", linejoin = "mitre") +
    geom_segment(aes(x = x2, xend = x2, y = y - tip_down, yend = y),
                 inherit.aes = FALSE, linewidth = BRACKET_LW,
                 lineend = "butt", linejoin = "mitre") +
    annotate("text", x = (x1 + x2) / 2, y = y + star_offset,
             label = label_text, size = STAR_SIZE)
}

# ---------- 2) Load GT_species ----------
GT_path <- file.path(BASE_INPUT, "config", "GT_species.tsv")
if (!file.exists(GT_path)) stop("Missing GT_species.tsv: ", GT_path)

GT_df <- read.table(GT_path, sep = "\t", header = TRUE,
                    check.names = FALSE, stringsAsFactors = FALSE)
GT_df$group <- as.character(GT_df$group)
GT_df$taxid <- as.character(GT_df$taxid)
GT_species <- split(GT_df$taxid, GT_df$group)

# ---------- 3) Load data ----------
mi_pair_df <- read_tsv_smart(mi_pair_path)

# ---------- 4) Extract ANI vectors ----------
S_constellatus_ani <- extract_ani_by_taxid(mi_pair_df, target_taxid = "76860",   min_count = 50)
M_radiodurans_ani  <- extract_ani_by_taxid(mi_pair_df, target_taxid = "2202828", min_count = 50)

if (length(S_constellatus_ani) == 0 || length(M_radiodurans_ani) == 0) {
  stop("One or more ANI vectors are empty. Check the input table and filters.")
}

ani_df <- data.frame(
  ANI = c(S_constellatus_ani, M_radiodurans_ani),
  Species = factor(
    c(rep("S. constellatus", length(S_constellatus_ani)),
      rep("M. radiodurans",  length(M_radiodurans_ani))),
    levels = c("S. constellatus", "M. radiodurans")
  )
)

# ---------- 5) Statistics (one-sided) ----------
x <- ani_df$ANI[ani_df$Species == "S. constellatus"]
y <- ani_df$ANI[ani_df$Species == "M. radiodurans"]

wilcox_res <- wilcox.test(x, y, alternative = ALT, exact = FALSE)
pval <- wilcox_res$p.value
cat(sprintf("[Wilcoxon %s]  p = %.6g (S.constellatus vs M.radiodurans)\n", ALT, pval))

p_signif <- if (pval <= 0.001) "***" else if (pval <= 0.01) "**" else if (pval <= 0.05) "*" else "ns"

# ---------- 6) Palette ----------
d3_cols <- pal_d3("category20")(20)
species_colors <- c(
  "S. constellatus" = d3_cols[10],
  "M. radiodurans"  = d3_cols[6]
)

# ---------- 7) Left plot (MEAN label only here) ----------
means_left <- ani_df %>%
  group_by(Species) %>%
  summarise(mean_val = mean(ANI, na.rm = TRUE), .groups = "drop") %>%
  mutate(x_pos = as.numeric(Species), y_pos = VALUE_LABEL_Y)

Sc_vs_Mr_ANI_compare <- ggboxplot(
  ani_df, x = "Species", y = "ANI",
  color = "Species", palette = species_colors,
  add = "jitter",
  add.params = list(size = DOT_SIZE, alpha = DOT_ALPHA),
  size = BOX_LW
) +
  geom_text(
    data = means_left,
    aes(x = x_pos, y = y_pos, label = round(mean_val, 3)),
    inherit.aes = FALSE, size = 5
  ) +
  scale_y_continuous(limits = Y_LIMITS, breaks = seq(0.95, 1.00, 0.05)) +
  theme_pubr() +
  theme(
    plot.title         = element_text(size = 15, hjust = 0.5, face = "bold"),
    axis.title.x       = element_blank(),
    axis.title.y       = element_blank(),
    axis.text.x        = element_text(size = 15, angle = 45, hjust = 1),
    axis.text.y        = element_text(size = 15),
    legend.position    = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "black", fill = NA, linewidth = BORDER_LW),
    aspect.ratio       = 1.5
  )

x1_left <- match("S. constellatus", levels(ani_df$Species))
x2_left <- match("M. radiodurans",  levels(ani_df$Species))
Sc_vs_Mr_ANI_compare <- add_fixed_bracket(Sc_vs_Mr_ANI_compare, x1_left, x2_left, p_signif)

print(Sc_vs_Mr_ANI_compare)

# =========================================================
# Right plot helper
#   - RIGHT plot stays as MEDIAN
# =========================================================
ani_boxplot <- function(x,
                        highlight_taxid,
                        title_text,
                        ani_col = NULL,
                        id_col  = NULL,
                        y_limits = Y_LIMITS,
                        tag = NULL,
                        jitter_size = DOT_SIZE,
                        jitter_alpha = DOT_ALPHA,
                        show_p_numeric = FALSE) {
  
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
  
  # ID column
  if (is.null(id_col)) {
    id_candidates <- c("Taxid","taxid","TaxID","tax_id","NCBI_taxid","NCBI TaxID",
                       "taxon_id","taxon","Taxon", cols[1])
    id_col <- id_candidates[id_candidates %in% cols][1]
  }
  if (is.na(id_col) || !(id_col %in% cols)) {
    stop(sprintf("ID column '%s' not found. Available: %s", id_col, paste(cols, collapse = ", ")))
  }
  
  # ANI column
  if (is.null(ani_col)) {
    ani_candidates <- c(
      "Mean.Pairwise.ANI","Mean-Pairwise-ANI","Mean_Pairwise_ANI","mean_pairwise_ani",
      "ANI_mean","ANI.Mean","ANI","meanANI","pairwise_ani_mean",
      "MeanANI","Mean.ANI","ANI.mean"
    )
    ani_col <- ani_candidates[ani_candidates %in% cols][1]
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
  if (is.na(ani_col) || !(ani_col %in% colnames(df))) {
    stop(sprintf("ANI column not found. Available: %s", paste(cols, collapse = ", ")))
  }
  
  df <- df %>%
    mutate(
      Taxid_chr   = as.character(.data[[id_col]]),
      ANI_value   = suppressWarnings(as.numeric(.data[[ani_col]])),
      groundtruth = ifelse(Taxid_chr %in% highlight_taxid, "Contaminant", "Non-contaminant"),
      groundtruth = factor(groundtruth, levels = c("Non-contaminant","Contaminant"))
    )
  
  # Wilcoxon (Contaminant < Non-contaminant)
  p_signif <- "ns"
  pval <- NA_real_
  ng <- nlevels(droplevels(df$groundtruth))
  if (ng == 2) {
    wil <- suppressWarnings(wilcox.test(ANI_value ~ groundtruth, data = df,
                                        alternative = "less", exact = FALSE))
    pval <- wil$p.value
    p_signif <- if (pval <= 0.001) "***" else if (pval <= 0.01) "**" else if (pval <= 0.05) "*" else "ns"
    
    cat(if (!is.null(tag)) paste0("[", tag, "] ") else "",
        "Wilcoxon p-value:", format(pval, digits = 3), "\tsignif:", p_signif, "\n", sep = "")
  } else {
    cat(if (!is.null(tag)) paste0("[", tag, "] ") else "",
        "Wilcoxon p-value: NA (only one group)\n", sep = "")
  }
  
  medians_right <- df %>%
    group_by(groundtruth) %>%
    summarise(median_val = median(ANI_value, na.rm = TRUE), .groups = "drop") %>%
    mutate(x_pos = as.numeric(groundtruth), y_pos = VALUE_LABEL_Y)
  
  nejm_cols <- pal_nejm("default")(8)
  cols_groundtruth <- c("Non-contaminant" = nejm_cols[5], "Contaminant" = nejm_cols[6])
  
  p <- ggboxplot(
    df, x = "groundtruth", y = "ANI_value",
    color = "groundtruth",
    add = "jitter",
    add.params = list(size = jitter_size, alpha = jitter_alpha),
    size = BOX_LW
  ) +
    scale_color_manual(values = cols_groundtruth, name = NULL) +
    geom_text(
      data = medians_right,
      aes(x = x_pos, y = y_pos, label = round(median_val, 3)),
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
      panel.border       = element_rect(color = "black", fill = NA, linewidth = BORDER_LW)
    )
  
  x1 <- match("Non-contaminant", levels(df$groundtruth))
  x2 <- match("Contaminant",     levels(df$groundtruth))
  p <- add_fixed_bracket(p, x1, x2, p_signif)
  
  if (isTRUE(show_p_numeric) && is.finite(pval)) {
    p <- p + annotate("text", x = (x1 + x2) / 2,
                      y = VALUE_LABEL_Y + 0.002,
                      label = paste0("p=", format(pval, digits = 2, scientific = TRUE)),
                      size = 4)
  }
  
  p
}

# =========================================================
# Build p_ANI2 (left + right) and align
# =========================================================
MI_output <- file.path(BASE_INPUT, "Metacontam_output_folder", "metacontam_MI_out")
PRED_MI   <- file.path(MI_output, "Final_prediction.txt")

p_ANI2 <- if (file.exists(PRED_MI)) {
  
  left_plot <- Sc_vs_Mr_ANI_compare + theme(aspect.ratio = 1.5)
  
  right_plot <- ani_boxplot(
    PRED_MI, GT_species$MI_permissive, "",
    ani_col = "Mean-Pairwise-ANI", id_col = "Taxid", tag = "MI",
    show_p_numeric = FALSE
  ) + theme(aspect.ratio = 1.5)
  
  left_clean  <- left_plot  + theme(axis.title = element_blank(), axis.text = element_blank())
  right_clean <- right_plot + theme(axis.title = element_blank(), axis.text = element_blank())
  
  aligned <- cowplot::align_plots(left_clean, right_clean, align = "hv", axis = "tblr")
  p_equal <- cowplot::plot_grid(aligned[[2]], aligned[[1]],
                                nrow = 1, rel_widths = c(1, 1), align = "h")
  print(p_equal)
  p_equal
}

Fig2b_c <- p_ANI2