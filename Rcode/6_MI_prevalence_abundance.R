# ============================================================
# MI — Abundance / Prevalence
# GT: histogram (green bars, count)
# Tools: binned count lines + filled markers
# Compare: Decontam, Squeegee, Metacontam
# ============================================================

# ===== Global style =====
BORDER_LW <- 0.8
CURVE_LW  <- 1.8
PT_SIZE   <- 3.6
PT_STROKE <- 1.0
GT_ALPHA  <- 0.65
GT_FILL   <- "#2E8B57"
GT_COLOR  <- "#1F5D3A"
HIST_BINS <- 7

# --- Libraries ---
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(ggsci)
  library(cowplot)
  library(tibble)
  library(scales)
  library(tidyr)
})

# --- Preconditions (expects the staging script already ran) ---
stopifnot(exists("MI_output"))
stopifnot(exists("GT_species"))
stopifnot("MI_permissive" %in% names(GT_species))
stopifnot(exists("decontam_MI_predict"))
stopifnot(exists("MI_matrix"))

if (!exists("squeegee_MI_predict")) {
  stopifnot(exists("Squeegee_output_folder"))
  squeegee_MI_predict <- read.table(
    file.path(Squeegee_output_folder, "squeegee_MI", "final_predictions.txt"),
    header = FALSE, sep = "\t", stringsAsFactors = FALSE
  )[[1]] %>% as.character()
}

# --- Paths (Metacontam outputs) ---
path_MI_Abundance  <- file.path(MI_output, "Abundance.txt")
path_MI_Prevalence <- file.path(MI_output, "Prevalence.txt")
path_MI_predict    <- file.path(MI_output, "Final_prediction.txt")

# ============================================================
# 1) Load predictions
# ============================================================

metacontam_pred <- read.table(
  path_MI_predict, sep = "\t", header = TRUE,
  check.names = FALSE, stringsAsFactors = FALSE
)
stopifnot(all(c("Taxid", "contamination_status") %in% names(metacontam_pred)))

MI_predict <- metacontam_pred %>%
  filter(contamination_status == "Contaminant") %>%
  pull(Taxid) %>%
  as.character()

# ============================================================
# 2) GT / true positives
# ============================================================

MI_gt_permissive <- as.character(GT_species$MI_permissive)

MI_permissive_true_positive          <- intersect(MI_predict,          MI_gt_permissive)
squeegee_MI_permissive_true_positive <- intersect(squeegee_MI_predict, MI_gt_permissive)
decontam_MI_permissive_true_positive <- intersect(decontam_MI_predict, MI_gt_permissive)

# ============================================================
# 3) Robust loaders
# ============================================================

read_metric_table <- function(path) {
  stopifnot(file.exists(path))
  df <- read.table(
    path, sep = "\t", header = TRUE,
    check.names = FALSE, stringsAsFactors = FALSE,
    quote = "", comment.char = "", fill = TRUE
  )
  
  nms <- names(df)
  cand <- c("taxid", "Taxid", "TaxID", "tax_id", "NCBI_taxid", "taxon_id")
  hit <- cand[cand %in% nms][1]
  
  if (!is.na(hit)) {
    names(df)[names(df) == hit] <- "taxid"
  } else {
    if (is.na(nms[1]) || nms[1] == "" ||
        all(grepl("^\\d+$", as.character(df[[1]]), perl = TRUE), na.rm = TRUE)) {
      names(df)[1] <- "taxid"
    } else {
      rn <- rownames(df)
      if (!is.null(rn) && any(rn != as.character(seq_len(nrow(df))))) {
        df <- tibble::rownames_to_column(df, "taxid")
      } else {
        stop("Cannot infer taxid column/rownames from: ", path)
      }
    }
  }
  
  df$taxid <- as.character(df$taxid)
  df
}

get_metric_df <- function(path, taxid_list, method_name, value_col) {
  df <- read_metric_table(path)
  stopifnot(value_col %in% names(df))
  
  df %>%
    filter(taxid %in% as.character(taxid_list)) %>%
    transmute(
      taxid,
      value  = as.numeric(.data[[value_col]]),
      Method = method_name
    ) %>%
    tidyr::drop_na(value)
}

get_abun_df_from_mat <- function(mat, taxid_list, method_name, as_percent = TRUE) {
  taxid_list <- as.character(taxid_list)
  present <- intersect(taxid_list, rownames(mat))
  
  if (length(present) == 0) {
    return(tibble::tibble())
  }
  
  vals <- rowMeans(mat[present, , drop = FALSE], na.rm = TRUE)
  if (as_percent) vals <- vals * 100
  
  tibble::tibble(
    taxid  = present,
    value  = as.numeric(vals),
    Method = method_name
  )
}

# ============================================================
# 4) Build tool data frames
# ============================================================

df_decontam_abun   <- get_abun_df_from_mat(MI_matrix, decontam_MI_permissive_true_positive, "Decontam",   as_percent = TRUE)
df_squeegee_abun   <- get_abun_df_from_mat(MI_matrix, squeegee_MI_permissive_true_positive, "Squeegee",   as_percent = TRUE)
df_metacontam_abun <- get_abun_df_from_mat(MI_matrix, MI_permissive_true_positive,          "Metacontam", as_percent = TRUE)

MI_abundance_df <- bind_rows(df_decontam_abun, df_squeegee_abun, df_metacontam_abun) %>%
  mutate(Method = factor(Method, levels = c("Decontam", "Squeegee", "Metacontam")))

df_metacontam_prev <- get_metric_df(path_MI_Prevalence, MI_permissive_true_positive,          "Metacontam", "preval")
df_squeegee_prev   <- get_metric_df(path_MI_Prevalence, squeegee_MI_permissive_true_positive, "Squeegee",   "preval")
df_decontam_prev   <- get_metric_df(path_MI_Prevalence, decontam_MI_permissive_true_positive, "Decontam",   "preval")

MI_prevalence_df <- bind_rows(df_decontam_prev, df_squeegee_prev, df_metacontam_prev) %>%
  mutate(Method = factor(Method, levels = c("Decontam", "Squeegee", "Metacontam")))

# ============================================================
# 5) Build GT data frames
# ============================================================

df_gt_abun <- get_abun_df_from_mat(
  MI_matrix,
  MI_gt_permissive,
  method_name = "GT",
  as_percent = TRUE
)

df_gt_prev <- get_metric_df(
  path_MI_Prevalence,
  MI_gt_permissive,
  method_name = "GT",
  value_col = "preval"
)

# ============================================================
# 6) Plot style
# ============================================================

nejm_cols <- pal_nejm("default")(8)
METHOD_COLORS <- c(
  "Decontam"   = nejm_cols[8],
  "Squeegee"   = nejm_cols[2],
  "Metacontam" = nejm_cols[3]
)

METHOD_SHAPES <- c(
  "Decontam"   = 16,
  "Squeegee"   = 15,
  "Metacontam" = 17
)

common_theme <- theme_pubr() +
  theme(
    plot.title         = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title.x       = element_text(size = 18),
    axis.title.y       = element_text(size = 18),
    axis.text.x        = element_text(size = 16),
    axis.text.y        = element_text(size = 16),
    legend.position    = "right",
    legend.title       = element_blank(),
    legend.text        = element_text(size = 14),
    legend.background  = element_blank(),
    legend.key         = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "black", fill = NA, linewidth = BORDER_LW),
    aspect.ratio       = 1,
    plot.margin        = margin(12, 12, 12, 12)
  )

# ============================================================
# 7) Helpers
# ============================================================

make_xlim <- function(x, pad_frac = 0.06) {
  x <- x[is.finite(x)]
  stopifnot(length(x) > 0)
  
  xr <- range(x, na.rm = TRUE)
  dx <- diff(xr)
  
  if (!is.finite(dx) || dx <= 0) {
    dx <- max(abs(xr[1]) * 0.1, 0.1)
  }
  
  c(xr[1] - dx * pad_frac, xr[2] + dx * pad_frac)
}

make_hist_breaks <- function(x, n_bins = 25) {
  x <- x[is.finite(x)]
  xr <- range(x, na.rm = TRUE)
  seq(xr[1], xr[2], length.out = n_bins + 1)
}

make_binned_count_df <- function(df, xvar, breaks) {
  df <- df %>% filter(is.finite(.data[[xvar]]), !is.na(Method))
  
  split_list <- split(df[[xvar]], df$Method)
  
  res <- lapply(names(split_list), function(m) {
    x <- split_list[[m]]
    h <- hist(x, breaks = breaks, plot = FALSE, include.lowest = TRUE, right = FALSE)
    data.frame(
      Method = m,
      x = h$mids,
      count = h$counts,
      stringsAsFactors = FALSE
    )
  })
  
  out <- bind_rows(res)
  out$Method <- factor(out$Method, levels = c("Decontam", "Squeegee", "Metacontam"))
  out
}

make_log10_power_breaks <- function(x, pad = 0) {
  x <- x[is.finite(x)]
  xr <- range(x, na.rm = TRUE)
  lo <- floor(xr[1]) - pad
  hi <- ceiling(xr[2]) + pad
  seq(lo, hi, by = 1)
}

make_log10_power_labels <- function(breaks) {
  parse(text = paste0("10^", breaks))
}

make_xlim_from_breaks <- function(breaks, expand_bins = 0.5) {
  breaks <- sort(unique(breaks))
  stopifnot(length(breaks) >= 2)
  
  d <- diff(breaks)
  
  if (length(unique(round(d, 10))) == 1) {
    step_left  <- d[1]
    step_right <- d[1]
  } else {
    step_left  <- d[1]
    step_right <- d[length(d)]
  }
  
  c(
    min(breaks) - expand_bins * step_left,
    max(breaks) + expand_bins * step_right
  )
}

overlay_hist_linecount_plot <- function(
    gt_df, tool_df, xvar, xlab_txt,
    x_limits = NULL, x_breaks = waiver(), x_labels = waiver(),
    hist_bins = 25, hist_breaks = NULL
) {
  gt_df   <- gt_df   %>% filter(is.finite(.data[[xvar]]))
  tool_df <- tool_df %>% filter(is.finite(.data[[xvar]]), !is.na(Method))
  
  all_x <- c(gt_df[[xvar]], tool_df[[xvar]])
  
  if (is.null(hist_breaks)) {
    hist_breaks <- make_hist_breaks(all_x, n_bins = hist_bins)
  }
  
  hist_breaks <- sort(unique(hist_breaks))
  stopifnot(length(hist_breaks) >= 2)
  
  line_df <- make_binned_count_df(tool_df, xvar = xvar, breaks = hist_breaks)
  
  ggplot() +
    geom_histogram(
      data = transform(gt_df, Group = "GT"),
      aes(x = .data[[xvar]], y = after_stat(count), fill = Group),
      breaks = hist_breaks,
      color = GT_COLOR,
      alpha = GT_ALPHA,
      linewidth = 0.35,
      inherit.aes = FALSE,
      show.legend = TRUE
    ) +
    geom_line(
      data = line_df,
      aes(x = x, y = count, color = Method, linetype = Method, group = Method),
      linewidth = CURVE_LW,
      lineend = "round",
      na.rm = TRUE,
      show.legend = TRUE,
      inherit.aes = FALSE
    ) +
    geom_point(
      data = line_df,
      aes(x = x, y = count, color = Method, shape = Method),
      size = PT_SIZE,
      stroke = PT_STROKE,
      na.rm = TRUE,
      show.legend = FALSE,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(
      values = c("GT" = GT_FILL),
      breaks = "GT"
    ) +
    scale_color_manual(values = METHOD_COLORS) +
    scale_shape_manual(values = METHOD_SHAPES) +
    scale_linetype_manual(
      values = c(
        "Decontam" = "solid",
        "Squeegee" = "solid",
        "Metacontam" = "solid"
      )
    ) +
    scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels,
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      breaks = scales::pretty_breaks(n = 5),
      labels = scales::label_number(accuracy = 1),
      expand = expansion(mult = c(0, 0.05))
    ) +
    coord_cartesian(xlim = x_limits, clip = "on") +
    labs(title = "", x = xlab_txt, y = "Count") +
    guides(
      fill = guide_legend(
        order = 1,
        override.aes = list(
          color = GT_COLOR,
          alpha = GT_ALPHA
        )
      ),
      linetype = "none",
      color = guide_legend(
        order = 2,
        override.aes = list(
          linewidth = CURVE_LW,
          shape = c(16, 15, 17),
          size = PT_SIZE,
          linetype = 1
        )
      ),
      shape = "none"
    ) +
    common_theme
}

# ============================================================
# 8) Abundance plot
# ============================================================

gt_abun_plot <- df_gt_abun %>%
  mutate(value = value / 100)

tool_abun_plot <- MI_abundance_df %>%
  mutate(value = value / 100)

all_abun_vals <- c(gt_abun_plot$value, tool_abun_plot$value)
min_pos_all_abun <- suppressWarnings(min(all_abun_vals[all_abun_vals > 0], na.rm = TRUE))
if (!is.finite(min_pos_all_abun)) min_pos_all_abun <- 1e-8
eps_abun <- min(min_pos_all_abun / 10, 1e-8)

gt_abun_plot <- gt_abun_plot %>%
  mutate(
    value_plot  = ifelse(value <= 0 | !is.finite(value), eps_abun, value),
    log10_value = log10(value_plot)
  )

tool_abun_plot <- tool_abun_plot %>%
  mutate(
    value_plot  = ifelse(value <= 0 | !is.finite(value), eps_abun, value),
    log10_value = log10(value_plot)
  )

abun_breaks <- make_log10_power_breaks(
  c(gt_abun_plot$log10_value, tool_abun_plot$log10_value)
)

abun_labels <- make_log10_power_labels(abun_breaks)

abun_xlim <- make_xlim_from_breaks(abun_breaks, expand_bins = 0.5)

p_MI_abundance_like_ANI <-
  overlay_hist_linecount_plot(
    gt_df = gt_abun_plot,
    tool_df = tool_abun_plot,
    xvar = "log10_value",
    xlab_txt = "Abundance",
    x_limits = abun_xlim,
    x_breaks = abun_breaks,
    x_labels = abun_labels,
    hist_breaks = abun_breaks
  )

# ============================================================
# 9) Prevalence plot
# ============================================================

gt_prev_plot <- df_gt_prev %>%
  filter(is.finite(value), !is.na(value))

tool_prev_plot <- MI_prevalence_df %>%
  filter(is.finite(value), !is.na(value), !is.na(Method))

prev_breaks <- pretty(c(gt_prev_plot$value, tool_prev_plot$value), n = 5)
prev_breaks <- sort(unique(prev_breaks))

prev_xlim <- make_xlim_from_breaks(prev_breaks, expand_bins = 0.5)

p_MI_prevalence_like_ANI <-
  overlay_hist_linecount_plot(
    gt_df = gt_prev_plot,
    tool_df = tool_prev_plot,
    xvar = "value",
    xlab_txt = "Prevalence",
    x_limits = prev_xlim,
    x_breaks = prev_breaks,
    x_labels = label_number(accuracy = 0.1)(prev_breaks),
    hist_breaks = prev_breaks
  )

# ============================================================
# 10) Side-by-side figure
# ============================================================

aligned <- cowplot::align_plots(
  p_MI_abundance_like_ANI,
  p_MI_prevalence_like_ANI,
  align = "hv", axis = "tblr"
)

p_methods_side_by_side <- cowplot::plot_grid(
  aligned[[1]], aligned[[2]],
  nrow = 1, rel_widths = c(1, 1), align = "h"
)

print(p_methods_side_by_side)

Figure2e_f <- p_methods_side_by_side


 
