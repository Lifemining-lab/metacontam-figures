suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(purrr)
  library(ggplot2)
  library(ggsci); library(patchwork)
})

# ============================================================
# Colors
# ============================================================
nejm_cols     <- pal_nejm("default")(8)
METHOD_COLORS <- c(
  "Original"   = "forestgreen",
  "Decontam"   = nejm_cols[8],
  "Squeegee"   = nejm_cols[2],
  "Metacontam" = nejm_cols[3]
)

if (!exists("squeegee_skin_standard_predict") && exists("squegee_skin_standard_predict")) {
  squeegee_skin_standard_predict <- squegee_skin_standard_predict
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ============================================================
# Utils
# ============================================================

# 입력 mat을 항상 relative abundance로 맞춤
to_rel <- function(mat) {
  cs <- colSums(mat, na.rm = TRUE)
  rel <- sweep(mat, 2, ifelse(cs == 0, 1, cs), "/")
  rel[is.na(rel)] <- 0
  rel
}

# contaminant 제거 후 alpha 계산 전에 다시 재정규화
renorm_rel <- function(mat) {
  to_rel(mat)
}

remove_contam_with_cutoff <- function(mat, contam_taxids_char, cutoff = 1.00) {
  if (length(contam_taxids_char) == 0) return(mat)
  
  rownames(mat) <- as.character(rownames(mat))
  contam_taxids_char <- unique(as.character(contam_taxids_char))
  hit <- intersect(rownames(mat), contam_taxids_char)
  if (length(hit) == 0) return(mat)
  
  rel <- to_rel(mat)
  
  mask <- matrix(FALSE, nrow = nrow(mat), ncol = ncol(mat),
                 dimnames = dimnames(mat))
  mask[hit, ] <- rel[hit, ] < cutoff
  
  out <- mat
  out[mask] <- 0
  out
}

alpha_from_mat <- function(mat, renormalize = TRUE) {
  rel <- if (renormalize) renorm_rel(mat) else to_rel(mat)
  
  apply(rel, 2, function(p) {
    p <- p[p > 0]
    if (length(p) == 0) return(c(shannon = 0, simpson = 0))
    sh <- -sum(p * log(p))
    si <- 1 - sum(p^2)
    c(shannon = sh, simpson = si)
  }) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("sample")
}

build_alpha_long <- function(mat,
                             contam_squeegee, contam_metacontam, contam_decontam,
                             cutoff = 1.00,
                             group_map = NULL,
                             dataset_label = NULL,
                             renormalize = TRUE) {
  stopifnot(!is.null(colnames(mat)))
  rownames(mat) <- as.character(rownames(mat))
  
  # 처음 input을 relative abundance로 한 번 맞춤
  mat <- to_rel(mat)
  
  contam_squeegee   <- unique(as.character(contam_squeegee))
  contam_metacontam <- unique(as.character(contam_metacontam))
  contam_decontam   <- unique(as.character(contam_decontam))
  
  mat_orig  <- mat
  mat_sq    <- remove_contam_with_cutoff(mat, contam_squeegee,   cutoff = cutoff)
  mat_meta  <- remove_contam_with_cutoff(mat, contam_metacontam, cutoff = cutoff)
  mat_decon <- remove_contam_with_cutoff(mat, contam_decontam,   cutoff = cutoff)
  
  d0 <- alpha_from_mat(mat_orig,  renormalize = renormalize) %>% mutate(method = "Original")
  d1 <- alpha_from_mat(mat_decon, renormalize = renormalize) %>% mutate(method = "Decontam")
  d2 <- alpha_from_mat(mat_sq,    renormalize = renormalize) %>% mutate(method = "Squeegee")
  d3 <- alpha_from_mat(mat_meta,  renormalize = renormalize) %>% mutate(method = "Metacontam")
  
  df <- bind_rows(d0, d1, d2, d3) %>%
    mutate(method = factor(method, levels = c("Original","Decontam","Squeegee","Metacontam")))
  
  if (is.null(group_map)) {
    df <- df %>% mutate(group = dataset_label %||% "All")
  } else {
    gm <- tibble(sample = names(group_map), group = unname(group_map))
    df <- df %>% left_join(gm, by = "sample") %>% filter(!is.na(group))
  }
  
  df %>% mutate(group = as.character(group))
}

# ============================================================
# Wilcoxon p-values (unpaired)
# ============================================================
fmt_p <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 1e-4) sprintf("p = %.2e", p) else sprintf("p = %.3f", p)
}

two_sided_wilcox_p <- function(df_long, metric = c("shannon","simpson")) {
  metric <- match.arg(metric)
  df_long <- df_long %>% mutate(group = as.character(group))
  
  p_to_stars <- function(p) {
    if (is.na(p)) return(NA_character_)
    as.character(stats::symnum(
      p, corr = FALSE, na = FALSE,
      cutpoints = c(0, .001, .01, .05, 1),
      symbols   = c("***", "**", "*", "n.s.")
    ))
  }
  
  df_long %>%
    group_by(group) %>%
    group_split() %>%
    purrr::map_dfr(function(g) {
      grp <- unique(g$group)
      n_samp <- dplyr::n_distinct(g$sample)
      
      wide <- g %>%
        select(sample, method, !!rlang::sym(metric)) %>%
        tidyr::pivot_wider(names_from = method, values_from = !!rlang::sym(metric))
      
      x <- wide$Original
      y <- wide$Squeegee
      z <- wide$Metacontam
      d <- wide$Decontam
      
      ok_OS <- is.finite(x) & is.finite(y)
      ok_OM <- is.finite(x) & is.finite(z)
      ok_OD <- is.finite(x) & is.finite(d)
      
      p_OS <- if (any(ok_OS)) stats::wilcox.test(
        x[ok_OS], y[ok_OS], paired = FALSE, exact = FALSE
      )$p.value else NA_real_
      
      p_OM <- if (any(ok_OM)) stats::wilcox.test(
        x[ok_OM], z[ok_OM], paired = FALSE, exact = FALSE
      )$p.value else NA_real_
      
      p_OD <- if (any(ok_OD)) stats::wilcox.test(
        x[ok_OD], d[ok_OD], paired = FALSE, exact = FALSE
      )$p.value else NA_real_
      
      tibble(
        group = grp, n = n_samp,
        p_Squeegee   = p_OS,
        p_Metacontam = p_OM,
        p_Decontam   = p_OD,
        star_S = p_to_stars(p_OS),
        star_M = p_to_stars(p_OM),
        star_D = p_to_stars(p_OD),
        txtp_S = fmt_p(p_OS),
        txtp_M = fmt_p(p_OM),
        txtp_D = fmt_p(p_OD)
      )
    })
}

# ============================================================
# Brackets
# ============================================================
.dodge_w <- 0.8
.offset <- function(k, n, dodge_w = .dodge_w) (k - (n + 1) / 2) * (dodge_w / n)

make_onepanel_brackets <- function(df_long,
                                   metric = c("shannon","simpson"),
                                   show_ns = TRUE,
                                   pval_shift_mult = 3.1,
                                   star_shift_mult = 0.1,
                                   ns_shift_mult   = 1.05,
                                   alpha_sig = 0.05) {
  metric <- match.arg(metric)
  df_long <- df_long %>% mutate(group = as.character(group))
  
  ptab <- two_sided_wilcox_p(df_long, metric = metric)
  
  g_levels <- df_long %>% distinct(group) %>% pull(group) %>% as.character()
  g_index  <- setNames(seq_along(g_levels), g_levels)
  
  method_levels <- levels(df_long$method)
  n_methods     <- length(method_levels)
  method_index  <- setNames(seq_along(method_levels), method_levels)
  
  y_tab <- df_long %>% group_by(group) %>%
    summarise(ymax = max(.data[[metric]], na.rm = TRUE), .groups = "drop")
  
  segs <- list()
  star_labs <- list()
  ns_labs <- list()
  pval_labs <- list()
  
  global_max <- suppressWarnings(max(df_long[[metric]], na.rm = TRUE))
  tip <- if (is.finite(global_max) && global_max > 0) global_max * 0.015 else 0.02
  
  # 작은 값 그룹에서도 최소 간격 보장
  min_gap <- if (is.finite(global_max) && global_max > 0) global_max * 0.11 else 0.1
  
  for (g in g_levels) {
    yi <- y_tab$ymax[y_tab$group == g]
    if (!is.finite(yi)) yi <- 0
    
    y1 <- yi + max(yi * 0.12, min_gap * 1.0)
    y2 <- yi + max(yi * 0.30, min_gap * 2.0)
    y3 <- yi + max(yi * 0.48, min_gap * 3.0)
    
    xg <- g_index[[g]]
    xO <- xg + .offset(method_index["Original"],   n_methods)
    xD <- xg + .offset(method_index["Decontam"],   n_methods)
    xS <- xg + .offset(method_index["Squeegee"],   n_methods)
    xM <- xg + .offset(method_index["Metacontam"], n_methods)
    
    row <- ptab %>% filter(group == g)
    
    star_S <- row$star_S; star_M <- row$star_M; star_D <- row$star_D
    txtp_S <- row$txtp_S; txtp_M <- row$txtp_M; txtp_D <- row$txtp_D
    p_S    <- row$p_Squeegee
    p_M    <- row$p_Metacontam
    p_D    <- row$p_Decontam
    
    draw_S <- length(star_S)==1 && !is.na(star_S) && (show_ns || star_S != "n.s.")
    draw_M <- length(star_M)==1 && !is.na(star_M) && (show_ns || star_M != "n.s.")
    draw_D <- length(star_D)==1 && !is.na(star_D) && (show_ns || star_D != "n.s.")
    
    # Original vs Squeegee
    if (draw_S) {
      segs[[length(segs)+1]] <- tibble(x=xO, xend=xS, y=y1, yend=y1, group=g)
      segs[[length(segs)+1]] <- tibble(x=xO, xend=xO, y=y1, yend=y1-tip, group=g)
      segs[[length(segs)+1]] <- tibble(x=xS, xend=xS, y=y1, yend=y1-tip, group=g)
      
      if (star_S == "n.s.") {
        ns_labs[[length(ns_labs)+1]] <- tibble(
          x = (xO+xS)/2, y = y1 + tip*ns_shift_mult, label = star_S, group = g
        )
      } else {
        star_labs[[length(star_labs)+1]] <- tibble(
          x = (xO+xS)/2, y = y1 + tip*star_shift_mult, label = star_S, group = g
        )
      }
      
      if (length(p_S)==1 && is.finite(p_S) && p_S < alpha_sig) {
        pval_labs[[length(pval_labs)+1]] <- tibble(
          x=(xO+xS)/2, y=y1+tip*pval_shift_mult, label=txtp_S, group=g
        )
      }
    }
    
    # Original vs Metacontam
    if (draw_M) {
      segs[[length(segs)+1]] <- tibble(x=xO, xend=xM, y=y2, yend=y2, group=g)
      segs[[length(segs)+1]] <- tibble(x=xO, xend=xO, y=y2, yend=y2-tip, group=g)
      segs[[length(segs)+1]] <- tibble(x=xM, xend=xM, y=y2, yend=y2-tip, group=g)
      
      if (star_M == "n.s.") {
        ns_labs[[length(ns_labs)+1]] <- tibble(
          x = (xO+xM)/2, y = y2 + tip*ns_shift_mult, label = star_M, group = g
        )
      } else {
        star_labs[[length(star_labs)+1]] <- tibble(
          x = (xO+xM)/2, y = y2 + tip*star_shift_mult, label = star_M, group = g
        )
      }
      
      if (length(p_M)==1 && is.finite(p_M) && p_M < alpha_sig) {
        pval_labs[[length(pval_labs)+1]] <- tibble(
          x=(xO+xM)/2, y=y2+tip*pval_shift_mult, label=txtp_M, group=g
        )
      }
    }
    
    # Original vs Decontam
    if (draw_D) {
      segs[[length(segs)+1]] <- tibble(x=xO, xend=xD, y=y3, yend=y3, group=g)
      segs[[length(segs)+1]] <- tibble(x=xO, xend=xO, y=y3, yend=y3-tip, group=g)
      segs[[length(segs)+1]] <- tibble(x=xD, xend=xD, y=y3, yend=y3-tip, group=g)
      
      if (star_D == "n.s.") {
        ns_labs[[length(ns_labs)+1]] <- tibble(
          x = (xO+xD)/2, y = y3 + tip*ns_shift_mult, label = star_D, group = g
        )
      } else {
        star_labs[[length(star_labs)+1]] <- tibble(
          x = (xO+xD)/2, y = y3 + tip*star_shift_mult, label = star_D, group = g
        )
      }
      
      if (length(p_D)==1 && is.finite(p_D) && p_D < alpha_sig) {
        pval_labs[[length(pval_labs)+1]] <- tibble(
          x=(xO+xD)/2, y=y3+tip*pval_shift_mult, label=txtp_D, group=g
        )
      }
    }
  }
  
  seg_df  <- if (length(segs)) bind_rows(segs) else tibble(x=numeric(0), xend=numeric(0), y=numeric(0), yend=numeric(0), group=character(0))
  star_df <- if (length(star_labs)) bind_rows(star_labs) else tibble(x=numeric(0), y=numeric(0), label=character(0), group=character(0))
  ns_df   <- if (length(ns_labs)) bind_rows(ns_labs) else tibble(x=numeric(0), y=numeric(0), label=character(0), group=character(0))
  pval_df <- if (length(pval_labs)) bind_rows(pval_labs) else tibble(x=numeric(0), y=numeric(0), label=character(0), group=character(0))
  
  list(seg_df=seg_df, star_df=star_df, ns_df=ns_df, pval_df=pval_df, g_levels=g_levels)
}

# ============================================================
# Plot
# ============================================================
plot_alpha_onepanel_stars <- function(df_long, metric = c("shannon","simpson"),
                                      show_ns = TRUE, title = NULL) {
  metric <- match.arg(metric)
  df_long <- df_long %>% mutate(group = as.character(group))
  
  ann <- make_onepanel_brackets(df_long, metric = metric, show_ns = show_ns)
  
  nlab <- tibble(
    group  = as.character(df_long$group),
    sample = as.character(df_long$sample)
  ) %>%
    distinct(group, sample) %>%
    dplyr::count(group, name = "n")
  
  xlabs <- setNames(paste0(nlab$group, " (n=", nlab$n, ")"), nlab$group)
  
  ymax_data <- suppressWarnings(max(df_long[[metric]], na.rm = TRUE))
  ymax_ann  <- suppressWarnings(max(bind_rows(ann$star_df, ann$ns_df, ann$pval_df)$y, na.rm = TRUE))
  ylim_top  <- max(ymax_data * 2.10, ymax_ann * 1.12, na.rm = TRUE)
  
  ggplot(df_long, aes(x = factor(group, levels = ann$g_levels),
                      y = .data[[metric]], fill = method)) +
    geom_boxplot(outlier.shape = 16,
                 width = 0.6,
                 position = position_dodge(width = .dodge_w)) +
    scale_fill_manual(values = METHOD_COLORS,
                      breaks = c("Original","Decontam","Squeegee","Metacontam")) +
    scale_x_discrete(labels = xlabs) +
    labs(
      x = NULL,
      y = ifelse(metric == "shannon", "Shannon index", "Simpson index"),
      fill = NULL,
      title = title
    ) +
    geom_segment(data = ann$seg_df,
                 aes(x = x, xend = xend, y = y, yend = yend),
                 inherit.aes = FALSE, linewidth = 0.4) +
    geom_text(data = ann$star_df,
              aes(x = x, y = y, label = label),
              inherit.aes = FALSE, size = 5, vjust = 0) +
    geom_text(data = ann$ns_df,
              aes(x = x, y = y, label = label),
              inherit.aes = FALSE, size = 4.5, vjust = 0) +
    geom_text(data = ann$pval_df,
              aes(x = x, y = y, label = label),
              inherit.aes = FALSE, size = 3.8, vjust = 0) +
    coord_cartesian(ylim = c(NA, ylim_top), clip = "off") +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "top",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# ============================================================
# MI group mapping
# ============================================================
path_MI_metadata <- "/media/junwoojo/18T/Metacontam_dataset/Metadata/MI_metadata"

make_MI_group_map <- function(path = path_MI_metadata, mi_samples = colnames(MI_matrix)) {
  meta <- read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                     col.names = c("sample","subtype","R1","R2")) %>%
    select(sample, subtype)
  
  agg <- meta %>%
    mutate(group = case_when(
      grepl("^placenta_", subtype, ignore.case = TRUE)                ~ "Placenta",
      grepl("^humanmilk_BreastMilk", subtype, ignore.case = TRUE)     ~ "BreastMilk",
      grepl("^humanoral_Gingiva_Buccal", subtype, ignore.case = TRUE) ~ "Oral",
      grepl("^humanvaginal_", subtype, ignore.case = TRUE)            ~ "Vaginal",
      grepl("^humangut_(Stool|Rectum|Meconium)_", subtype, TRUE)      ~ "Fecal",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(group)) %>%
    semi_join(tibble(sample = mi_samples), by = "sample")
  
  setNames(agg$group, agg$sample)
}

# ============================================================
# p-value print utilities
# ============================================================
pval_table <- function(df_long, metric = c("shannon","simpson"),
                       alpha_sig = 0.05, show_only_sig = TRUE) {
  metric <- match.arg(metric)
  
  tab <- two_sided_wilcox_p(df_long, metric = metric) %>%
    transmute(
      group, n,
      `Orig vs Decontam (p)`   = p_Decontam,
      `Orig vs Squeegee (p)`   = p_Squeegee,
      `Orig vs Metacontam (p)` = p_Metacontam
    )
  
  if (show_only_sig) {
    tab <- tab %>%
      filter(
        (is.finite(`Orig vs Decontam (p)`)   & `Orig vs Decontam (p)`   < alpha_sig) |
          (is.finite(`Orig vs Squeegee (p)`)   & `Orig vs Squeegee (p)`   < alpha_sig) |
          (is.finite(`Orig vs Metacontam (p)`) & `Orig vs Metacontam (p)` < alpha_sig)
      )
  }
  
  fmt_num <- function(p) ifelse(is.na(p), NA_character_,
                                ifelse(p < 1e-4, formatC(p, format = "e", digits = 2),
                                       formatC(p, format = "f", digits = 3)))
  
  tab %>%
    mutate(
      `Orig vs Decontam (p)`   = fmt_num(`Orig vs Decontam (p)`),
      `Orig vs Squeegee (p)`   = fmt_num(`Orig vs Squeegee (p)`),
      `Orig vs Metacontam (p)` = fmt_num(`Orig vs Metacontam (p)`)
    ) %>%
    arrange(group)
}

print_pvals <- function(df_long, metric, cutoff_label = "cutoff 100%",
                        alpha_sig = 0.05, show_only_sig = TRUE) {
  cat("\n[", ifelse(metric == "shannon", "Shannon", "Simpson"),
      ", ", cutoff_label, "] p-values\n", sep = "")
  tt <- pval_table(df_long, metric = metric, alpha_sig = alpha_sig, show_only_sig = show_only_sig)
  if (nrow(tt) == 0) cat("  - No significant results.\n") else print(tt, n = Inf)
}

# ============================================================
# RUN: cutoff = 100% only
# relative abundance input + renormalize after removal + unpaired Wilcoxon
# ============================================================
MI_group_map <- make_MI_group_map()

mi_long_100p <- build_alpha_long(
  mat               = MI_matrix,
  contam_squeegee   = squeegee_MI_predict,
  contam_metacontam = MI_predict,
  contam_decontam   = decontam_MI_predict,
  cutoff            = 1.00,
  group_map         = MI_group_map,
  renormalize       = TRUE
)

skin_micro_long_100p <- build_alpha_long(
  mat               = microbiome_matrix,
  contam_squeegee   = squeegee_skin_microbiome_predict,
  contam_metacontam = Skin_microbiome_predict,
  contam_decontam   = decontam_skin_microbiome_predict,
  cutoff            = 1.00,
  dataset_label     = "Skin (microbiome)",
  renormalize       = TRUE
)

skin_std_long_100p <- build_alpha_long(
  mat               = standard_matrix,
  contam_squeegee   = squeegee_skin_standard_predict,
  contam_metacontam = Skin_standard_predict,
  contam_decontam   = decontam_skin_standard_predict,
  cutoff            = 1.00,
  dataset_label     = "Skin (standard)",
  renormalize       = TRUE
)

nasal_long_100p <- build_alpha_long(
  mat               = Nasal_matrix,
  contam_squeegee   = squeegee_nasal_predict,
  contam_metacontam = Nasal_predict,
  contam_decontam   = decontam_nasal_predict,
  cutoff            = 1.00,
  dataset_label     = "Nasal",
  renormalize       = TRUE
)

# ============================================================
# Combine all datasets
# ============================================================
mi_prefix <- function(df) df %>% mutate(group = paste0("MI: ", group))

mi_levels <- mi_long_100p %>% distinct(group) %>% pull(group) %>% as.character()
all_levels <- c(paste0("MI: ", mi_levels), "Skin (microbiome)", "Skin (standard)", "Nasal")

relevel_groups <- function(df, levels) {
  df %>% mutate(group = factor(group, levels = levels))
}

all_long_100p <- bind_rows(
  mi_prefix(mi_long_100p),
  skin_micro_long_100p,
  skin_std_long_100p,
  nasal_long_100p
) %>% relevel_groups(all_levels)

# ============================================================
# Plot
# ============================================================
p_all_shan_100p <- plot_alpha_onepanel_stars(
  all_long_100p, metric = "shannon", show_ns = TRUE,
  title = "All datasets (cutoff 100%): Shannon"
)

p_all_simp_100p <- plot_alpha_onepanel_stars(
  all_long_100p, metric = "simpson", show_ns = TRUE,
  title = "All datasets (cutoff 100%): Simpson"
)


# ============================================================
# Output
# ============================================================
print(p_all_shan_100p)
print(p_all_simp_100p)
print_pvals(all_long_100p, metric = "shannon", cutoff_label = "cutoff 100%")
print_pvals(all_long_100p, metric = "simpson", cutoff_label = "cutoff 100%")


Supp_Fig3 <- p_all_shan_100p+p_all_simp_100p

