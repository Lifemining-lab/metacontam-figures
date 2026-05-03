# =========================
# R script: PCoA with group ellipses, labels, and legend control
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(vegan)
  library(tibble)
  library(ape)
  library(scales)
  library(ggsci)
})

set.seed(42)

# ---------- Color palettes ----------
nejm_cols <- pal_nejm("default")(8)
aaas_cols <- pal_aaas("default")(8)
bmj_cols  <- pal_bmj("default")(8)

# ---------- Paths ----------
metadata_path <- "/media/junwoojo/18T/Submission/Rcode/inputs/Skin_negative/Skin_negative_metadata"
input_dir     <- "/media/junwoojo/18T/Submission/Rcode/inputs/Skin_negative/Skin_negative_bracken_files"
output_dir    <- "/home/junwoojo/18T/Submission/Rcode/inputs/Skin_negative"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- Analysis options ----------
group_col      <- "extraction_kit_id"
kits_main      <- c("standard", "microbiome")
permutations_n <- 999

min_rel_all   <- 1e-4
min_prev_all  <- 0.00
min_rel_main  <- 1e-4
min_prev_main <- 0.00

# ---------- Plot aesthetics ----------
base_size  <- 18
axis_txt   <- 16
axis_title <- 18
legend_txt <- 16
title_size <- 20

pal_main <- c(
  standard   = "#7570b3",
  microbiome = bmj_cols[2]
)

# =========================================================
# Helper functions
# =========================================================

fmt_p <- function(p) {
  ifelse(is.na(p), "p=NA", ifelse(p < 1e-3, "p<0.001", sprintf("p=%.3g", p)))
}

fmt_r2 <- function(r) {
  ifelse(is.na(r), "R²=NA", sprintf("R²=%.3f", r))
}

assert_cols <- function(df, cols, name = "data") {
  missing_cols <- setdiff(cols, colnames(df))
  if (length(missing_cols)) {
    stop(sprintf("[%s] missing columns: %s", name, paste(missing_cols, collapse = ", ")))
  }
}

read_kraken2_report <- function(file_path) {
  df <- suppressWarnings(read_tsv(file_path, col_names = FALSE, show_col_types = FALSE))
  
  if (ncol(df) < 5) {
    stop(sprintf("Invalid file format: %s", file_path))
  }
  
  df <- df[, 1:5]
  colnames(df) <- c("percentage", "reads_count", "reads_assigned", "rank", "name")
  
  species_df <- df %>%
    filter(rank == "S") %>%
    select(name, reads_count) %>%
    as.data.frame()
  
  rownames(species_df) <- species_df$name
  species_df$name <- NULL
  
  species_df
}

merge_species_counts <- function(files) {
  merged <- NULL
  sample_names <- character(0)
  
  for (file_path in files) {
    sample_id <- str_replace(basename(file_path), "_kraken_report_bracken_species$", "")
    sample_names <- c(sample_names, sample_id)
    
    species_df <- read_kraken2_report(file_path)
    colnames(species_df) <- sample_id
    
    merged <- if (is.null(merged)) {
      species_df
    } else {
      merged %>%
        rownames_to_column("name") %>%
        full_join(species_df %>% rownames_to_column("name"), by = "name") %>%
        column_to_rownames("name")
    }
  }
  
  if (any(duplicated(sample_names))) {
    stop("Duplicated sample names were parsed from filenames.")
  }
  
  merged[is.na(merged)] <- 0
  merged
}

tss <- function(mat) {
  col_sums <- colSums(mat)
  keep <- col_sums > 0
  
  if (!all(keep)) {
    warning(sprintf("Dropping %d zero-sum samples.", sum(!keep)))
    mat <- mat[, keep, drop = FALSE]
    col_sums <- col_sums[keep]
  }
  
  sweep(mat, 2, col_sums, "/")
}

prevalence_filter <- function(mat_rel, min_prev = 0, min_rel = 1e-4) {
  if (min_prev <= 0) {
    return(mat_rel)
  }
  
  keep <- rowMeans(mat_rel >= min_rel) >= min_prev
  kept_n <- sum(keep)
  
  message(sprintf(
    "[prevalence] kept %d / %d taxa (min_prev=%.3f, min_rel=%g)",
    kept_n, nrow(mat_rel), min_prev, min_rel
  ))
  
  if (kept_n == 0) {
    stop("No taxa passed the prevalence thresholds. Please relax the filters.")
  }
  
  mat_rel[keep, , drop = FALSE]
}

bray_dist <- function(mat_rel) {
  vegan::vegdist(t(mat_rel), method = "bray")
}

get_pcoa_coords <- function(dist_mat) {
  pcoa_res <- ape::pcoa(dist_mat)
  
  coords <- data.frame(
    PC1 = pcoa_res$vectors[, 1],
    PC2 = pcoa_res$vectors[, 2],
    row.names = rownames(pcoa_res$vectors)
  )
  
  var_exp <- 100 * pcoa_res$values$Relative_eig[1:2]
  
  list(df = coords, var_exp = var_exp)
}

run_permanova <- function(dist_mat, labels, perms = 999) {
  group_factor <- factor(labels)
  
  if (nlevels(group_factor) < 2) {
    warning("PERMANOVA skipped: factor has fewer than 2 levels.")
    return(list(R2 = NA_real_, p = NA_real_))
  }
  
  meta_df <- data.frame(group = group_factor)
  adonis_res <- vegan::adonis2(dist_mat ~ group, data = meta_df, permutations = perms)
  
  list(
    R2 = adonis_res$R2[1],
    p  = adonis_res$`Pr(>F)`[1]
  )
}

normalize_kit <- function(x) {
  lx <- tolower(x)
  
  dplyr::case_when(
    # Explicitly exclude these labels
    stringr::str_detect(lx, "^homebrew$|^microbiome beta$") ~ "other",
    
    # microbiome: exact label or MagMAX Microbiome Ultra keywords
    stringr::str_detect(lx, "^microbiome$|magmax|microbiome ultra|a42357") ~ "microbiome",
    
    # pro: exact label or PowerSoil Pro keywords
    stringr::str_detect(lx, "^pro$|power\\s*soil\\s*pro|powersoil\\s*pro|47109") ~ "pro",
    
    # standard: exact label or general PowerSoil / standardized keywords
    stringr::str_detect(
      lx,
      "^standard$|power\\s*soil(?!\\s*pro)|powersoil(?!\\s*pro)|27000-4-kf|standardized"
    ) ~ "standard",
    
    # norgen
    stringr::str_detect(lx, "^norgen$|norgen|65600") ~ "norgen",
    
    # nucleomag
    stringr::str_detect(lx, "^nucleomag$|nucleomag|744945") ~ "nucleomag",
    
    # zymo
    stringr::str_detect(lx, "^zymo$|zymo|magbead|d4302") ~ "zymo",
    
    TRUE ~ "other"
  )
}

theme_clean_box <- function(show_legend = FALSE) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title      = element_text(size = title_size, face = "bold"),
      axis.title      = element_text(size = axis_title),
      axis.text       = element_text(size = axis_txt),
      legend.position = if (show_legend) "right" else "none",
      legend.text     = element_text(size = legend_txt),
      legend.title    = element_text(size = legend_txt),
      panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      panel.grid      = element_blank(),
      plot.margin     = margin(12, 12, 12, 12)
    )
}

# =========================================================
# Load input data
# =========================================================

metadata <- suppressWarnings(read_tsv(metadata_path, col_types = cols(.default = "c")))
assert_cols(metadata, c("Run", group_col), "metadata")
metadata <- metadata %>% distinct(Run, .keep_all = TRUE)

files <- list.files(
  input_dir,
  pattern = "_kraken_report_bracken_species$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No *_kraken_report_bracken_species files were found.")
}

counts <- merge_species_counts(files)

# Keep only samples present in metadata
counts <- counts[, colnames(counts) %in% metadata$Run, drop = FALSE]
if (ncol(counts) < 2) {
  stop("Not enough samples remain after matching metadata$Run.")
}

# Total sum scaling
rel_all <- tss(counts)

labs_all <- metadata[[group_col]][match(colnames(rel_all), metadata$Run)]
names(labs_all) <- colnames(rel_all)

# =========================================================
# Main analysis: standard vs microbiome
# PCoA + group ellipses + centroid labels
# =========================================================

keep_main <- labs_all %in% kits_main
rel_main  <- rel_all[, keep_main, drop = FALSE]

labs_main <- droplevels(factor(labs_all[keep_main], levels = kits_main))
names(labs_main) <- colnames(rel_main)

rel_main_f <- prevalence_filter(
  rel_main,
  min_prev = min_prev_main,
  min_rel  = min_rel_main
)

D_main <- bray_dist(rel_main_f)
per_m  <- run_permanova(D_main, labs_main, perms = permutations_n)
pco_m  <- get_pcoa_coords(D_main)

df_main <- pco_m$df %>%
  tibble::rownames_to_column("Sample") %>%
  mutate(Group = labs_main[Sample])

vx_main <- pco_m$var_exp
xlab_main <- sprintf("PCoA1 (%.1f%%)", vx_main[1])
ylab_main <- sprintf("PCoA2 (%.1f%%)", vx_main[2])

stat_label_main <- paste(fmt_p(per_m$p), fmt_r2(per_m$R2), sep = ", ")
cat(sprintf("[MAIN PERMANOVA] %s\n", stat_label_main))

centroids_main <- df_main %>%
  group_by(Group) %>%
  summarise(
    PC1 = mean(PC1),
    PC2 = mean(PC2),
    .groups = "drop"
  )

p_main <- ggplot(df_main, aes(PC1, PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.95) +
  stat_ellipse(
    aes(group = Group, color = Group),
    type = "t",
    level = 0.80,
    linetype = "solid",
    linewidth = 1.0,
    show.legend = FALSE
  ) +
  geom_text(
    data = centroids_main,
    aes(PC1, PC2, label = Group, color = Group),
    fontface = "bold",
    size = 6.5,
    show.legend = FALSE
  ) +
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = stat_label_main,
    hjust = 1.05,
    vjust = -0.8,
    size = 5.5
  ) +
  scale_color_manual(values = pal_main, drop = FALSE) +
  labs(
    title = "PCoA: standard vs microbiome",
    x = xlab_main,
    y = ylab_main
  ) +
  theme_clean_box(show_legend = FALSE) +
  coord_cartesian(clip = "off")

print(p_main)

# =========================================================
# Supplementary analysis: selected 6 kits only
# PCoA with legend showing one-vs-other PERMANOVA stats + n
# Keep: microbiome, standard, pro, norgen, nucleomag, zymo
# Exclude: homebrew, microbiome beta, and all other labels
# =========================================================

rel_all_f <- prevalence_filter(
  rel_all,
  min_prev = min_prev_all,
  min_rel  = min_rel_all
)

labs_simpl_all <- normalize_kit(labs_all)
names(labs_simpl_all) <- names(labs_all)

keep_kits <- c("microbiome", "standard", "pro", "norgen", "nucleomag", "zymo")
keep_idx  <- labs_simpl_all %in% keep_kits

if (!any(keep_idx)) {
  stop("No samples matched the selected six kits.")
}

# Report how many samples were excluded by kit
excluded_tab <- sort(table(labs_simpl_all[!keep_idx]), decreasing = TRUE)
if (length(excluded_tab)) {
  message(sprintf(
    "[supplementary] excluded %d samples | by kit: %s",
    sum(excluded_tab),
    paste(names(excluded_tab), as.integer(excluded_tab), sep = "=", collapse = ", ")
  ))
}

rel_six <- rel_all_f[, keep_idx, drop = FALSE]

labs_six <- factor(labs_simpl_all[keep_idx], levels = keep_kits)
labs_six <- droplevels(labs_six)
names(labs_six) <- colnames(rel_six)

D_six   <- bray_dist(rel_six)
pco_six <- get_pcoa_coords(D_six)

df_six <- pco_six$df %>%
  tibble::rownames_to_column("Sample") %>%
  mutate(Kit = labs_six[Sample])

df_six$Kit <- droplevels(df_six$Kit)

vx_six <- pco_six$var_exp
xlab_six <- sprintf("PCoA1 (%.1f%%)", vx_six[1])
ylab_six <- sprintf("PCoA2 (%.1f%%)", vx_six[2])

# Sample counts for the legend
n_by_kit <- table(df_six$Kit)

# One-vs-other PERMANOVA for each kit
kits_present <- levels(df_six$Kit)

ova <- lapply(kits_present, function(k) {
  labels_k <- ifelse(labs_six == k, k, "Other")
  per <- run_permanova(D_six, labels_k, perms = permutations_n)
  
  data.frame(
    kit = k,
    p = per$p,
    R2 = per$R2,
    stringsAsFactors = FALSE
  )
}) %>%
  dplyr::bind_rows()

legend_labels <- sapply(ova$kit, function(k) {
  n_here <- unname(as.integer(n_by_kit[k]))
  paste0(
    k,
    " (n=", n_here,
    "; ", fmt_p(ova$p[ova$kit == k]),
    ", ", fmt_r2(ova$R2[ova$kit == k]), ")"
  )
}, USE.NAMES = FALSE)

label_map <- setNames(legend_labels, ova$kit)

p_supp <- ggplot(df_six, aes(PC1, PC2, color = Kit)) +
  geom_point(size = 6.0, alpha = 0.95) +
  scale_color_discrete(
    name = "Kit",
    labels = label_map[levels(df_six$Kit)],
    drop = FALSE
  ) +
  labs(
    title = "PCoA: Selected 6 kits (legend = one-vs-other stats + n)",
    x = xlab_six,
    y = ylab_six
  ) +
  theme_clean_box(show_legend = TRUE)

print(p_supp)

cat(
  "\nDone. PCoA with Bray-Curtis distance. Main figure: legend off. ",
  "Supplementary figure: only 6 selected kits included, with legend showing n and PERMANOVA statistics (p, R²).\n",
  sep = ""
)
Figure3a <- p_main
Supp_Fig4 <- p_supp
