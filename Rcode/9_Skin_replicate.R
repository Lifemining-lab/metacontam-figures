# ============================================================
# Skin WGS replicate analysis
# - Jaccard distance across filtering methods
# - Decontam taxids loaded from files
# - Plots kept as close as possible to original
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggsci)
  library(ggpubr)     # ggboxplot, stat_pvalue_manual
  library(scales)
  library(grid)
  library(openxlsx)
})

# ------------------------------------------------------------
# User paths (edit if needed)
# ------------------------------------------------------------

decon_std_path <- file.path(BASE_INPUT, "Decontam_output_folder", "skin_standard_contaminant_taxids.txt")

decon_mic_path <- file.path(BASE_INPUT, "Decontam_output_folder", "skin_microbiome_contaminant_taxids.txt")

metadata_path  <- file.path(BASE_INPUT, "Skin_sra_metadata")

# ------------------------------------------------------------
# Global palette & ordering
# ------------------------------------------------------------
nejm_cols <- pal_nejm("default")(8)
aaas_cols <- pal_aaas("default")(8)
bmj_cols  <- pal_bmj("default")(8)

METHOD_LEVELS <- c("Original","Decontam","Squeegee","Metacontam")
METHOD_COLORS <- c(
  "Original"   = "forestgreen",
  "Decontam"   = nejm_cols[8],
  "Squeegee"   = nejm_cols[2],
  "Metacontam" = nejm_cols[3]
)

# Spaghetti plot colors
SLOPE_COLORS <- c(
  "Increased" = aaas_cols[2],
  "Decreased" = bmj_cols[7],
  "Unchanged" = "#636363"
)

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

# Canonicalize sample/run IDs (keep consistent with your original behavior)
canon <- function(x) sub("_.*$|\\..*$", "", as.character(x))

# Clean taxid lists from vectors (remove empty, comments, headers)
cleanup_lines_to_taxids <- function(v) {
  v <- as.character(v)
  v <- v[!grepl("^\\s*$", v)]
  v <- v[!grepl("^\\s*#", v)]
  v <- sub("\\t.*$", "", v)
  v <- trimws(v)
  unique(v[v != "" & tolower(v) != "taxid"])
}

# Extract sample id from Bracken report filename
sample_id_from_filename <- function(fname) sub("_.*$", "", fname)

# Build abundance matrix from Bracken *.report files
make_matrix <- function(input_path,
                        pattern      = "report$",
                        exclude_taxa = c("9606"),
                        use_rank     = "S") {
  
  files <- list.files(input_path, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) stop(sprintf("No files matching '%s' in %s", pattern, input_path))
  
  vecs <- list()
  
  for (f in files) {
    sid <- sample_id_from_filename(basename(f))
    
    df <- tryCatch(
      readr::read_tsv(
        f,
        col_names = c("percent","reads","assigned","rank","taxid","name"),
        col_types = readr::cols(
          percent  = readr::col_double(),
          reads    = readr::col_double(),
          assigned = readr::col_double(),
          rank     = readr::col_character(),
          taxid    = readr::col_character(),
          name     = readr::col_character()
        ),
        progress = FALSE,
        show_col_types = FALSE
      ),
      error = function(e) {
        warning(sprintf("Read error: %s (%s)", f, e$message))
        return(NULL)
      }
    )
    
    if (is.null(df) || nrow(df) == 0) next
    
    df$rank  <- trimws(df$rank)
    df$taxid <- trimws(df$taxid)
    
    df <- df[df$rank == use_rank & !(df$taxid %in% exclude_taxa), c("taxid","percent")]
    if (nrow(df) == 0) next
    
    df <- stats::aggregate(percent ~ taxid, data = df, sum)
    
    v <- df$percent
    names(v) <- df$taxid
    vecs[[sid]] <- v
  }
  
  if (length(vecs) == 0) stop("No usable species rows after filtering")
  
  all_taxids <- sort(unique(unlist(lapply(vecs, names))), na.last = NA)
  sample_ids <- sort(names(vecs))
  
  mat_list <- lapply(sample_ids, function(sid) {
    v <- vecs[[sid]][all_taxids]
    v[is.na(v)] <- 0
    as.numeric(v)
  })
  
  mat <- do.call(cbind, mat_list)
  rownames(mat) <- all_taxids
  colnames(mat) <- sample_ids
  storage.mode(mat) <- "double"
  mat
}

# Parse metadata sample name to (body_site, subject, replicate)
extract_body_subject_replicate <- function(x) {
  x <- as.character(x)
  m <- stringr::str_match(
    x,
    "\\.skin\\.([a-z0-9._-]+(?:\\.[a-z0-9._-]+)?)\\.(?:fresh\\.)?([A-Z])\\.(\\d+)"
  )
  if (is.na(m[1, 1])) c(NA_character_, NA_character_, NA_character_) else m[1, 2:4]
}

add_parsed_cols <- function(df) {
  tmp <- t(vapply(as.character(df$Sample_name), extract_body_subject_replicate,
                  FUN.VALUE = rep(NA_character_, 3L)))
  tmp <- as.data.frame(tmp, stringsAsFactors = FALSE)
  names(tmp) <- c("body_site","subject","replicate")
  dplyr::bind_cols(df, tmp)
}

# Jaccard distance (binary)
compute_jaccard <- function(p, q, threshold = 0) {
  p_bin <- as.integer(p > threshold)
  q_bin <- as.integer(q > threshold)
  inter <- sum(p_bin & q_bin)
  uni   <- sum(p_bin | q_bin)
  if (uni == 0) 0 else 1 - inter / uni
}

# Pairwise Jaccard calculation with prevalence/abundance filtering
calculate_pairwise_jaccard <- function(repl_pairs, std_mat, mic_mat,
                                       prevalence_threshold = 0.30,
                                       abundance_threshold  = 0.01,
                                       binarize_threshold   = 0) {
  
  common_taxids <- intersect(rownames(std_mat), rownames(mic_mat))
  std_c <- std_mat[common_taxids, , drop = FALSE]
  mic_c <- mic_mat[common_taxids, , drop = FALSE]
  
  prev  <- (rowMeans(std_c > 0) + rowMeans(mic_c > 0)) / 2
  meanA <- rowMeans(cbind(std_c, mic_c))
  
  keep <- names(which(prev >= prevalence_threshold & meanA >= abundance_threshold))
  if (!length(keep)) {
    return(list(
      jaccard_by_pair = numeric(0),
      median_jaccard  = NA_real_,
      kept_taxa       = character(0),
      kept_count      = 0
    ))
  }
  
  std_f <- std_c[keep, , drop = FALSE]
  mic_f <- mic_c[keep, , drop = FALSE]
  
  d <- c()
  for (k in seq_len(nrow(repl_pairs))) {
    s <- repl_pairs$std[k]
    m <- repl_pairs$mic[k]
    if (!(s %in% colnames(std_f)) || !(m %in% colnames(mic_f))) next
    d[paste0(s, " | ", m)] <- compute_jaccard(std_f[, s], mic_f[, m], binarize_threshold)
  }
  
  list(
    jaccard_by_pair = d,
    median_jaccard  = if (length(d)) median(d) else NA_real_,
    kept_taxa       = keep,
    kept_count      = length(keep)
  )
}

safe_median <- function(x) if (length(x)) median(x, na.rm = TRUE) else NA_real_

mw_p_greater <- function(x, y) {
  if (length(x) < 1 || length(y) < 1) return(NA_real_)
  tryCatch(wilcox.test(x, y, alternative = "greater")$p.value, error = function(e) NA_real_)
}

p_to_signif <- function(p) {
  if (is.na(p)) "ns" else if (p <= 0.01) "***" else if (p <= 0.05) "*" else "ns"
}

# Spaghetti plot data builder
# Spaghetti plot data builder (fixed order + direction)
make_spaghetti_df <- function(res_before, res_after, label_after = "Metacontam") {
  common_pairs <- intersect(names(res_before$jaccard_by_pair), names(res_after$jaccard_by_pair))
  if (!length(common_pairs)) {
    return(tibble(
      Pair = character(0),
      State = factor(character(0), levels = c("Original", label_after)),
      Distance = numeric(0),
      Direction = factor(character(0), levels = c("Increased","Decreased","Unchanged"))
    ))
  }
  
  orig <- as.numeric(res_before$jaccard_by_pair[common_pairs])
  aft  <- as.numeric(res_after$jaccard_by_pair[common_pairs])
  names(orig) <- names(aft) <- common_pairs
  
  # Direction based on (after - before) = (Metacontam - Original)
  delta <- aft - orig
  dir_lab <- ifelse(delta >  1e-12, "Increased",
                    ifelse(delta < -1e-12, "Decreased", "Unchanged"))
  
  tibble(
    Pair      = rep(common_pairs, each = 2),
    State     = rep(c("Original", label_after), times = length(common_pairs)),
    Distance  = c(rbind(orig, aft)),
    Direction = rep(dir_lab, each = 2)
  ) %>%
    mutate(
      # Force x-axis order: left Original, right Metacontam
      State = factor(State, levels = c("Original", label_after)),
      # Optional: force legend order
      Direction = factor(Direction, levels = c("Increased","Decreased","Unchanged"))
    )
}


# ------------------------------------------------------------
# A) Build replicate pairs if not already provided
# ------------------------------------------------------------
if (!exists("replicate_pairs")) {
  
  if (!exists("metadata_tbl") || !inherits(metadata_tbl, c("data.frame","tbl","tbl_df"))) {
    stopifnot(file.exists(metadata_path))
    metadata_tbl <- readr::read_csv(metadata_path, show_col_types = FALSE)
  }
  
  required_cols <- c(
    "sample_type","Assay Type","extraction_kit_id","Instrument",
    "Sample_name","collection_timestamp","Run"
  )
  missing_cols <- setdiff(required_cols, names(metadata_tbl))
  if (length(missing_cols) > 0) stop("metadata_tbl missing: ", paste(missing_cols, collapse = ", "))
  
  ms <- metadata_tbl %>%
    dplyr::filter(
      .data$sample_type == "skin",
      .data$`Assay Type` == "WGS",
      .data$extraction_kit_id == "standard",
      .data$Instrument == "Illumina NovaSeq 6000"
    )
  
  mm <- metadata_tbl %>%
    dplyr::filter(
      .data$sample_type == "skin",
      .data$`Assay Type` == "WGS",
      .data$extraction_kit_id == "microbiome",
      .data$Instrument == "Illumina NovaSeq 6000"
    )
  
  ms <- add_parsed_cols(ms)
  mm <- add_parsed_cols(mm)
  
  matched_df <- dplyr::inner_join(
    ms, mm,
    by = c("body_site","subject","replicate","collection_timestamp"),
    suffix = c("_standard","_microbiome")
  )
  
  if (!all(c("Run_standard","Run_microbiome") %in% names(matched_df))) {
    stop("Run_standard/Run_microbiome missing in joined metadata")
  }
  
  replicate_pairs <- matched_df %>%
    dplyr::transmute(std = .data$Run_standard, mic = .data$Run_microbiome) %>%
    dplyr::distinct()
  
  message(sprintf("Rebuilt replicate_pairs: %d pairs", nrow(replicate_pairs)))
}

# ------------------------------------------------------------
# B) Bracken reports -> abundance matrices (exclude 9606)
# ------------------------------------------------------------
if (!exists("Skin_standard_output"))   stop("Skin_standard_output is not defined in the environment.")
if (!exists("Skin_microbiome_output")) stop("Skin_microbiome_output is not defined in the environment.")

standard_path   <- file.path(Skin_standard_output,   "Bracken_dir")
microbiome_path <- file.path(Skin_microbiome_output, "Bracken_dir")

standard_matrix   <- make_matrix(standard_path)
microbiome_matrix <- make_matrix(microbiome_path)

# Normalize IDs (columns + replicate pairs)
colnames(standard_matrix)   <- canon(colnames(standard_matrix))
colnames(microbiome_matrix) <- canon(colnames(microbiome_matrix))
replicate_pairs <- replicate_pairs %>% mutate(std = canon(std), mic = canon(mic))

# ------------------------------------------------------------
# C) Apply filters: Metacontam / Squeegee / Decontam
# ------------------------------------------------------------

# Provide safe defaults if external objects are missing
if (!exists("squeegee_skin_microbiome_predict")) squeegee_skin_microbiome_predict <- character(0)
# Keep backward-compatibility for the typo variable name in the original script
if (!exists("squegee_skin_standard_predict") && exists("squeegee_skin_standard_predict")) {
  squegee_skin_standard_predict <- squeegee_skin_standard_predict
}
if (!exists("squegee_skin_standard_predict")) squegee_skin_standard_predict <- character(0)

if (!exists("Skin_standard_predict"))   Skin_standard_predict   <- character(0)
if (!exists("Skin_microbiome_predict")) Skin_microbiome_predict <- character(0)

# Squeegee / Metacontam taxid vectors
sq_mic_tax <- cleanup_lines_to_taxids(squeegee_skin_microbiome_predict)
sq_std_tax <- cleanup_lines_to_taxids(squegee_skin_standard_predict)

# Apply Metacontam & Squeegee filters
standard_matrix_meta   <- standard_matrix  [!(rownames(standard_matrix)   %in% Skin_standard_predict),   , drop = FALSE]
microbiome_matrix_meta <- microbiome_matrix[!(rownames(microbiome_matrix) %in% Skin_microbiome_predict), , drop = FALSE]

standard_matrix_sq     <- standard_matrix  [!(rownames(standard_matrix)   %in% sq_std_tax),              , drop = FALSE]
microbiome_matrix_sq   <- microbiome_matrix[!(rownames(microbiome_matrix) %in% sq_mic_tax),              , drop = FALSE]

# Decontam taxids from files
if (!file.exists(decon_std_path)) stop("Decontam file not found: ", decon_std_path)
if (!file.exists(decon_mic_path)) stop("Decontam file not found: ", decon_mic_path)

decon_std_tax <- cleanup_lines_to_taxids(readr::read_lines(decon_std_path))
decon_mic_tax <- cleanup_lines_to_taxids(readr::read_lines(decon_mic_path))

standard_matrix_decon   <- standard_matrix  [!(rownames(standard_matrix)   %in% decon_std_tax), , drop = FALSE]
microbiome_matrix_decon <- microbiome_matrix[!(rownames(microbiome_matrix) %in% decon_mic_tax), , drop = FALSE]

# ------------------------------------------------------------
# D) Main comparison (Boxplot) - Decontam included
# ------------------------------------------------------------
prev_thr <- 0.35
abun_thr <- 0.01
bin_thr  <- 0.00

res_orig <- calculate_pairwise_jaccard(replicate_pairs, standard_matrix,         microbiome_matrix,
                                       prevalence_threshold = prev_thr, abundance_threshold = abun_thr, binarize_threshold = bin_thr)
res_meta <- calculate_pairwise_jaccard(replicate_pairs, standard_matrix_meta,    microbiome_matrix_meta,
                                       prevalence_threshold = prev_thr, abundance_threshold = abun_thr, binarize_threshold = bin_thr)
res_sq   <- calculate_pairwise_jaccard(replicate_pairs, standard_matrix_sq,      microbiome_matrix_sq,
                                       prevalence_threshold = prev_thr, abundance_threshold = abun_thr, binarize_threshold = bin_thr)
res_dec  <- calculate_pairwise_jaccard(replicate_pairs, standard_matrix_decon,   microbiome_matrix_decon,
                                       prevalence_threshold = prev_thr, abundance_threshold = abun_thr, binarize_threshold = bin_thr)

cat(sprintf("[Kept species at prev ≥ %.2f & mean abun ≥ %.3f]\n  Original:   %4d\n  Decontam:   %4d\n  Squeegee:   %4d\n  Metacontam: %4d\n\n",
            prev_thr, abun_thr, res_orig$kept_count, res_dec$kept_count, res_sq$kept_count, res_meta$kept_count))

jaccard_orig <- unname(res_orig$jaccard_by_pair)
jaccard_meta <- unname(res_meta$jaccard_by_pair)
jaccard_sq   <- unname(res_sq$jaccard_by_pair)
jaccard_dec  <- unname(res_dec$jaccard_by_pair)

ani_df <- tibble(
  Method  = c(rep("Original",   length(jaccard_orig)),
              rep("Metacontam", length(jaccard_meta)),
              rep("Squeegee",   length(jaccard_sq)),
              rep("Decontam",   length(jaccard_dec))),
  Jaccard = c(jaccard_orig, jaccard_meta, jaccard_sq, jaccard_dec)
)
ani_df$Method <- factor(ani_df$Method, levels = METHOD_LEVELS)

# One-sided p-values: Original vs others (alternative = "greater")
p_o_meta <- tryCatch(wilcox.test(jaccard_orig, jaccard_meta, alternative = "greater")$p.value, error = function(e) NA_real_)
p_o_sq   <- tryCatch(wilcox.test(jaccard_orig, jaccard_sq,   alternative = "greater")$p.value, error = function(e) NA_real_)
p_o_dec  <- tryCatch(wilcox.test(jaccard_orig, jaccard_dec,  alternative = "greater")$p.value, error = function(e) NA_real_)

y_upper <- 1.10
pval_df <- tibble::tibble(
  group1     = c("Original","Original","Original"),
  group2     = c("Decontam","Squeegee","Metacontam"),
  y.position = c(y_upper - 0.03, y_upper - 0.06, y_upper - 0.09),
  p.signif   = c(p_to_signif(p_o_dec), p_to_signif(p_o_sq), p_to_signif(p_o_meta))
)
median_line_y <- min(pval_df$y.position) - 0.02

med_df <- ani_df %>%
  dplyr::group_by(Method) %>%
  dplyr::summarise(median_val = median(Jaccard), .groups = "drop") %>%
  dplyr::mutate(x_pos = as.numeric(Method), y_pos = median_line_y)

# --- Boxplot (kept as-is as much as possible) ---
Skin_replicate <- ggboxplot(
  ani_df, x = "Method", y = "Jaccard",
  color = "Method", palette = METHOD_COLORS,
  add = "jitter", add.params = list(size = 2),
  size = 0.9
) +
  geom_text(data = med_df,
            aes(x = x_pos, y = y_pos, label = sprintf("%.3f", median_val)),
            inherit.aes = FALSE, size = 5, vjust = 1) +
  stat_pvalue_manual(
    pval_df, label = "p.signif",
    tip.length  = 0.012,
    size        = 6,
    bracket.size= 0.6
  ) +
  scale_y_continuous(
    limits = c(0.1, 1.1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(title = "Pairwise Jaccard Across Filtering Methods") +
  theme_pubr(base_size = 15) +
  theme(
    plot.title         = element_text(size = 15, hjust = 0.5, face = "bold"),
    axis.title.x       = element_blank(),
    axis.title.y       = element_blank(),
    axis.text.x        = element_text(size = 15),
    axis.text.y        = element_text(size = 15),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
    axis.ticks.y       = element_line(color = "grey70", linewidth = 0.3),
    axis.ticks.length  = unit(3, "pt"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "none",
    panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.margin        = margin(t = 10, r = 5, b = 5, l = 5),
    aspect.ratio       = 1
  ) +
  coord_cartesian(clip = "off")

print(Skin_replicate)

# ------------------------------------------------------------
# E) Prevalence sweep (Decontam included)
# - Stars: red = Orig > Meta, blue = Orig > Decon
# ------------------------------------------------------------
prevalence_range <- seq(0, 0.50, by = 0.01)
min_species      <- 10
ab_thresh        <- 0.01

# Species count (based on Original matrices)
species_counts <- numeric(length(prevalence_range))
{
  common    <- intersect(rownames(standard_matrix), rownames(microbiome_matrix))
  std_c     <- standard_matrix [common, , drop = FALSE]
  mic_c     <- microbiome_matrix[common, , drop = FALSE]
  prev_std  <- rowMeans(std_c > 0)
  prev_mic  <- rowMeans(mic_c > 0)
  combined  <- (prev_std + prev_mic) / 2
  mean_abun <- rowMeans(cbind(std_c, mic_c))
  
  for (i in seq_along(prevalence_range)) {
    prev <- prevalence_range[i]
    species_counts[i] <- sum(combined >= prev & mean_abun >= ab_thresh)
  }
}

cat("=== Species kept by prevalence threshold (Original; mean abun ≥ ", ab_thresh, ") ===\n", sep = "")
for (i in seq_along(prevalence_range)) {
  cat(sprintf("  prev ≥ %4.2f : %4d species\n", prevalence_range[i], species_counts[i]))
}
cat("\n")

res_list <- vector("list", length(prevalence_range))
pvals_meta  <- numeric(length(prevalence_range))
pvals_decon <- numeric(length(prevalence_range))

for (i in seq_along(prevalence_range)) {
  prev <- prevalence_range[i]
  
  j_orig <- unname(calculate_pairwise_jaccard(
    replicate_pairs, standard_matrix, microbiome_matrix,
    prevalence_threshold = prev, abundance_threshold = ab_thresh, binarize_threshold = 0
  )$jaccard_by_pair)
  
  j_meta <- unname(calculate_pairwise_jaccard(
    replicate_pairs, standard_matrix_meta, microbiome_matrix_meta,
    prevalence_threshold = prev, abundance_threshold = ab_thresh, binarize_threshold = 0
  )$jaccard_by_pair)
  
  j_sque <- unname(calculate_pairwise_jaccard(
    replicate_pairs, standard_matrix_sq, microbiome_matrix_sq,
    prevalence_threshold = prev, abundance_threshold = ab_thresh, binarize_threshold = 0
  )$jaccard_by_pair)
  
  j_deco <- unname(calculate_pairwise_jaccard(
    replicate_pairs, standard_matrix_decon, microbiome_matrix_decon,
    prevalence_threshold = prev, abundance_threshold = ab_thresh, binarize_threshold = 0
  )$jaccard_by_pair)
  
  res_list[[i]] <- tibble(
    Prevalence = prev,
    Original   = safe_median(j_orig),
    Decontam   = safe_median(j_deco),
    Squeegee   = safe_median(j_sque),
    Metacontam = safe_median(j_meta)
  )
  
  pvals_meta[i]  <- mw_p_greater(j_orig, j_meta)  # Orig > Meta
  pvals_decon[i] <- mw_p_greater(j_orig, j_deco)  # Orig > Decon
}

df <- bind_rows(res_list)
df$SpeciesCount        <- species_counts
df$Pval_Orig_vs_Meta   <- pvals_meta
df$Pval_Orig_vs_Decon  <- pvals_decon

df_long <- df %>%
  pivot_longer(cols = c(Original, Decontam, Squeegee, Metacontam),
               names_to = "Method", values_to = "MedianJaccard") %>%
  mutate(MedianJaccard_plot = if_else(Method %in% c("Metacontam","Decontam") & SpeciesCount < min_species,
                                      NA_real_, MedianJaccard))

df_long$Method <- factor(df_long$Method, levels = METHOD_LEVELS)

star_df_meta <- df %>%
  filter(!is.na(Metacontam), SpeciesCount >= min_species, Pval_Orig_vs_Meta  <= 0.05) %>%
  transmute(Prevalence, y = Metacontam)

star_df_deco <- df %>%
  filter(!is.na(Decontam),   SpeciesCount >= min_species, Pval_Orig_vs_Decon <= 0.05) %>%
  transmute(Prevalence, y = Decontam)

y_min <- min(df_long$MedianJaccard, na.rm = TRUE)
y_max <- max(df_long$MedianJaccard, na.rm = TRUE)
offset <- (y_max - y_min) * 0.001

# --- Prevalence sweep plot (kept as-is) ---
p_prev <- ggplot(df_long, aes(x = Prevalence, y = MedianJaccard_plot, color = Method)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  geom_text(data = star_df_meta, aes(x = Prevalence, y = y + offset),
            label = "*", inherit.aes = FALSE, vjust = 0, size = 6, color = "red3") +
  geom_text(data = star_df_deco, aes(x = Prevalence, y = y + offset),
            label = "*", inherit.aes = FALSE, vjust = 0, size = 6, color = "blue3") +
  scale_color_manual(values = METHOD_COLORS, breaks = METHOD_LEVELS) +
  scale_x_continuous(breaks = seq(0, 0.5, by = 0.1), limits = c(0, 0.5)) +
  labs(
    title = sprintf("Median Jaccard vs. Prevalence Threshold\n(Meta/Decon shown only when ≥%d taxa; * p≤0.05)", min_species),
    x = "Prevalence Threshold", y = "Median Jaccard Distance", color = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5),
    panel.grid.major = element_line(linetype = "dashed", linewidth = 0.3, color = "grey80"),
    panel.grid.minor = element_blank(),
    legend.position  = c(0.18, 0.2),
    legend.background     = element_rect(fill = scales::alpha("white", 0.9), color = "black"),
    legend.key            = element_rect(fill = scales::alpha("white", 0), color = NA),
    legend.box.background = element_rect(color = "black", fill = scales::alpha("white", 0.9)),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.background = element_rect(fill = "white", color = NA)
  )

print(p_prev)

# ------------------------------------------------------------
# F) Replicate pairs check + Excel export
# ------------------------------------------------------------
print(as.data.frame(replicate_pairs))

out <- replicate_pairs
colnames(out) <- c("Standard kit", "Microbiome kit")

write.xlsx(
  x = list("Supplementary — Replicate pairs" = out),
  file = "/media/junwoojo/18T/Metacontam_dataset/Supplementary_information/Supplementary_replicate_pairs.xlsx",
  asTable = TRUE
)

# ------------------------------------------------------------
# G) Spaghetti plot (Original -> Metacontam change)
# ------------------------------------------------------------
df_slope_long <- make_spaghetti_df(res_orig, res_meta, label_after = "Metacontam")

ymin <- suppressWarnings(min(df_slope_long$Distance, na.rm = TRUE))
ymax <- suppressWarnings(max(df_slope_long$Distance, na.rm = TRUE))
yrng <- ymax - ymin

pad_base   <- if (is.finite(yrng) && yrng > 0) yrng * 0.02 else 0.01
pad_legend <- if (is.finite(yrng) && yrng > 0) yrng * 0.20 else 0.20
ylims <- c(max(0, ymin - pad_base), ymax + pad_base + pad_legend)

# --- Spaghetti plot (kept as-is; removed redundant re-scale duplication) ---
p_spaghetti_replicate <- ggplot(df_slope_long, aes(x = State, y = Distance, group = Pair, color = Direction)) +
  geom_line(linewidth = 0.9, alpha = 0.9) +
  geom_point(size = 2.1) +
  scale_color_manual(values = SLOPE_COLORS,
                     breaks = c("Increased","Decreased","Unchanged")) +
  scale_x_discrete(expand = expansion(mult = 0, add = 0.09)) +
  scale_y_continuous(limits = ylims, breaks = scales::pretty_breaks(5),
                     expand = expansion(mult = 0, add = 0)) +
  labs(
    title = "Pairwise Jaccard Change by Metacontam Filtering",
    x = NULL, y = "Jaccard Distance", color = NULL
  ) +
  theme_pubr(base_size = 14) +
  theme(
    plot.margin        = margin(t = 16, r = 8, b = 6, l = 8),
    plot.title         = element_text(size = 14, hjust = 0.5, face = "bold", margin = margin(b = 4)),
    axis.text.x        = element_text(size = 13, margin = margin(t = 10)),
    axis.text.y        = element_text(size = 12, margin = margin(r = 10)),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.ticks.length  = unit(2, "pt"),
    legend.position      = c(0.5, 0.95),
    legend.justification = c(0.5, 1),
    legend.direction     = "horizontal",
    legend.margin        = margin(t = 0, r = 2, b = 2, l = 2),
    legend.box.margin    = margin(0, 0, 0, 0),
    legend.key.height    = unit(10, "pt"),
    legend.key.width     = unit(14, "pt"),
    panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.9),
    panel.background   = element_rect(fill = "white", color = NA)
  ) +
  coord_cartesian(clip = "on")

print(p_spaghetti_replicate)
 
# ------------------------------------------------------------
# H) (Optional) Print one-sided Wilcoxon direction + raw p-values
# ------------------------------------------------------------
cat("=== One-sided Wilcoxon (Boxplot: Jaccard) ===\n")
cat("wilcox.test(x, y, alternative='greater') was used.\n")
cat("H1: x(Original) > y(Method)  => tests whether the method yields smaller Jaccard distances.\n\n")

cat(sprintf("H1: Original > Decontam   (Decontam smaller)   p = %.6g   med(O)=%.3f  med(D)=%.3f\n",
            p_o_dec,  safe_median(jaccard_orig), safe_median(jaccard_dec)))
cat(sprintf("H1: Original > Squeegee   (Squeegee smaller)   p = %.6g   med(O)=%.3f  med(S)=%.3f\n",
            p_o_sq,   safe_median(jaccard_orig), safe_median(jaccard_sq)))
cat(sprintf("H1: Original > Metacontam (Metacontam smaller) p = %.6g   med(O)=%.3f  med(M)=%.3f\n\n",
            p_o_meta, safe_median(jaccard_orig), safe_median(jaccard_meta)))




Figure3b <- Skin_replicate
Figure3c <- p_spaghetti_replicate

