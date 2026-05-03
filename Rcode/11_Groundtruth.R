suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(tibble)
  library(ggsci)  # <- pal_aaas
})


# Palette (place after ggsci is loaded)
nejm_colors <- pal_aaas("default")(8)

# ---------------------------
# 0) Utility: whitelist summary & print
# ---------------------------
summarize_whitelist <- function(label, whitelist, blacklist) {
  wl <- unique(as.character(whitelist))
  bl <- unique(as.character(blacklist))
  wl <- wl[!is.na(wl) & nzchar(wl)]
  bl <- bl[!is.na(bl) & nzchar(bl)]
  
  n_total <- length(wl)
  n_removed <- sum(wl %in% bl)
  wl_used <- setdiff(wl, bl)
  n_used <- length(wl_used)
  
  cat(sprintf("[Whitelist - %s] total=%d, removed_by_blacklist=%d, used=%d\n",
              label, n_total, n_removed, n_used))
  return(list(filtered = wl_used,
              total = n_total,
              removed = n_removed,
              used = n_used))
}

# ---------------------------
# 1) Prevalence calculation (optionally exclude 9606, etc.)
# ---------------------------
calculate_relative_prevalence <- function(report_paths, output_filename,
                                          abundance_threshold = 0.001,
                                          drop_taxids = c("9606")) {
  df_list <- lapply(report_paths, function(file) {
    df <- readr::read_tsv(file, col_names = FALSE, show_col_types = FALSE) %>%
      dplyr::select(X1, X4, X5) %>%
      dplyr::rename(abundance = X1, rank = X4, taxid = X5) %>%
      dplyr::filter(rank == "S") %>%
      dplyr::select(taxid, abundance) %>%
      dplyr::mutate(taxid = as.character(taxid))
    
    df <- df %>%
      dplyr::group_by(taxid) %>%
      dplyr::summarise(abundance = sum(abundance), .groups = "drop") %>%
      dplyr::mutate(rel_abun = abundance / sum(abundance)) %>%
      dplyr::select(taxid, rel_abun)
    
    colnames(df)[2] <- file
    df
  })
  
  merged_df <- purrr::reduce(df_list, dplyr::full_join, by = "taxid") %>% 
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ tidyr::replace_na(.x, 0)))
  
  prevalence <- merged_df %>%
    tibble::column_to_rownames("taxid") %>%
    dplyr::mutate(prevalence = rowSums(. >= abundance_threshold) / ncol(.))
  
  prevalence_df <- prevalence %>%
    dplyr::select(prevalence) %>%
    tibble::rownames_to_column("taxid")
  
  if (!is.null(drop_taxids) && length(drop_taxids) > 0) {
    prevalence_df <- prevalence_df %>% dplyr::filter(!taxid %in% drop_taxids)
  }
  
  readr::write_tsv(prevalence_df, output_filename)
  prevalence_df
}
# ---------------------------
# 2) Plot function
# ---------------------------
plot_threshold_ratio <- function(prevalence_df, whitelist, blacklist, title, save_path) {
  # Use whitelist after excluding blacklist entries
  whitelist <- setdiff(as.character(whitelist), as.character(blacklist))
  thresholds <- round(seq(0.1, 1.0, by = 0.1), 1)
  
  stats <- purrr::map_dfr(thresholds, function(t) {
    final_contam <- prevalence_df %>%
      dplyr::filter(prevalence >= t) %>%
      dplyr::pull(taxid) %>%
      as.character()
    
    total <- length(final_contam)
    white_count <- sum(final_contam %in% whitelist)
    
    tibble::tibble(
      threshold     = t,
      total_species = total,
      white_ratio   = ifelse(total > 0, white_count / total, 0)
    )
  })
  
  valid <- dplyr::filter(stats, total_species > 0)
  min_ratio <- if (nrow(valid) > 0) min(valid$white_ratio, na.rm = TRUE) else 0
  max_ratio <- if (nrow(valid) > 0) max(valid$white_ratio, na.rm = TRUE) else 1
  if (!is.finite(min_ratio)) min_ratio <- 0
  if (!is.finite(max_ratio)) max_ratio <- 1
  
  pad <- if (max_ratio > min_ratio) max(0.02, 0.10 * (max_ratio - min_ratio)) else 0.05
  y_lower <- max(0, min_ratio - pad)
  y_upper <- min(1, max_ratio + pad)
  if (y_upper - y_lower < 0.05) y_upper <- min(1, y_lower + 0.05)
  
  max_n <- max(stats$total_species, na.rm = TRUE)
  lift <- (y_upper - y_lower) * 0.015
  species_y_true <- if (max_n > 0) ((stats$total_species / max_n) * (y_upper - y_lower)) + y_lower else rep(y_lower, nrow(stats))
  species_y_plot <- ifelse(stats$total_species > 0, species_y_true, y_lower + lift)
  
  min_thr <- if (nrow(valid) > 0) valid$threshold[which.min(valid$white_ratio)] else NA_real_
  
  p <- ggplot(stats, aes(x = threshold)) +
    geom_line(aes(y = white_ratio, color = "Whitelist ratio"), size = 1.3) +
    geom_point(aes(y = white_ratio, color = "Whitelist ratio"), size = 2.4) +
    geom_line(aes(y = species_y_plot, color = "Species count"), size = 1.3) +
    geom_point(aes(y = species_y_plot, color = "Species count"), shape = 17, size = 2.4) +
    geom_text(aes(y = species_y_plot, label = total_species),
              vjust = -0.6, color = nejm_colors[2], size = 4) +
    { if (!is.na(min_thr)) geom_vline(xintercept = min_thr, linetype = "dashed") } +
    scale_x_continuous(breaks = thresholds, limits = c(0.1, 1.0)) +
    scale_y_continuous(
      name   = "Whitelist Ratio",
      limits = c(y_lower, y_upper),
      breaks = scales::pretty_breaks(n = 5),
      sec.axis = sec_axis(~ (.- y_lower) * (if (y_upper > y_lower) max_n / (y_upper - y_lower) else 0),
                          name = "Species Count")
    ) +
    scale_color_manual(values = c("Whitelist ratio" = nejm_colors[1], "Species count" = nejm_colors[2])) +
    theme_bw(base_size = 16) +
    theme(
      axis.text.x        = element_text(size = 14),
      axis.text.y        = element_text(size = 14),
      axis.title.x       = element_text(size = 16),
      axis.title.y       = element_text(size = 16),
      axis.title.y.right = element_text(size = 16),
      legend.text        = element_text(size = 13),
      legend.title       = element_blank()
    ) +
    labs(title = title, x = "Prevalence Threshold")
  
  # ggsave(save_path, p, width = 8, height = 5, dpi = 600)  # removed
  p
}

# ---------------------------
# 3) Run per dataset (Skin: fix set by intersecting file lists)
# ---------------------------
file.path(BASE_INPUT, "Skin_sra_metadata")

# Skin
skin_path   <- file.path(BASE_INPUT, "Skin_negative", "Skin_negative_bracken_files")
metadata    <- read_csv(file.path(BASE_INPUT, "Skin_sra_metadata"), show_col_types = FALSE)
report_path <- list.files(skin_path, full.names = TRUE, pattern = "_kraken_report_bracken_species$")

metadata2 <- metadata %>%
  filter(sample_type == "control blank",
         Instrument  == "Illumina NovaSeq 6000")

standard_runs   <- metadata2 %>% filter(extraction_kit_id == "standard")   %>% pull(Run)
microbiome_runs <- metadata2 %>% filter(extraction_kit_id == "microbiome") %>% pull(Run)

standard_kit_paths   <- file.path(skin_path, paste0(standard_runs,   "_kraken_report_bracken_species"))
microbiome_kit_paths <- file.path(skin_path, paste0(microbiome_runs, "_kraken_report_bracken_species"))

standard_reports   <- intersect(report_path, standard_kit_paths)
microbiome_reports <- intersect(report_path, microbiome_kit_paths)

standard_prevalence_df <- calculate_relative_prevalence(
  standard_reports, "standard_prevalence_relative.tsv",
  abundance_threshold = 0.001, drop_taxids = c("9606")
)

microbiome_prevalence_df <- calculate_relative_prevalence(
  microbiome_reports, "microbiome_prevalence_relative.tsv",
  abundance_threshold = 0.001, drop_taxids = c("9606")
)

# Nasal
nasal_path <- file.path(BASE_INPUT, "Nasal_negative_bracken")
nasal_reports <- list.files(nasal_path, full.names = TRUE, pattern = "_bracken_species\\.report$")
nasal_prevalence_df <- calculate_relative_prevalence(
  nasal_reports, "nasal_prevalence_relative.tsv",
  abundance_threshold = 0.001, drop_taxids = c("9606")
)

# Oral
oral_path <- file.path(BASE_INPUT, "Oral_negative_bracken")
oral_reports <- list.files(oral_path, full.names = TRUE, pattern = "_bracken_species\\.report$")
oral_prevalence_df <- calculate_relative_prevalence(
  oral_reports, "oral_prevalence_relative.tsv",
  abundance_threshold = 0.001, drop_taxids = c("9606")
)

# ---------------------------
# 4) Whitelist / Blacklist + count prints
# ---------------------------
file.path(BASE_INPUT, "BlackWhite_list" ,"blacklist.csv")

skin_whitelist_raw  <- read_csv(file.path(BASE_INPUT, "BlackWhite_list" ,"Healthy_Skin_Associate"),
                                col_names = FALSE, show_col_types = FALSE)[[1]] %>% as.character()
nasal_whitelist_raw <- read_csv(file.path(BASE_INPUT, "BlackWhite_list" ,"Healthy_Nasal_associated"),
                                col_names = FALSE, show_col_types = FALSE)[[1]] %>% as.character()
oral_whitelist_raw  <- read_csv(file.path(BASE_INPUT, "BlackWhite_list" ,"healthy_Associate_oral_taxid"),
                                col_names = FALSE, show_col_types = FALSE)[[1]] %>% as.character()

blacklist <- read_tsv(file.path(BASE_INPUT, "BlackWhite_list" ,"blacklist.csv"),
                      show_col_types = FALSE)$blacklist %>% as.character()

wl_skin  <- summarize_whitelist("Skin",  skin_whitelist_raw,  blacklist)
wl_nasal <- summarize_whitelist("Nasal", nasal_whitelist_raw, blacklist)
wl_oral  <- summarize_whitelist("Oral",  oral_whitelist_raw,  blacklist)

# ---------------------------
# 5) Plot objects (uses whitelist after blacklist exclusion)
# ---------------------------
plot_standard_gt <- plot_threshold_ratio(
  standard_prevalence_df, wl_skin$filtered, blacklist,
  "Skin-Associated species ratio - Standard kit",
  "/home/skin_associate_ratio_standard.pdf"
)

plot_microbiome_gt <- plot_threshold_ratio(
  microbiome_prevalence_df, wl_skin$filtered, blacklist,
  "Skin-Associated species ratio - Microbiome kit",
  "/home/skin_associate_ratio_microbiome.pdf"
)

plot_nasal_gt <- plot_threshold_ratio(
  nasal_prevalence_df, wl_nasal$filtered, blacklist,
  "Nasal-Associated species ratio",
  "/home/nasal_associate_ratio.pdf"
)

plot_oral_gt <- plot_threshold_ratio(
  oral_prevalence_df, wl_oral$filtered, blacklist,
  "Oral-Associated species ratio",
  "/home/oral_associate_ratio.pdf"
)

# Note: to combine multiple plots, patchwork/cowplot is recommended
library(patchwork)
p_ground_truth <- plot_standard_gt + plot_microbiome_gt + plot_nasal_gt + plot_oral_gt

p_ground_truth

# ---------------------------
# (Extra) Helper: print sample counts
# ---------------------------
print_sample_count <- function(label, report_paths, expected_runs = NULL) {
  n_files <- length(report_paths)
  msg <- sprintf("[Sample count] %s : n = %d", label, n_files)
  
  if (!is.null(expected_runs)) {
    msg <- sprintf("%s (expected from metadata = %d)", msg, length(expected_runs))
    
    # Check missing Run IDs (Skin only)
    missing <- setdiff(expected_runs,
                       sub("_kraken_report_bracken_species$", "", basename(report_paths)))
    if (length(missing) > 0) {
      msg <- paste0(msg, sprintf(" | missing=%d", length(missing)))
      cat(msg, "\n  Missing Run IDs: ",
          paste(head(missing, 20), collapse = ", "),
          if (length(missing) > 20) " ...", "\n", sep = "")
      return(invisible(n_files))
    }
  }
  
  cat(msg, "\n")
  invisible(n_files)
}

# ---------------------------
# (Extra) Print sample counts per panel
# ---------------------------
print_sample_count("Skin - Standard (negative controls)",   standard_reports,   expected_runs = standard_runs)
print_sample_count("Skin - Microbiome (negative controls)", microbiome_reports, expected_runs = microbiome_runs)
print_sample_count("Nasal (negative controls)",             nasal_reports)
print_sample_count("Oral  (negative controls)",             oral_reports)

###################################################################################################

Supp_Fig9 <- p_ground_truth



