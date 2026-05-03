# ============================================================
# Project Paths & Env Setup (Submission-ready)
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(reshape2)
  library(ggplot2)
  library(ggsci)
})


if (!requireNamespace("readr", quietly = TRUE)) {
  stop("Package 'readr' is required. install.packages('readr')")
}


# ------------------------------------------------------------
# 0) BASE_INPUT (staged inputs root)
# ------------------------------------------------------------
BASE_INPUT <- Sys.getenv("BASE_INPUT", unset = "/media/junwoojo/18T/Submission/Rcode/inputs")

assert_exists <- function(paths) {
  miss <- paths[!file.exists(paths)]
  if (length(miss)) stop("Missing path(s):\n- ", paste(miss, collapse = "\n- "))
}

# ------------------------------------------------------------
# 1) GT_species load from TSV
#    /.../inputs/config/GT_species.tsv
# ------------------------------------------------------------
GT_path <- file.path(BASE_INPUT, "config", "GT_species.tsv")
assert_exists(GT_path)

GT_df <- read.table(GT_path, sep = "\t", header = TRUE,
                    stringsAsFactors = FALSE, check.names = FALSE)
GT_df$group <- as.character(GT_df$group)
GT_df$taxid <- as.character(GT_df$taxid)

GT_species <- split(GT_df$taxid, GT_df$group)

# (선택) 기존 코드 호환용 alias
# highlight_list <- GT_species
 
print(sapply(GT_species, length))

# ------------------------------------------------------------
# 2) Output folders (Metacontam) - staged
# ------------------------------------------------------------
MI_output              <- file.path(BASE_INPUT, "Metacontam_output_folder", "metacontam_MI_out")
Skin_microbiome_output <- file.path(BASE_INPUT, "Metacontam_output_folder", "metacontam_microbiome_out")
Skin_standard_output   <- file.path(BASE_INPUT, "Metacontam_output_folder", "metacontam_standard_out")
Nasal_output           <- file.path(BASE_INPUT, "Metacontam_output_folder", "metacontam_nasal_out")
Oral_output            <- file.path(BASE_INPUT, "Metacontam_output_folder", "metacontam_oral_out")

# Simulation outputs (staged)
path_Simul_1_50_output <- file.path(BASE_INPUT, "Metacontam_output_folder", "Simulation", "metacontam_1_50_out")
path_Simul_1_30_output <- file.path(BASE_INPUT, "Metacontam_output_folder", "Simulation", "metacontam_1_30_out")
path_Simul_1_10_output <- file.path(BASE_INPUT, "Metacontam_output_folder", "Simulation", "metacontam_1_10_out")

path_Simul_0.5_50_output <- file.path(BASE_INPUT, "Metacontam_output_folder", "Simulation", "metacontam_0.5_50_out")
path_Simul_0.5_30_output <- file.path(BASE_INPUT, "Metacontam_output_folder", "Simulation", "metacontam_0.5_30_out")
path_Simul_0.5_10_output <- file.path(BASE_INPUT, "Metacontam_output_folder", "Simulation", "metacontam_0.5_10_out")

path_Simul_0.1_50_output <- file.path(BASE_INPUT, "Metacontam_output_folder", "Simulation", "metacontam_0.1_50_out")
path_Simul_0.1_30_output <- file.path(BASE_INPUT, "Metacontam_output_folder", "Simulation", "metacontam_0.1_30_out")
path_Simul_0.1_10_output <- file.path(BASE_INPUT, "Metacontam_output_folder", "Simulation", "metacontam_0.1_10_out")

# 기본 output 폴더 존재 체크
assert_exists(c(
  MI_output, Skin_microbiome_output, Skin_standard_output, Nasal_output, Oral_output
))

# ------------------------------------------------------------
# 3) Squeegee / Decontam (staged)
# ------------------------------------------------------------
Squeegee_output_folder <- file.path(BASE_INPUT, "Squeegee_output_folder")
Decontam_output_folder <- file.path(BASE_INPUT, "Decontam_output_folder")

read_vec1 <- function(path, sep = "", header = FALSE, col = 1) {
  if (!file.exists(path)) stop("Missing file: ", path)
  x <- read.table(path, header = header, sep = sep, stringsAsFactors = FALSE, check.names = FALSE)
  as.character(x[[col]])
}

# Squeegee
squeegee_skin_microbiome_predict <- read_vec1(
  file.path(Squeegee_output_folder, "squeegee_skin_microbiome", "final_predictions.txt")
)

squegee_skin_standard_predict <- read_vec1(
  file.path(Squeegee_output_folder, "squeegee_skin_standard", "final_predictions.txt")
)

squeegee_MI_predict <- read_vec1(
  file.path(Squeegee_output_folder, "squeegee_MI", "final_predictions.txt"),
  sep = "\t"
)

squeegee_nasal_predict <- read_vec1(
  file.path(Squeegee_output_folder, "squeegee_nasal", "final_predictions.txt")
)

squeegee_nasal_oral_predict <- read_vec1(
  file.path(Squeegee_output_folder, "squeegee_nasal_oral_out", "final_predictions.txt"),
  sep = "\t"
)

# Decontam (주의: MI decontam 파일은 Squeegee_output_folder 쪽에 있음)
decontam_MI_predict <- read_vec1(
  file.path(Squeegee_output_folder, "Decontam_MI", "decontam_final.txt"),
  sep = "\t"
)

decontam_skin_microbiome_predict <- read_vec1(
  file.path(Decontam_output_folder, "skin_microbiome_contaminant_taxids.txt"),
  sep = "\t"
)

decontam_skin_standard_predict <- read_vec1(
  file.path(Decontam_output_folder, "skin_standard_contaminant_taxids.txt"),
  sep = "\t"
)

decontam_nasal_predict <- read_vec1(
  file.path(Decontam_output_folder, "nasal_contaminant_taxids.txt"),
  sep = "\t"
)

decontam_oral_predict <- read_vec1(
  file.path(Decontam_output_folder, "oral_contaminant_taxids.txt"),
  sep = "\t"
)





# partitions
path_MI_partition                 <- file.path(MI_output, "Network_Output", "MI_partition")
path_Skin_microbiome_partition    <- file.path(Skin_microbiome_output, "Network_Output", "Skin_microbiome_partition")
path_Skin_standard_partition      <- file.path(Skin_standard_output, "Network_Output", "Skin_standard_partition")
path_Nasal_partition              <- file.path(Nasal_output, "Network_Output", "Nasal_partition")
path_Oral_partition              <- file.path(Oral_output, "Network_Output", "Oral_partition")

path_Simulation_1_50_partition    <- file.path(path_Simul_1_50_output,   "Network_Output", "Simulation_1percent_50_partition")
path_Simulation_1_30_partition    <- file.path(path_Simul_1_30_output,   "Network_Output", "Simulation_1percent_30_partition")
path_Simulation_1_10_partition    <- file.path(path_Simul_1_10_output,   "Network_Output", "Simulation_1percent_10_partition")

path_Simulation_0.5_50_partition  <- file.path(path_Simul_0.5_50_output, "Network_Output", "Simulation_0.5percent_50_partition")
path_Simulation_0.5_30_partition  <- file.path(path_Simul_0.5_30_output, "Network_Output", "Simulation_0.5percent_30_partition")
path_Simulation_0.5_10_partition  <- file.path(path_Simul_0.5_10_output, "Network_Output", "Simulation_0.5percent_10_partition")

path_Simulation_0.1_50_partition  <- file.path(path_Simul_0.1_50_output, "Network_Output", "Simulation_0.1percent_50_partition")
path_Simulation_0.1_30_partition  <- file.path(path_Simul_0.1_30_output, "Network_Output", "Simulation_0.1percent_30_partition")
path_Simulation_0.1_10_partition  <- file.path(path_Simul_0.1_10_output, "Network_Output", "Simulation_0.1percent_10_partition")













# edges
path_MI_edge              <- file.path(MI_output, "Network_Output", "edge.tsv")
path_Skin_microbiome_edge <- file.path(Skin_microbiome_output, "Network_Output", "edge.tsv")
path_Skin_standard_edge   <- file.path(Skin_standard_output, "Network_Output", "edge.tsv")
path_Nasal_edge           <- file.path(Nasal_output, "Network_Output", "edge.tsv")
path_Oral_edge            <- file.path(Oral_output, "Network_Output", "edge.tsv")

path_Simulation_1_50_edge <- file.path(path_Simul_1_50_output, "Network_Output", "edge.tsv")
path_Simulation_1_30_edge <- file.path(path_Simul_1_30_output, "Network_Output", "edge.tsv")
path_Simulation_1_10_edge <- file.path(path_Simul_1_10_output, "Network_Output", "edge.tsv")

path_Simulation_0.5_50_edge <- file.path(path_Simul_0.5_50_output, "Network_Output", "edge.tsv")
path_Simulation_0.5_30_edge <- file.path(path_Simul_0.5_30_output, "Network_Output", "edge.tsv")
path_Simulation_0.5_10_edge <- file.path(path_Simul_0.5_10_output, "Network_Output", "edge.tsv")

path_Simulation_0.1_50_edge <- file.path(path_Simul_0.1_50_output, "Network_Output", "edge.tsv")
path_Simulation_0.1_30_edge <- file.path(path_Simul_0.1_30_output, "Network_Output", "edge.tsv")
path_Simulation_0.1_10_edge <- file.path(path_Simul_0.1_10_output, "Network_Output", "edge.tsv")

# community
path_MI_community              <- file.path(MI_output, "Network_Output", "Contaminant_community_taxid")
path_Skin_microbiome_community <- file.path(Skin_microbiome_output, "Network_Output", "Contaminant_community_taxid")
path_Skin_standard_community   <- file.path(Skin_standard_output, "Network_Output", "Contaminant_community_taxid")
path_Nasal_community           <- file.path(Nasal_output, "Network_Output", "Contaminant_community_taxid")
path_Oral_community            <- file.path(Oral_output, "Network_Output", "Contaminant_community_taxid")

path_Simulation_1_50_community <- file.path(path_Simul_1_50_output, "Network_Output", "Contaminant_community_taxid")
path_Simulation_1_30_community <- file.path(path_Simul_1_30_output, "Network_Output", "Contaminant_community_taxid")
path_Simulation_1_10_community <- file.path(path_Simul_1_10_output, "Network_Output", "Contaminant_community_taxid")

path_Simulation_0.5_50_community <- file.path(path_Simul_0.5_50_output, "Network_Output", "Contaminant_community_taxid")
path_Simulation_0.5_30_community <- file.path(path_Simul_0.5_30_output, "Network_Output", "Contaminant_community_taxid")
path_Simulation_0.5_10_community <- file.path(path_Simul_0.5_10_output, "Network_Output", "Contaminant_community_taxid")

path_Simulation_0.1_50_community <- file.path(path_Simul_0.1_50_output, "Network_Output", "Contaminant_community_taxid")
path_Simulation_0.1_30_community <- file.path(path_Simul_0.1_30_output, "Network_Output", "Contaminant_community_taxid")
path_Simulation_0.1_10_community <- file.path(path_Simul_0.1_10_output, "Network_Output", "Contaminant_community_taxid")

read_community <- function(path, col = 1, sep = "", header = FALSE,
                           trim = TRUE, drop_empty = TRUE) {
  if (!file.exists(path)) {
    warning(sprintf("File not found: %s", path)); return(character())
  }
  df <- tryCatch(
    read.table(path, sep = sep, header = header, quote = "",
               comment.char = "", check.names = FALSE,
               stringsAsFactors = FALSE, fill = TRUE),
    error = function(e) { warning(sprintf("Read error for %s: %s", path, e$message)); data.frame() }
  )
  if (nrow(df) == 0 || ncol(df) < col) return(character())
  x <- as.character(df[[col]])
  if (trim) x <- trimws(x)
  if (drop_empty) x <- x[nzchar(x)]
  x
}

MI_community              <- read_community(path_MI_community)
skin_microbiome_community <- read_community(path_Skin_microbiome_community)
skin_standard_community   <- read_community(path_Skin_standard_community)
nasal_community           <- read_community(path_Nasal_community)
oral_community            <- read_community(path_Oral_community)

simul_1_50_community <- read_community(path_Simulation_1_50_community)
simul_1_30_community <- read_community(path_Simulation_1_30_community)
simul_1_10_community <- read_community(path_Simulation_1_10_community)

simul_0.5_50_community <- read_community(path_Simulation_0.5_50_community)
simul_0.5_30_community <- read_community(path_Simulation_0.5_30_community)
simul_0.5_10_community <- read_community(path_Simulation_0.5_10_community)

simul_0.1_50_community <- read_community(path_Simulation_0.1_50_community)
simul_0.1_30_community <- read_community(path_Simulation_0.1_30_community)
simul_0.1_10_community <- read_community(path_Simulation_0.1_10_community)

# ------------------------------------------------------------
# 5) Prediction (Final_prediction.txt) -> Contaminant Taxids
# ------------------------------------------------------------
path_MI_predict              <- file.path(MI_output, "Final_prediction.txt")
path_Skin_microbiome_predict <- file.path(Skin_microbiome_output, "Final_prediction.txt")
path_Skin_standard_predict   <- file.path(Skin_standard_output, "Final_prediction.txt")
path_Nasal_predict           <- file.path(Nasal_output, "Final_prediction.txt")
path_Oral_predict            <- file.path(Oral_output, "Final_prediction.txt")

path_Simul_1_50_predict <- file.path(path_Simul_1_50_output, "Final_prediction.txt")
path_Simul_1_30_predict <- file.path(path_Simul_1_30_output, "Final_prediction.txt")
path_Simul_1_10_predict <- file.path(path_Simul_1_10_output, "Final_prediction.txt")

path_Simul_0.5_50_predict <- file.path(path_Simul_0.5_50_output, "Final_prediction.txt")
path_Simul_0.5_30_predict <- file.path(path_Simul_0.5_30_output, "Final_prediction.txt")
path_Simul_0.5_10_predict <- file.path(path_Simul_0.5_10_output, "Final_prediction.txt")

path_Simul_0.1_50_predict <- file.path(path_Simul_0.1_50_output, "Final_prediction.txt")
path_Simul_0.1_30_predict <- file.path(path_Simul_0.1_30_output, "Final_prediction.txt")
path_Simul_0.1_10_predict <- file.path(path_Simul_0.1_10_output, "Final_prediction.txt")

read_predict <- function(path,
                         status = "Contaminant",
                         status_col = "contamination_status",
                         id_col = "Taxid") {
  x <- read.table(path, sep = "\t", header = TRUE,
                  check.names = FALSE, stringsAsFactors = FALSE)
  as.character(x[x[[status_col]] == status, id_col])
}

MI_predict              <- read_predict(path_MI_predict)
Skin_microbiome_predict <- read_predict(path_Skin_microbiome_predict)
Skin_standard_predict   <- read_predict(path_Skin_standard_predict)
Nasal_predict           <- read_predict(path_Nasal_predict)
Oral_predict            <- read_predict(path_Oral_predict)

Simul_1_50_predict <- read_predict(path_Simul_1_50_predict)
Simul_1_30_predict <- read_predict(path_Simul_1_30_predict)
Simul_1_10_predict <- read_predict(path_Simul_1_10_predict)

Simul_0.5_50_predict <- read_predict(path_Simul_0.5_50_predict)
Simul_0.5_30_predict <- read_predict(path_Simul_0.5_30_predict)
Simul_0.5_10_predict <- read_predict(path_Simul_0.5_10_predict)

Simul_0.1_50_predict <- read_predict(path_Simul_0.1_50_predict)
Simul_0.1_30_predict <- read_predict(path_Simul_0.1_30_predict)
Simul_0.1_10_predict <- read_predict(path_Simul_0.1_10_predict)

# ------------------------------------------------------------
# 6) Blacklist predict
# ------------------------------------------------------------
blacklist_path <- file.path(BASE_INPUT, "BlackWhite_list", "blacklist.csv")
assert_exists(blacklist_path)

blacklist <- read.csv(blacklist_path, sep = "\t")$blacklist %>%
  as.character() %>% as.vector()

path_MI_mt              <- file.path(MI_output, "kraken_filtered_matrix.txt")
path_Skin_microbiome_mt <- file.path(Skin_microbiome_output, "kraken_filtered_matrix.txt")
path_Skin_standard_mt   <- file.path(Skin_standard_output, "kraken_filtered_matrix.txt")
path_Nasal_mt           <- file.path(Nasal_output, "kraken_filtered_matrix.txt")
path_Oral_mt            <- file.path(Oral_output, "kraken_filtered_matrix.txt")

path_Simul_1_50_mt <- file.path(path_Simul_1_50_output, "kraken_filtered_matrix.txt")
path_Simul_1_30_mt <- file.path(path_Simul_1_30_output, "kraken_filtered_matrix.txt")
path_Simul_1_10_mt <- file.path(path_Simul_1_10_output, "kraken_filtered_matrix.txt")

path_Simul_0.5_50_mt <- file.path(path_Simul_0.5_50_output, "kraken_filtered_matrix.txt")
path_Simul_0.5_30_mt <- file.path(path_Simul_0.5_30_output, "kraken_filtered_matrix.txt")
path_Simul_0.5_10_mt <- file.path(path_Simul_0.5_10_output, "kraken_filtered_matrix.txt")

path_Simul_0.1_50_mt <- file.path(path_Simul_0.1_50_output, "kraken_filtered_matrix.txt")
path_Simul_0.1_30_mt <- file.path(path_Simul_0.1_30_output, "kraken_filtered_matrix.txt")
path_Simul_0.1_10_mt <- file.path(path_Simul_0.1_10_output, "kraken_filtered_matrix.txt")

MI_all_taxid <- read.table(path_MI_mt) %>% rownames()

black_predict <- function(path, blacklist, id_col = "Taxid") {
  if (!file.exists(path)) {
    warning(sprintf("Missing file: %s", path))
    return(character())
  }
  x <- read.table(path, sep = "\t", header = TRUE,
                  check.names = FALSE, stringsAsFactors = FALSE)
  
  ids <-
    if (id_col %in% names(x)) x[[id_col]]
  else if (!is.null(rownames(x)) &&
           !identical(rownames(x), as.character(seq_len(nrow(x))))) rownames(x)
  else x[[1]]
  
  intersect(blacklist, as.character(ids))
}

MI_black_predict              <- black_predict(path_MI_mt, blacklist)
skin_microbiome_black_predict <- black_predict(path_Skin_microbiome_mt, blacklist)
skin_standard_black_predict   <- black_predict(path_Skin_standard_mt, blacklist)
nasal_black_predict           <- black_predict(path_Nasal_mt, blacklist)
oral_black_predict            <- black_predict(path_Oral_mt, blacklist)

simul_1_50_black_predict <- black_predict(path_Simul_1_50_mt, blacklist)
simul_1_30_black_predict <- black_predict(path_Simul_1_30_mt, blacklist)
simul_1_10_black_predict <- black_predict(path_Simul_1_10_mt, blacklist)

simul_0.5_50_black_predict <- black_predict(path_Simul_0.5_50_mt, blacklist)
simul_0.5_30_black_predict <- black_predict(path_Simul_0.5_30_mt, blacklist)
simul_0.5_10_black_predict <- black_predict(path_Simul_0.5_10_mt, blacklist)

simul_0.1_50_black_predict <- black_predict(path_Simul_0.1_50_mt, blacklist)
simul_0.1_30_black_predict <- black_predict(path_Simul_0.1_30_mt, blacklist)
simul_0.1_10_black_predict <- black_predict(path_Simul_0.1_10_mt, blacklist)

# ------------------------------------------------------------
# 7) Bracken -> make_matrix
# ------------------------------------------------------------
standard_bracken_path   <- file.path(Skin_standard_output, "Bracken_dir")
microbiome_bracken_path <- file.path(Skin_microbiome_output, "Bracken_dir")
MI_bracken_path         <- file.path(MI_output, "Bracken_dir")
nasal_bracken_path      <- file.path(Nasal_output, "Bracken_dir")
oral_bracken_path       <- file.path(Oral_output, "Bracken_dir")

# negative bracken (staged)
MI_negative_bracken_path         <- file.path(BASE_INPUT, "Negative", "MI_blank", "Bracken_dir")
microbiome_negative_bracken_path <- file.path(BASE_INPUT, "Negative", "Skin_microbiome", "Bracken_dir")
standard_negative_bracken_path   <- file.path(BASE_INPUT, "Negative", "Skin_standard", "Bracken_dir")
Oral_negative_bracken_path       <- file.path(BASE_INPUT, "Negative", "Oral", "Bracken_dir")
Nasal_negative_bracken_path      <- file.path(BASE_INPUT, "Negative", "Nasal", "Bracken_dir")

# 파일명 -> 샘플ID (예: "ERR4701466_..." -> "ERR4701466")
sample_id_from_filename <- function(fname) sub("_.*$", "", fname)

make_matrix <- function(input_path,
                        pattern      = c("report$", "_kraken_report_bracken_species$"),
                        exclude_taxa = c("9606"),
                        use_rank     = "S") {
  
  pattern_regex <- if (length(pattern) > 1) paste0("(", paste(pattern, collapse = "|"), ")") else pattern
  
  files <- list.files(input_path, pattern = pattern_regex, full.names = TRUE)
  if (length(files) == 0) {
    stop(sprintf("No files matching '%s' in %s", pattern_regex, input_path))
  }
  
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
        progress = FALSE, show_col_types = FALSE
      ),
      error = function(e) { warning(sprintf("Read error: %s (%s)", f, e$message)); return(NULL) }
    )
    if (is.null(df) || nrow(df) == 0) next
    
    df$rank  <- trimws(df$rank)
    df$taxid <- trimws(df$taxid)
    
    df <- df[df$rank == use_rank & !(df$taxid %in% exclude_taxa), c("taxid","reads")]
    if (nrow(df) == 0) next
    
    names(df)[2] <- "count"
    df <- stats::aggregate(count ~ taxid, data = df, sum, na.rm = TRUE)
    
    v <- df$count
    names(v) <- df$taxid
    vecs[[sid]] <- v
  }
  
  if (length(vecs) == 0) {
    stop("No usable species rows after filtering (check rank / exclude_taxa).")
  }
  
  all_taxids <- sort(unique(unlist(lapply(vecs, names))), na.last = NA)
  sample_ids <- sort(names(vecs))
  
  mat_list <- lapply(sample_ids, function(sid) {
    v <- vecs[[sid]][all_taxids]
    v[is.na(v)] <- 0
    as.numeric(v)
  })
  
  count_mat <- do.call(cbind, mat_list)
  rownames(count_mat) <- all_taxids
  colnames(count_mat) <- sample_ids
  storage.mode(count_mat) <- "double"
  
  col_sums <- colSums(count_mat, na.rm = TRUE)
  denom <- ifelse(col_sums > 0, col_sums, NA_real_)
  rel_mat <- sweep(count_mat, 2, denom, "/")
  rel_mat[is.na(rel_mat)] <- 0
  
  rel_mat
}

# matrices
MI_matrix         <- make_matrix(MI_bracken_path)
oral_matrix       <- make_matrix(oral_bracken_path)
standard_matrix   <- make_matrix(standard_bracken_path)
microbiome_matrix <- make_matrix(microbiome_bracken_path)
Nasal_matrix      <- make_matrix(nasal_bracken_path)

MI_negative_matrix <- make_matrix(MI_negative_bracken_path)


cat("\n[OK] Script finished successfully.\n")
