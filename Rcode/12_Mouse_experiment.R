############################################################
## 0) Libraries
############################################################
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(ggsci)
  library(ggpubr)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(NetCoMi)
})

############################################################
## 1) TaxID -> scientific name (use prebuilt CSV)
############################################################
library(dplyr)
library(readr)

tax_name_file <- file.path(BASE_INPUT, "Taxa_name", "names_dmp_scientific.csv")

if (!file.exists(tax_name_file)) {
  stop("Not found: ", tax_name_file)
}

sci_name_df <- readr::read_csv(
  tax_name_file,
  show_col_types = FALSE,
  progress = FALSE
)

# Your CSV has an extra first column (row index) like ",taxid,scientific_name"
# Drop it safely if present.
first_col <- names(sci_name_df)[1]
if (first_col %in% c("", "...1", "X1")) {
  sci_name_df <- sci_name_df %>% select(-1)
}

# Validate required columns
if (!all(c("taxid", "scientific_name") %in% names(sci_name_df))) {
  stop("CSV must contain columns: taxid, scientific_name")
}

# Clean types/whitespace and keep unique mapping
sci_name_df <- sci_name_df %>%
  mutate(
    taxid = trimws(as.character(taxid)),
    scientific_name = trimws(as.character(scientific_name))
  ) %>%
  filter(taxid != "", !is.na(taxid)) %>%
  distinct(taxid, .keep_all = TRUE)


############################################################
## 2) make_raw_rel_matrix
############################################################
sample_id_from_filename <- function(fname) sub("_.*$", "", fname)

make_raw_rel_matrix <- function(input_path,
                                pattern      = c("report$", "_kraken_report_bracken_species$"),
                                exclude_taxa = c("9606"),
                                use_rank     = "S") {
  pattern_regex <- if (length(pattern) > 1) paste0("(", paste(pattern, collapse="|"), ")") else pattern
  files <- list.files(input_path, pattern = pattern_regex, full.names = TRUE)
  if (length(files) == 0) stop("No matching files found")
  
  vecs <- list()
  for (f in files) {
    sid <- sample_id_from_filename(basename(f))
    
    df <- readr::read_tsv(
      f,
      col_names = c("percent", "reads", "assigned", "rank", "taxid", "name"),
      col_types = cols(
        percent  = col_double(),
        reads    = col_double(),
        assigned = col_double(),
        rank     = col_character(),
        taxid    = col_character(),
        name     = col_character()
      ),
      progress = FALSE,
      show_col_types = FALSE
    )
    
    df$rank  <- trimws(df$rank)
    df$taxid <- trimws(df$taxid)
    
    df <- df[df$rank == use_rank & !(df$taxid %in% exclude_taxa), c("taxid", "reads")]
    if (nrow(df) == 0) next
    
    names(df)[2] <- "count"
    df <- aggregate(count ~ taxid, data = df, sum)
    
    v <- df$count
    names(v) <- df$taxid
    vecs[[sid]] <- v
  }
  
  all_taxids <- sort(unique(unlist(lapply(vecs, names))))
  sample_ids <- sort(names(vecs))
  
  raw_list <- lapply(sample_ids, function(sid) {
    v <- vecs[[sid]][all_taxids]
    v[is.na(v)] <- 0
    as.numeric(v)
  })
  
  raw_mat <- do.call(cbind, raw_list)
  rownames(raw_mat) <- all_taxids
  colnames(raw_mat) <- sample_ids
  
  rel_mat <- sweep(raw_mat, 2, colSums(raw_mat), "/")
  rel_mat[is.na(rel_mat)] <- 0
  
  list(raw = raw_mat, rel = rel_mat)
}

############################################################
## 3) Load fecal/kit/BALF
##    NOTE: BASE_INPUT must exist in your environment.
############################################################
fecal_kit_balf_mt <- make_raw_rel_matrix(
  file.path(BASE_INPUT, "Mouse_experiment", "Fecal_BALF_kit_Bracken_files")
)
rel3 <- fecal_kit_balf_mt$rel

############################################################
## 4) Define enriched taxa
############################################################
fecal_only_taxa <- rownames(rel3)[rel3[, "kit"] < 1e-4 & rel3[, "fecal"] >= 1e-3]
lung_only_taxa  <- rownames(rel3)[rel3[, "BALF"] >= 1e-3]
kit_only_taxa   <- rownames(rel3)[rel3[, "kit"] >= 1e-4]

############################################################
## 5) Load lung mixed dataset + filter (mean >= 0.001)
############################################################
lung_mt <- make_raw_rel_matrix(
  file.path(BASE_INPUT, "Mouse_experiment", "BALF_FECAL_Mixed_Bracken_files")
)
rel_mat <- lung_mt$rel
raw_mat <- lung_mt$raw

row_means    <- rowMeans(rel_mat)
rel_filtered <- rel_mat[row_means >= 0.001, ]
raw_filtered <- raw_mat[row_means >= 0.001, ]

############################################################
## 6) NetCoMi network + save edge list
############################################################
raw_for_net <- t(raw_filtered)

net_obj <- netConstruct(
  raw_for_net,
  measure     = "pearson",
  normMethod  = "clr",
  zeroMethod  = "multRepl",
  sparsMethod = "threshold",
  thresh      = 0.3,
  verbose     = 3
)

net <- netAnalyze(
  net_obj,
  clustMethod = "cluster_louvain",
  hubPar      = "degree"
)

output_edge_file <- file.path(BASE_INPUT, "Mouse_experiment", "Network", "edge.tsv")
dir.create(dirname(output_edge_file), recursive = TRUE, showWarnings = FALSE)
write.table(net_obj$edgelist1, file = output_edge_file, sep = "\t",
            row.names = FALSE, quote = FALSE)

############################################################
## 7) Safe readers (partition/edge)
############################################################
.clean_names <- function(nm) {
  nm <- as.character(nm)
  nm[is.na(nm)] <- ""
  nm <- trimws(nm)
  nm <- sub("^\ufeff", "", nm)
  nm <- gsub('^["\']+|["\']+$', "", nm)
  nm
}
.unquote <- function(x) gsub('^["\']+|["\']+$', "", as.character(x))

.read_table_raw <- function(path) {
  df <- tryCatch(
    read.table(path, sep = "\t", header = TRUE, check.names = FALSE,
               stringsAsFactors = FALSE, quote = "", comment.char = "", fill = TRUE),
    error = function(e) data.frame()
  )
  if (ncol(df) == 0) {
    df <- read.table(path, sep = "\t", header = FALSE, check.names = FALSE,
                     stringsAsFactors = FALSE, quote = "", comment.char = "", fill = TRUE)
  }
  df
}

read_partition_safe <- function(path) {
  df <- .read_table_raw(path)
  names(df) <- .clean_names(names(df))
  ti <- if ("taxid" %in% names(df)) which(names(df) == "taxid")[1] else 1
  ci <- if ("cluster" %in% names(df)) which(names(df) == "cluster")[1] else 2
  data.frame(
    taxid   = .unquote(df[[ti]]),
    cluster = .unquote(df[[ci]]),
    stringsAsFactors = FALSE
  )
}

read_edge_safe <- function(path) {
  df <- .read_table_raw(path)
  names(df) <- .clean_names(names(df))
  
  pick <- function(cands, default) {
    x <- which(tolower(names(df)) %in% tolower(cands))
    if (length(x) == 0) default else x[1]
  }
  
  v1i   <- pick(c("v1", "from", "source"), 1)
  v2i   <- pick(c("v2", "to", "target"),   2)
  assoi <- pick(c("asso", "association", "corr", "weight"), NA)
  adjai <- pick(c("adja", "weight", "score"), assoi)
  
  data.frame(
    v1   = .unquote(df[[v1i]]),
    v2   = .unquote(df[[v2i]]),
    asso = if (!is.na(assoi)) suppressWarnings(as.numeric(df[[assoi]])) else 1,
    adja = if (!is.na(adjai)) suppressWarnings(as.numeric(df[[adjai]])) else 1,
    stringsAsFactors = FALSE
  )
}

get_comm_palette <- function(partition_path) {
  pt  <- read_partition_safe(partition_path)
  cls <- pt %>%
    count(cluster) %>%
    filter(n >= 10) %>%
    pull(cluster) %>%
    unique() %>%
    sort()
  
  n <- max(length(cls), 1)
  cols <- pal_nejm("default")(n)
  names(cols) <- cls
  list(colors = cols, classes = cls, labels = paste0("Community ", seq_len(n)))
}

############################################################
## 8) ANI map (mean_conANI_table.tsv)
############################################################
clean_taxid <- function(x) {
  x <- gsub('^["\']+|["\']+$', "", x)
  trimws(x)
}

mean_ani_file <- file.path(BASE_INPUT, "Mouse_experiment", "Metacontam_out", "mean_conANI_table.tsv")
mean_ani_df <- read.table(mean_ani_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
mean_ani_df$taxid <- clean_taxid(mean_ani_df$taxid)

ani_map <- mean_ani_df$mean_conANI
names(ani_map) <- mean_ani_df$taxid

############################################################
## 9) Fecal vs BALF stacked bar (topN + Other) : SAFE ALL-IN-ONE
############################################################
if (!exists("topN_map")) topN_map <- c(fecal = 5, BALF = 5)
stopifnot(exists("rel3"))
stopifnot(all(c("fecal", "BALF") %in% colnames(rel3)))

bf_mat <- rel3[, c("fecal", "BALF"), drop = FALSE]

df_long <- as.data.frame(bf_mat) %>%
  tibble::rownames_to_column("taxid") %>%
  pivot_longer(-taxid, names_to = "sample", values_to = "rel") %>%
  mutate(taxid = as.character(taxid), sample = as.character(sample))

# Join scientific names if available; otherwise fall back to taxid
if (exists("sci_name_df")) {
  df_long <- df_long %>%
    left_join(
      sci_name_df %>%
        mutate(taxid = as.character(taxid)) %>%
        distinct(taxid, scientific_name),
      by = "taxid"
    )
} else {
  df_long$scientific_name <- NA_character_
}

df_long <- df_long %>%
  mutate(species = ifelse(is.na(scientific_name) | scientific_name == "", taxid, scientific_name))

# If kit_only_taxa exists, remove them from fecal only
if (exists("kit_only_taxa")) {
  df_long <- df_long %>%
    filter(!(sample == "fecal" & taxid %in% as.character(kit_only_taxa)))
}

df_long$sample <- factor(df_long$sample, levels = c("fecal", "BALF"))

df_sum <- df_long %>%
  group_by(sample, species) %>%
  summarise(rel = sum(rel, na.rm = TRUE), .groups = "drop") %>%
  filter(rel > 0)

df_top <- df_sum %>%
  group_by(sample) %>%
  arrange(desc(rel), .by_group = TRUE) %>%
  mutate(
    topN = unname(topN_map[as.character(sample)]),
    keep = row_number() <= dplyr::first(topN)
  ) %>%
  ungroup()

df_plot <- df_top %>%
  mutate(species2 = ifelse(keep, species, "Other")) %>%
  group_by(sample, species2) %>%
  summarise(rel = sum(rel), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(rel_prop = rel / sum(rel)) %>%
  ungroup()

species_levels <- df_plot %>% distinct(species2) %>% pull(species2)
species_levels <- c(setdiff(species_levels, "Other"), "Other")
df_plot$species2 <- factor(df_plot$species2, levels = species_levels)

non_other <- setdiff(species_levels, "Other")
n_non <- length(non_other)
primer_cols <- pal_primer("mark17")(17)
col_vec <- c(setNames(primer_cols[seq_len(n_non)], non_other), Other = "grey70")

topN_fecal <- topN_map[["fecal"]]
topN_balf  <- topN_map[["BALF"]]

# --- Plot (kept as-is to match your original output as closely as possible) ---
Mouse_fecal_balf_bar <- ggplot(df_plot, aes(x = sample, y = rel_prop, fill = species2, order = -rel_prop)) +
  geom_col(width = 0.70, color = "white", linewidth = 0.25) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1.01),
    expand = expansion(mult = c(0, 0))   # remove extra space at bottom/top of y-axis
  ) +
  scale_x_discrete(
    labels = c(fecal = "Fecal", BALF = "BALF"),
    expand = expansion(add = c(0.65, 0.65))
  ) +
  scale_fill_manual(values = col_vec, drop = FALSE) +
  labs(x = NULL, y = "Proportion (100%)", fill = "Species",
       title = paste0("Fecal (Top ", topN_fecal, ") vs BALF (Top ", topN_balf, ") + Other")) +
  theme_pubr() +
  theme(
    plot.title  = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 15),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "right",
    plot.margin        = margin(4, 4, 4, 4)
  ) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
  coord_cartesian(clip = "off")

Mouse_fecal_balf_bar

############################################################
## 10) Colors + ANI box (2 groups)
############################################################
bmj_cols    <- pal_bmj("default")(8)
cosmic_cols <- pal_cosmic("hallmarks_light")(8)

fecal_col <- bmj_cols[6]
lung_col  <- cosmic_cols[2]

plot_ani_two_groups <- function(ani_map, fecal_taxa, balf_taxa,
                                title = "ANI Comparison: BALF vs Fecal",
                                jitter_size = 3, jitter_alpha = 1.0) {
  ani_df <- data.frame(taxid = names(ani_map), ANI = as.numeric(ani_map), stringsAsFactors = FALSE)
  ani_df$group <- dplyr::case_when(
    ani_df$taxid %in% balf_taxa  ~ "BALF",
    ani_df$taxid %in% fecal_taxa ~ "Fecal",
    TRUE ~ NA_character_
  )
  ani_df <- ani_df[!is.na(ani_df$group) & !is.na(ani_df$ANI), ]
  ani_df$group <- factor(ani_df$group, levels = c("Fecal", "BALF"))
  
  p_FB <- wilcox.test(
    ani_df$ANI[ani_df$group == "Fecal"],
    ani_df$ANI[ani_df$group == "BALF"],
    alternative = "greater", exact = FALSE
  )$p.value
  message("\n[Wilcoxon one-sided: Fecal > BALF] p = ", signif(p_FB, 5))
  
  col_map <- c(BALF = lung_col, Fecal = fecal_col)
  set.seed(3)
  
  ymax <- max(ani_df$ANI); ymin <- min(ani_df$ANI); span <- ymax - ymin
  h1 <- ymax + span * 0.04; h1_top <- h1 + span * 0.02
  label_offset <- span * 0.015
  ylim_top <- h1_top + span * 0.03; ylim_bottom <- ymin - span * 0.01
  
  ggplot(ani_df, aes(x = group, y = ANI, color = group)) +
    geom_boxplot(fill = NA, size = 1.4, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = jitter_size, alpha = jitter_alpha) +
    scale_color_manual(values = col_map) +
    ggtitle(title) +
    geom_segment(aes(x = 1, xend = 1, y = h1, yend = h1_top), inherit.aes = FALSE, color = "black", size = 0.8) +
    geom_segment(aes(x = 2, xend = 2, y = h1, yend = h1_top), inherit.aes = FALSE, color = "black", size = 0.8) +
    geom_segment(aes(x = 1, xend = 2, y = h1_top, yend = h1_top), inherit.aes = FALSE, color = "black", size = 0.8) +
    annotate("text", x = 1.5, y = h1_top + label_offset, label = paste0("p = ", signif(p_FB, 3)), size = 4) +
    scale_y_continuous(limits = c(ylim_bottom, ylim_top), breaks = seq(0.96, 1.00, 0.01)) +
    theme_pubr() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      axis.title = element_blank(),
      axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 15),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      plot.margin = margin(6, 6, 6, 6),
      aspect.ratio = 1.5
    )
}

set.seed(1234)
Mouse_ANI_2groups <- plot_ani_two_groups(
  ani_map,
  fecal_taxa = fecal_only_taxa,
  balf_taxa  = lung_only_taxa,
  title      = "ANI: BALF vs Fecal"
)
Mouse_ANI_2groups

############################################################
## 11) Network plot: dual highlight + species label
############################################################
plot_network_dual_highlight_specieslabel <- function(edgelist,
                                                     group1_taxa, group2_taxa,
                                                     group1_name, group2_name,
                                                     group1_color, group2_color,
                                                     sci_name_df,
                                                     min_weight = 0, min_degree = 0,
                                                     title = "",
                                                     label_mode = c("both", "group1", "group2", "none"),
                                                     label_size = 2.6) {
  label_mode <- match.arg(label_mode)
  
  ef <- edgelist %>%
    filter(!is.na(adja), adja > 0, adja >= min_weight) %>%
    select(from = v1, to = v2)
  
  ef$from <- as.character(ef$from)
  ef$to   <- as.character(ef$to)
  group1_taxa <- as.character(group1_taxa)
  group2_taxa <- as.character(group2_taxa)
  
  g_tmp <- igraph::graph_from_data_frame(ef, directed = FALSE)
  keep_nodes <- names(igraph::degree(g_tmp)[igraph::degree(g_tmp) >= min_degree])
  
  ef_filt <- ef %>% filter(from %in% keep_nodes, to %in% keep_nodes)
  
  g <- igraph::graph_from_data_frame(
    ef_filt, directed = FALSE,
    vertices = data.frame(name = unique(c(ef_filt$from, ef_filt$to)), stringsAsFactors = FALSE)
  )
  
  set.seed(1234)
  lay <- ggraph::create_layout(g, layout = "stress")
  
  lay$group <- case_when(
    lay$name %in% group1_taxa ~ group1_name,
    lay$name %in% group2_taxa ~ group2_name,
    TRUE ~ "Other"
  )
  lay$group <- factor(lay$group, levels = c(group1_name, group2_name, "Other"))
  
  sci_map <- sci_name_df %>%
    mutate(taxid = as.character(taxid)) %>%
    distinct(taxid, scientific_name)
  
  lay <- lay %>%
    left_join(sci_map, by = c("name" = "taxid")) %>%
    mutate(
      label = if_else(is.na(scientific_name) | scientific_name == "", name, scientific_name),
      label_show = case_when(
        label_mode == "both"   ~ group %in% c(group1_name, group2_name),
        label_mode == "group1" ~ group == group1_name,
        label_mode == "group2" ~ group == group2_name,
        TRUE ~ FALSE
      ),
      label_plot = if_else(label_show, label, "")
    )
  
  pal <- c(
    setNames(as.character(group1_color)[1], group1_name),
    setNames(as.character(group2_color)[1], group2_name),
    Other = "grey90"
  )
  
  ggraph(lay) +
    geom_edge_link(color = "grey60", alpha = 0.7) +
    geom_node_point(aes(fill = group), shape = 21, size = 7, colour = "grey20", stroke = 0.6) +
    scale_fill_manual(values = pal, drop = FALSE) +
    geom_node_text(aes(label = label_plot), size = label_size, repel = TRUE) +
    ggtitle(title) +
    theme_void() +
    theme(
      aspect.ratio = 1,
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

min_degree <- 2
weight     <- 0.45

Mouse_fecal_balf_network <- plot_network_dual_highlight_specieslabel(
  net_obj$edgelist1,
  fecal_only_taxa, lung_only_taxa,
  "Fecal enriched", "BALF enriched",
  fecal_col, lung_col,
  sci_name_df,
  min_weight = weight,
  min_degree = min_degree,
  title      = "",
  label_mode = "group2",
  label_size = 2.6
)

Mouse_fecal_balf_network


Supp_Fig10b <- Mouse_fecal_balf_bar
Supp_Fig10c <- Mouse_fecal_balf_network
Supp_Fig10d <- Mouse_ANI_2groups
