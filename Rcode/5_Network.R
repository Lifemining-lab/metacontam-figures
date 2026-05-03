# ================================
# Network plotting
#   - Keep existing path_* variables
#   - Replace highlight_list -> GT_species
#   - Restore original figure behavior
# ================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(ggplot2)
  library(ggsci)
  library(patchwork)
})

# ----------------
# Paths expected somewhere above or here:
# path_MI_partition              <- "..."
# path_MI_edge                   <- "..."
# path_Skin_microbiome_partition <- "..."
# path_Skin_microbiome_edge      <- "..."
# path_Skin_standard_partition   <- "..."
# path_Skin_standard_edge        <- "..."
# path_Nasal_partition           <- "..."
# path_Nasal_edge                <- "..."
# path_Oral_partition            <- "..."
# path_Oral_edge                 <- "..."
# ----------------

# ================================
# 1) Load GT_species
# ================================
BASE_INPUT <- Sys.getenv("BASE_INPUT", unset = "/media/junwoojo/18T/Submission/Rcode/inputs")
GT_path <- file.path(BASE_INPUT, "config", "GT_species.tsv")

if (!file.exists(GT_path)) stop("Missing GT_species.tsv: ", GT_path)

GT_df <- read.table(
  GT_path,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

GT_df$group <- as.character(GT_df$group)
GT_df$taxid <- as.character(GT_df$taxid)

as_chr_vec <- function(x) {
  if (is.null(x)) return(character(0))
  x <- unlist(x, recursive = TRUE, use.names = FALSE)
  x <- as.character(x)
  x[!is.na(x)]
}

GT_species <- lapply(split(GT_df$taxid, GT_df$group), as_chr_vec)

.valid_clusters <- function(cluster_vec, min_n = 10) {
  cluster_vec <- as_chr_vec(cluster_vec)
  tb <- table(cluster_vec)
  sort(names(tb)[tb >= min_n])
}









 
# ================================
# 2) Blacklist
# ================================
blacklist <- read.csv(
  file.path(BASE_INPUT, "BlackWhite_list", "blacklist.csv"),
  sep = "\t",
  header = TRUE
)$blacklist %>% as_chr_vec()

# ================================
# 3) Safe readers
# ================================
.clean_names <- function(nm) {
  nm <- as.character(nm)
  nm[is.na(nm)] <- ""
  nm <- trimws(nm)
  nm <- sub("^\ufeff", "", nm)
  nm <- gsub('^["\']+|["\']+$', "", nm)
  nm
}

.unquote <- function(x) gsub('^["\']+|["\']+$', "", as_chr_vec(x))

.read_table_raw <- function(path) {
  df <- tryCatch(
    read.table(
      path,
      sep = "\t",
      header = TRUE,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = "",
      fill = TRUE
    ),
    error = function(e) data.frame()
  )
  
  if (ncol(df) == 0) {
    df <- read.table(
      path,
      sep = "\t",
      header = FALSE,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = "",
      fill = TRUE
    )
  }
  
  df
}

read_partition_safe <- function(path) {
  df <- .read_table_raw(path)
  names(df) <- .clean_names(names(df))
  
  ti <- if ("taxid" %in% names(df)) which(names(df) == "taxid")[1] else {
    x <- which(tolower(names(df)) %in% c("taxid", "tax_id", "id"))
    if (length(x) == 0) 1 else x[1]
  }
  
  ci <- if ("cluster" %in% names(df)) which(names(df) == "cluster")[1] else {
    x <- which(grepl("cluster|community", tolower(names(df))))
    if (length(x) == 0) {
      if (ncol(df) >= 2) 2 else NA_integer_
    } else x[1]
  }
  
  if (is.na(ci)) stop(sprintf("'%s': cannot find a cluster column.", path))
  
  out <- data.frame(
    taxid   = .unquote(df[[ti]]),
    cluster = .unquote(df[[ci]]),
    stringsAsFactors = FALSE
  )
  
  out$taxid   <- as_chr_vec(out$taxid)
  out$cluster <- as_chr_vec(out$cluster)
  out
}

read_edge_safe <- function(path) {
  df <- .read_table_raw(path)
  names(df) <- .clean_names(names(df))
  
  has_index6 <- ncol(df) >= 6 &&
    suppressWarnings(all(!is.na(as.numeric(df[[1]])))) &&
    suppressWarnings(all(!is.na(as.numeric(df[[4]])))) &&
    mean(suppressWarnings(as.numeric(df[[4]]) >= -1 & as.numeric(df[[4]]) <= 1), na.rm = TRUE) > 0.9
  
  if (has_index6) {
    v1   <- .unquote(df[[2]])
    v2   <- .unquote(df[[3]])
    asso <- suppressWarnings(as.numeric(df[[4]]))
    adja <- suppressWarnings(as.numeric(df[[6]]))
    return(data.frame(v1 = v1, v2 = v2, asso = asso, adja = adja, stringsAsFactors = FALSE))
  }
  
  pick <- function(cands, default) {
    x <- which(tolower(names(df)) %in% tolower(cands))
    if (length(x) == 0) default else x[1]
  }
  
  v1i   <- pick(c("v1", "from", "source", "taxid1", "node1"), 1)
  v2i   <- pick(c("v2", "to", "target", "taxid2", "node2"), if (ncol(df) >= 2) 2 else 1)
  assoi <- pick(c("asso", "association", "corr", "weight", "score", "w"), NA)
  adjai <- pick(c("adja", "weight", "score", "w"), if (is.na(assoi)) NA else assoi)
  
  out <- data.frame(
    v1   = .unquote(df[[v1i]]),
    v2   = .unquote(df[[v2i]]),
    asso = if (!is.na(assoi)) suppressWarnings(as.numeric(df[[assoi]])) else 1,
    adja = if (!is.na(adjai)) suppressWarnings(as.numeric(df[[adjai]])) else 1,
    stringsAsFactors = FALSE
  )
  
  out$v1 <- as_chr_vec(out$v1)
  out$v2 <- as_chr_vec(out$v2)
  out
}

# ================================
# 4) Helpers
# ================================
.placeholder_plot <- function(title_txt, subtitle_txt = "No edges among communities with ≥10 nodes") {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = subtitle_txt, size = 4) +
    ggtitle(title_txt) +
    theme_void() +
    theme(
      plot.title   = element_text(size = 10, face = "bold", hjust = 0.5),
      plot.margin  = margin(5, 5, 5, 5),
      aspect.ratio = 1
    )
}

get_comm_palette <- function(partition_path,
                             swap_comm_indices = NULL,
                             override_colors_by_comm = NULL) {
  pt <- read_partition_safe(partition_path)
  cls <- .valid_clusters(pt$cluster, min_n = 10)
  
  n <- length(cls)
  if (n < 1) n <- 1
  
  cols <- pal_nejm("default")(n)
  
  if (!is.null(swap_comm_indices) && length(swap_comm_indices) == 2) {
    i <- swap_comm_indices[1]
    j <- swap_comm_indices[2]
    if (i >= 1 && i <= n && j >= 1 && j <= n) {
      tmp <- cols[i]
      cols[i] <- cols[j]
      cols[j] <- tmp
    }
  }
  
  if (!is.null(override_colors_by_comm)) {
    for (k in names(override_colors_by_comm)) {
      idx <- suppressWarnings(as.integer(k))
      if (!is.na(idx) && idx >= 1 && idx <= n) cols[idx] <- override_colors_by_comm[[k]]
    }
  }
  
  names(cols) <- as.character(seq_len(n))
  labels <- paste0("Community ", seq_len(n))
  
  list(
    colors  = setNames(cols, as.character(cls)),
    labels  = setNames(labels, as.character(cls)),
    classes = cls
  )
}

print_comm_index_map <- function(partition_path) {
  pt  <- read_partition_safe(partition_path)
  cls <- .valid_clusters(pt$cluster, min_n = 10)
  
  cat("Community index -> cluster value\n")
  for (i in seq_along(cls)) cat(sprintf("%2d -> %s\n", i, cls[i]))
}

# ================================
# 5) Contaminant network plot
# ================================
plot_network_pair <- function(part_p, edge_p, highlight_ids, title,
                              swap_comm_indices = NULL, override_colors_by_comm = NULL) {
  highlight_ids <- as_chr_vec(highlight_ids)
  
  pt <- read_partition_safe(part_p)
  valid <- .valid_clusters(pt$cluster, min_n = 10)
  
  if (length(valid) == 0) {
    return(.placeholder_plot(title, "No communities with \u226510 nodes"))
  }
  
  Pt <- pt[pt$cluster %in% valid, , drop = FALSE]
  
  et <- read_edge_safe(edge_p)
  et <- et[is.na(et$asso) | et$asso > 0, , drop = FALSE]
  
  ef <- et[et$v1 %in% Pt$taxid & et$v2 %in% Pt$taxid, , drop = FALSE]
  ef <- data.frame(
    from = ef$v1,
    to   = ef$v2,
    weight = ef$adja,
    stringsAsFactors = FALSE
  )
  
  nds <- data.frame(
    name = Pt$taxid,
    cluster = factor(Pt$cluster),
    highlight = ifelse(Pt$taxid %in% highlight_ids, "yes", "no"),
    stringsAsFactors = FALSE
  )
  
  g <- as_tbl_graph(graph_from_data_frame(ef, vertices = nds, directed = FALSE))
  if (ecount(g) == 0) return(.placeholder_plot(title))
  
  set.seed(1234)
  coords <- layout_with_graphopt(g, charge = 0.1, niter = 1000)
  lay <- create_layout(g, layout = "manual", x = coords[, 1], y = coords[, 2])
  
  pal <- get_comm_palette(
    part_p,
    swap_comm_indices = swap_comm_indices,
    override_colors_by_comm = override_colors_by_comm
  )
  
  ggraph(lay) +
    geom_edge_link(color = "grey50", width = 0.3, alpha = 0.7, show.legend = FALSE) +
    geom_node_point(
      aes(fill = cluster, shape = highlight, alpha = highlight),
      size = 5, stroke = 0.8, color = "black"
    ) +
    scale_shape_manual(
      name   = "Ground truth",
      breaks = c("yes"),
      values = c(yes = 24, no = 21),
      labels = c("Ground-truth contaminant")
    ) +
    scale_alpha_manual(values = c(no = 0.3, yes = 1), guide = "none") +
    scale_fill_manual(
      name   = "Community",
      values = pal$colors,
      breaks = pal$classes,
      labels = pal$labels[as.character(pal$classes)],
      guide  = guide_legend(
        direction = "vertical",
        title.position = "top",
        override.aes = list(shape = 21, size = 5)
      )
    ) +
    ggtitle(title) +
    theme_void() +
    theme(
      plot.title      = element_text(size = 10, face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.box      = "vertical",
      plot.margin     = margin(0, 0, 0, 0),
      aspect.ratio    = 1
    ) +
    guides(shape = guide_legend(direction = "horizontal", title.position = "top", order = 1))
}

# ================================
# 6) Blacklist network plot
# ================================
plot_network_blacklist <- function(part_p, edge_p, blacklist_ids, title,
                                   swap_comm_indices = NULL, override_colors_by_comm = NULL) {
  blacklist_ids <- as_chr_vec(blacklist_ids)
  
  pt <- read_partition_safe(part_p)
  valid <- .valid_clusters(pt$cluster, min_n = 10)
  
  if (length(valid) == 0) {
    return(.placeholder_plot(title, "No communities with \u226510 nodes"))
  }
  
  Pt <- pt[pt$cluster %in% valid, , drop = FALSE]
  
  et <- read_edge_safe(edge_p)
  et <- et[is.na(et$asso) | et$asso > 0, , drop = FALSE]
  
  ef <- et[et$v1 %in% Pt$taxid & et$v2 %in% Pt$taxid, , drop = FALSE]
  ef <- data.frame(
    from = ef$v1,
    to   = ef$v2,
    weight = ef$adja,
    stringsAsFactors = FALSE
  )
  
  nds <- data.frame(
    name = Pt$taxid,
    cluster = factor(Pt$cluster),
    black = ifelse(Pt$taxid %in% blacklist_ids, "yes", "no"),
    stringsAsFactors = FALSE
  )
  
  g <- as_tbl_graph(graph_from_data_frame(ef, vertices = nds, directed = FALSE))
  if (ecount(g) == 0) return(.placeholder_plot(title))
  
  set.seed(1234)
  coords <- layout_with_graphopt(g, charge = 0.1, niter = 1000)
  lay <- create_layout(g, layout = "manual", x = coords[, 1], y = coords[, 2])
  
  pal <- get_comm_palette(
    part_p,
    swap_comm_indices = swap_comm_indices,
    override_colors_by_comm = override_colors_by_comm
  )
  
  ggraph(lay) +
    geom_edge_link(color = "grey50", width = 0.3, alpha = 0.7, show.legend = FALSE) +
    geom_node_point(
      aes(fill = cluster, shape = black, alpha = black),
      size = 5, stroke = 0.8, color = "black"
    ) +
    scale_shape_manual(
      name   = "Blacklist",
      breaks = c("yes"),
      values = c(yes = 22, no = 21),
      labels = c("Blacklist")
    ) +
    scale_alpha_manual(values = c(no = 0.3, yes = 1), guide = "none") +
    scale_fill_manual(
      name   = "Community",
      values = pal$colors,
      breaks = pal$classes,
      labels = pal$labels[as.character(pal$classes)],
      guide  = guide_legend(
        direction = "vertical",
        title.position = "top",
        override.aes = list(shape = 21, size = 5)
      )
    ) +
    ggtitle(title) +
    theme_void() +
    theme(
      plot.title      = element_text(size = 10, face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.box      = "vertical",
      plot.margin     = margin(0, 0, 0, 0),
      aspect.ratio    = 1
    ) +
    guides(shape = guide_legend(direction = "horizontal", title.position = "top", order = 1))
}

# ================================
# 7) Console report
# ================================
report_comm_summary <- function(part_p, h_ids, core) {
  h_ids <- as_chr_vec(h_ids)
  
  pt <- read_partition_safe(part_p)
  valid <- pt %>%
    group_by(cluster) %>%
    tally() %>%
    filter(n >= 10) %>%
    pull(cluster)
  
  if (length(valid) == 0) {
    cat(sprintf(
      "\n[%s] highlight total = %d (present in data = 0)\n - No communities with >=10 nodes.\n",
      core, length(unique(h_ids))
    ))
    return(invisible(NULL))
  }
  
  P <- pt %>%
    filter(cluster %in% valid) %>%
    transmute(taxid = .unquote(taxid), cluster = .unquote(cluster))
  
  P$taxid   <- as_chr_vec(P$taxid)
  P$cluster <- as_chr_vec(P$cluster)
  
  P$is_black <- P$taxid %in% blacklist
  P$is_high  <- P$taxid %in% h_ids
  
  cls_sorted <- sort(unique(as.character(P$cluster)))
  idx_map <- setNames(seq_along(cls_sorted), as.character(cls_sorted))
  
  total_high         <- length(unique(h_ids))
  present_high_total <- sum(P$is_high)
  
  summaries <- P %>%
    group_by(cluster) %>%
    summarise(
      nodes_in_comm = n(),
      black_in_comm = sum(is_black),
      high_in_comm  = sum(is_high),
      .groups = "drop"
    ) %>%
    filter(black_in_comm > 0) %>%
    mutate(comm_no = idx_map[as.character(cluster)]) %>%
    arrange(comm_no)
  
  cat(sprintf(
    "\n[%s] highlight total = %d (present in data = %d)\n",
    core, total_high, present_high_total
  ))
  
  if (nrow(summaries) == 0) {
    cat(" - No communities (>=10 nodes) contain blacklist taxa.\n")
    return(invisible(summaries))
  }
  
  for (i in seq_len(nrow(summaries))) {
    cl <- summaries$cluster[i]
    
    a <- summaries$high_in_comm[i]
    b <- summaries$nodes_in_comm[i] - a
    
    out <- P %>% filter(cluster != cl)
    c2 <- sum(out$is_high)
    d  <- nrow(out) - c2
    
    pval <- tryCatch(
      fisher.test(matrix(c(a, b, c2, d), nrow = 2), alternative = "greater")$p.value,
      error = function(e) NA_real_
    )
    
    cat(sprintf(
      " - Community %d: nodes=%d; blacklist_in_comm=%d; highlight_in_comm=%d/%d; p=%.3g\n",
      idx_map[as.character(cl)],
      summaries$nodes_in_comm[i],
      summaries$black_in_comm[i],
      a, total_high, pval
    ))
  }
  
  cat("\n")
  invisible(summaries)
}

# ================================
# 8) Pair wrapper
# ================================
draw_pair <- function(part_p, edge_p, h_ids, core,
                      swap_comm_indices = NULL, override_colors_by_comm = NULL) {
  h_ids <- as_chr_vec(h_ids)
  
  report_comm_summary(part_p, h_ids, core)
  
  p1 <- plot_network_blacklist(
    part_p, edge_p, blacklist, paste0(core, " (blacklist)"),
    swap_comm_indices = swap_comm_indices,
    override_colors_by_comm = override_colors_by_comm
  ) + guides(fill = "none")
  
  p2 <- plot_network_pair(
    part_p, edge_p, h_ids, paste0(core, " (contaminant)"),
    swap_comm_indices = swap_comm_indices,
    override_colors_by_comm = override_colors_by_comm
  )
  
  wrap_plots(p1, p2, guides = "keep") + plot_layout(guides = "collect")
}

# ================================
# 9) Build pairs
# ================================
cols_aaas <- pal_aaas()(8)
cols_jco  <- pal_jco()(10)

MI_pair <- draw_pair(
  path_MI_partition,
  path_MI_edge,
  GT_species[["MI_permissive"]],
  "MI",
  override_colors_by_comm = c(
    "1" = "#4DBBD5FF",
    "2" = cols_jco[1],
    "4" = cols_aaas[5],
    "5" = cols_aaas[2]
  )
)
print(MI_pair)

Skin_microbiome_pair <- draw_pair(
  path_Skin_microbiome_partition,
  path_Skin_microbiome_edge,
  GT_species[["Skin_microbiome"]],
  "Skin (microbiome)"
)

Skin_standard_pair <- draw_pair(
  path_Skin_standard_partition,
  path_Skin_standard_edge,
  GT_species[["skin_standard"]],
  "Skin (standard)"
)

Nasal_pair <- draw_pair(
  path_Nasal_partition,
  path_Nasal_edge,
  GT_species[["nasal"]],
  "Nasal"
)

Oral_pair <- draw_pair(
  path_Oral_partition,
  path_Oral_edge,
  GT_species[["oral"]],
  "Oral"
)

# ================================
# 10) Stack & print
# ================================
all_pairs_list <- list(
  Skin_microbiome_pair,
  Nasal_pair,
  Skin_standard_pair,
  Oral_pair
)

ncol_plots <- 2

all_pairs_plot <- wrap_plots(all_pairs_list, ncol = ncol_plots, guides = "keep") +
  plot_layout(guides = "collect")

print(all_pairs_plot)

# ================================
# 11) Save (optional)
# ================================
# nrows_plots <- ceiling(length(all_pairs_list) / ncol_plots)
# width_in  <- 8 * ncol_plots
# height_in <- 6 * nrows_plots
# ggsave("/media/junwoojo/18T/Metacontam_Figure/Network_all_pair_with_simulation.pdf",
#        plot = all_pairs_plot, width = width_in, height = height_in,
#        units = "in", dpi = 600, useDingbats = TRUE)
 
Figure2a <- MI_pair
#Supp_Fig1 <- all_pairs_plot