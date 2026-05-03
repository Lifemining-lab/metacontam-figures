suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggpubr)
  library(ggsci)
  library(scales)
  library(grid)
  library(stringr)
  library(patchwork)
})

# ============================================================
# 0) Required objects assumed to exist
# ============================================================
required_objs <- c(
  "Nasal_predict",
  "squeegee_nasal_predict",
  "squeegee_nasal_oral_predict",
  "Nasal_output",
  "path_Nasal_partition"
)
missing_objs <- required_objs[!vapply(required_objs, exists, logical(1))]
if (length(missing_objs) > 0) {
  stop("Missing required objects in environment: ", paste(missing_objs, collapse = ", "))
}

# ============================================================
# 1) Paths and list loading
# ============================================================
path_blacklist <- file.path(BASE_INPUT, "BlackWhite_list", "blacklist.csv")
path_whitelist <- file.path(BASE_INPUT, "BlackWhite_list", "Healthy_Nasal_associated")
path_names <- file.path(BASE_INPUT, "Taxa_name", "names_dmp_scientific.csv")

blacklist <- read.csv(path_blacklist, sep = "\t")$blacklist %>% as.character()
whitelist <- readLines(path_whitelist) %>% as.character()

taxid_sciname <- read.csv(path_names, stringsAsFactors = FALSE) %>%
  transmute(taxid = as.character(taxid), scientific_name) %>%
  distinct(taxid, .keep_all = TRUE)

# ============================================================
# 2) Helper utilities
# ============================================================
clean_colnames <- function(nm) {
  nm <- as.character(nm)
  nm[is.na(nm)] <- ""
  nm <- trimws(nm)
  nm <- sub("^\ufeff", "", nm)
  nm <- gsub('^["\']+|["\']+$', "", nm)
  nm
}

strip_quotes <- function(x) {
  gsub('^["\']+|["\']+$', "", as.character(x))
}

read_tsv_loose <- function(path) {
  df <- tryCatch(
    read.table(
      path, sep = "\t", header = TRUE, check.names = FALSE,
      stringsAsFactors = FALSE, quote = "", comment.char = "", fill = TRUE
    ),
    error = function(e) data.frame()
  )
  
  if (ncol(df) == 0) {
    df <- read.table(
      path, sep = "\t", header = FALSE, check.names = FALSE,
      stringsAsFactors = FALSE, quote = "", comment.char = "", fill = TRUE
    )
  }
  df
}

read_partition_safe <- function(path) {
  df <- read_tsv_loose(path)
  names(df) <- clean_colnames(names(df))
  
  ti <- if ("taxid" %in% names(df)) {
    which(names(df) == "taxid")[1]
  } else {
    x <- which(tolower(names(df)) %in% c("taxid", "tax_id", "id"))
    if (length(x) == 0) 1 else x[1]
  }
  
  ci <- if ("cluster" %in% names(df)) {
    which(names(df) == "cluster")[1]
  } else {
    x <- which(grepl("cluster|community", tolower(names(df))))
    if (length(x) == 0) {
      if (ncol(df) >= 2) 2 else NA_integer_
    } else x[1]
  }
  
  if (is.na(ci)) stop(sprintf("'%s': cannot find a cluster/community column.", path))
  
  data.frame(
    taxid   = strip_quotes(df[[ti]]),
    cluster = strip_quotes(df[[ci]]),
    stringsAsFactors = FALSE
  )
}

norm_cluster <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_replace_all("^Community\\s*", "") %>%
    str_replace_all("^C\\s*", "") %>%
    str_replace_all("\\s+", "") %>%
    str_replace_all("[^0-9A-Za-z_-]", "")
}

get_gt_nasal <- function(GT_species) {
  if (!exists("GT_species")) stop("GT_species is not defined in the environment.")
  if (is.data.frame(GT_species) && "nasal" %in% names(GT_species)) {
    return(as.character(GT_species$nasal))
  }
  if (is.list(GT_species) && "nasal" %in% names(GT_species)) {
    return(as.character(GT_species[["nasal"]]))
  }
  as.character(GT_species)
}

# ============================================================
# 3) Tile plot function (Metacontam vs Squeegee)
# ============================================================
make_nasal_tile <- function(squeegee_vec, squeegee_label,
                            nasal_predict, blacklist, whitelist, taxid_sciname_tbl) {
  
  squeegee_vec  <- as.character(squeegee_vec)
  nasal_predict <- as.character(nasal_predict)
  
  all_taxid <- unique(c(nasal_predict, squeegee_vec))
  df <- tibble(
    taxid    = all_taxid,
    Category = case_when(
      taxid %in% blacklist & taxid %in% whitelist ~ "Mixed",
      taxid %in% blacklist                        ~ "Black List",
      taxid %in% whitelist                        ~ "White List",
      TRUE                                        ~ "Other"
    ),
    Metacontam = taxid %in% nasal_predict,
    !!squeegee_label := taxid %in% squeegee_vec
  )
  
  df_long <- df %>%
    pivot_longer(
      cols = c("Metacontam", all_of(squeegee_label)),
      names_to = "Tool",
      values_to = "Predicted"
    ) %>%
    filter(Predicted) %>%
    select(-Predicted)
  
  df_long <- df_long %>%
    left_join(taxid_sciname_tbl, by = "taxid") %>%
    mutate(label = ifelse(is.na(scientific_name) | scientific_name == "", taxid, scientific_name)) %>%
    select(-scientific_name)
  
  keep_cats <- c("White List", "Black List", "Mixed")
  category_colors <- c(
    "White List" = "#CCCCCC",
    "Black List" = "#000000",
    "Mixed"      = "#FDAE61"
  )
  
  df_long_core <- df_long %>%
    filter(Category %in% keep_cats) %>%
    mutate(
      Tool     = factor(Tool, levels = c(squeegee_label, "Metacontam")),
      Category = factor(Category, levels = keep_cats)
    )
  
  labels_sorted <- sort(unique(as.character(df_long_core$label)))
  df_long_core  <- df_long_core %>%
    mutate(label = factor(label, levels = rev(labels_sorted)))
  
  ggplot(df_long_core, aes(x = Tool, y = label, fill = Category)) +
    geom_tile(color = "white", width = 0.98, height = 0.9) +
    scale_fill_manual(values = category_colors, drop = FALSE) +
    scale_x_discrete(expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_discrete(
      limits = levels(df_long_core$label),
      expand = expansion(mult = c(0.003, 0.003)),
      position = "right"
    ) +
    labs(title = "", x = "", y = "", fill = "List Type") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title        = element_text(hjust = 0.5),
      axis.text.x       = element_text(face = "bold", size = 12),
      axis.text.y.left  = element_blank(),
      axis.text.y.right = element_text(size = 20),
      axis.ticks        = element_blank(),
      panel.grid        = element_blank(),
      panel.spacing.x   = unit(1, "pt"),
      plot.margin       = margin(t = 0, r = 0, b = 0, l = 5),
      legend.position   = "right"
    ) +
    coord_cartesian(clip = "off")
}

p_nasal_summary_core_sq_nasal <- make_nasal_tile(
  squeegee_vec      = squeegee_nasal_predict,
  squeegee_label    = "Squeegee (nasal-only)",
  nasal_predict     = Nasal_predict,
  blacklist         = blacklist,
  whitelist         = whitelist,
  taxid_sciname_tbl = taxid_sciname
)

p_nasal_summary_core_sq_no <- make_nasal_tile(
  squeegee_vec      = squeegee_nasal_oral_predict,
  squeegee_label    = "Squeegee (nasal+oral)",
  nasal_predict     = Nasal_predict,
  blacklist         = blacklist,
  whitelist         = whitelist,
  taxid_sciname_tbl = taxid_sciname
)

p_nasal_summary_core_sq_nasal
# p_nasal_summary_core_sq_no

# ============================================================
# 4) Prevalence vs relative abundance plot (Moraxella + GT)
# ============================================================
path_abun <- file.path(Nasal_output, "Abundance.txt")
path_prev <- file.path(Nasal_output, "Prevalence.txt")

abun <- read.table(path_abun, header = TRUE, check.names = FALSE) %>%
  rownames_to_column("taxid") %>%
  mutate(taxid = as.character(taxid))

prev <- read.table(path_prev, header = TRUE, check.names = FALSE) %>%
  rownames_to_column("taxid") %>%
  mutate(taxid = as.character(taxid))

gt_nasal <- get_gt_nasal(GT_species)

df_all <- abun %>%
  select(taxid, mean_abundance) %>%
  left_join(prev %>% select(taxid, preval), by = "taxid") %>%
  left_join(taxid_sciname %>% transmute(taxid, Species = scientific_name), by = "taxid") %>%
  mutate(
    mean_abundance = as.numeric(mean_abundance) / 100,
    mean_abundance = pmax(mean_abundance, 1e-10),
    Species = coalesce(Species, taxid),
    Genus = sub(" .*", "", Species),
    is_GT = taxid %in% gt_nasal,
    is_Moraxella = Genus == "Moraxella"
  ) %>%
  filter(is.finite(mean_abundance), is.finite(preval))

df_others <- df_all %>% filter(!is_Moraxella & !is_GT)
df_morax  <- df_all %>% filter(is_Moraxella)
df_gt     <- df_all %>% filter(is_GT)

size_base <- 2.0
size_big  <- size_base * 2

cosmic <- pal_cosmic()(8)
col_morax <- "#00FFFF"
col_gt    <- cosmic[2]

ymin <- floor(log10(min(df_all$mean_abundance, na.rm = TRUE)))
ymax <- ceiling(log10(max(df_all$mean_abundance, na.rm = TRUE)))
break_exps   <- seq(ymin, ymax, by = 2)
log10_breaks <- 10^break_exps

sup10_labels <- function(breaks) {
  exps <- as.integer(round(log10(breaks)))
  parse(text = paste0("10^", exps))
}

TXT <- list(title = 20, axis_title = 18, axis_text = 16, legend_text = 16)

p_nasal_prev_abun <- ggplot() +
  geom_bin2d(
    data = df_others,
    aes(x = preval, y = mean_abundance),
    bins = 60, fill = "#9e9e9e", alpha = 0.5, color = NA, show.legend = FALSE
  ) +
  geom_point(
    data = df_morax,
    aes(x = preval, y = mean_abundance, color = "Moraxella", shape = "Moraxella"),
    size = size_big, alpha = 0.95
  ) +
  geom_point(
    data = df_gt,
    aes(x = preval, y = mean_abundance, color = "GT", shape = "GT"),
    size = size_big, alpha = 0.95
  ) +
  scale_color_manual(
    values = c(Moraxella = col_morax, GT = col_gt),
    name = NULL
  ) +
  scale_shape_manual(
    values = c(Moraxella = 16, GT = 17),
    name = NULL,
    guide = "none"
  ) +
  scale_y_log10(breaks = log10_breaks, labels = sup10_labels) +
  labs(
    title = "Nasal — Prevalence vs Relative abundance",
    x = "Prevalence",
    y = "Relative abundance"
  ) +
  theme_pubr() +
  theme(
    plot.title      = element_text(size = TXT$title, face = "bold", hjust = 0.5),
    axis.title.x    = element_text(size = TXT$axis_title, face = "bold"),
    axis.title.y    = element_text(size = TXT$axis_title, face = "bold"),
    axis.text.x     = element_text(size = TXT$axis_text),
    axis.text.y     = element_text(size = TXT$axis_text),
    legend.text     = element_text(size = TXT$legend_text),
    legend.key.size = unit(14, "pt"),
    panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.9),
    aspect.ratio    = 1
  ) +
  guides(
    color = guide_legend(
      override.aes = list(shape = c(16, 17), size = 5)
    )
  )

p_nasal_prev_abun

# ============================================================
# 5) Community bar plot using partition (>=10 nodes only)
# ============================================================
comm_raw <- read_partition_safe(path_Nasal_partition) %>%
  transmute(
    taxid       = as.character(taxid),
    cluster_key = norm_cluster(cluster)
  )

big_keys <- comm_raw %>%
  dplyr::count(cluster_key, name = "n") %>%
  dplyr::filter(n >= 10) %>%
  dplyr::arrange(as.numeric(cluster_key)) %>%
  dplyr::pull(cluster_key)

key_to_idx <- setNames(seq_along(big_keys), big_keys)

comm_map_big <- comm_raw %>%
  filter(cluster_key %in% big_keys) %>%
  select(taxid, cluster_key) %>%
  distinct()

focus_long <- bind_rows(
  df_morax  %>% mutate(Group = "Moraxella"),
  df_gt     %>% mutate(Group = "GT"),
  df_others %>% mutate(Group = "Other")
) %>%
  select(taxid, Group) %>%
  left_join(comm_map_big, by = "taxid") %>%
  filter(!is.na(cluster_key)) %>%
  mutate(
    Group       = factor(Group, levels = c("Moraxella", "GT", "Other")),
    cluster_idx = as.character(key_to_idx[cluster_key])
  )

idx_levels <- sort(as.numeric(unique(focus_long$cluster_idx)), na.last = NA)
idx_levels <- as.character(idx_levels)

focus_long <- focus_long %>%
  mutate(cluster_idx = factor(cluster_idx, levels = idx_levels))

nejm_colors <- pal_nejm("default")(max(5, length(idx_levels)))
comm_colors <- setNames(nejm_colors[seq_along(idx_levels)], idx_levels)

# C1과 C2의 색만 서로 교체
if (all(c("1", "2") %in% names(comm_colors))) {
  tmp_col <- comm_colors["1"]
  comm_colors["1"] <- comm_colors["2"]
  comm_colors["2"] <- tmp_col
}

# 라벨은 C1, C2, C3, C4 ...
comm_labels <- setNames(paste0("C", idx_levels), idx_levels)
if ("1" %in% names(comm_labels)) comm_labels["1"] <- "C2"
if ("2" %in% names(comm_labels)) comm_labels["2"] <- "C1"

legend_breaks <- idx_levels
if (all(c("1", "2") %in% legend_breaks)) {
  legend_breaks <- c("2", "1", setdiff(legend_breaks, c("1", "2")))
}


count_df <- focus_long %>%
  dplyr::count(Group, cluster_idx, name = "n") %>%
  tidyr::complete(
    Group = levels(focus_long$Group),
    cluster_idx = idx_levels,
    fill = list(n = 0)
  ) %>%
  dplyr::group_by(Group) %>%
  dplyr::mutate(
    total = sum(n),
    pct   = ifelse(total > 0, 100 * n / total, 0)
  ) %>%
  dplyr::ungroup()

p_nasal_comm_counts <- ggplot(count_df, aes(x = Group, y = n, fill = cluster_idx)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.2) +
  scale_fill_manual(
    name   = "Community",
    values = comm_colors,
    breaks = legend_breaks,
    labels = unname(comm_labels[legend_breaks]),
    drop   = FALSE
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  labs(
    x = NULL,
    y = "Number of species",
    title = "Communities of Moraxella vs GT vs Other (Nasal; ≥10 nodes only)"
  ) +
  theme_pubr() +
  theme(
    plot.title  = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(face = "bold")
  )

p_nasal_comm_percent <- ggplot(count_df, aes(x = Group, y = pct, fill = cluster_idx)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.2) +
  scale_fill_manual(
    name   = "Community",
    values = comm_colors,
    breaks = idx_levels,
    labels = unname(comm_labels),
    drop   = FALSE
  ) +
  scale_y_continuous(
    labels = function(x) sprintf("%.0f%%", x),
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = NULL,
    y = "Percentage of species",
    title = "Communities of Moraxella vs GT vs Other (%, Nasal; ≥10 nodes only)"
  ) +
  theme_pubr() +
  theme(
    plot.title  = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(face = "bold")
  )

# ============================================================
# 6) Optional: patchwork layout
# ============================================================
p_nasal_prev_abun_noleg <- p_nasal_prev_abun +
  theme(legend.position = "none")

p_nasal_comm_counts_tight <- p_nasal_comm_counts +
  theme(
    plot.margin = margin(0, 0, 0, 4),
    legend.position = "right",
    legend.direction = "vertical"
  )

p_nasal_comm_counts_tight

Figure4b <- p_nasal_summary_core_sq_nasal 
Figure4c <- p_nasal_prev_abun
Figure4d <- p_nasal_comm_counts_tight
