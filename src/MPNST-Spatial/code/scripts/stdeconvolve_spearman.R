library(tidyverse)
library(ComplexHeatmap)
library(future.apply)
library(data.table)
library(circlize)
plan(multisession, workers = 30)

beta_paths <- list.files(
  path      = "/projects/yegnalab-hpc/Marvin_Li/STDeconvolve results",
  pattern   = "beta\\.csv$",
  recursive = TRUE,
  full.names= TRUE
)

# 2) for each file, read it in and turn it into a matrix of topics×genes
topic_mats <- future_lapply(beta_paths, function(path) {
  # read with data.table
  dt <- fread(path)
  
  # first column is the topic name
  topic_ids <- dt[[1]]
  
  # the rest are gene columns, so coerce to numeric matrix
  mat <- as.matrix(dt[, -1, with = FALSE])
  
  # set rownames = topic IDs
  rownames(mat) <- topic_ids
  
  # figure out sample and K from the path
  samp <- basename(dirname(dirname(path)))  # e.g. "mpnst_11502"
  krun <- basename(dirname(path))          # e.g. "10"
  
  # rename each topic (row) to be unique
  new_rn <- paste0(samp, "_k=", krun, "_topic", rownames(mat))
  rownames(mat) <- new_rn
  
  return(mat)
})

all_genes <- sort(unique(unlist(lapply(topic_mats, colnames))))

# 3) Join them all by gene name
topic_mats_union <- lapply(topic_mats, function(mat) {
  missing_genes <- setdiff(all_genes, colnames(mat))
  if (length(missing_genes) > 0) {
    # create a zero‐filled padding matrix
    pad <- matrix(
      0,
      nrow = nrow(mat),
      ncol = length(missing_genes),
      dimnames = list(rownames(mat), missing_genes)
    )
    mat <- cbind(mat, pad)
  }
  # reorder to the global gene order
  mat[, all_genes, drop = FALSE]
})

big_mat <- do.call(rbind, topic_mats_union)

corr_topics <- cor(
  t(big_mat),
  method = "spearman",
  use    = "pairwise.complete.obs"
)





# 1) big_mat[t, g] is beta for topic t and gene g.

# 2) Compute the marginal (mean over topics) for each gene:
marginal <- colMeans(big_mat, na.rm = TRUE)   # length = number of genes

# 3) Divide each column of big_mat by its marginal:
lift_mat <- sweep(big_mat, 2, marginal, FUN = "/")
# now lift_mat[t, g] = beta[t, g] / marginal[g]

# 4) (Optional) log-transform to symmetrize:
log_lift <- log2(lift_mat + 1e-6)  # add small pseudocount to avoid log(0)

# 5) Use log_lift for downstream correlation & visualization:
corr_lift <- cor(
  t(log_lift),
  method = "spearman",
  use    = "pairwise.complete.obs"
)









# 1) Extract the sample name from each topic name
topic_names <- rownames(corr_topics)
# split on the first underscore
samples <- sub("^([^_]*_[^_]*)_.*$", "\\1", topic_names)
# 2) Define a unique color per sample
sample_levels <- unique(samples)
# e.g. use a qualitative palette
sample_cols <- structure(
  RColorBrewer::brewer.pal(length(sample_levels), "Set3"),
  names = sample_levels
)

ha_col <- HeatmapAnnotation(
  Sample = samples,
  col    = list(Sample = sample_cols),
  annotation_legend_param = list(
    Sample = list(title = "Sample", at = sample_levels)
  )
)

# 1) Define a three‐point color ramp that puts 0 in the middle
col_fun <- colorRamp2(
  breaks = c(-1, 0, 1),
  colors = c("blue", "white", "red")
)

# 2) Pass that into your Heatmap call via the `col` argument
Heatmap(
  corr_lift,
  name                   = "Spearman\nρ",
  col                    = col_fun,
  top_annotation         = ha_col,
  show_row_names         = FALSE,
  show_column_names      = FALSE,
  cluster_rows           = TRUE,
  cluster_columns        = TRUE,
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  heatmap_legend_param   = list(
    at    = c(-1, -0.5, 0, 0.5, 1),
    labels = c("-1", "-0.5", "0", "0.5", "1"),
    title  = "ρ"
  )
)


library(igraph)

# 1a) pick a correlation cutoff
cutoff <- 0.20

# 1b) build adjacency: keep only edges ≥ cutoff
adj   <- corr_topics
adj[adj < cutoff] <- 0
diag(adj)        <- 0  # drop self‐loops

# 1c) graph
g     <- graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE)

# 1d) detect communities
comms <- cluster_louvain(g)  

# 1e) extract membership and sample labels
df <- data.frame(
  topic   = names(membership(comms)),
  cluster = membership(comms),
  sample  = sub("^([^_]*_[^_]*)_.*$", "\\1", names(membership(comms)))
)

# 1f) count how many samples per community
table(df$cluster, df$sample)

# 1g) identify cross‐sample communities
cross_sample <- names(which(rowSums(table(df$cluster, df$sample) > 0) > 1))
df[df$cluster %in% cross_sample, ]
