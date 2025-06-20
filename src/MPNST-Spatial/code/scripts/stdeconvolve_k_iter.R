library(methods)
library(STdeconvolve)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(tools)

custom_lib <- "/projects/yegnalab-hpc/Marvin_Li/Rlibs"
.libPaths(custom_lib)

required_ver <- "5.1.0"
have_seurat <- requireNamespace("Seurat", quietly = TRUE) &&
  as.character(packageVersion("Seurat")) == required_ver

if (!have_seurat) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", lib = custom_lib)
  }
  remotes::install_github(
    "satijalab/seurat@v5.1.0",
    lib     = custom_lib,
    upgrade = "never"
  )
}
library(Seurat)
message("Seurat ", packageVersion("Seurat"), " loaded from ", .libPaths()[1])


input_dir <- "/projects/yegnalab-hpc/Marvin_Li/MPNST rds objects"
output_dir <- "/projects/yegnalab-hpc/Marvin_Li/STDeconvolve results 2"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

rds_files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)
k_range <- 6:15

lapply(rds_files, function(rds_file) {
  sample_id <- file_path_sans_ext(basename(rds_file))
  message("Processing sample: ", sample_id)
  
  seuratobj <- readRDS(rds_file)
  DefaultAssay(seuratobj) <- "Spatial"
  counts_mat <- GetAssayData(seuratobj, assay = "Spatial", layer = "counts")
  
  message("  class(counts_mat) = ", class(counts_mat)[1])
  
  counts_clean <- cleanCounts(counts_mat, min.lib.size = 100)
  
  message("  class(counts_clean) = ", class(counts_clean)[1])
  
  corpus <- restrictCorpus(counts_clean, plot = FALSE, removeBelow = 0, nTopOD = 2000, alpha = 0.1)
  corpus_mat <- t(as.matrix(corpus))
  message("Made it to GetTissueCoordinates")
  pos <- setNames(
    as.data.frame(GetTissueCoordinates(seuratobj))[, 1:2],
    c("x","y")
  )
  rownames(pos) <- rownames(seuratobj@meta.data)  
  morphology <- seuratobj$Morphology
  group_colors <- c("MPNST" = "black", "NF" = "lightgray", "Muscle" = "lightgreen", "Capsule" = "lightblue", "Necrosis" = "lightyellow")
  
  # Run LDA for Ks = 6 to 20
  model_list <- fitLDA(
    corpus_mat,
    Ks = k_range,
    perc.rare.thresh = 0.05,
    plot = TRUE,
    ncores = 15,
    verbose = TRUE
  )

    # Loop through each K's model
  for (k in names(model_list$models)) {
    message("Working on K = ", k)
    lda_fit <- model_list$models[[k]]
    
    results <- getBetaTheta(lda_fit, perc.filt = 0.05, betaScale = 1000)
    beta <- results$beta
    theta <- results$theta
    
    # Make K-specific folder
    k_dir <- file.path(output_dir, sample_id, k)
    dir.create(k_dir, recursive = TRUE, showWarnings = FALSE)
    
    saveRDS(results, file.path(k_dir, "results.rds"))
    write.csv(beta, file.path(k_dir, "beta.csv"))
    write.csv(theta, file.path(k_dir, "theta.csv"))
    
    # Align metadata for plotting
    common_spots <- intersect(rownames(theta), rownames(pos))
    pos_subset <- pos[common_spots, ]
    morph_subset <- morphology[common_spots]
    
    pos_rotated <- pos_subset
    pos_rotated[, "x"] <- pos_subset[, "y"]
    pos_rotated[, "y"] <- -pos_subset[, "x"]
    
    # Plot and save each topic
    for (i in seq_along(colnames(theta))) {
      topic_name <- colnames(theta)[i]
      p <- vizTopic(
        theta = theta,
        pos = pos_rotated,
        topic = topic_name,
        plotTitle = paste0(sample_id, "_", k, "_", topic_name),
        size = 1.5,
        stroke = 0.2,
        alpha = 1,
        groups = morph_subset,
        group_cols = group_colors,
        low = "white",
        high = "red"
      )
      ggsave(
        filename = file.path(k_dir, paste0(sample_id, "_", "k=", k, "_", "topic", topic_name, ".png")),
        plot = p,
        width = 6, height = 6, dpi = 300
      )
    }
  }
})
