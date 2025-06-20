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
output_dir <- "/projects/yegnalab-hpc/Marvin_Li/STDeconvolve test"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

rds_files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)

lapply(rds_files, function(rds_file) {
  sample_id <- file_path_sans_ext(basename(rds_file))
  message("Processing sample: ", sample_id)
  
  seuratobj <- readRDS(rds_file)
  DefaultAssay(seuratobj) <- "Spatial"
  counts_mat <- GetAssayData(seuratobj, assay = "Spatial", layer = "counts")
  
  message("  class(counts_mat) = ", class(counts_mat)[1])
  
  counts_clean <- cleanCounts(counts_mat, min.lib.size = 100)
  
  message("  class(counts_clean) = ", class(counts_clean)[1])
  
  corpus <- restrictCorpus(counts_clean, plot = FALSE)
  corpus_mat <- t(as.matrix(corpus))
  
  message("Made it to GetTissueCoordinates")
  pos <- GetTissueCoordinates(seuratobj)
  morphology <- seuratobj$Morphology
  group_colors <- c("MPNST" = "black", "NF" = "lightgray", "Muscle" = "lightgreen", "Capsule" = "lightblue", "Necrosis" = "lightyellow")

  
  message("  Saving pos for ", sample_id, head(pos))
 

  
})
