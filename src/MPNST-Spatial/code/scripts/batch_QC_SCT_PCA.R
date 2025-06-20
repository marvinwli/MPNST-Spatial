# ────────────────────────────────────────────────────────────────────────
# Seurat Workflow through PCA & Elbow Plot, with PDF Saving for Plots
# ────────────────────────────────────────────────────────────────────────

library(Seurat)
library(ggplot2)

# Increase allowable future size if needed
options(future.globals.maxSize = 700 * 1024^3)

# Load the pre‐made list of Seurat objects
sample_list <- readRDS("/home/mli110/MPNST-Spatial/rds objects/validation_sample_list.RDS")

# Directory to save figures
fig_dir <- "/home/mli110/MPNST-Spatial/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ────────────────────────────────────────────────────────────────────────
# 1) QC: Calculate % mitochondrial & generate pre‐filter violin plots.
#    Save all pre‐filter violin plots into one PDF.
# ────────────────────────────────────────────────────────────────────────

pdf(file = file.path(fig_dir, "pre_filter_QC.pdf"), width = 10, height = 5)
for (sample_name in names(sample_list)) {
  seu <- sample_list[[sample_name]]
  
  # Compute % mitochondrial reads
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  
  # Violin plot of nFeature_RNA, nCount_RNA, percent.mt
  p_pre <- VlnPlot(
    object = seu,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3,
    pt.size = 0.1
  ) + ggtitle(paste0(sample_name, " (pre‐filter QC)"))
  print(p_pre)
  
  # Save back the updated Seurat object (with percent.mt added)
  sample_list[[sample_name]] <- seu
}
dev.off()

# ────────────────────────────────────────────────────────────────────────
# 2) FILTER CELLS
# ────────────────────────────────────────────────────────────────────────

for (sample_name in names(sample_list)) {
  seu <- sample_list[[sample_name]]
  
  # Keep cells with 200 ≤ nFeature_RNA ≤ 6000 and percent.mt < 10
  seu <- subset(
    x = seu,
    subset = nFeature_RNA >= 200 &
      nFeature_RNA <= 10000 &
      percent.mt   <  15
  )
  
  sample_list[[sample_name]] <- seu
}

# ────────────────────────────────────────────────────────────────────────
# 3) Post‐filter QC: recompute % mitochondrial & generate post‐filter violin plots.
#    Save all post‐filter violin plots into one PDF.
# ────────────────────────────────────────────────────────────────────────

pdf(file = file.path(fig_dir, "post_filter_QC.pdf"), width = 10, height = 5)
for (sample_name in names(sample_list)) {
  seu <- sample_list[[sample_name]]
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  
  p_post <- VlnPlot(
    object = seu,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3,
    pt.size = 0.1
  ) + ggtitle(paste0(sample_name, " (post‐filter QC)"))
  print(p_post)
  
  sample_list[[sample_name]] <- seu
}
dev.off()

# ────────────────────────────────────────────────────────────────────────
# 4) NORMALIZE via SCTransform (regressing out percent.mt)
# ────────────────────────────────────────────────────────────────────────

for (sample_name in names(sample_list)) {
  seu <- sample_list[[sample_name]]
  
  seu <- SCTransform(
    object          = seu,
    assay           = "RNA",
    vars.to.regress = "percent.mt",
    verbose         = FALSE
  )
  
  sample_list[[sample_name]] <- seu
}

# ────────────────────────────────────────────────────────────────────────
# 5) PCA & Elbow Plot: run PCA and save all elbow plots into one PDF.
# ────────────────────────────────────────────────────────────────────────

pdf(file = file.path(fig_dir, "pca_elbow_plots.pdf"), width = 6, height = 5)
for (sample_name in names(sample_list)) {
  seu <- sample_list[[sample_name]]
  
  # Ensure SCT assay is active
  DefaultAssay(seu) <- "SCT"
  
  # Run PCA
  seu <- RunPCA(seu, verbose = FALSE)
  
  # Elbow plot (up to 50 PCs) for manual inspection
  e <- ElbowPlot(seu, ndims = 50) + ggtitle(paste0(sample_name, " PCA Elbow"))
  print(e)
  
  # Save back the updated object (with PCA slot)
  sample_list[[sample_name]] <- seu
}
dev.off()

# ────────────────────────────────────────────────────────────────────────
# End of Script: 
# • All pre‐filter QC violins saved in pre_filter_QC.pdf 
# • All post‐filter QC violins saved in post_filter_QC.pdf 
# • All elbow plots saved in pca_elbow_plots.pdf 
# • sample_list now contains Seurat objects with percent.mt, filtered cells, SCT assay, and PCA slots. 
# ────────────────────────────────────────────────────────────────────────
