library(Seurat)
library(ggplot2)

base_dir <- "/home/mli110/MPNST-Spatial/validation datasets"
sample_dirs <- list.dirs(path = base_dir, recursive = F)

sample_list <- list()

for(dir_path in sample_dirs){
  # 1) Extract the folder name to use as a sample identifier
  sample_name <- basename(dir_path)
  
  # 2) Read the 10X data (expects barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz inside)
  counts_matrix <- Read10X(data.dir = dir_path)
  
  # 3) Create a Seurat object from the counts
  seu <- CreateSeuratObject(
    counts  = counts_matrix,
    project = sample_name,
    min.cells    = 3,   # adjust as desired
    min.features = 200  # adjust as desired
  )
  
  # 4) (Optional) Add the sample name into metadata for tracking
  seu$sample_id <- sample_name
  
  # 5) Store the Seurat object in your list, named by folder
  sample_list[[sample_name]] <- seu
}

saveRDS(sample_list, file = "/home/mli110/MPNST-Spatial/rds objects/validation_sample_list.RDS")
