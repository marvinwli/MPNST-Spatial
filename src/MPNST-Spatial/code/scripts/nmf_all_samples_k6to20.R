#!/usr/bin/env Rscript
# nmf_all_samples_k6to20.R

library(Seurat)
library(STdeconvolve)
library(NMF)
library(tools)

# params
input_dir  <- "/projects/yegnalab-hpc/Marvin_Li/MPNST rds objects"
out_dir    <- "/projects/yegnalab-hpc/Marvin_Li/MPNST_nmf_results"    
k_range    <- 6:20
nrun       <- 2
method     <- "brunet"
seed       <- 42

# make output directory
if (!dir.exists(out_dir)) dir.create(out_dir)

# find all samples
rds_files <- list.files(input_dir, "\\.rds$", full.names=TRUE, recursive=TRUE)

for (rds in rds_files) {
  sample_id <- file_path_sans_ext(basename(rds))
  samp_dir  <- file.path(out_dir, sample_id)
  if (!dir.exists(samp_dir)) dir.create(samp_dir)
  
  # load & preprocess once
  so            <- readRDS(rds)
  DefaultAssay(so) <- "Spatial"
  counts_clean  <- cleanCounts(
    GetAssayData(so, assay="Spatial", slot="counts"),
    min.lib.size = 100
  )
  
  nmf_input     <- as.matrix(log1p(counts_clean))
  
  for (k in k_range) {
    # output file for this W
    W_file <- file.path(samp_dir, paste0("W_k", k, ".rds"))
    H_file <- file.path(samp_dir, paste0("H_k", k, ".rds"))
    if (file.exists(out_file)) {
      message("Skipping ", sample_id, " k=", k, " (already exists)")
      next
    }
    
    message("Running NMF for ", sample_id, " k=", k)
    
    res <- nmf(
      nmf_input, 
      rank   = k, 
      nrun   = nrun, 
      method = method, 
      seed   = seed
    )
    
    # extract & name W
    W <- basis(res)    # genes × k
    colnames(W) <- paste0(sample_id, "_k=", k, "_Comp", seq_len(ncol(W)))
    saveRDS(W, W_file)
    
    # extract & name H
    H <- coef(res)     # components × spots
    rownames(H) <- colnames(W)   # same component names
    saveRDS(H, H_file)
  }
}

message("Congrats")
# Once all W matrices are on disk, you can reload them in bulk:
#beta_files <- list.files(out_dir, pattern="W_k.*\\.rds$", recursive=TRUE, full.names=TRUE)
#W_list     <- lapply(beta_files, readRDS)
