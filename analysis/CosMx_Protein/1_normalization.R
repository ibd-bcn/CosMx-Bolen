library(Seurat)
library(TidyDensity)

## Read object -----------------------------------------------------------------
seu <-
  readRDS(
    "~/SPATIAL/Mackensy_analysis/CosMx_Prot/0_Curation/Objects/qc_seurats_filt.RDS"
  )

## Normalize -------------------------------------------------------------------

#Total intensity Scaling + Arcsinh Normalization
total_intensity <- Matrix::colSums(seu@assays$Prot$counts)
average_total_intensity <- mean(total_intensity)
normalized_data <-
  sweep(seu@assays$Prot$counts, 2, total_intensity, FUN = "/")
scaled_data <- normalized_data * average_total_intensity
cofactor <- 50
arcsinh_data <- asinh(scaled_data / cofactor)

seu@assays$Prot$data <- arcsinh_data
saveRDS(seu,
        "~/SPATIAL/Mackensy_analysis/CosMx_Prot/01_Normalization/norm.RDS")
