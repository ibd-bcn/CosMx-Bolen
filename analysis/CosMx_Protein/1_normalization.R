library(Seurat)
library(TidyDensity)

## Read object -----------------------------------------------------------------
seu <-
  readRDS(
    "analysis/CosMx_Protein/Objects/qc_seurats_filt.RDS"
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
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu,npcs = 40)
ElbowPlot(seu,ndims = 40)
seu <- RunUMAP(seu,dims = 1:25)
saveRDS(seu,
        "analysis/CosMx_Protein/Objects/norm.RDS")

