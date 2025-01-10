library(Seurat)

#Normalization -- Scale --------------------------------------------------------
seurats_all <- readRDS("~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Objects/seurats_all.RDS")
#Normalizaton and PCA seurats90
seurats_all <- SCTransform(seurats_all)
seurats_all <- RunPCA(seurats_all, npcs = 100)
#Elbowplot
ElbowPlot(seurats_all, ndims = 50)
#Find neighboors
seurats_all <- FindNeighbors(seurats_all, dims = 1:30)
seurats_all <- RunUMAP(seurats_all, dims=1:30)
saveRDS(seurats_all, "~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Objects/seurats_all_norm.RDS")

