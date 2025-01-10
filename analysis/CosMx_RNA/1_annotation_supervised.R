library(Seurat)
library(readr)
library(BiocParallel)
library(InSituType)
library(plyr)
library(paletteer)
library(harmony)
library(stringr)
library(SeuratWrappers)
options(stringsAsFactors = FALSE,bitmapType = "cairo")
cols<- c(
  paletteer_d("ggsci::default_igv"),
  paletteer_d("ggsci::category20_d3"),
  paletteer_d("ggsci::default_ucscgb")
)
#Subset classification----------------------------------------------------------
#Read seurats
seurats <- readRDS( "~/SPATIAL/Mackensy_analysis/CosMx_RNA/0_Curation/Objects/qc_seurats_pass.RDS")

#Reference for cell typing
todas <- readRDS("/home/acorraliza/data_Albas/All_Together_31082021/piezas_finales/tog_numeric.RDS")
todas@active.ident <- as.factor(todas$subset)
#Average EXPRESSION of SUBSET
data <- AverageExpression(object = todas)
#Reference Markers of SUBSET
Idents(todas) <- "subset"
todas <- todas[rownames(seurats),]
ref.markers <- RunPrestoAll(todas)
ref.markers.ss <- ref.markers[
  ref.markers$pct.1 > 0.20 & ref.markers$pct.2 < 0.05,
]

#Compare with AMOUNT of count of COSMX object
counts <- as.data.frame(colSums(t(as.matrix(seurats@assays$RNA$counts))))
counts$genes <- rownames(counts)
colnames(counts) <- c("Count","Gene")
#Keep only markers from Reference
counts <- counts[ref.markers.ss$gene,]
#Assign gene to SUBSET
counts$Subset <- mapvalues(x = counts$Gene, from = ref.markers.ss$gene,to = ref.markers.ss$cluster)
ggplot(counts, aes(x = Gene, y = Count, fill = Subset)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  theme_minimal() +  labs(x = "Gene", y = "Count", title = "Cluster by SC reference") + scale_x_discrete(labels = NULL)

#InSituType: Annotate CELLS-----------------------------------------------------
meta <- seurats@meta.data
counts <- seurats[["RNA"]]$counts
negpb_c <- seurats@assays$Negprob$data
counts_c <- seurats@assays$RNA$counts
ifdata = as.matrix(meta[, c("Area","AspectRatio","Mean.Membrane","Mean.PanCK","Mean.CD45","Mean.DAPI","Mean.CD68")])

#Prepare input
cohort = fastCohorting(mat = ifdata,gaussian_transform = TRUE)
columns_to_keep <- setdiff(rownames(counts_c), c("IGKC","IGHA1"))
counts_c <- t(as.matrix(counts_c))[,columns_to_keep]
negpb_c <-  Matrix::rowMeans(t(as.matrix(negpb_c)))

sup <- insitutypeML(x = counts_c,
                    neg = negpb_c,
                    cohort = cohort,
                    reference_profiles = as.matrix(data$RNA)) 

clust_names <- as.data.frame(sup$clust)
clust_names$cell_names <- rownames(clust_names)
clust_names$prob <- sup$prob
seurats@meta.data$subset <- mapvalues(x = seurats@meta.data$cell, from = clust_names$cell_names, to = clust_names$`sup$clust`)
seurats@meta.data$subset_prob <- mapvalues(x = seurats@meta.data$cell, from = clust_names$cell_names, to= clust_names$prob)
#SaveRDS
saveRDS(seurats,"~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Objects/seurats.RDS")
#Threshold of probability >0.99
seurats.99 <- seurats[,rownames(seurats@meta.data[seurats@meta.data$subset_prob > 0.99,])]
saveRDS(seurats.99,"~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Objects/seurats99.RDS")
#Threshold of probability >0.90
seurats.90 <- seurats[,rownames(seurats@meta.data[seurats@meta.data$subset_prob > 0.90,])]
saveRDS(seurats.90,"~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Objects/seurats90.RDS")
#Threshold of probability >0.85
seurats.85 <- seurats[,rownames(seurats@meta.data[seurats@meta.data$subset_prob > 0.85,])]
saveRDS(seurats.85,"~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Objects/seurats85.RDS")
#Threshold of probability >0.75
seurats.75 <- seurats[,rownames(seurats@meta.data[seurats@meta.data$subset_prob > 0.75,])]
saveRDS(seurats.75,"~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Objects/seurats75.RDS")

#Plots
cols <- cols[seq_along(unique(sup$clust))]
names(cols) <- unique(sup$clust)
# make the flightpath plot
fp <- flightpath_plot(flightpath_result = NULL, insitutype_result = sup, col = cols[sup$clust])
print(fp)

#Write meta file
write.csv(seurats.75@meta.data,"~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Objects/meta.csv")


#Refined annotation-------------------------------------------------------------

#Epi----------------------------------------------------------------------------

#Read Reference
todas_epi <- readRDS('/home/acorraliza/000_GitHub/ibd-bcn_single_cell/Analysis of our data/02_Samples_Together/SUBSETS/ON_THEIR_OWN/epi_annotated.RDS')

#Delete cell types that DO NOT match with CosMx genes
grab_annot <- rownames(todas_epi@meta.data[!(todas_epi@meta.data$annotation_intermediate %in% c("Inflammatory colonocyte", "BEST4 OTOP2","Epithelium Ribhi","Paneth-like","Tuft cells", "Secretory progenitor")),])
todas_epi <- todas_epi[,grab_annot]

todas_epi@active.ident <- as.factor(todas_epi$annotation_intermediate)
data <- AverageExpression(object = todas_epi)
#Read Subset Classification
seu <- readRDS("~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Objects/seurats75.RDS")
#Grab only epithilium
epi <- seu[,rownames(seu@meta.data[seu@meta.data$subset == "epi",])]

#Insitutype
meta <- epi@meta.data
counts_c <- epi[["RNA"]]$counts
negpb_c <- epi@assays$Negprob$data
ifdata = as.matrix(meta[,c("Area","AspectRatio","Mean.Membrane","Mean.PanCK","Mean.CD45","Mean.DAPI","Mean.CD68")])
cohort = fastCohorting(mat = ifdata,gaussian_transform = TRUE)
sup <- insitutypeML(x = t(counts_c),
                    neg = Matrix::rowMeans(t(negpb_c)),
                    cohort = cohort,
                    reference_profiles = as.matrix(data$RNA)) 

clust_names <- as.data.frame(sup$clust)
clust_names$cell_names <- rownames(clust_names)
clust_names$prob <- sup$prob
epi@meta.data$refined <- mapvalues(x = epi@meta.data$cell, from = clust_names$cell_names, to= clust_names$`sup$clust`)
epi@meta.data$refined_prob <- mapvalues(x = epi@meta.data$cell, from = clust_names$cell_names, to= clust_names$prob)
saveRDS(epi,"~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Objects/epi.RDS")

#FlightPlot
cols_epi <- cols[seq_along(unique(sup$clust))]
names(cols_epi) <- unique(sup$clust)
fp <- flightpath_plot_black(flightpath_result = NULL, insitutype_result = sup, col = cols_epi[sup$clust])
png(
  "~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Markers/epi.png",
  width = 15,
  height = 10,
  units = "in",
  res = 600
)
print(fp) 
dev.off()

# Get Markers
markers <- get_markers_ref(ref = todas_epi, object = epi, ss= FALSE, anot = "annotation_intermediate")
write.csv(markers,"~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Markers/epi.csv")

#Myeloids-----------------------------------------------------------------------

#Read Reference
todas_myeloids <- readRDS('/home/acorraliza/000_GitHub/ibd-bcn_single_cell/Analysis of our data/02_Samples_Together/SUBSETS/ON_THEIR_OWN/myeloids_annotated.RDS')

#Delete clusters that have no markers for it
grab_annot <- rownames(todas_myeloids@meta.data[!(todas_myeloids@meta.data$annotation_intermediate %in% c("Cycling myeloid", "IDA macrophage")),])
todas_myeloids <- todas_myeloids[,grab_annot]
#Grab clusters to annotate
todas_myeloids@active.ident <- as.factor(todas_myeloids$annotation_intermediate)
data <- AverageExpression(object = todas_myeloids)

#Grab myeloids
seu <- readRDS("~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Objects/seurats75.RDS")
myeloids <- seu[,rownames(seu@meta.data[seu@meta.data$subset == "myeloids",])]

#Insitutype Supervised Clustering
meta <- myeloids@meta.data
counts_c <- myeloids[["RNA"]]$counts
negpb_c <- myeloids@assays$Negprob$data
ifdata = as.matrix(meta[, c("Area","AspectRatio","Mean.Membrane","Mean.PanCK","Mean.CD45","Mean.DAPI","Mean.CD68")])
cohort = fastCohorting(mat = ifdata,gaussian_transform = TRUE)
sup <- insitutypeML(x = t(counts_c),
                    neg = Matrix::rowMeans(t(negpb_c)),
                    cohort = cohort,
                    reference_profiles = as.matrix(data$RNA)) 

clust_names <- as.data.frame(sup$clust)
clust_names$cell_names <- rownames(clust_names)
clust_names$prob <- sup$prob
myeloids@meta.data$refined <- mapvalues(x = myeloids@meta.data$cell, from = clust_names$cell_names, to= clust_names$`sup$clust`)
myeloids@meta.data$refined_prob <- mapvalues(x = myeloids@meta.data$cell, from = clust_names$cell_names, to= clust_names$prob)
saveRDS(myeloids,"~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Objects/myeloids.RDS")

#Flight plot
cols_my <- cols[seq_along(unique(sup$clust))]
names(cols_my) <- unique(sup$clust)
fp <- flightpath_plot_black(flightpath_result = NULL, insitutype_result = sup, col = cols_my[sup$clust])
png(
  "~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Markers/myeloids.png",
  width = 15,
  height = 10,
  units = "in",
  res = 600
)
print(fp) 
dev.off()

#Markers
markers <- get_markers_ref(ref = todas_myeloids, object = myeloids, ss= FALSE, anot = "annotation_intermediate")
write.csv(markers,"~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Markers/myeloids.csv")

#Tcells-------------------------------------------------------------------------

#Read Reference
todas_tcells <- readRDS('/home/acorraliza/000_GitHub/ibd-bcn_single_cell/Analysis of our data/02_Samples_Together/SUBSETS/ON_THEIR_OWN/tcells_annotated.RDS')

#Tcells
grab_annot <- rownames(todas_tcells@meta.data[!(todas_tcells@meta.data$annotation_intermediate %in% c("gd IEL", "DN","T cells CCL20","Ribhi T cells","Cycling T cells","Tregs","ILC4","MT T cells","MAIT")),])
todas_tcells <- todas_tcells[,grab_annot]

#Grab clusters to annotate
todas_tcells@active.ident <- as.factor(todas_tcells$annotation_intermediate)
data <- AverageExpression(object = todas_tcells)

#Grab Tcells
seu <- readRDS("~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Objects/seurats75.RDS")
tcells <- seu[,rownames(seu@meta.data[seu@meta.data$subset == "tcells",])]

#Insitutype Supervised Clustering
meta <- tcells@meta.data
counts_c <- tcells[["RNA"]]$counts
negpb_c <- tcells@assays$Negprob$data
ifdata = as.matrix(meta[, c("Area","AspectRatio","Mean.Membrane","Mean.PanCK","Mean.CD45","Mean.DAPI","Mean.CD68")])
cohort = fastCohorting(mat = ifdata,gaussian_transform = TRUE)
sup <- insitutypeML(x = t(counts_c),
                    neg = Matrix::rowMeans(t(negpb_c)),
                    cohort = cohort,
                    reference_profiles = as.matrix(data$RNA)) 

clust_names <- as.data.frame(sup$clust)
clust_names$cell_names <- rownames(clust_names)
clust_names$prob <- sup$prob
tcells@meta.data$refined <- mapvalues(x = tcells@meta.data$cell, from = clust_names$cell_names, to= clust_names$`sup$clust`)
tcells@meta.data$refined_prob <- mapvalues(x = tcells@meta.data$cell, from = clust_names$cell_names, to= clust_names$prob)
saveRDS(tcells,"~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Objects/tcells.RDS")

#Flightpath Plot
cols_tc <- cols[seq_along(unique(sup$clust))]
names(cols_tc) <- unique(sup$clust)
fp <- flightpath_plot_black(flightpath_result = NULL, insitutype_result = sup, col = cols_tc[sup$clust])
png(
  "~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Markers/tcells.png",
  width = 15,
  height = 10,
  units = "in",
  res = 600
)
print(fp) 
dev.off()

#Tcells Markers
markers <- get_markers_ref(ref = todas_tcells, object = tcells, ss= FALSE,anot = "annotation_intermediate")
write.csv(markers,"~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Markers/tcells.csv")


#Plasmas------------------------------------------------------------------------

#Read Reference Plasmas
todas_plasmas <- readRDS('/home/acorraliza/000_GitHub/ibd-bcn_single_cell/Analysis of our data/02_Samples_Together/SUBSETS/ON_THEIR_OWN/plasmas_annotated.RDS')

#Plasmas grab
grab_annot <- rownames(todas_plasmas@meta.data[!(todas_plasmas@meta.data$annotation_intermediate %in% c("PC IER", "PC IgA heat shock","Naïve B cell","Memory B cell","GC B cell","Cycling cells")),])
todas_plasmas <- todas_plasmas[,grab_annot]

#Grab Clusters for annotation
todas_plasmas@active.ident <- as.factor(todas_plasmas$annotation_intermediate)
data <- AverageExpression(object = todas_plasmas)

#Grab plasmas
seu <- readRDS("~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Objects/seurats75.RDS")
plasmas <- seu[,rownames(seu@meta.data[seu@meta.data$subset == "plasmas",])]

#Insitutype Supervised Clustering
meta <- plasmas@meta.data
counts_c <- plasmas[["RNA"]]$counts
negpb_c <- plasmas@assays$Negprob$data
ifdata = as.matrix(meta[, c("Area","AspectRatio","Mean.Membrane","Mean.PanCK","Mean.CD45","Mean.DAPI","Mean.CD68")])
cohort = fastCohorting(mat = ifdata,gaussian_transform = TRUE)
sup <- insitutypeML(x = t(counts_c),
                    neg = Matrix::rowMeans(t(negpb_c)),
                    cohort = cohort,
                    reference_profiles = as.matrix(data$RNA)) 

clust_names <- as.data.frame(sup$clust)
clust_names$cell_names <- rownames(clust_names)
clust_names$prob <- sup$prob
plasmas@meta.data$refined <- mapvalues(x = plasmas@meta.data$cell, from = clust_names$cell_names, to= clust_names$`sup$clust`)
plasmas@meta.data$refined_prob <- mapvalues(x = plasmas@meta.data$cell, from = clust_names$cell_names, to= clust_names$prob)
saveRDS(plasmas,"~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Objects/plasmas.RDS")

#Flight plot
cols_pl <- cols[seq_along(unique(sup$clust))]
names(cols_pl) <- unique(sup$clust)
fp <- flightpath_plot_black(flightpath_result = NULL, insitutype_result = sup, col = cols_pl[sup$clust])
png(
  "~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Markers/plasmas.png",
  width = 15,
  height = 10,
  units = "in",
  res = 600
)
print(fp) 
dev.off()

#Plasmas Markers
markers <- get_markers_ref(ref = todas_plasmas, object = plasmas, ss= FALSE,anot = "annotation_intermediate")
write.csv(markers,"~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Markers/plasmas.csv")

#Stroma-------------------------------------------------------------------------

#Read Reference Stroma
todas_stroma <- readRDS('/home/acorraliza/000_GitHub/ibd-bcn_single_cell/Analysis of our data/02_Samples_Together/SUBSETS/ON_THEIR_OWN/stroma_annotated.RDS')

#Delete clusters that have no markers for it
grab_annot <- rownames(todas_stroma@meta.data[!(todas_stroma@meta.data$annotation_intermediate %in% c("Inflammatory fibroblasts")),])
todas_stroma <- todas_stroma[,grab_annot]

#Annotation change 
todas_stroma$annotation_intermediate <- gsub("S1","Fibroblasts",todas_stroma$annotation_intermediate)
todas_stroma$annotation_intermediate <- gsub("S3","Fibroblasts",todas_stroma$annotation_intermediate)
#Grab clusters to annotate
todas_stroma@active.ident <- as.factor(todas_stroma$annotation_intermediate)
data <- AverageExpression(object = todas_stroma)

#Grab Stroma
seu <- readRDS("~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Objects/seurats75.RDS")
stroma <- seu[,rownames(seu@meta.data[seu@meta.data$subset == "stroma",])]

#Insitutype Supervised Clustering
meta <- stroma@meta.data
counts_c <- stroma[["RNA"]]$counts
negpb_c <- stroma@assays$Negprob$data
ifdata = as.matrix(meta[, c("Area","AspectRatio","Mean.Membrane","Mean.PanCK","Mean.CD45","Mean.DAPI","Mean.CD68")])
cohort = fastCohorting(mat = ifdata,gaussian_transform = TRUE)
sup <- insitutypeML(x = t(counts_c),
                    neg = Matrix::rowMeans(t(negpb_c)),
                    cohort = cohort,
                    reference_profiles = as.matrix(data$RNA)) 

clust_names <- as.data.frame(sup$clust)
clust_names$cell_names <- rownames(clust_names)
clust_names$prob <- sup$prob
stroma@meta.data$refined <- mapvalues(x = stroma@meta.data$cell, from = clust_names$cell_names, to= clust_names$`sup$clust`)
stroma@meta.data$refined_prob <- mapvalues(x = stroma@meta.data$cell, from = clust_names$cell_names, to= clust_names$prob)
saveRDS(stroma,"~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Objects/stroma.RDS")

#Flightpath Plot
cols_st <- cols[seq_along(unique(sup$clust))]
names(cols_st) <- unique(sup$clust)
fp <- flightpath_plot_black(flightpath_result = NULL, insitutype_result = sup, col = cols_st[sup$clust])
png(
  "~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Markers/stroma.png",
  width = 15,
  height = 10,
  units = "in",
  res = 600
)
print(fp) 
dev.off()

#Stroma Markers
markers <- get_markers_ref(ref = todas_stroma, object = stroma, ss= FALSE,anot = "annotation_intermediate")
write.csv(markers,"~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Markers/stroma.csv")

#Fuse all together--------------------------------------------------------------

tcells <- readRDS("~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Objects/tcells.RDS") 
plasmas <- readRDS("~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Objects/plasmas.RDS")
epi <- readRDS("~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Objects/epi.RDS")
myeloids <- readRDS("~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Objects/myeloids.RDS")
stroma <- readRDS("~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Objects/stroma.RDS")
seurats_r <- list(tcells,plasmas,epi,myeloids,stroma)

seurats_all <- seurats_r[[1]]
for(i in 2:length(seurats_r)) {
  seurats_all <- merge(seurats_all, y = seurats_r[[i]])
}
seurats_all <- JoinLayers(seurats_all)

saveRDS(seurats_all,"~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Objects/seurats_all.RDS")
write.csv(seurats_all@meta.data,"~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Objects/meta.csv")


