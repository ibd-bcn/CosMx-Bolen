library(Seurat)
library(readr)
library(BiocParallel)
library(InSituType)
library(plyr)
library(paletteer)
library(harmony)
library(stringr)
library(outliers)

## Read DATA and create Seurat object ------------------------------------------
fov <- 1:66
tissue <-
  c(
    "IBD",
    "IBD",
    "PD",
    "PD",
    "PD",
    "PD",
    "IBD",
    "IBD",
    "PD",
    "PD",
    "IBD",
    "IBD",
    "NHC",
    "NHC",
    "PD",
    "PD",
    "IBD",
    "IBD",
    "IBD",
    "IBD",
    "NHC",
    "NHC",
    "PD",
    "PD",
    "IBD",
    "IBD",
    "IBD",
    "IBD",
    "NHC",
    "NHC",
    "IBD",
    "IBD",
    "IBD",
    "IBD",
    "NHC",
    "NHC",
    "PD",
    "PD",
    "PD",
    "PD",
    "IBD",
    "IBD",
    "NHC",
    "NHC",
    "PD",
    "PD",
    "IBD",
    "IBD",
    "NHC",
    "NHC",
    "PD",
    "PD",
    "PD",
    "PD",
    "IBD",
    "IBD",
    "NHC",
    "NHC",
    "PD",
    "PD",
    "PD",
    "PD",
    "IBD",
    "IBD",
    "NHC",
    "NHC"
  )
num_patients <- ceiling(length(fov) / 2)
patients <-
  rep(paste0("patient", 1:num_patients),
      each = 2,
      length.out = length(fov))
# Creating the dictionary
fov_tissue_dict <- setNames(tissue, fov)
fov_patient_dict <- setNames(patients, fov)
#Read RAW FILES-----------------------------------------------------------------
c <-
  list.dirs(path = "~/SPATIAL/Mackensy_analysis/CosMx_Prot/Raw_data", recursive = F)[1]
#Obtain list for each cell
setwd(c)
#Raw
raw <- read_csv(list.files(pattern = "exprMat_file.csv"))
cell_names <- raw$cell
rownames(raw) <- cell_names
raw <- raw[, 4:ncol(raw)]
#Count Matrix
count <- raw[, !(colnames(raw) %in% c("Rb IgG", "Ms IgG1"))]
#Negprob Matrix
neg <- raw[, c("Rb IgG", "Ms IgG1")]
#Metadata
meta <- read_csv(list.files(pattern = "metadata_file.csv"))
rownames(meta) <- meta$cell
meta$patient <- fov_patient_dict[meta$fov]
meta$tissue <- fov_tissue_dict[meta$fov]
#Pols
pols <- read_csv(list.files(pattern = "polygons.csv"))
write.csv(
  pols,
  paste0(
    "~/SPATIAL/Mackensy_analysis/CosMx_Prot/0_Curation/Polygons/",
    strsplit(c, "/")[[1]][8],
    ".csv",
    sep = ""
  )
)
#Create seurat
seurats <-
  CreateSeuratObject(counts = t(count),
                     meta.data = meta,
                     assay = "Prot")
seurats <- RenameCells(seurats, new.names = cell_names)
Negprob <- CreateSeuratObject(counts = t(neg), assay = "Negprob")
Negprob <- RenameCells(Negprob, new.names = cell_names)
seurats[["Negprob"]] <-
  CreateAssayObject(data = Negprob@assays$Negprob$counts)

#SaveRDS
saveRDS(
  seurats,
  "~/SPATIAL/Mackensy_analysis/CosMx_Prot/0_Curation/Objects/seurats.RDS"
)


## Quality Control--------------------------------------------------------------

#FLAG1
#High expression proteins' proportion and threshold: flag cells where 50% or
#more (0.5, default; range 0-1) of proteins are in the 90th percentile or
#higher (0.9, default; range 0-1). Low expression proteins' number and
# threshold: flag cells where fewer than 10 (default; range 0-N where
#N is total number of proteins in panel) proteins are in the 50th percentile
#or higher (0.5, default; range 0-1)
counts <- t(as.data.frame(seurats@assays$Prot$counts))
bp <- MulticoreParam(workers = 20, progressbar = TRUE)
#Flagged cells high prot
cells_high_prot <-
  BiocParallel::bplapply(colnames(counts), BPPARAM = , \(c) {
    lower_threshold <- quantile(count[[c]], 0.50) # 50st
    upper_threshold <- quantile(count[[c]], 0.90) # 90th
    cl <- as.vector(counts[, c]) > as.numeric(lower_threshold)
    cu <- as.vector(counts[, c]) < as.numeric(upper_threshold)
    return(list(cl, cu))
  })
#Grab CL == LOW THRESSHOLD & CU == UPPER THRESHOLD
cl <-
  as.vector(rowSums(as.data.frame(lapply(cells_high_prot, function(x)
    x[[1]]))))
cu <-
  as.vector(rowSums(as.data.frame(lapply(cells_high_prot, function(x)
    x[[2]]))))
seurats$flag1 <- (cu > 34) * (cl > 10)

#Message of total of cells
message(paste0("A total of ", sum(seurats$flag1 == 1)),
        " cells have passed the LOW count and HIGH count filter.",
        sep = "")


#FLAG2
#Negative probe range: flag cells with negative probe mean below the lower
#threshold (2) or above the upper threshold (15).

inicial_num <- ncol(seurats)
negp <- seurats@assays$Negprob$data
avg_exp <- as.data.frame(colMeans(negp))
avg_exp$cells <- rownames(avg_exp)
colnames(avg_exp) <- c("mean", "cells")
seurats$flag2 <- (avg_exp$mean > 2) & (avg_exp$mean < 15)
message(paste0("A total of ", sum(seurats$flag2 == 1)),
        " cells have passed the LOW count and HIGH count filter.",
        sep = "")

#FLAG3
#Flags area outliers
data_vector <- seurats$Area
repeat {
  grubbs_test <- grubbs.test(data_vector, type = 10)
  p_value <- grubbs_test$p.value
  if (p_value > 0.01) {
    break
  }
  # Determine the outlier
  outliers <-
    regmatches(grubbs_test$alternative,
               gregexpr("\\d+", grubbs_test$alternative))
  # Convert the extracted numbers to numeric
  outliers <- as.numeric(unlist(outliers))
  # Remove the outlier from the data vector
  data_vector <- data_vector[!(data_vector %in% outliers)]
}
seurats$flag3 <- seurats$Area <= max(data_vector)
message(paste0("A total of ", sum(seurats$flag3 == 1)),
        " cells have passed the Area Outlier filterfilter.",
        sep = "")

#Do CELLS pass 3 FLAGS?
#Grab cells that pass the flagging
seurats$pass <-
  seurats$flag1 * seurats$flag2 * seurats$flag3
message(paste0("A total of ", sum(seurats$pass == 1)), " cells have passed the ENTIRE qc.", sep = "")


## Save QC --- -----------------------------------------------------------------
saveRDS(
  seurats,
  "~/SPATIAL/Mackensy_analysis/CosMx_Prot/0_Curation/Objects/qc_seurats.RDS"
)
#SaveRDS of filtered seurat
seurats_c <-
  seurats[, rownames(seurats@meta.data[seurats@meta.data$pass == 1, ])]
saveRDS(
  seurats_c,
  "~/SPATIAL/Mackensy_analysis/CosMx_Prot/0_Curation/Objects/qc_seurats_filt.RDS"
)

#Save objects for further analysis
write.csv(
  seurats_c@meta.data,
  "~/SPATIAL/Mackensy_analysis/CosMx_Prot/0_Curation/Objects/meta.csv"
)
saveRDS(
  seurats_c@assays$Prot$counts,
  "~/SPATIAL/Mackensy_analysis/CosMx_Prot/0_Curation/Objects/counts.RDS"
)


