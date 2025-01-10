library(Seurat)
library(readr)
library(BiocParallel)
library(InSituType)
library(plyr)
library(paletteer)
library(harmony)
library(stringr)
library(ggplot2)
library(outliers)

#IDs patient
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
  list.dirs(path = "~/SPATIAL/Mackensy_analysis/CosMx_RNA/raw_data", recursive = F)[1]
#Obtain list for each cell
setwd(c)
#Raw
raw <- read_csv(list.files(pattern = "exprMat_file.csv"))
cell_names <- raw$cell
rownames(raw) <- cell_names
raw <- raw[, 4:ncol(raw)]
#NegativeProbes
neg <- raw[, grep(pattern = "Negative", x = colnames(raw))]
#SystemControl
system <- raw[, grep(pattern = "SystemControl", x = colnames(raw))]
#CountMatrix
diff_cols <-
  setdiff(colnames(raw), c(colnames(neg), colnames(system)))
count <- raw[, diff_cols]
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
    "~/SPATIAL/Mackensy_analysis/CosMx_RNA/0_Curation/Polygons/",
    strsplit(c, "/")[[1]][8],
    ".csv",
    sep = ""
  )
)
#Mols
mols <- read_csv(list.files(pattern = "tx_file.csv"))
mols <- mols[mols$cell_ID != 0,]
write.csv(
  mols,
  paste0(
    "~/SPATIAL/Mackensy_analysis/CosMx_RNA/0_Curation/Molecules/",
    strsplit(c, "/")[[1]][8],
    ".csv",
    sep = ""
  )
)
#Create seurat
seurats <-
  CreateSeuratObject(counts = t(count),
                     meta.data = meta,
                     assay = "RNA")
seurats <- RenameCells(seurats, new.names = cell_names)
Negprob <- CreateSeuratObject(counts = t(neg), assay = "Negprob")
Negprob <- RenameCells(Negprob, new.names = cell_names)
SystemControl <-
  CreateSeuratObject(counts = t(system), assay = "System")
SystemControl <- RenameCells(SystemControl, new.names = cell_names)
seurats[["Negprob"]] <-
  CreateAssayObject(data = Negprob@assays$Negprob$counts)
seurats[["SystemControl"]] <-
  CreateAssayObject(data = SystemControl@assays$System$counts)

#SaveRDS
saveRDS(seurats,
        "~/SPATIAL/Mackensy_analysis/CosMx_RNA/0_Curation/Objects/seurats.RDS")

#Quality Control ---------------------------------------------------------------

#Outlier test for Negative Probes: % Pass is % of Negative Probes which that
#pass Grubb's test (i.e. are not outliers) based
#on the p-value threshold set in module parameters.
#Outlier P-value cutoff: (default: 0.01; range: 0-1) to
#flag outlier negative probes.
negp <- rowSums(seurats@assays$Negprob$data)
grubbs_test <- grubbs.test(negp, type = 10)
grubbs_test$p.value
message("No outliers, p val > 0.01.")

#Cell QC Summary
#Minimal counts per cell: recommend 50 or 100 for 6K panel; 20 for
#1000-plex panel; 5 for 100-plex
#panel; must be >1. Increase the threshold to make QC more conservative.
#Delete cells with less than 20 counts
#FLAG1
seurats$flag1 <- seurats$nCount_RNA > 20
message(paste0(
  "A total of ",
  sum(seurats$flag1 == 1),
  " cells have passed the LOW count filter.",
  sep = ""
))

#Proportion of negative counts: flag cells where >10% (0.1, the default value)
#of the counts per cell are
#negative probes. Decrease the threshold to make QC more conservative
#FLAG2
seurats$flag2 <-
  ifelse(
    seurats$nCount_negprobes == 0 & seurats$nCount_RNA == 0,
    0,
    seurats$nCount_negprobes / seurats$nCount_RNA < 0.1
  )
message(
  paste0(
    "A total of ",
    sum(seurats$flag2 == 1),
    " cells have passed the proportion of negative counts filter.",
    sep = ""
  )
)

#Complexity of cells
#Count distribution: flag cells where
#(total counts) / (number of detected genes) ≤1 (default value;
#range 1-200). In other words, total counts must exceed the number of detected
#genes in the cell. Increase threshold to make QC more conservative.
#Also, cells must have more than 10 feat
#FLAG3
seurats$flag3 <-
  ifelse(
    seurats$nCount_RNA == 0 & seurats$nFeature_RNA == 0,
    0,
    (seurats$nCount_RNA / seurats$nFeature_RNA) > 1
  )
seurats$flag3 <- seurats$flag3 * (seurats$nFeature_RNA > 10)
message(
  paste0(
    "A total of ",
    sum(seurats$flag3 == 1),
    " cells have passed the proportion of negative counts filter.",
    sep = ""
  )
)

#Area outlier: Grubb's test p-value (default: 0.01, range 0-1) to flag outlier
# cells based on cell area.
#FLAG4
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
seurats$flag4 <- seurats$Area <= max(data_vector)
message(paste0("A total of ", sum(seurats$flag4 == 1)),
        " cells have passed the Area Outlier filterfilter.",
        sep = "")

#Target level QC:
#• Negative control probe quantile cutoff: set the threshold at which to
#flag probes. A value of 0.5 (default) will flag probes with lower total counts
#than the median (50th percentile) of the negative control probes' counts.
#Range 0-1.
#FLAG5
negp <- seurats@assays$Negprob$data
seurats$flag5 <- colMeans(negp) < 0.5
message(paste0("A total of ", sum(seurats$flag5 == 1)),
        " cells have passed the target level QC.",
        sep = "")


#Do CELLS pass 5 FLAGS?
#Grab cells that pass the flagging
seurats$pass <-
  seurats$flag1 * seurats$flag2 * seurats$flag3 * seurats$flag4 * seurats$flag5
message(paste0("A total of ", sum(seurats$pass == 1)), " cells have passed the ENTIRE qc.", sep = "")


#SaveRDS with ALL cells
saveRDS(
  seurats,
  "~/SPATIAL/Mackensy_analysis/CosMx_RNA/0_Curation/Objects/qc_seurats.RDS"
)

#SaveRDS only cells that have passed QC
saveRDS(
  seurats[, rownames(seurats@meta.data[seurats@meta.data$pass == 1,])],
  "~/SPATIAL/Mackensy_analysis/CosMx_RNA/0_Curation/Objects/qc_seurats_pass.RDS"
)

