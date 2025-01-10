## Libraries--------------------------------------------------------------------
library(readr)
library(plyr)
library(dplyr)


## Read connections ------------------------------------------------------------
connections <- read_csv("CosMx_Prot/connections.csv")
meta <-
  read_csv(
    "/home/mmoro/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Objects/meta_SC_cut.csv"
  )
subset <- as.vector(meta$anot_maxfuse)
connections$subset <- subset[connections$mod1_indx + 1]
norm <-
  readRDS("~/SPATIAL/Mackensy_analysis/CosMx_Prot/01_Normalization/norm.RDS")
cell_names <- as.vector(norm$cell_id)
connections$cell_name <- cell_names[connections$mod2_indx + 1]

meta_p <-
  read.csv("/home/mmoro/SPATIAL/Mackensy_analysis/CosMx_Prot/0_Curation/Objects/meta.csv")

# Perform a left join to map cell_id to subset, leaving unmatched as NA
meta_p$subset <- plyr::mapvalues(x = meta_p$cell_id,
                                 from = connections$cell_name,
                                 to = connections$subset)
meta_p$subset <-
  ifelse(
    test = meta_p$subset %in% unique(connections$subset),
    yes = meta_p$subset,
    no = NA
  )

write.csv(
  meta_p,
  "/home/mmoro/SPATIAL/Mackensy_analysis/CosMx_Prot/02_celltyping/meta.csv"
)