## Libraries--------------------------------------------------------------------
library(readr)
library(plyr)
library(dplyr)


## Read connections ------------------------------------------------------------
connections <- read_csv("analysis/CosMx_Protein/Objects/connections.csv")
meta <-
  read_csv(
    "analysis/CosMx_Protein/Objects/meta_SC_cut.csv"
  )
subset <- as.vector(meta$anot_maxfuse)
connections$subset <- subset[connections$mod1_indx + 1]
norm <-
  readRDS("analysis/CosMx_Protein/Objects/norm.RDS")
cell_names <- as.vector(norm$cell_id)
connections$cell_name <- cell_names[connections$mod2_indx + 1]

meta_p <-
  read.csv("analysis/CosMx_Protein/Objects/meta.csv")

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

norm@meta.data$subset <- mapvalues(x = norm@meta.data$cell_id,from = meta_p$cell_id,to = meta_p$subset)

write.csv(
  meta_p,
  "analysis/CosMx_Protein/Results/meta.csv"
)

saveRDS(norm,
        "analysis/CosMx_Protein/Objects/seurats_norm_all.RDS")

