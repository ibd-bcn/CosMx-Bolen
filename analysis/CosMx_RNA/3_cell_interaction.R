#Read libraries ----------------------------------------------------------------
library(reticulate)
library(Seurat)
library(readr)
library(ggplot2)
library(plyr)
library(dplyr)
library(ggrepel)

#Read object -------------------------------------------------------------------
seu <-
  readRDS(
    "~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Objects/seurats_all_norm.RDS"
  )

#Obtain files for each FOV to perform SCOTIA
genes <- rownames(seu)
meta_seu <- as.data.frame(seu@meta.data)
write_csv(meta_seu,
          "~/SPATIAL/Mackensy_analysis/CosMx_RNA/4_SCOTIA/Files/meta_seu.csv")

patients <- unique(seu$patient)

for (patient in patients) {
  #Obtain
  cell_id <-
    seu@meta.data[seu@meta.data$patient == patient,]$cell_id
  #Order
  meta <- meta_seu[cell_id,]
  fov <- meta$fov
  annotation <- meta$subset
  cords <- meta[, c("CenterX_global_px", "CenterY_global_px")]
  x_pos <- cords$CenterX_global_px
  y_pos <- cords$CenterY_global_px
  
  #Patient
  df <- data.frame(
    cell_id = cell_id,
    fov = fov,
    annotation = annotation,
    x_positions = x_pos,
    y_positions = y_pos
  )
  
  
  
  write_csv(
    df,
    paste0(
      "~/SPATIAL/Mackensy_analysis/CosMx_RNA/4_SCOTIA/Files/",
      patient,
      "_meta_def.csv",
      sep = ""
    )
  )
  
  #Exp
  data_ex <- as.data.frame(seu@assays$RNA$counts)
  data_ex <- as.data.frame(t(data_ex[, cell_id]))
  genes <- colnames(data_ex)
  data_ex$cell_id <- rownames(data_ex)
  data_ex$fov <- df$fov
  data_ex <- data_ex[c("cell_id", "fov", genes)]
  
  write_csv(
    data_ex,
    paste0(
      "~/SPATIAL/Mackensy_analysis/CosMx_RNA/4_SCOTIA/Files/",
      patient,
      "_exp_def.csv",
      sep = ""
    )
  )
  
}

## LR database -----------------------------------------------------------------
## Obtained from this database: https://github.com/ZJUFanLab/CellTalkDB
lr_pair <-
  readRDS("~/SPATIAL/Cell_neigh/Neigh_v2/SCOTIA/human_lr_pair.rds")
lr_pair <-
  lr_pair[lr_pair$ligand_gene_symbol %in% genes &
            lr_pair$receptor_gene_symbol %in% genes,]

lr_pair <-
  lr_pair[, c("ligand_gene_symbol", "receptor_gene_symbol")]
colnames(lr_pair) <- c("l_gene", "r_gene")

write_csv(lr_pair,
          "~/SPATIAL/Mackensy_analysis/CosMx_RNA/4_SCOTIA/Files/lr_pair.csv")

## Run SCOTIA ------------------------------------------------------------------
system(
  "taskset -c 0,30 python3 ~/SPATIAL/Mackensy_analysis/CosMx_RNA/4_SCOTIA/scotia_run.py"
)


## Create a file for final interactions ----------------------------------------

meta_seu <-
  read_csv("/home/mmoro/SPATIAL/Mackensy_analysis/CosMx_RNA/4_SCOTIA/Files/meta_seu.csv")
col_names <-
  c(
    "source_cell_idx",
    "receptor_cell_idx",
    "likelihood",
    "ligand_recptor",
    "source_cell_type",
    "target_cell_type",
    "cell_pairs",
    "id_source",
    "id_receptor",
    "subset_source",
    "subset_receptor"
  )

df <- data.frame(matrix(ncol = length(col_names), nrow = 0))
colnames(df) <- col_names
patients <- unique(meta_seu$patient)

for (patient in patients) {
  fovs <- unique(meta_seu[meta_seu$patient == patient, ]$fov)
  
  for (fov in fovs) {
    if (file.exists(
      paste0(
        "/home/mmoro/SPATIAL/Mackensy_analysis/CosMx_RNA/4_SCOTIA/Files/",
        patient,
        "_fov_",
        fov,
        ".ot.csv",
        sep = ""
      )
    )) {
      meta <-
        read_delim(
          paste0(
            "/home/mmoro/SPATIAL/Mackensy_analysis/CosMx_RNA/4_SCOTIA/Files/",
            patient,
            "_fov_",
            fov,
            ".csv",
            sep = ""
          ),
          delim = "\t",
          escape_double = FALSE,
          trim_ws = TRUE
        )
      
      ot <-
        read_delim(
          paste0(
            "/home/mmoro/SPATIAL/Mackensy_analysis/CosMx_RNA/4_SCOTIA/Files/",
            patient,
            "_fov_",
            fov,
            ".ot.csv",
            sep = ""
          ),
          delim = "\t",
          escape_double = FALSE,
          trim_ws = TRUE
        )
      
      ot$id_source <-
        mapvalues(
          x = ot$source_cell_idx,
          from = meta$index ,
          to = meta$cell_id,
          warn_missing = F
        )
      ot$id_receptor <-
        mapvalues(
          x = ot$receptor_cell_idx,
          from = meta$index ,
          to = meta$cell_id,
          warn_missing = F
        )
      ot$refined_source <-
        mapvalues(
          x = ot$id_source,
          from = meta_seu$cell ,
          to = meta_seu$refined,
          warn_missing = F
        )
      ot$refined_receptor <-
        mapvalues(
          x = ot$id_receptor,
          from = meta_seu$cell ,
          to = meta_seu$refined,
          warn_missing = F
        )
      ot$x_receptor <-
        as.numeric(
          mapvalues(
            x = ot$id_receptor,
            from = meta$cell_id ,
            to = meta$x_positions,
            warn_missing = F
          )
        )
      ot$y_receptor <-
        as.numeric(
          mapvalues(
            x = ot$id_receptor,
            from = meta$cell_id ,
            to = meta$y_positions,
            warn_missing = F
          )
        )
      ot$x_source <-
        as.numeric(
          mapvalues(
            x = ot$id_source,
            from = meta$cell_id ,
            to = meta$x_positions,
            warn_missing = F
          )
        )
      ot$y_source <-
        as.numeric(
          mapvalues(
            x = ot$id_source,
            from = meta$cell_id ,
            to = meta$y_positions,
            warn_missing = F
          )
        )
      
      ot$fov <-
        as.numeric(
          mapvalues(
            x = ot$id_source,
            from = meta_seu$cell ,
            to = meta_seu$fov,
            warn_missing = F
          )
        )
      
      
      
      
      df <- rbind(df, ot)
    }
    
    
  }
  
}

#Save
write_csv(
  df,
  "/home/mmoro/SPATIAL/Mackensy_analysis/CosMx_RNA/4_SCOTIA/Results/all_int.csv"
)


