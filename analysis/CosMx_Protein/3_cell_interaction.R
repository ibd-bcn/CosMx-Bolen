##Libraries --------------------------------------------------------------------
library(dplyr)
library(biomaRt)
library(reticulate)
library(Seurat)
library(readr)
library(ggplot2)
library(plyr)
library(dplyr)
library(ggrepel)
library(readxl)

## Prepare protein data --------------------------------------------------------
cos <- read_excel("CosMx_Prot/3_SCOTIA/Ensembl Cosmx protein.xlsx")
dic <- cos$`Protein Name`
names(dic) <- cos$`Ensembl Prot. code`
df <-
  read.table(
    "~/SPATIAL/Mackensy_analysis/CosMx_Prot/3_SCOTIA/9606.protein.links.v12.0.txt",
    header = TRUE,
    sep = " ",
    stringsAsFactors = FALSE
  )
df$protein1 <- sub("9606\\.", "", df$protein1)
df$protein2 <- sub("9606\\.", "", df$protein2)

df <-
  df[df$protein1 %in% cos$`Ensembl Prot. code` &
       df$protein2 %in% cos$`Ensembl Prot. code`, ]
df$protein1 <- dic[df$protein1]
df$protein2 <- dic[df$protein2]

norm <-
  readRDS("~/SPATIAL/Mackensy_analysis/CosMx_Prot/01_Normalization/norm.RDS")

norm <- norm[unique(c(df$protein1, df$protein2)), ]
genes <- rownames(norm)

df_unique <- df[, c(1, 2)]
write.csv(df_unique,
          "~/SPATIAL/Mackensy_analysis/CosMx_Prot/3_SCOTIA/db.csv")

## Create intermediate files ---------------------------------------------------
meta <- read_csv("CosMx_Prot/02_celltyping/meta.csv")
seu <- norm
seu$subset <-
  mapvalues(x = seu$cell_id ,
            from = meta$cell_id,
            to = meta$subset)
genes <- rownames(seu)
meta_seu <- as.data.frame(seu@meta.data)
write_csv(meta_seu,
          "~/SPATIAL/Mackensy_analysis/CosMx_Prot/3_SCOTIA/Files/meta_seu.csv")
patients <- unique(seu$patient)

for (patient in patients) {
  #Obtain
  cell_id <- seu@meta.data[seu@meta.data$patient == patient, ]$cell_id
  #Order
  meta <- meta_seu[cell_id, ]
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
      "~/SPATIAL/Mackensy_analysis/CosMx_Prot/3_SCOTIA/Files/",
      patient,
      "_meta_def.csv",
      sep = ""
    )
  )
  
  #Exp
  data_ex <- as.data.frame(seu@assays$Prot$data)
  data_ex <- as.data.frame(t(data_ex[, cell_id]))
  genes <- colnames(data_ex)
  data_ex$cell_id <- rownames(data_ex)
  data_ex$fov <- df$fov
  data_ex <- data_ex[c("cell_id", "fov", genes)]
  
  write_csv(
    data_ex,
    paste0(
      "~/SPATIAL/Mackensy_analysis/CosMx_Prot/3_SCOTIA/Files/",
      patient,
      "_exp_def.csv",
      sep = ""
    )
  )
  
}

## LR database -----------------------------------------------------------------
## Obtained from this database: https://github.com/ZJUFanLab/CellTalkDB
colnames(df_unique) <- c("l_gene", "r_gene")

write_csv(df_unique,
          "~/SPATIAL/Mackensy_analysis/CosMx_Prot/3_SCOTIA/Files/lr_pair.csv")

## Run SCOTIA ------------------------------------------------------------------
system(
  "taskset -c 0-20 python3 ~/SPATIAL/Mackensy_analysis/CosMx_Prot/3_SCOTIA/scotia_cell_int.py"
)

## Create a file for final interactions ----------------------------------------

meta_seu <-
  read_csv("/home/mmoro/SPATIAL/Mackensy_analysis/CosMx_Prot/3_SCOTIA/Files/meta_seu.csv")
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
        "/home/mmoro/SPATIAL/Mackensy_analysis/CosMx_Prot/3_SCOTIA/Files/",
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
            "/home/mmoro/SPATIAL/Mackensy_analysis/CosMx_Prot/3_SCOTIA/Files/",
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
            "/home/mmoro/SPATIAL/Mackensy_analysis/CosMx_Prot/3_SCOTIA/Files/",
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

#Save final result
write_csv(
  df,
  "/home/mmoro/SPATIAL/Mackensy_analysis/CosMx_Prot/3_SCOTIA/Results/all_int.csv"
)

