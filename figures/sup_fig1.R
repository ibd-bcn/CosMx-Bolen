# ============================================================================ #
#   CosMx Analysis - Sup Figure 1 Script
# ============================================================================ #

# Libraries -------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(plyr)

# Read objects ----------------------------------------------------------------
seurat <- readRDS("analysis/CosMx_RNA/Objects/seurats_all_norm.RDS")
meta   <- seurat@meta.data

# ============================================================================ #
#   Palette
# ============================================================================ #

refined_col <- c(
  # --- T cells ---
  "MT T cells"       = "#5050FFFF",
  "CD8"              = "#CE3D32FF",
  "Cycling T cells"  = "#802268FF",
  "ILC4"             = "#749B58FF",
  "NK"               = "#466983FF",
  "T cells CCL20"    = "#BA6338FF",
  "Tregs"            = "#F0E685FF",
  "CD4"              = "#5DB1DDFF",
  "MAIT"             = "#6BD76BFF",
  "Ribhi T cells"    = "#D595A7FF",
  "gd IEL"           = "#00FFFFFF",
  "DN"               = "#7A65A5FF",

  # --- Plasmas ---
  "PC IgG"           = "#CC9900FF",
  "Memory B cell"    = "#99CC00FF",
  "Cycling cells"    = "#FF1463FF",
  "PC IgA"           = "#0000CCFF",
  "PC IgA heat shock" = "#3B1B53FF",
  "PC IER"           = "#CCCC99FF",
  "Naïve B cell"     = "#FF0000FF",
  "B cell"           = "#F7B6D2FF",
  "GC B cell"        = "#990080FF",

  # --- Epithelium ---
  "Secretory progenitor"    = "#FFFF00FF",
  "Epithelium Ribhi"        = "#FF7F0EFF",
  "Cycling TA"              = "#C75127FF",
  "Colonocytes"             = "#9EDAE5FF",
  "Inflammatory colonocyte" = "#9467BDFF",
  "BEST4 OTOP2"             = "#33CC00FF",
  "Goblet"                  = "#CC0000FF",
  "Enteroendocrine"         = "#003399FF",
  "Tuft cells"              = "#FFC20AFF",
  "Paneth-like"             = "#FF00CCFF",

  # --- Myeloids ---
  "M2"                 = "#ECFF00",
  "M1"                 = "#FFCCCCFF",
  "M0"                 = "#924822FF",
  "DCs"                = "#489C97",
  "Mast"               = "#FF5733FF",
  "Cycling myeloid"    = "#2E8B57FF",
  "IDA macrophage"     = "#8A2BE2FF",
  "Inflammatory monocytes" = "#FF1493FF",
  "Neutrophil"         = "#00FA9AFF",
  "Eosinophils"        = "#FFD700FF",

  # --- Stroma ---
  "Endothelium"             = "#4682B4FF",
  "Myofibroblasts"          = "#7FFF00FF",
  "Pericytes"               = "#00CED1FF",
  "S3"                      = "#BDB76BFF",
  "Fibroblasts"             = "#556B2FFF",
  "Glia"                    = "#32CD32FF",
  "Inflammatory fibroblasts" = "#924822FF",
  "S1"                      = "#FF5733FF",
  "FRCs"                    = "#8A2BE2FF"
)

# ============================================================================ #
#   Sup Figure 1
# ============================================================================ #

df <-
  data.frame(matrix(0, nrow = length(unique(meta[["patient"]])) *
                      length(unique(meta[["refined"]])),
                    ncol = 3))
colnames(df) <- c("patient", "refined", "value")

val1 <- unique(meta[["patient"]])
val2 <- unique(meta[["refined"]])

df[["patient"]]  <- rep(val1, each  = length(val2))
df[["refined"]] <- rep(val2, times = length(val1))
df[["refined"]] <- as.factor(df[["refined"]])
df[["patient"]]  <- as.factor(df[["patient"]])

for (x in val1) {
  for (p in val2) {
    num_of_cells <-
      nrow(meta[meta[["patient"]] == x & meta[["refined"]] == p,])
    df[which(df[["patient"]] == x &
               df[["refined"]] == p), "value"] <- num_of_cells
  }
}
df$patient <- factor(df$patient,
                     levels = paste0("patient", 1:33)  # or use your vector of unique patients)
                     p <-
                       ggplot(df, aes(fill = .data[["refined"]], y = value, x = .data[["patient"]])) +
                       geom_bar(position = "fill", stat = "identity") +
                       scale_fill_manual(values = refined_col) +
                       theme_bw() +
                       theme(axis.text.x = element_text(
                         angle = 90,
                         hjust = 1,
                         vjust = 0.5
                       ),
                       text = element_text(size = 10))

                     ggsave(
                       "figures/plots/supfig1.png",
                       plot = p,
                       width = 8,
                       height = 6,
                       dpi = 300
                     )
