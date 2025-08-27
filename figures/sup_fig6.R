# ============================================================================ #
#   CosMx Analysis - Figure 1 Script
# ============================================================================ #

# Libraries -------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(plyr)
library(ggdark)

# Read objects ----------------------------------------------------------------
seurat <- readRDS("analysis/CosMx_RNA/Objects/seurats_all_norm.RDS")
meta_rna   <- seurat@meta.data
pols_rna   <- read.csv("analysis/CosMx_RNA/Polygons/pols.csv")

seurats_prot <-
  readRDS("analysis/CosMx_Protein/Objects/seurats_norm_all.RDS")
meta_prot <- seurats_prot@meta.data
pols_prot    <- read.csv("analysis/CosMx_Protein/Polygons/pols.csv")

# ============================================================================ #
#   Helper functions
# ============================================================================ #

# Plot polygons
plot_pol <- function(object,
                     fov,
                     poly,
                     mols,
                     pt_size,
                     annotation,
                     pal,
                     mols_c = FALSE,
                     genes = FALSE) {
  if (fov != "all") {
    object <- object[object$fov == fov, ]
  }

  cells <- object$cell
  poly  <- poly[poly$cell %in% cells, ]
  poly[[annotation]] <- mapvalues(x = poly$cell,
                                  from = object$cell,
                                  to   = object[[annotation]])

  p <- ggplot(poly, aes(x = x_global_px, y = y_global_px)) +
    geom_polygon(aes(group = cell, fill = .data[[annotation]]),
                 color = 'black')

  if (mols_c) {
    mols <- mols[mols$cell %in% cells, ]
    mols <- mols[mols$target %in% genes, ]
    p <- p + geom_point(
      data = mols,
      aes(x = x_global_px, y = y_global_px, color = target),
      size = pt_size
    )
  }

  p <- p +
    dark_theme_gray(base_family = "Fira Sans Condensed Light", base_size = 20) +
    scale_fill_manual(values = pal) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title       = element_blank(),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      panel.background = element_blank()
    ) +
    labs(x = "x", y = "y")

  return(p)
}

plot_spatial_dis_gene <-
  function(object,
           meta ,
           sample,
           fov,
           gene,
           pol_path,
           cell_name_col = "cell_names",
           mult = FALSE) {
    meta[[gene]] <- as.vector(object[gene,])
    #Modify Data
    if (sample == FALSE) {
      meta <- meta[meta$fov %in% fov,]



      poly <- pol_path

    }

    cells <- meta[[cell_name_col]]
    poly <- poly[poly[[cell_name_col]] %in% cells, ]
    poly[[gene]] <- as.numeric(mapvalues(x = poly[[cell_name_col]],
                                         from = meta[[cell_name_col]],
                                         to = meta[[gene]]))
    poly[[gene]] <- log(poly[[gene]] + 0.1)

    #Plot
    p <- ggplot(poly, aes(x = x_local_px, y = y_local_px)) +
      geom_polygon(aes(group = .data[[cell_name_col]], fill = .data[[gene]]),
                   color = 'black') +
      dark_theme_gray(base_family = "Fira Sans Condensed Light", base_size = 20) +  scale_fill_gradient(low = "white", high = "darkred") +
      facet_wrap(~ fov) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        legend.position = "right"
      ) +
      labs(x = "x",
           y = "y")

    return(p)
  }

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


proti <-
  c(
    "NA" = "grey",
    "Plasma" = "#FF5733FF",
    "Bcell" = "#FF00CCFF",
    "Epithelium" = "#33CC00FF",
    "Tcells" = "#ECFF00",
    "Fibroblasts" =  "#8A2BE2FF"
  )

# ============================================================================ #
#   Figure 1
# ============================================================================ #

## C - RNA ---------------------------------------------------------------------
pC1 <-
  plot_pol(
    meta_rna,
    fov = 10,
    poly = pols_rna,
    annotation = "subset",
    pal = refined_col
  )
ggsave(
  "figures/plots/sup_fig6c_sub.png",
  plot = pC1,
  width = 10,
  height = 6,
  dpi = 300
)

pC2 <-
  plot_pol(
    meta_rna,
    fov = 10,
    poly = pols_rna,
    annotation = "refined",
    pal = refined_col
  )
ggsave(
  "figures/plots/sup_fig6c_ref.png",
  plot = pC2,
  width = 10,
  height = 6,
  dpi = 300
)

## D - PROT---------------------------------------------------------------------
pD1 <-
  plot_pol(
    meta_prot,
    fov = 10,
    poly = pols_prot,
    annotation = "subset",
    pal = proti
  )
ggsave(
  "figures/plots/sup_fig6d_sub.png",
  plot = pD1,
  width = 10,
  height = 6,
  dpi = 300
)

pD2 <- plot_spatial_dis_gene(
  object = seurats_prot ,
  fov = 10,
  sample = FALSE,
  gene = "4-1BB",
  pol_path = pols_prot,
  cell_name_col = "cell",
  meta = meta_prot
)
ggsave(
  "figures/plots/sup_fig6d_exp.png",
  plot = pD2,
  width = 10,
  height = 6,
  dpi = 300
)
