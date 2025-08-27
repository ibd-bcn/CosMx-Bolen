# ============================================================================ #
#   CosMx Analysis - Figure 2 Script
# ============================================================================ #

# Libraries -------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(plyr)
library(ggdark)

# Read objects ----------------------------------------------------------------
seurat <- readRDS("analysis/CosMx_RNA/Objects/seurats_all_norm.RDS")
meta   <- seurat@meta.data
pols   <- read.csv("analysis/CosMx_RNA/Polygons/pols.csv")
mols   <- read.csv("analysis/CosMx_RNA/Molecules/mols.csv")

# ============================================================================ #
#   Helper functions
# ============================================================================ #

# Plot polygons ---------------------------------------------------------------
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
    object <- object[object$fov == fov,]
  }

  cells <- object$cell
  poly  <- poly[poly$cell %in% cells,]
  poly[[annotation]] <- mapvalues(x = poly$cell,
                                  from = object$cell,
                                  to   = object[[annotation]])

  p <- ggplot(poly, aes(x = x_global_px, y = y_global_px)) +
    geom_polygon(aes(group = cell, fill = .data[[annotation]]),
                 color = 'black')

  if (mols_c) {
    mols <- mols[mols$cell %in% cells,]
    mols <- mols[mols$target %in% genes,]
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

# Volcano plots ---------------------------------------------------------------
volcano <- function(anot = "subset",
                    ct,
                    dif_col = "tissue",
                    seu,
                    id1,
                    id2) {
  metadata <- seu@meta.data
  myeloid_metadata <- metadata[metadata[[anot]] == ct,]
  myeloid_counts_table <- table(myeloid_metadata[[dif_col]])

  object <-
    seu[, rownames(seu@meta.data[seu@meta.data[[anot]] == ct,])]
  object <- NormalizeData(object)
  object <- ScaleData(object)
  object@meta.data[[dif_col]] <-
    as.factor(object@meta.data[[dif_col]])
  object <- SetIdent(object, value = object@meta.data[[dif_col]])

  deg_results <- FindMarkers(object, ident.1 = id1, ident.2 = id2)
  deg_results <- na.omit(deg_results)
  deg_results$genes <- rownames(deg_results)

  deg_results$diffexpressed <- "NO"
  deg_results$diffexpressed[deg_results$avg_log2FC > log2(1.2) &
                              deg_results$p_val < 0.05] <- "p.val<0.05 & FC>1.2"
  deg_results$diffexpressed[deg_results$avg_log2FC < -log2(1.2) &
                              deg_results$p_val < 0.05] <- "p.val<0.05 & FC<0.83"
  deg_results$diffexpressed[deg_results$avg_log2FC > log2(1.2) &
                              deg_results$p_val_adj < 0.05] <- "p.adj<0.05 & FC>1.2"
  deg_results$diffexpressed[deg_results$avg_log2FC < -log2(1.2) &
                              deg_results$p_val_adj < 0.05] <- "p.adj<0.05 & FC<0.83"

  deg_results$delabel <- NA
  deg_results$delabel[deg_results$diffexpressed != "NO"] <-
    deg_results$genes[deg_results$diffexpressed != "NO"]

  deg_results$p_val <-
    ifelse(deg_results$p_val < 1e-300, 1e-300, deg_results$p_val)

  p <- ggplot(data = deg_results,
              aes(
                x = avg_log2FC,
                y = -log10(p_val),
                col = diffexpressed,
                label = delabel
              )) +
    geom_point(size = 3) +
    theme_bw() +
    geom_text_repel(size = 8, max.overlaps = 8) +
    scale_color_manual(
      values = c(
        "p.val<0.05 & FC<0.83"  = "#5cbde7",
        "p.adj<0.05 & FC<0.83" = "darkblue",
        "p.val<0.05 & FC>1.2"  = "red",
        "p.adj<0.05 & FC>1.2"  = "darkred"
      )
    ) +
    geom_vline(
      xintercept = c(-log2(1.2), log2(1.2)),
      col = "black",
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      col = "black",
      linetype = "dashed"
    ) +
    labs(color = "Legend") +
    theme(text = element_text(size = 25)) +
    guides(color = guide_legend(override.aes = list(size = 5)))

  return(list(p, deg_results))
}

# Density plots ---------------------------------------------------------------
density_plots <- function(gene, cells, seurat) {
  ccl18 <- seurat@assays$RNA$counts[gene,]
  meta <- seurat@meta.data
  meta$ccl18 <- ccl18
  meta <- meta[meta$refined %in% cells,]

  ccl18 <- meta %>%
    group_by(tissue, refined) %>%
    summarize(mean_ccl18 = mean(ccl18, na.rm = TRUE))

  p <- ggplot(ccl18,
              aes(
                x = refined,
                y = mean_ccl18,
                group = tissue,
                fill = tissue,
                color = tissue
              )) +
    geom_area(alpha = 0.5, position = 'identity') +
    geom_line(size = 1) +
    theme_bw() +
    labs(title = gene,
         y = paste0(gene, ": Average raw gene count")) +
    theme(
      text = element_text(size = 30),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank()
    )

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

# ============================================================================ #
#   Figure 2
# ============================================================================ #

## A --------------------------------------------------------------------------
p <- volcano(
  anot = "subset",
  ct   = "epi",
  seu  = seurat,
  id1  = "IBD",
  id2  = "NHC",
  dif_col = "tissue"
)[[1]]

ggsave(
  "figures/plots/fig2A.png",
  plot = p,
  width = 12,
  height = 9,
  dpi = 300
)

## Ai -------------------------------------------------------------------------
mmtt <- meta
mmtt[mmtt$subset != "epi",]$subset <- "other"
p <- plot_pol(
  object = mmtt,
  fov = "26",
  annotation = "subset",
  pal = subset_col,
  mols_c = TRUE,
  mols = mols,
  poly = pols,
  pt_size = 2,
  genes = "FERT"
)

ggsave(
  "figures/plots/fig2Ai.png",
  plot = p,
  width = 13,
  height = 10,
  dpi = 300
)

## Aii ------------------------------------------------------------------------
mmtt <- meta
mmtt[mmtt$subset != "epi",]$subset <- "other"
p <- plot_pol(
  object = mmtt,
  fov = "13",
  annotation = "subset",
  pal = subset_col,
  mols_c = TRUE,
  mols = mols,
  poly = pols,
  pt_size = 2,
  genes = "FERT"
)

ggsave(
  "figures/plots/fig2Aii.png",
  plot = p,
  width = 13,
  height = 10,
  dpi = 300
)

## B --------------------------------------------------------------------------
p <- volcano(
  anot   = "subset",
  ct     = "epi",
  seu    = seurat,
  id1    = "PD",
  id2    = "NHC",
  dif_col = "tissue"
)[[1]]

ggsave(
  "figures/plots/fig2B.png",
  plot = p,
  width = 12,
  height = 9,
  dpi = 300
)

## Bi -------------------------------------------------------------------------
mmtt <- meta
mmtt[mmtt$subset != "epi",]$subset <- "other"
p <- plot_pol(
  object = mmtt,
  fov = "10",
  annotation = "subset",
  pal = subset_col,
  mols_c = TRUE,
  mols = mols,
  poly = pols,
  pt_size = 2,
  genes = "FERT"
)

ggsave(
  "figures/plots/fig2Bi.png",
  plot = p,
  width = 13,
  height = 10,
  dpi = 300
)

## Bii ------------------------------------------------------------------------
mmtt <- meta
mmtt[mmtt$subset != "epi",]$subset <- "other"
p <- plot_pol(
  object = mmtt,
  fov = "30",
  annotation = "subset",
  pal = subset_col,
  mols_c = TRUE,
  mols = mols,
  poly = pols,
  pt_size = 2,
  genes = "FERT"
)

ggsave(
  "figures/plots/fig2Bii.png",
  plot = p,
  width = 13,
  height = 10,
  dpi = 300
)

## C --------------------------------------------------------------------------
p <- volcano(
  anot   = "refined",
  ct     = "Colonocytes",
  seu    = seurat,
  id1    = "IBD",
  id2    = "NHC",
  dif_col = "tissue"
)[[1]]

ggsave(
  "figures/plots/fig2C.png",
  plot = p,
  width = 12,
  height = 9,
  dpi = 300
)

## Ci -------------------------------------------------------------------------
mmtt <- meta
mmtt[mmtt$refined != "Colonocytes",]$refined <- "other"
p <- plot_pol(
  object = mmtt,
  fov = "26",
  annotation = "refined",
  pal = refined_col,
  mols_c = TRUE,
  mols = mols,
  poly = pols,
  pt_size = 2,
  genes = "FERT"
)

ggsave(
  "figures/plots/fig2Ci.png",
  plot = p,
  width = 13,
  height = 10,
  dpi = 300
)

## Cii ------------------------------------------------------------------------
mmtt <- meta
mmtt[mmtt$refined != "Colonocytes",]$refined <- "other"
p <- plot_pol(
  object = mmtt,
  fov = "13",
  annotation = "refined",
  pal = refined_col,
  mols_c = TRUE,
  mols = mols,
  poly = pols,
  pt_size = 2,
  genes = "FERT"
)

ggsave(
  "figures/plots/fig2Cii.png",
  plot = p,
  width = 13,
  height = 10,
  dpi = 300
)

## D --------------------------------------------------------------------------
p <- volcano(
  anot   = "refined",
  ct     = "Colonocytes",
  seu    = seurat,
  id1    = "PD",
  id2    = "NHC",
  dif_col = "tissue"
)[[1]]

ggsave(
  "figures/plots/fig2D.png",
  plot = p,
  width = 12,
  height = 9,
  dpi = 300
)

## Di -------------------------------------------------------------------------
mmtt <- meta
mmtt[mmtt$refined != "Colonocytes",]$refined <- "other"
p <- plot_pol(
  object = mmtt,
  fov = "10",
  annotation = "refined",
  pal = refined_col,
  mols_c = TRUE,
  mols = mols,
  poly = pols,
  pt_size = 2,
  genes = "FERT"
)

ggsave(
  "figures/plots/fig2Di.png",
  plot = p,
  width = 13,
  height = 10,
  dpi = 300
)

## Dii ------------------------------------------------------------------------
mmtt <- meta
mmtt[mmtt$refined != "Colonocytes",]$refined <- "other"
p <- plot_pol(
  object = mmtt,
  fov = "30",
  annotation = "refined",
  pal = refined_col,
  mols_c = TRUE,
  mols = mols,
  poly = pols,
  pt_size = 2,
  genes = "FERT"
)

ggsave(
  "figures/plots/fig2Dii.png",
  plot = p,
  width = 13,
  height = 10,
  dpi = 300
)

## E --------------------------------------------------------------------------
p <- density_plots(
  gene = "FERT",
  cells = c(
    "B cell",
    "CD4",
    "CD8",
    "Colonocytes",
    "Cycling TA",
    "DCs",
    "Endothelium",
    "Enteroendocrine",
    "Eosinophils",
    "Fibroblasts",
    "FRCs",
    "Glia",
    "Goblet",
    "Inflammatory monocytes",
    "M0",
    "M1",
    "M2",
    "Mast",
    "Myofibroblasts",
    "Neutrophil",
    "NK",
    "PC IgA",
    "PC IgG",
    "Pericytes",
    "S2"
  ),
  seurat = seurat
)

ggsave(
  "figures/plots/fig2E.png",
  plot = p,
  width = 15,
  height = 10,
  dpi = 300
)

## F --------------------------------------------------------------------------
p <- density_plots(
  gene = "SLC40A1",
  cells = c(
    "B cell",
    "CD4",
    "CD8",
    "Colonocytes",
    "Cycling TA",
    "DCs",
    "Endothelium",
    "Enteroendocrine",
    "Eosinophils",
    "Fibroblasts",
    "FRCs",
    "Glia",
    "Goblet",
    "Inflammatory monocytes",
    "M0",
    "M1",
    "M2",
    "Mast",
    "Myofibroblasts",
    "Neutrophil",
    "NK",
    "PC IgA",
    "PC IgG",
    "Pericytes",
    "S2"
  ),
  seurat = seurat
)

ggsave(
  "figures/plots/fig2F.png",
  plot = p,
  width = 15,
  height = 10,
  dpi = 300
)
