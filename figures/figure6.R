# ============================================================================ #
#   CosMx Analysis - Figure 6 Script
# ============================================================================ #

# Libraries -------------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(ggdark)
library(plyr)
library(ggpubr)
library(dplyr)
library(tidyr)
library(readr)
library(circlize)
library(paletteer)

# Read objects ----------------------------------------------------------------
akoya        <- read.csv("data/Akoya/slide40/D-6_measurements.csv")
seurats_rna  <-
  readRDS("analysis/CosMx_RNA/Objects/seurats_all_norm.RDS")
pols_rna     <- read.csv("analysis/CosMx_RNA/Polygons/pols.csv")
seurats_prot <-
  readRDS("analysis/CosMx_Protein/Objects/seurats_norm_all.RDS")
pols_prot    <- read.csv("analysis/CosMx_Protein/Polygons/pols.csv")
enrichment   <-
  read.csv("analysis/CosMx_Protein/Results/enrichment.csv")
all_int      <-
  read.csv("analysis/CosMx_Protein/Results/all_int.csv")
# ============================================================================ #
#   Helper functions
# ============================================================================ #

# Plot polygons ---------------------------------------------------------------
plot_pol <- function(object,
                     fov,
                     poly,
                     mols = NULL,
                     pt_size = 1,
                     annotation,
                     pal,
                     mols_c = FALSE,
                     genes = NULL) {
  if (fov != "all")
    object <- object[object$fov == fov,]
  cells <- object$cell
  poly  <- poly[poly$cell %in% cells,]
  poly[[annotation]] <- mapvalues(x = poly$cell,
                                  from = object$cell,
                                  to   = object[[annotation]])

  p <- ggplot(poly, aes(x = x_global_px, y = y_global_px)) +
    geom_polygon(aes(group = cell, fill = .data[[annotation]]), color = 'black')

  if (mols_c && !is.null(mols)) {
    mols <- mols[mols$cell %in% cells,]
    mols <- mols[mols$target %in% genes,]
    p <- p + geom_point(
      data = mols,
      aes(x = x_global_px, y = y_global_px, color = target),
      size = pt_size
    )
  }

  p + dark_theme_gray(base_family = "Fira Sans Condensed Light", base_size = 20) +
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
}

# Plot points for FOV ---------------------------------------------------------
plot_fov <- function(object,
                     fov,
                     pt_size = 1,
                     x = "Centroid X µm",
                     y = "Centroid Y µm",
                     annotation) {
  if (fov != "all")
    object <- object[object$fov == fov,]

  ggplot(object, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(aes(color = .data[[annotation]]), size = pt_size) +
    dark_theme_gray(base_family = "Fira Sans Condensed Light", base_size = 20) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_blank()
    )
}

# Chord diagram matrix --------------------------------------------------------
chordiagram_prot <-
  function(chord,
           meta_seu,
           fov,
           subset_source,
           subset_target,
           freq,
           down = FALSE) {
    if ("all" %in%  subset_source) {
      subset_source <-
        c("NA"  ,
          "Plasma"  ,
          "Epithelium"     ,
          "Fibroblasts" ,
          "Bcell" ,
          "Tcells")
    }
    if ("all" %in%  subset_target) {
      subset_target <-
        c("NA"  ,
          "Plasma"  ,
          "Epithelium"     ,
          "Fibroblasts" ,
          "Bcell" ,
          "Tcells")
    }
    if ("all" %in% fov) {
      fov <- c(1:66)

    }
    #Select fov
    pat_cells <- meta_seu[meta_seu$fov %in% fov, ]$cell

    chord <-
      chord[chord$id_source %in% pat_cells &
              chord$source_cell_type %in% subset_source &
              chord$target_cell_type %in% subset_target , ]

    chord <-
      as.data.frame(table(chord$source_cell_type, chord$target_cell_type))
    chord <- chord[chord$Freq > freq, ]
    chord_1 <- chord %>% spread(key = Var2, value = Freq)
    rownames(chord_1) <- chord_1$Var1
    chord_1 <- chord_1[, 2:ncol(chord_1)]

    if (down == TRUE) {
      return(chord)
    } else{
      return(chord_1)
    }


  }

# Chord plot ------------------------------------------------------------------
chord_plot_prot <- function(meta_seu , mat, cell_type, color_chord) {
  cells <- unique(meta_seu$subset)
  names(color_chord) <- cells

  col_mat <- matrix("#F0F0F0", nrow = nrow(mat), ncol = ncol(mat))
  colnames(col_mat) <- colnames(mat)
  rownames(col_mat) <- rownames(mat)

  for (z in rownames(mat)) {
    col_mat[z, ] <- color_chord[z]

  }

  if (!("all" %in% cell_type)) {
    cols <- setdiff(colnames(col_mat), c(cell_type))

    # Get all row names except "hola23" and "hola333"
    rows <- setdiff(rownames(col_mat), c(cell_type))

    col_mat[rows, cols] <- "#F0F0F0"


  }

  chordDiagram(
    mat,
    annotationTrack = "grid",
    preAllocateTracks = 1,
    grid.col = color_chord,
    col = col_mat
  )


  circos.trackPlotRegion(
    track.index = 1,
    panel.fun = function(x, y) {
      xlim = get.cell.meta.data("xlim")
      ylim = get.cell.meta.data("ylim")
      sector.name = get.cell.meta.data("sector.index")
      circos.text(
        mean(xlim),
        ylim[1] + .1,
        sector.name,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5)
      )
      circos.axis(
        h = "top",
        labels.cex = 0.5,
        major.tick.length = 0.2,
        sector.index = sector.name,
        track.index = 2
      )
    },
    bg.border = NA
  )




}

# ============================================================================ #
#   Palettes
# ============================================================================ #

refined_col <- c(
  # T cells
  "MT T cells" = "#5050FFFF",
  "CD8" = "#CE3D32FF",
  "Cycling T cells" = "#802268FF",
  "ILC4" = "#749B58FF",
  "NK" = "#466983FF",
  "T cells CCL20" = "#BA6338FF",
  "Tregs" = "#F0E685FF",
  "CD4" = "#5DB1DDFF",
  "MAIT" = "#6BD76BFF",
  "Ribhi T cells" = "#D595A7FF",
  "gd IEL" = "#00FFFFFF",
  "DN" = "#7A65A5FF",
  # Plasmas
  "PC IgG" = "#CC9900FF",
  "Memory B cell" = "#99CC00FF",
  "Cycling cells" = "#FF1463FF",
  "PC IgA" = "#0000CCFF",
  "PC IgA heat shock" = "#3B1B53FF",
  "PC IER" = "#CCCC99FF",
  "Naïve B cell" = "#FF0000FF",
  "B cell" = "#F7B6D2FF",
  "GC B cell" = "#990080FF",
  # Epithelium
  "Secretory progenitor" = "#FFFF00FF",
  "Epithelium Ribhi" = "#FF7F0EFF",
  "Cycling TA" = "#C75127FF",
  "Colonocytes" = "#9EDAE5FF",
  "Inflammatory colonocyte" = "#9467BDFF",
  "BEST4 OTOP2" = "#33CC00FF",
  "Goblet" = "#CC0000FF",
  "Enteroendocrine" = "#003399FF",
  "Tuft cells" = "#FFC20AFF",
  "Paneth-like" = "#FF00CCFF",
  # Myeloids
  "M2" = "#ECFF00",
  "M1" = "#FFCCCCFF",
  "M0" = "#924822FF",
  "DCs" = "#489C97",
  "Mast" = "#FF5733FF",
  "Cycling myeloid" = "#2E8B57FF",
  "IDA macrophage" = "#8A2BE2FF",
  "Inflammatory monocytes" = "#FF1493FF",
  "Neutrophil" = "#00FA9AFF",
  "Eosinophils" = "#FFD700FF",
  # Stroma
  "Endothelium" = "#4682B4FF",
  "Myofibroblasts" = "#7FFF00FF",
  "Pericytes" = "#00CED1FF",
  "S3" = "#BDB76BFF",
  "Fibroblasts" = "#556B2FFF",
  "Glia" = "#32CD32FF",
  "Inflammatory fibroblasts" = "#924822FF",
  "S1" = "#FF5733FF",
  "FRCs" = "#8A2BE2FF"
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
#   Figures
# ============================================================================ #

## A - Akoya + RNA + Protein polygons ----------------------------------------
# Akoya
pA1 <-
  plot_fov(object = akoya,
           fov = "all",
           annotation = "Classification")
ggsave(
  "figures/plots/fig6A_Akoya.png",
  plot = pA1,
  width = 12,
  height = 9,
  dpi = 300
)

# RNA polygons
pA2 <-
  plot_pol(
    object = seurats_rna@meta.data,
    fov = 36,
    poly = pols_rna,
    annotation = "subset",
    pal = refined_col
  )
ggsave(
  "figures/plots/fig6A_RNA.png",
  plot = pA2,
  width = 12,
  height = 9,
  dpi = 300
)

# Protein polygons
pA3 <-
  plot_pol(
    object = seurats_prot@meta.data,
    fov = 36,
    poly = pols_prot,
    annotation = "subset",
    pal = proti
  )
ggsave(
  "figures/plots/fig6A_Prot.png",
  plot = pA3,
  width = 12,
  height = 9,
  dpi = 300
)

## B - Protein UMAP ------------------------------------------------------------
#B UMAP obtained in Python following code of MaxFuse
# Code of MaxFuse present in celltyping: analysis/CosMx_Protein
# https://github.com/shuxiaoc/maxfuse

#Seurat UMAP Bi
pBi <- DimPlot(seurats_prot, group.by = "subset", cols = proti) +
  ggtitle("Protein UMAP (Seurat) - Python UMAP not included")
ggsave(
  "figures/plots/fig6B_UMAP.png",
  plot = pBi,
  width = 12,
  height = 9,
  dpi = 300
)

## C - Enrichment heatmap ----------------------------------------------------
pC <- ggplot(enrichment, aes(Comparison, Cluster, fill = e.score)) +
  geom_raster() +
  geom_text(aes(label = ast_Chisq)) +
  scale_fill_gradient2(
    low = 'blue',
    mid = 'white',
    high = 'red',
    limits = c(min(enrichment$e.score), max(enrichment$e.score)),
    midpoint = 0
  ) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank()) +
  labs(fill = 'Enrichment\nscore\nvs\nNHC')
ggsave(
  "figures/plots/fig6C_Enrichment.png",
  plot = pC,
  width = 12,
  height = 9,
  dpi = 300
)

##D - Chordplots ---------------------------------------------------------------
color_palette <- c(
  paletteer_d("ggsci::default_igv"),
  paletteer_d("ggsci::category20_d3"),
  paletteer_d("ggsci::default_ucscgb")
)
color_chord <- color_palette[1:52]
all_int$tissue <-
  mapvalues(all_int$id_source,
            seurats_prot@meta.data$cell_id,
            seurats_prot@meta.data$tissue)
conditions <- c("IBD", "PD", "NHC")
meta <- seurats_prot@meta.data

for (cond in conditions) {
  mat <- as.matrix(
    chordiagram_prot(
      chord = all_int[all_int$tissue == cond,],
      meta_seu = meta,
      fov = "all",
      subset_source = "Tcells",
      subset_target = "all",
      freq = 1
    )
  )

  png(
    filename = paste0("figures/plots/fig6D_", cond, "_Chord.png"),
    width = 10,
    height = 8,
    units = "in",
    res = 300
  )
  chord_plot_prot(
    meta_seu = meta,
    mat = mat,
    cell_type = "Epithelium",
    color_chord = color_chord
  )
  dev.off()
}
##E - Violin plots -------------------------------------------------------------
meta_p <- seurats_prot@meta.data
gene   <- c("CTLA4", "GITR", "LAG3", "PD-1", "Tim-3")
counts_p <- seurats_prot@assays$Prot$counts

for (g in gene) {
  # Prepare data
  ccl18 <- counts_p[g, meta_p$cell_id]
  meta_tt <- meta_p
  meta_tt$ccl18 <- log2(ccl18 + 1)
  meta_tt <- meta_tt[meta_tt$subset == "Tcells",]

  # File path for saving
  out_file <- paste0("figures/plots/fig6E_", g, ".png")

  # Open PNG device
  png(
    filename = out_file,
    width = 4,
    height = 6,
    units = "in",
    res = 600
  )

  # Plot
  p <- ggplot(meta_tt, aes(x = tissue, y = ccl18, fill = tissue)) +
    geom_violin(trim = TRUE, alpha = 0.9) +
    #geom_jitter(width = 0.1, size = 0.01, alpha = 0.2) +
    labs(
      x = 'Health',
      y = paste0('log2 ', g, ' fluorescence'),
      title = g
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 23,
      size = 3,
      fill = "black"
    ) +
    stat_summary(
      fun = mean,
      geom = "crossbar",
      width = 0.5,
      color = "black",
      fatten = 2
    ) +
    theme_minimal() +
    theme(legend.position = "none") +
    # Statistical tests
    stat_compare_means(method = "kruskal.test", label.y = max(meta_tt$ccl18) * 1.35) +
    stat_compare_means(
      comparisons = list(c("IBD", "NHC"),
                         c("NHC", "PD"),
                         c("IBD", "PD")),
      method = "wilcox.test",
      p.adjust.method = 'bonferroni'
    )
  print(p)

  # Close device
  dev.off()
}
