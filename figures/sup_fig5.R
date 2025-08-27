# ============================================================================ #
#   CosMx Analysis - Sup Figure 5 Script
# ============================================================================ #

# Libraries -------------------------------------------------------------------
library(Seurat)
library(ggplot2)

# Read objects ----------------------------------------------------------------
todas     <- readRDS("data/single_cell_ref.RDS")
todas_cut <-
  readRDS("analysis/CosMx_Protein/Objects/sc_ref_cut.RDS")

# ============================================================================ #
#   Sup Figure 5
# ============================================================================ #

## A ---------------------------------------------------------------------------
pA <- DimPlot(todas,
              group.by = "annotation_intermediate",
              label = TRUE) + theme_bw()

ggsave(
  "figures/plots/supfig5A.png",
  plot = pA,
  width = 7,
  height = 6,
  dpi = 300
)

## B ---------------------------------------------------------------------------
pB <- DimPlot(todas_cut,
              group.by = "anot_maxfuse",
              label = TRUE) + theme_bw()

ggsave(
  "figures/plots/supfig5B.png",
  plot = pB,
  width = 7,
  height = 6,
  dpi = 300
)

## C ---------------------------------------------------------------------------
genes <- c(
  "TNFRSF9",
  "CTNNB1",
  "CCR7",
  "PECAM1",
  "CD38",
  "ENTPD1",
  "IL7R",
  "SDC1",
  "CD14",
  "PTPRC",
  "CD68",
  "CD8A",
  "MS4A1",
  "CD27",
  "CD3E",
  "EPCAM",
  "FABP4",
  "FN1",
  "TNFRSF18",
  "GZMA",
  "GZMB",
  "LAMP1",
  "ACTA2",
  "HLA-DRA",
  "ICAM1",
  "ICOS",
  "IL1B",
  "IGHD",
  "MKI67"
)

for (g in genes) {
  pC <- FeaturePlot(todas, features = g) + theme_bw()
  ggsave(
    paste0("figures/plots/supfig5C_", g, ".png"),
    plot = pC,
    width = 6,
    height = 5,
    dpi = 300
  )
}
