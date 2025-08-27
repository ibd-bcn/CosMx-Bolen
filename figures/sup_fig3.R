# ============================================================================ #
#   CosMx Analysis - Sup Figure 3 Script
# ============================================================================ #

# Libraries -------------------------------------------------------------------
library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(ggplot2)

# Read objects ----------------------------------------------------------------
seurat <- readRDS("analysis/CosMx_RNA/Objects/seurats_all_norm.RDS")

# ============================================================================ #
#   Sup Figure 3
# ============================================================================ #

#LogNormalize
seurat <- NormalizeData(seurat)
#Extract Markers
markers <- RunPrestoAll(seurat)
#Extact top markers
ref.markers.ss <- markers[
  markers$pct.1 > 0.20 & markers$pct.2 < 0.05,
]
# Define the number of top markers you want per cluster
top_n <- 5
# Function to get top markers for each cluster
get_top_markers <- function(df, top_n) {
  df %>%
    group_by(cluster) %>%
    arrange(cluster, desc(avg_log2FC), desc(pct.1), desc(pct.2), p_val_adj) %>%
    slice_head(n = top_n) %>%
    ungroup()
}
# Get the top markers for each cluster
top_markers <- get_top_markers(markers, top_n)
dot_plot <- DotPlot(object = seu,features = top_markers$gene,group.by = "subset",assay = "RNA")
dot_plot <- dot_plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("figures/plots/supfig3.png", plot = dot_plot, width = 8, height = 6, dpi = 300)
