# ============================================================================ #
#   CosMx Analysis - Figure 4 Script
# ============================================================================ #

# Libraries -------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(plyr)
library(readr)
library(circlize)
library(paletteer)
library(BiocParallel)

# Read objects ----------------------------------------------------------------
seurat <- readRDS("analysis/CosMx_RNA/Objects/seurats_all_norm.RDS")
meta   <- seurat@meta.data
all_int <- read_csv("analysis/CosMx_RNA/Results/all_int.csv")

# Add refined annotations
all_int$refined_receptor <- mapvalues(all_int$id_receptor, from = meta$cell_id, to = meta$refined)
all_int$refined_source   <- mapvalues(all_int$id_source,   from = meta$cell_id, to = meta$refined)
all_int$tissue           <- mapvalues(all_int$id_receptor, from = meta$cell_id, to = meta$tissue)

# ============================================================================ #
#   Helper functions
# ============================================================================ #

# Chord diagram matrix --------------------------------------------------------
chordiagram <- function(chord, meta_seu, fov = "all",
                        subset_source = "all", subset_target = "all",
                        refined_source = "all", refined_receptor = "all",
                        freq = 1, down = FALSE) {
  if ("all" %in% subset_source) subset_source <- unique(meta_seu$subset)
  if ("all" %in% subset_target) subset_target <- unique(meta_seu$subset)
  if ("all" %in% refined_source) refined_source <- unique(meta_seu$refined)
  if ("all" %in% refined_receptor) refined_receptor <- unique(meta_seu$refined)
  if ("all" %in% fov) fov <- unique(meta_seu$fov)

  # Select fov
  pat_cells <- meta_seu[meta_seu$fov %in% fov, ]$cell

  chord <- chord[chord$id_source %in% pat_cells &
                   chord$source_cell_type %in% subset_source &
                   chord$target_cell_type %in% subset_target &
                   chord$refined_source %in% refined_source &
                   chord$refined_receptor %in% refined_receptor, ]

  chord <- as.data.frame(table(chord$refined_source, chord$refined_receptor))
  chord <- chord[chord$Freq > freq, ]
  chord_1 <- chord %>% spread(key = Var2, value = Freq)
  rownames(chord_1) <- chord_1$Var1
  chord_1 <- chord_1[, -1]

  if (down) return(chord)
  return(chord_1)
}

# Chord plot ------------------------------------------------------------------
chord_plot <- function(meta_seu, mat, cell_type, color_chord) {
  cells <- unique(meta_seu$refined)
  names(color_chord) <- cells

  col_mat <- matrix("#F0F0F0", nrow = nrow(mat), ncol = ncol(mat),
                    dimnames = list(rownames(mat), colnames(mat)))

  for (z in rownames(mat)) col_mat[z, ] <- color_chord[z]

  if (!("all" %in% cell_type)) {
    cols <- setdiff(colnames(col_mat), cell_type)
    rows <- setdiff(rownames(col_mat), cell_type)
    col_mat[rows, cols] <- "#F0F0F0"
  }

  chordDiagram(mat,
               annotationTrack = "grid", preAllocateTracks = 1,
               grid.col = color_chord, col = col_mat)

  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim  <- get.cell.meta.data("xlim")
    ylim  <- get.cell.meta.data("ylim")
    sector <- get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1] + .1, sector, facing = "clockwise",
                niceFacing = TRUE, adj = c(0, 0.5), cex = 0.6)
    circos.axis(h = "top", labels.cex = 0.5,
                major.tick.length = 0.2, sector.index = sector, track.index = 2)
  }, bg.border = NA)
}

#Connection map ----------------------------------------------------------------
plot_spatial_connections <- function(seurat_obj, fov_id, radius = 139, workers = 35, output_file) {
   # Extract meta and subset by FOV
  meta <- seurat_obj@meta.data
  meta <- meta[meta$fov == fov_id, ]

  # Add FERT info
  fert <- seurat_obj@assays$RNA$counts["FERT", rownames(meta)]
  meta$fert <- ifelse(fert == 0, "neg", "pos")

  # Setup parallel
  bp <- MulticoreParam(workers = workers, progressbar = TRUE)
  cells <- meta$cell_id

  # Build network
  network <- BiocParallel::bplapply(cells, BPPARAM = bp, \(c) {
    x <- meta[c, "CenterX_global_px"]
    y <- meta[c, "CenterY_global_px"]
    tt <- sqrt((meta$CenterX_global_px - x)^2 + (meta$CenterY_global_px - y)^2)
    cn <- rownames(meta)
    your_dataframe <- data.frame(Column1 = cn, Column2 = as.vector(tt))
    ordered_df <- your_dataframe[order(your_dataframe$Column2), ]
    keep <- ordered_df[ordered_df$Column2 < radius, ]
    return(keep$Column1)
  })

  # Connection matrix
  lines_df <- BiocParallel::bplapply(1:length(network), BPPARAM = bp, \(c) {
    conn <- network[[c]]
    start_cell <- conn[1]
    connected_cells <- conn[-1]
    start_coords <- meta %>% filter(cell_id == start_cell)
    end_coords <- meta %>% filter(cell_id %in% connected_cells)

    if (nrow(end_coords) != 0) {
      df <- data.frame(
        x_start = start_coords$CenterX_global_px,
        y_start = start_coords$CenterY_global_px,
        cell_start = start_coords$refined,
        x_end = end_coords$CenterX_global_px,
        y_end = end_coords$CenterY_global_px,
        cell_end = end_coords$refined,
        fert_end = end_coords$fert,
        fert_start = start_coords$fert
      )
      return(df)
    }
  })

  big_df <- do.call(rbind, lines_df)

  # Adjust cell_start for plotting
  big_df$cell_start <- ifelse(big_df$cell_start == "Colonocytes", "Colonocytes", "Other")
  big_df <- big_df %>%
    mutate(cell_start = case_when(
      cell_start == "Colonocytes" & fert_end == "pos" ~ "Colonocytes+",
      cell_start == "Colonocytes" & fert_end == "neg" ~ "Colonocytes-",
      TRUE ~ cell_start
    ))

  # Plot
  p <- ggplot() +
    geom_point(data = meta, aes(x = CenterX_global_px, y = CenterY_global_px),
               color = "black", size = 0.01) +
    geom_segment(data = big_df[big_df$cell_start == "Colonocytes+", ],
                 aes(x = x_start, y = y_start, xend = x_end, yend = y_end,
                     color = "blue", linewidth = 0.2, alpha = 1)) +
    geom_segment(data = big_df[big_df$cell_start == "Colonocytes-", ],
                 aes(x = x_start, y = y_start, xend = x_end, yend = y_end,
                     color = "red", linewidth = 0.2, alpha = 1)) +
    geom_segment(data = big_df[big_df$cell_start == "Other", ],
                 aes(x = x_start, y = y_start, xend = x_end, yend = y_end,
                     color = "white", linewidth = 0.05, alpha = 0.8)) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "black", color = NA),
      plot.background = element_rect(fill = "black", color = NA),
      strip.background = element_rect(fill = "black"),
      strip.text = element_blank(),
      panel.grid = element_line(color = "gray40"),
      axis.text = element_text(color = "white"),
      axis.title = element_text(color = "white", face = "bold"),
      plot.title = element_text(color = "white", face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.text = element_text(color = "white")
    ) +
    labs(title = paste0("Spatial Connection Map - FOV ", fov_id),
         x = "X Coordinate", y = "Y Coordinate") +
    scale_color_identity() +
    scale_linewidth_identity() +
    scale_alpha_identity()

  # Save plot
  ggsave(filename = output_file, plot = p, width = 8, height = 8, dpi = 300)
}

# ============================================================================ #
#   Figure 4
# ============================================================================ #

## A - Chord plots ------------------------------------------------------------
color_palette <- c(
  paletteer_d("ggsci::default_igv"),
  paletteer_d("ggsci::category20_d3"),
  paletteer_d("ggsci::default_ucscgb")
)
color_chord <- color_palette[1:52]

conditions <- c("IBD", "PD", "NHC")

for (cond in conditions) {
  mat <- as.matrix(
    chordiagram(
      chord = all_int[all_int$tissue == cond, ],
      meta_seu = meta,
      fov = "all",
      subset_source = "all",
      subset_target = "all",
      refined_source = c("Colonocytes_FERT+", "Colonocytes_FERT-"),
      refined_receptor = c("B cell", "CD4", "CD8", "Eosinophils", "Glia",
                           "Inflammatory monocytes", "NK", "Neutrophil",
                           "Mast", "PC IgA", "PC IgG", "Myofibroblasts",
                           "M0", "M1", "M2"),
      freq = 1
    )
  )

  png(
    filename = paste0("figures/plots/fig4A_", cond, "_Chord.png"),
    width = 10, height = 8, units = "in", res = 300
  )
  chord_plot(meta_seu = meta, mat = mat, cell_type = "all", color_chord = color_chord)
  dev.off()
}

## B - Statistics between interactions ----------------------------------------
colonocyte_types <- c("Colonocytes_FERT-", "Colonocytes_FERT+")

df_filtered <- all_int %>%
  filter(refined_source %in% colonocyte_types | refined_receptor %in% colonocyte_types)

interaction_counts <- df_filtered %>%
  mutate(cell_type = case_when(
    refined_source %in% colonocyte_types ~ refined_source,
    refined_receptor %in% colonocyte_types ~ refined_receptor
  )) %>%
  group_by(tissue, cell_type, ligand_recptor) %>%
  summarise(count = n(), .groups = 'drop') %>%
  filter(cell_type == "Colonocytes_FERT+" & tissue %in% c("IBD", "NHC", "PD"))

long_df <- interaction_counts %>%
  pivot_wider(names_from = tissue, values_from = count, values_fill = list(count = 0)) %>%
  pivot_longer(cols = c("IBD", "PD", "NHC"), names_to = "Condition", values_to = "Value")

custom_colors <- c("IBD" = "#E57373", "NHC" = "#66BB6A", "PD" = "#42A5F5")

pB <- ggplot(long_df, aes(x = Condition, y = Value, fill = Condition)) +
  geom_boxplot(outlier.shape = NA, color = "black", size = 0.8, width = 0.6) +
  geom_point(shape = 16, position = position_jitter(0.1), size = 1, alpha = 0.5) +
  scale_fill_manual(values = custom_colors) +
  theme_classic() +
  theme(text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10)) +
  labs(title = "Number of interactions",
       y = "Number of interactions", x = "")

ggsave("figures/plots/fig4B_boxplot.png", plot = pB, width = 6, height = 4, dpi = 300)

#C - Spatial Connection maps ---------------------------------------------------
plot_spatial_connections(seurat_obj = seurat, fov_id = 26, output_file = "figures/plots/fig4C_fov26.png")
plot_spatial_connections(seurat_obj = seurat, fov_id = 10, output_file = "figures/plots/fig4Ci_fov10.png")
plot_spatial_connections(seurat_obj = seurat, fov_id = 13, output_file = "figures/plots/fig4Cii_fov13.png")
#D - Volcano plots -------------------------------------------------------------
pC <- volcano(anot = "refined", ct = "Colonocytes_FERT+", seu = seurat,
              id1 = "PD", id2 = "NHC", dif_col = "tissue")[[1]]
ggsave("figures/plots/fig4C_PDvsNHC.png", plot = pC, width = 6, height = 4, dpi = 300)

pD <- volcano(anot = "refined", ct = "Colonocytes_FERT+", seu = seurat,
              id1 = "IBD", id2 = "NHC", dif_col = "tissue")[[1]]
ggsave("figures/plots/fig4D_IBDvsNHC.png", plot = pD, width = 6, height = 4, dpi = 300)
