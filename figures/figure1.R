library(Seurats)
library(dplyr)
library(plyr)
library(ggdark)
source("~/SPATIAL/Mackensy_analysis/CosMx-Bolen_clean/figures/functions.R")


# Read objects------------------------------------------------------------------

#Seurat
seurat <-
  readRDS(
    "analysis/CosMx_RNA/Objects/seurats_all_norm.RDS"
  )
#Meta
meta <- seurat@meta.data
#Read polygons
pols <- read.csv("analysis/CosMx_RNA/Polygons/bolen_colon_RNA.csv")
#Enrichment
enrich <- read.csv("analysis/CosMx_RNA/Results/enrichment.csv")

#Functions----------------------------------------------------------------------
#Plot polygons
plot_pol <- function(object,
                     fov,
                     #patient,
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
  poly <- poly[poly$cell %in% cells, ]
  poly[[annotation]] <-
    mapvalues(x = poly$cell,
              from = object$cell,
              to = object[[annotation]])


  p <- ggplot(poly, aes(x = x_global_px, y = y_global_px)) +
    geom_polygon(aes(group = cell, fill = .data[[annotation]]),
                 color = 'black')

  if(mols_c == TRUE){
    mols <- mols[mols$cell %in% cells,]
    mols <- mols[mols$target %in% genes,]
    p <- p + geom_point(data = mols, aes(x = x_global_px, y =y_global_px, color = target), size = pt_size)

  }

  p <- p +
    dark_theme_gray(base_family = "Fira Sans Condensed Light", base_size = 20) +  scale_fill_manual(values = pal) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_blank()
    ) +
    labs(x = "x",
         y = "y")
  return(p)
}
#Palette of colours ------------------------------------------------------------
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


## Figure 1 --------------------------------------------------------------------

##A ----------------------------------------------------------------------------
p <- DimPlot(seurat, group.by = "refined",cols = refined_col)
ggsave(
  "figures/plots/fig1A.png",
  plot = p,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 10,
  height = 10,
  dpi = 300,

)

##B-----------------------------------------------------------------------------
df <-  data.frame(
  matrix(0, nrow = length(unique(meta[["subset"]]))*length(unique(meta[["refined"]])) , ncol = 3)
)

colnames(df) <- c("subset","refined", "value")


val1 <- unique(meta[["subset"]])
val2 <- unique(meta[["refined"]])

df[["subset"]]<- rep(val1, each = length(unique(meta[["refined"]])))
df[["refined"]] <- rep(val2, times= length(unique(meta[["subset"]])))

df[["refined"]] <- as.factor(df[["refined"]] )
df[["subset"]] <- as.factor(df[["subset"]])

for(x in val1){
  for(p in val2){
    num_of_cells <- nrow(meta[meta[["subset"]] == x & meta[["refined"]] == p , ])

    df[which(df[["subset"]] == x & df[["refined"]] == p),3] <- num_of_cells

  }


}

p1 <- ggplot(df, aes(fill = .data[["refined"]], y = value , x = .data[["subset"]])) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = refined_col) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), text = element_text(size = 10))  # Add legend title

ggsave(
  "figures/plots/fig1B.png",
  plot = p1,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 6,
  height = 6,
  dpi = 300,

)

##C ----------------------------------------------------------------------------
#Fov 10
p3 <- plot_pol(object = meta,fov =10,poly = pols,annotation = "refined",pal = refined_col )
ggsave(
  "figures/plots/fig1C_fov10.png",
  plot = p3,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 10,
  height = 6,
  dpi = 300,
)
#Fov 12
p4 <- plot_pol(object = meta,fov =12,poly = pols,annotation = "refined",pal = refined_col )
ggsave(
  "figures/plots/fig1C_fov12.png",
  plot = p4,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 10,
  height = 6,
  dpi = 300,
)
#Fov 30
p5 <- plot_pol(object = meta,fov =10,poly = pols,annotation = "refined",pal = refined_col )
ggsave(
  "figures/plots/fig1C_fov30.png",
  plot = p5,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 10,
  height = 6,
  dpi = 300,
)


##D-----------------------------------------------------------------------------
fused <- fused[,c(
  "subset",
  "refined",
  "tissue",
  "patient")]

## Refined ---------------------------------------------------------------------
df <-  data.frame(
  matrix(0, nrow = length(unique(meta[["tissue"]]))*length(unique(meta[["refined"]])) , ncol = 3)
)

colnames(df) <- c("tissue","refined", "value")


val1 <- unique(meta[["tissue"]])
val2 <- unique(meta[["refined"]])

df[["tissue"]]<- rep(val1, each = length(unique(meta[["refined"]])))
df[["refined"]] <- rep(val2, times= length(unique(meta[["tissue"]])))

df[["refined"]] <- as.factor(df[["refined"]] )
df[["tissue"]] <- as.factor(df[["tissue"]])

for(x in val1){
  for(p in val2){
    num_of_cells <- nrow(meta[meta[["tissue"]] == x & meta[["refined"]] == p , ])

    df[which(df[["tissue"]] == x & df[["refined"]] == p),3] <- num_of_cells

  }


}

p6 <- ggplot(df, aes(fill = .data[["refined"]], y = value , x = .data[["tissue"]])) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = refined_col) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), text = element_text(size = 10))  # Add legend title

ggsave(
  "figures/plots/fig1B.png",
  plot = p6,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 6,
  height = 6,
  dpi = 300,

)

##E-----------------------------------------------------------------------------
p7 <- ggplot(enrich, aes(Comparison, Cluster, fill= e.score)) +
  geom_raster() +
  geom_text(mapping = aes(label = ast_Chisq)) +
  scale_fill_gradient2(low = 'blue',
                       mid = 'white',
                       high = 'red',
                       limits = c(min(df_final$e.score),max(df_final$e.score)),
                       midpoint = 0) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank()) +
  labs(fill = 'Enrichment\nscore\nvs\nNHC')

ggsave(
  "figures/plots/fig1E.png",
  plot = p6,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 6,
  height = 6,
  dpi = 300,

)

