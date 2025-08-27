# ============================================================================ #
#   CosMx Analysis - Sup Figure 7 Script
# ============================================================================ #

# Libraries -------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(plyr)

# Read objects ----------------------------------------------------------------
seurats_prot <-
  readRDS("analysis/CosMx_Protein/Objects/seurats_norm_all.RDS")
meta <- seurats_prot@meta.data

# ============================================================================ #
#   Palettes
# ============================================================================ #

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
#   Sup Figure 7
# ============================================================================ #

df <- data.frame(matrix(0, nrow = length(unique(meta[["tissue"]])) *
                          length(unique(meta[["subset"]])),
                        ncol = 3))
colnames(df) <- c("tissue", "subset", "value")

val1 <- unique(meta[["tissue"]])
val2 <- unique(meta[["subset"]])

df[["tissue"]]  <- rep(val1, each  = length(val2))
df[["subset"]] <- rep(val2, times = length(val1))
df[["subset"]] <- as.factor(df[["subset"]])
df[["tissue"]]  <- as.factor(df[["tissue"]])

for (x in val1) {
  for (p in val2) {
    num_of_cells <-
      nrow(meta[meta[["tissue"]] == x & meta[["subset"]] == p,])
    df[which(df[["tissue"]] == x &
               df[["subset"]] == p), "value"] <- num_of_cells
  }
}

p1 <-
  ggplot(df, aes(fill = .data[["subset"]], y = value, x = .data[["tissue"]])) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = proti) +
  theme_bw() +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0.5
  ),
  text = element_text(size = 10))

ggsave(
  "figures/plots/sup_fig7.png",
  plot = p1,
  width = 5,
  height = 6,
  dpi = 300
)
