library(Seurat)
library(dplyr)
#Read single cell object of REFERENCE Trigo. et al -----------------------------
todas <-
  readRDS(
    "/home/acorraliza/data_Albas/All_Together_31082021/piezas_finales/tog_numeric.RDS"
  )

cells <-
  rownames(todas@meta.data[todas@meta.data$annotation_general_ok %in% c(
    "CD4",
    "CD8",
    "PC IgA",
    "B cell",
    "Undifferentiated epithelium",
    "Macrophages",
    "Colonocyte",
    "OLFM4 epithelium",
    "Fibroblast",
    "Myofibroblasts",
    "Inflammatory macrophages",
    "PC IgG"
  ),])

todas_cut <- todas[, cells]

#Create annotation for Maxfuse -------------------------------------------------
todas_cut@meta.data$anot_maxfuse <-
  recode(
    todas_cut@meta.data$annotation_general_ok,
    "CD4" = "Tcells",
    "CD8" = "Tcells",
    "PC IgA" = "Plasma",
    "B cell" = "Bcell",
    "Undifferentiated epithelium" = "Epithelium",
    "Macrophages" = "Macrophages",
    "Colonocyte" = "Epithelium",
    "OLFM4 epithelium" = "Epithelium",
    "Fibroblast" = "Fibroblasts",
    "Myofibroblasts" = "Fibroblasts",
    "Inflammatory macrophages" = "Macrophages",
    "PC IgG" = "Plasma"
  )

todas_cut_sub@meta.data$nCount_RNA <-
  Matrix::colSums(todas_cut_sub@assays$RNA@counts)
todas_cut_sub <- todas_cut_sub[, todas_cut_sub$nCount_RNA > 5]

cell_names <- todas_cut_sub$id
cell_type <- todas_cut_sub$anot_maxfuse

# Create a data frame with cell names and cell types
cell_data <-
  data.frame(cell_names = cell_names, cell_type = cell_type)

# Sample 2,000 cells per cell type (if available)
sampled_cells <- cell_data %>%
  group_by(cell_type) %>%
  slice_sample(n = 2000, replace = FALSE) %>%
  pull(cell_names)

todas_cut <- todas_cut[, sampled_cells]
DimPlot(todas_cut, group.by = "anot_maxfuse", label = TRUE)

#Save CSV files Maxfuse --------------------------------------------------------
write.csv(
  todas_cut@assays$RNA$counts,
  "~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Objects/counts_SC_cut.csv"
)
write.csv(
  todas_cut@meta.data,
  "~/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Objects/meta_SC_cut.csv"
)

#Create correspondence matrix

proteins <-
  c(
    "4-1BB",
    "B7-H3",
    "Bcl-2",
    "Beta-catenin",
    "CCR7",
    "CD11b",
    "CD11c",
    "CD123",
    "CD127",
    "CD138",
    "CD14",
    "CD15",
    "CD16",
    "CD163",
    "CD19",
    "CD20",
    "CD27",
    "CD3",
    "CD31",
    "CD34",
    "CD38",
    "CD39",
    "CD4",
    "CD40",
    "CD45",
    "CD45RA",
    "CD56",
    "CD68",
    "CD8",
    "CTLA4",
    "Channel-CD45",
    "Channel-CD68",
    "Channel-DNA",
    "Channel-Membrane",
    "Channel-PanCK",
    "EGFR",
    "EpCAM",
    "FABP4",
    "FOXP3",
    "Fibronectin",
    "GITR",
    "GZMA",
    "GZMB",
    "HLA-DR",
    "Her2",
    "ICAM1",
    "ICOS",
    "IDO1",
    "IL-18",
    "IL-1b",
    "IgD",
    "Ki-67",
    "LAG3",
    "LAMP1",
    "NF-kB p65",
    "PD-1",
    "PD-L1",
    "PD-L2",
    "SMA",
    "STING",
    "TCF7",
    "Tim-3",
    "VISTA",
    "Vimentin",
    "iNOS",
    "p53",
    "pan-RAS"
  )

# Corresponding gene names (manual conversion for known proteins)
genes <-
  c(
    "TNFRSF9",
    "CD276",
    "BCL2",
    "CTNNB1",
    "CCR7",
    "ITGAM",
    "ITGAX",
    "IL3RA",
    "IL7R",
    "SDC1",
    "CD14",
    "FUT4",
    "FCGR3A",
    "CD163",
    "CD19",
    "MS4A1",
    "CD27",
    "CD3E",
    "PECAM1",
    "CD34",
    "CD38",
    "ENTPD1",
    "CD4",
    "CD40",
    "PTPRC",
    "PTPRC",
    "NCAM1",
    "CD68",
    "CD8A",
    "CTLA4",
    NA,
    NA,
    NA,
    NA,
    NA,
    "EGFR",
    "EPCAM",
    "FABP4",
    "FOXP3",
    "FN1",
    "TNFRSF18",
    "GZMA",
    "GZMB",
    "HLA-DRA",
    "ERBB2",
    "ICAM1",
    "ICOS",
    "IDO1",
    "IL18",
    "IL1B",
    "IGHD",
    "MKI67",
    "LAG3",
    "LAMP1",
    "RELA",
    "PDCD1",
    "CD274",
    "PDCD1LG2",
    "ACTA2",
    "TMEM173",
    "TCF7",
    "HAVCR2",
    "VISTA",
    "VIM",
    "NOS2",
    "TP53",
    "HRAS"
  )

# Create a data frame
protein_gene_df <-
  data.frame("Protein_name" = proteins, "RNA_name" = genes)
protein_gene_df <- na.omit(protein_gene_df)
write.csv(
  protein_gene_df,
  "~/SPATIAL/Mackensy_analysis/CosMx_Prot/protein_gene_conversion.csv",
  row.names = FALSE
)

