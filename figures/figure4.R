library(Seurat)
library(ggplot2)
library(ggrepel)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(plyr)
library(ggdark)
library(dplyr)
options(bitmapType = "cairo")

## Data ------------------------------------------------------------------------
seu <-
  readRDS(
    "/home/mmoro/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Objects/seurats_all_norm.RDS"
  )
poly <-
  read.csv(
    "/home/mmoro/SPATIAL/Mackensy_analysis/CosMx_RNA/0_Curation/Polygons/bolen_colon_RNA.csv"
  )
mols <-
  read.csv(
    "/home/mmoro/SPATIAL/Mackensy_analysis/CosMx_RNA/0_Curation/Molecules/bolen_colon_RNA.csv"
  )
meta <-
  read.csv(
    "/home/mmoro/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Objects/meta.csv"
  )

## Palette of colours ----------------------------------------------------------
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
  "FRCs" = "#8A2BE2FF",

  #Other
  "other" =  "#515151"
)

subset_col <- c(
  "epi" = "#BA6338FF",
  "stroma" = "#5DB1DDFF",
  "tcells" = "#802268FF",
  "plasmas" = "#6BD76BFF",
  "myeloids" = "#D595A7FF",

  #Other
  "other" =  "#4b4b4b"
)

##Functions --------------------------------------------------------------------

#Volcanos
volcano <-
  function(anot = "subset",
           ct,
           dif_col = "tissue",
           seu,
           id1,
           id2) {
    metadata <- seu@meta.data
    myeloid_metadata <- metadata[metadata[[anot]] == ct,]
    myeloid_counts_table <- table(myeloid_metadata[[dif_col]])
    object <-
      seu[, rownames(seu@meta.data[seu@meta.data[[anot]] == ct, ])]
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
    p <-
      ggplot(data = deg_results,
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
          "p.val<0.05 & FC<0.83" = "#5cbde7",
          "p.adj<0.05 & FC<0.83" = "darkblue",
          "p.val<0.05 & FC>1.2" = "red",
          "p.adj<0.05 & FC>1.2" = "darkred"
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

#Pathway analysis
pathway_anal <-
  function(ct,
           id1,
           id2,
           deg_results,
           cat,
           max = 10,
           text_size = 1) {
    msigdb_hallmark <- msigdbr(species = "Homo sapiens", category = cat)
    up <-
      deg_results[deg_results$diffexpressed == "p.val<0.05 & FC>1.2" | deg_results$diffexpressed == "p.adj<0.05 & FC>1.2", ]$genes
    gene_list_up <- bitr(up,
                         fromType = "SYMBOL",
                         toType = "ENTREZID",
                         OrgDb = org.Hs.eg.db)
    down <-
      deg_results[deg_results$diffexpressed == "p.val<0.05 & FC<0.83" | deg_results$diffexpressed == "p.adj<0.05 & FC>1.2", ]$genes
    gene_list_down <- bitr(down,
                           fromType = "SYMBOL",
                           toType = "ENTREZID",
                           OrgDb = org.Hs.eg.db)


    gene_sets <- msigdb_hallmark %>%
      dplyr::select(gs_name, entrez_gene)

    up_enrichment_results <-
      as.data.frame(enricher(gene = gene_list_up$ENTREZID,
                             TERM2GENE = gene_sets))

    up_enrichment_results <- up_enrichment_results %>%
      mutate(GeneRatio = as.numeric(sapply(strsplit(GeneRatio, "/"), function(x)
        as.numeric(x[1]) / as.numeric(x[2]))))

    down_enrichment_results <-
      as.data.frame(enricher(gene = gene_list_down$ENTREZID,
                             TERM2GENE = gene_sets))

    down_enrichment_results <- down_enrichment_results %>%
      mutate(GeneRatio = as.numeric(sapply(strsplit(GeneRatio, "/"), function(x)
        as.numeric(x[1]) / as.numeric(x[2]))))

    up_enrichment_results <-
      up_enrichment_results[order(up_enrichment_results$qvalue),]

    up_enrichment_results <- head(up_enrichment_results, max)

    down_enrichment_results <-
      down_enrichment_results[order(down_enrichment_results$qvalue),]

    down_enrichment_results <- head(down_enrichment_results, max)

    up_enrichment_results$s1 <- "Upregulated"
    down_enrichment_results$s1 <- "Downregulated"

    final_enrichment_result <-
      rbind(up_enrichment_results, down_enrichment_results)

    final_enrichment_result$log10pval <-
      -log10(final_enrichment_result$pvalue)

    final_enrichment_result <- final_enrichment_result %>%
      arrange(GeneRatio)

    format_pathway <- function(pathway) {
      pathway <- sub("_", ": ", pathway)
      pathway <- gsub("_", " ", pathway)
      return(pathway)
    }

    final_enrichment_result$Description <- sapply(final_enrichment_result$Description, format_pathway)

    p <-
      ggplot(final_enrichment_result,
             aes(
               x = s1 ,
               y = Description,
               size = GeneRatio,
               color = s1
             )) +
      geom_point() +
      scale_color_manual(values = c(
        "Upregulated" = "darkred",
        "Downregulated" = "darkblue"
      )) +
      scale_size_continuous(
        name = "Gene ratio",
        range = c(3, 10),
        guide = guide_legend(override.aes = list(
          color = "black", fill = "black"
        ))
      ) +
      theme_bw() +
      labs(
        title = "",
        x = "",
        y = "Pathway description",
        size = "Gene ratio",
        color = "Pathway"
      ) +
      theme(
        text = element_text(family = "Helvetica"),
        plot.title = element_text(face = "bold", size = text_size),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = text_size),
        legend.text = element_text(size = text_size),
        legend.title = element_text(face = "bold", size = text_size)
      ) + guides(color = guide_legend(override.aes = list(size = 5)))
    return(p)
  }

## Figure 3A  ------------------------------------------------------------------

df <-
  volcano(
    anot = "subset",
    ct = "epi" ,
    seu = seu ,
    id1 = "IBD",
    id2 = "NHC" ,
    dif_col = "tissue"
  )
df <- df[[2]]
p <-
  pathway_anal(
    ct = "epi",
    id1 = "IBD" ,
    id2 = "NHC",
    deg_results = df,
    max = 10,
    cat = "C5",
    text_size =  20
  )
png(
  filename = "~/SPATIAL/Mackensy_analysis/Figures/plots/figure3B.png",
  width = 16,
  height = 12,
  units = "in",
  res = 800
)
p
dev.off()


## Figure 3B  ------------------------------------------------------------------
df <-
  volcano(
    anot = "subset",
    ct = "epi" ,
    seu = seu ,
    id1 = "PD",
    id2 = "NHC" ,
    dif_col = "tissue"
  )
df <- df[[2]]
p <-
  pathway_anal(
    ct = "epi",
    id1 = "PD" ,
    id2 = "NHC",
    deg_results = df,
    max = 10,
    cat = "C5",
    text_size =  20
  )
png(
  filename = "~/SPATIAL/Mackensy_analysis/CosMx-Bolen/figures/plots/figure3D.png",
  width = 22,
  height = 12,
  units = "in",
  res = 800
)
p
dev.off()

## Figure 3C  ------------------------------------------------------------------
df <-
  volcano(
    anot = "refined",
    ct = "Colonocytes" ,
    seu = seu ,
    id1 = "IBD",
    id2 = "NHC" ,
    dif_col = "tissue"
  )
df <- df[[2]]
p <-
  pathway_anal(
    ct = "Colonocytes",
    id1 = "IBD" ,
    id2 = "NHC",
    deg_results = df,
    max = 10,
    cat = "C5",
    text_size =  20
  )
png(
  filename = "~/SPATIAL/Mackensy_analysis/Figures/plots/figure3F.png",
  width = 16,
  height = 12,
  units = "in",
  res = 800
)
p
dev.off()

## Figure 3D  ------------------------------------------------------------------

df <-
  volcano(
    anot = "refined",
    ct = "Colonocytes" ,
    seu = seu ,
    id1 = "PD",
    id2 = "NHC" ,
    dif_col = "tissue"
  )
df <- df[[2]]
p <-
  pathway_anal(
    ct = "Colonocytes",
    id1 = "PD" ,
    id2 = "NHC",
    deg_results = df,
    max = 10,
    cat = "C5",
    text_size =  20
  )
png(
  filename = "~/SPATIAL/Mackensy_analysis/Figures/plots/figure3H.png",
  width = 18,
  height = 12,
  units = "in",
  res = 800
)
p
dev.off()
