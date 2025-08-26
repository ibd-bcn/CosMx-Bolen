# ============================================================================ #
#   CosMx Analysis - Figure 3 Script
# ============================================================================ #

# Libraries -------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(plyr)
library(ggplot2)
library(ggrepel)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)

# Read objects ----------------------------------------------------------------
seurat <- readRDS("analysis/CosMx_RNA/Objects/seurats_all_norm.RDS")
meta   <- seurat@meta.data

# ============================================================================ #
#   Helper functions
# ============================================================================ #

# Volcano plot ----------------------------------------------------------------
volcano <- function(anot = "subset", ct, dif_col = "tissue",
                    seu, id1, id2) {

  metadata <- seu@meta.data
  subset_meta <- metadata[metadata[[anot]] == ct, ]
  object <- seu[, rownames(subset_meta)]
  object <- NormalizeData(object)
  object <- ScaleData(object)
  object@meta.data[[dif_col]] <- as.factor(object@meta.data[[dif_col]])
  object <- SetIdent(object, value = object@meta.data[[dif_col]])

  deg_results <- FindMarkers(object, ident.1 = id1, ident.2 = id2)
  deg_results <- na.omit(deg_results)
  deg_results$genes <- rownames(deg_results)

  # Classification
  deg_results$diffexpressed <- "NO"
  deg_results$diffexpressed[deg_results$avg_log2FC > log2(1.2) & deg_results$p_val < 0.05] <- "p.val<0.05 & FC>1.2"
  deg_results$diffexpressed[deg_results$avg_log2FC < -log2(1.2) & deg_results$p_val < 0.05] <- "p.val<0.05 & FC<0.83"
  deg_results$diffexpressed[deg_results$avg_log2FC > log2(1.2) & deg_results$p_val_adj < 0.05] <- "p.adj<0.05 & FC>1.2"
  deg_results$diffexpressed[deg_results$avg_log2FC < -log2(1.2) & deg_results$p_val_adj < 0.05] <- "p.adj<0.05 & FC<0.83"

  # Labels
  deg_results$delabel <- NA
  deg_results$delabel[deg_results$diffexpressed != "NO"] <- deg_results$genes[deg_results$diffexpressed != "NO"]
  deg_results$p_val <- ifelse(deg_results$p_val < 1e-300, 1e-300, deg_results$p_val)

  p <- ggplot(data = deg_results,
              aes(x = avg_log2FC, y = -log10(p_val),
                  col = diffexpressed, label = delabel)) +
    geom_point(size = 2) +
    geom_text_repel(size = 3) +
    theme_bw() +
    scale_color_manual(values = c(
      "p.val<0.05 & FC<0.83" = "green",
      "p.adj<0.05 & FC<0.83" = "darkgreen",
      "p.val<0.05 & FC>1.2"  = "red",
      "p.adj<0.05 & FC>1.2"  = "darkred"
    )) +
    geom_vline(xintercept = c(-log2(1.2), log2(1.2)), col = "black", linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dashed") +
    labs(color = "Legend") +
    theme(text = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.text  = element_text(size = 10))

  return(list(p, deg_results))
}

# Pathway analysis ------------------------------------------------------------
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
# ============================================================================ #
#   Figure 3
# ============================================================================ #

## A - Epi IBD vs NHC ---------------------------------------------------------
df <- volcano(anot = "subset", ct = "epi", seu = seurat,
              id1 = "IBD", id2 = "NHC", dif_col = "tissue")
deg <- df[[2]]
pA <- pathway_anal(ct = "epi", id1 = "IBD", id2 = "NHC", deg_results = deg, max = 10, cat = "C5",text_size = 10)
ggsave("figures/plots/fig3A.png", plot = pA, width = 10, height = 4, dpi = 300)

## B - Epi PD vs NHC ----------------------------------------------------------
df <- volcano(anot = "subset", ct = "epi", seu = seurat,
              id1 = "PD", id2 = "NHC", dif_col = "tissue")
deg <- df[[2]]
pB <- pathway_anal(ct = "epi", id1 = "PD", id2 = "NHC", deg_results = deg, max = 10, cat = "C5",text_size = 10)
ggsave("figures/plots/fig3B.png", plot = pB, width = 10, height = 4, dpi = 300)

## C - Colonocytes IBD vs NHC -------------------------------------------------
df <- volcano(anot = "refined", ct = "Colonocytes", seu = seurat,
              id1 = "IBD", id2 = "NHC", dif_col = "tissue")
deg <- df[[2]]
pC <- pathway_anal(ct = "Colonocytes", id1 = "IBD", id2 = "NHC", deg_results = deg, max = 10, cat = "C5",text_size = 10)
ggsave("figures/plots/fig3C.png", plot = pC, width = 10, height = 4, dpi = 300)

## D - Colonocytes PD vs NHC --------------------------------------------------
df <- volcano(anot = "refined", ct = "Colonocytes", seu = seurat,
              id1 = "PD", id2 = "NHC", dif_col = "tissue")
deg <- df[[2]]
pD <- pathway_anal(ct = "Colonocytes", id1 = "PD", id2 = "NHC", deg_results = deg, max = 10, cat = "C5",text_size = 10)
ggsave("figures/plots/fig3D.png", plot = pD, width = 10, height = 4, dpi = 300)
