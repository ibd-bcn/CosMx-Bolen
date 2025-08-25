density_plots <- function(gene, cells, seurat){
  ccl18 <- seurat@assays$RNA$counts[gene,]
  meta <- seurat@meta.data
  meta$ccl18 <- ccl18
  meta <- meta[meta$refined %in% cells,]

  ccl18 <- meta  %>%
    group_by(tissue,refined) %>%
    summarize(mean_ccl18 = mean(ccl18, na.rm = TRUE))

  # Create the plot with geom_line and geom_area
  p <- ggplot(ccl18, aes(x = refined, y = mean_ccl18, group = tissue, fill = tissue, color = tissue)) +
    geom_area(alpha = 0.5, position = 'identity') +  # Fill the area under the lines
    geom_line(size = 1) +  # Connect the dots
    theme_bw() +
    labs(title = gene,
         x = "Variable",
         y = "Value") +
    theme(text = element_text(size = 16),
          axis.text.x = element_text(angle = 45, hjust = 1))

  # Print the plot
  return(p)

}

volcano <- function(anot = "subset", ct, dif_col = "tissue", seu,id1,id2){
  metadata <- seu@meta.data
  # Filter for myeloid cells
  myeloid_metadata <- metadata[metadata[[anot]] == ct, ]
  # Use table to count cells in each tissue
  myeloid_counts_table <- table(myeloid_metadata[[dif_col]])
  #seu
  object <- seu[,rownames(seu@meta.data[seu@meta.data[[anot]]== ct,])]
  #object <- object[list_of_genes,]
  object <- NormalizeData(object)
  object <- ScaleData(object)
  #object@active.ident <- as.factor(object@meta.data[[dif_col]])

  # Convert the specified column to a factor
  object@meta.data[[dif_col]] <- as.factor(object@meta.data[[dif_col]])
  # Update the active identity with the factor column
  object <- SetIdent(object, value = object@meta.data[[dif_col]])

  #IBD vs NHC
  deg_results <- FindMarkers(object, ident.1 = id1, ident.2 = id2)
  deg_results <- na.omit(deg_results)
  deg_results$genes <- rownames(deg_results)
  #Diff expressed
  deg_results$diffexpressed <- "NO"
  deg_results$diffexpressed[deg_results$avg_log2FC > log2(1.2) & deg_results$p_val < 0.05] <- "UP"
  deg_results$diffexpressed[deg_results$avg_log2FC < -log2(1.2) & deg_results$p_val < 0.05] <- "DOWN"
  deg_results$diffexpressed[deg_results$avg_log2FC > log2(1.2) & deg_results$p_val_adj < 0.05] <- "UPP"
  deg_results$diffexpressed[deg_results$avg_log2FC < -log2(1.2) & deg_results$p_val_adj < 0.05] <- "DOWNN"
  #Label of most down-up regulated genes
  deg_results$delabel <- NA
  deg_results$delabel[deg_results$diffexpressed != "NO"] <- deg_results$genes[deg_results$diffexpressed != "NO"]
  deg_results$p_val <- ifelse(deg_results$p_val < 1e-300, 1e-300, deg_results$p_val)

  #Plot
  p <- ggplot(data = deg_results, aes(x = avg_log2FC, y = -log10(p_val), col = diffexpressed, label = delabel)) +
    geom_point() +
    theme_bw() +
    geom_text_repel() +
    scale_color_manual(values = c("UP" = "red", "UPP" = "darkred", "DOWN" = "green", "DOWNN" = "darkgreen")) +
    #geom_vline(xintercept=c(-(0.5),0.5 ), col="red") +
    #geom_hline(yintercept=-log10(0.05), col="red") +
    theme(text = element_text(size = 18)) +
    ggtitle(paste0(ct, " : ",id1, " vs ", id2, sep = ""))

  # Print the plot
  return(list(p,deg_results))

}


scotia <-  function(meta ,ot,fov,l_gene,r_gene,labels){
  ot <- ot[ot$fov == fov,]

  if(l_gene != "all"){
    cut_ot <- ot[ot$ligand_recptor == paste(l_gene,"_",r_gene,sep = ""),]
  }else{
    cut_ot <- ot
  }

  meta_cells <- meta[meta$fov == fov,]

  p <- ggplot(data = meta_cells, aes(x = CenterX_global_px, y =CenterY_global_px)) +
    geom_point(size = 0.5, color = "grey") +
    geom_point(data = cut_ot, aes(x = x_receptor, y = y_receptor,size = likelihood), color = "red") +
    geom_point(data = cut_ot, aes(x = x_source, y = y_source,size = likelihood), color = "purple")+
    theme_bw()

  if(labels == TRUE){
    p <- p + geom_label_repel(
      data = cut_ot,
      mapping = aes(x = x_receptor, y = y_receptor, label = ligand_recptor),
      size = 3,
      box.padding = unit(1, "lines"),
      fill = "white",        # Set the background color of the label to white
      color = "black",       # Set the text color to black
      max.overlaps = Inf
    )
  }

  return(p)
}

#Plot diff expression in FOV
plot_spatial_dis_gene <- function(object, meta , sample, fov, gene, pol_path, cell_name_col = "cell_names", mult = FALSE){
  meta[[gene]] <- as.vector(object[gene,])
  #Modify Data
  if(sample == FALSE){

    meta <- meta[meta$fov %in% fov,]



    poly <- pol_path

  }

  cells <- meta[[cell_name_col]]
  poly <- poly[poly[[cell_name_col]] %in% cells, ]
  poly[[gene]] <-as.numeric(
    mapvalues(x = poly[[cell_name_col]],
              from = meta[[cell_name_col]],
              to = meta[[gene]]))
  poly[[gene]] <- log(poly[[gene]]+0.1)

  #Plot
  p <- ggplot(poly, aes(x = x_local_px, y = y_local_px)) +
    geom_polygon(aes(group = .data[[cell_name_col]], fill = .data[[gene]]),
                 color = 'black') +
    dark_theme_gray(base_family = "Fira Sans Condensed Light", base_size = 20) +  scale_fill_gradient(low = "white", high = "darkred") +
    facet_wrap(~fov) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_blank(),
      legend.position = "right"
    ) +
    labs(x = "x",
         y = "y")

  return(p)
}

#Plot polygons
plot_pol <- function(object,
                     fov,
                     #patient,
                     poly,
                     mols,
                     pt_size,
                     annotation,
                     pal,
                     mols_c,
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

#Plot dots
plot_fov <- function(object,
                     fov,
                     pt_size,
                     x = "CenterX_global_px",
                     y ="CenterY_global_px",
                     pal,
                     annotation) {
  if (fov != "all") {
    object <- object[object$fov == fov, ]
  }
  p <-
    ggplot(data = object, aes(x = .data[[x]] , y = .data[[y]]))  + geom_point(aes(color = .data[[annotation]]),size = pt_size) +
    dark_theme_gray(base_family = "Fira Sans Condensed Light", base_size = 20) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_blank()
    ) + scale_color_manual(values = pal)
  return(p)
}

#Plot Pathway analysis
pathway_anal <- function(anot = "subset", ct, dif_col = "tissue", seu,id1,id2, cat, max = 10, text_size = 1){
  metadata <- seu@meta.data
  # Filter for myeloid cells
  myeloid_metadata <- metadata[metadata[[anot]] == ct, ]
  # Use table to count cells in each tissue
  myeloid_counts_table <- table(myeloid_metadata[[dif_col]])
  #seu
  object <- seu[,rownames(seu@meta.data[seu@meta.data[[anot]]== ct,])]
  #object <- object[list_of_genes,]
  object <- NormalizeData(object)
  object <- ScaleData(object)
  object@active.ident <- as.factor(object$tissue)
  #IBD vs NHC
  deg_results <- FindMarkers(object, ident.1 = id1, ident.2 = id2)
  deg_results <- na.omit(deg_results)
  deg_results$genes <- rownames(deg_results)
  #Diff expressed
  deg_results$diffexpressed <- "NO"
  deg_results$diffexpressed[deg_results$avg_log2FC > log2(1.2) & deg_results$p_val < 0.05] <- "UP"
  deg_results$diffexpressed[deg_results$avg_log2FC < -log2(1.2) & deg_results$p_val < 0.05] <- "DOWN"
  #Label of most down-up regulated genes+
  msigdb_hallmark <- msigdbr(species = "Homo sapiens", category = cat)
  up <- deg_results[deg_results$diffexpressed == "UP",]$genes
  gene_list_up <- bitr(up, fromType = "SYMBOL",
                       toType = "ENTREZID",
                       OrgDb = org.Hs.eg.db)
  down <- deg_results[deg_results$diffexpressed == "DOWN",]$genes
  gene_list_down <- bitr(down, fromType = "SYMBOL",
                         toType = "ENTREZID",
                         OrgDb = org.Hs.eg.db)
  gene_sets <- msigdb_hallmark %>%
    dplyr::select(gs_name, entrez_gene)
  # Perform enrichment analysis
  up_enrichment_results <- as.data.frame(enricher(gene = gene_list_up$ENTREZID,
                                                  TERM2GENE = gene_sets))
  if (nrow(up_enrichment_results) >= max) {
    # Order the dataframe by the q_value column in ascending order
    up_enrichment_results <- up_enrichment_results[order(up_enrichment_results$qvalue), ]

    # Select the first 10 rows
    up_enrichment_results <- head(up_enrichment_results, max)
  }
  up_enrichment_results <- up_enrichment_results %>%
    mutate(GeneRatio = as.numeric(sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))))


  down_enrichment_results <- as.data.frame(enricher(gene = gene_list_down$ENTREZID,
                                                    TERM2GENE = gene_sets))

  if (nrow(down_enrichment_results) >= max) {
    # Order the dataframe by the q_value column in ascending order
    down_enrichment_results <- down_enrichment_results[order(down_enrichment_results$qvalue), ]

    # Select the first 10 rows
    down_enrichment_results <- head(down_enrichment_results, max)
  }

  down_enrichment_results <- down_enrichment_results %>%
    mutate(GeneRatio = as.numeric(sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))))

  if(nrow(up_enrichment_results) != 0 & nrow(down_enrichment_results) != 0){

    up_enrichment_results$Direction <- "Upregulated"
    up_enrichment_results$AdjustedGeneRatio <- up_enrichment_results$GeneRatio
    down_enrichment_results$Direction <- "Downregulated"
    down_enrichment_results$AdjustedGeneRatio <- -down_enrichment_results$GeneRatio

    combined_results <- rbind(as.data.frame(up_enrichment_results),
                              as.data.frame(down_enrichment_results))

    p <- ggplot(combined_results, aes(x = AdjustedGeneRatio, y = reorder(Description, GeneRatio), color = Direction)) +
      geom_point(aes(size = Count)) +
      theme_minimal() +
      ggtitle(paste0("Pathway Enrichment Analysis ",id1, " vs ",id2, ": ",ct, sep = "")) +
      xlab("Gene GeneRatio") +
      ylab("Pathway") +
      scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
      theme(legend.position = "bottom",text = element_text(size = text_size)) + geom_vline(xintercept = 0, linewidth= 1)
  }

  if(nrow(up_enrichment_results) == 0 & nrow(down_enrichment_results) == 0){
    message("no enriched pathway analysis, sorry sors")
  }

  if(nrow(up_enrichment_results) == 0 & nrow(down_enrichment_results) != 0){
    down_enrichment_results$Direction <- "Downregulated"
    down_enrichment_results$AdjustedGeneRatio <- -down_enrichment_results$GeneRatio

    p <- ggplot(down_enrichment_results, aes(x = AdjustedGeneRatio, y = reorder(Description, AdjustedGeneRatio), color = Direction)) +
      geom_point(aes(size = Count)) +
      theme_minimal() +
      ggtitle(paste0("Pathway Enrichment Analysis ",id1, " vs ",id2, ": ",ct, sep = "")) +
      xlab("Gene GeneRatio") +
      ylab("Pathway") +
      scale_color_manual(values = c("Downregulated" = "blue")) +
      theme(legend.position = "bottom",text = element_text(size = text_size)) + geom_vline(xintercept = 0, linewidth= 1)


  }

  if(nrow(up_enrichment_results) != 0 & nrow(down_enrichment_results) == 0){

    up_enrichment_results$Direction <- "Upregulated"

    p <- ggplot(up_enrichment_results, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = Direction)) +
      geom_point(aes(size = Count)) +
      theme_minimal() +
      ggtitle(paste0("Pathway Enrichment Analysis ",id1, " vs ",id2, ": ",ct, sep = "")) +
      xlab("Gene GeneRatio") +
      ylab("Pathway") +
      scale_color_manual(values = c("Upregulated" = "red")) +
      theme(legend.position = "bottom",text = element_text(size = text_size)) + geom_vline(xintercept = 0, linewidth= 1)
  }

  return(p)
}


chordiagram <- function(chord, meta_seu, fov, subset_source, subset_target,refined_source, refined_receptor, freq, down = FALSE){

  if("all" %in%  subset_source){
    subset_source <- c("plasmas"  ,"epi"     , "stroma" ,  "myeloids" ,"tcells"  )
  }
  if("all" %in%  subset_target){
    subset_target <- c("plasmas"  ,"epi"     , "stroma" ,  "myeloids" ,"tcells"  )
  }
  if("all" %in%  refined_source){
    refined_source <- unique(meta_seu$refined)

  }
  if("all" %in%  refined_receptor){
    refined_receptor <- unique(meta_seu$refined)
  }
  if("all" %in% fov){
    fov <- c(1:66)

  }
  #Select fov
  pat_cells <- meta_seu[meta_seu$fov %in% fov,]$cell

  chord <- chord[chord$id_source %in% pat_cells & chord$source_cell_type %in% subset_source & chord$target_cell_type %in% subset_target & chord$refined_source %in% refined_source & chord$refined_receptor %in% refined_receptor,]

  chord <- as.data.frame(table(chord$refined_source, chord$refined_receptor))
  chord <- chord[chord$Freq > freq,]
  chord_1 <- chord %>% spread(key = Var2, value = Freq)
  rownames(chord_1) <- chord_1$Var1
  chord_1 <- chord_1[,2:ncol(chord_1)]

  if(down == TRUE){
    return(chord)
  }else{
    return(chord_1)
  }


}

## Chordiagram prot
chordiagram_prot <- function(chord, meta_seu, fov, subset_source, subset_target, freq, down = FALSE){

  if("all" %in%  subset_source){
    subset_source <- c("NA"  ,"Plasma"  ,"Epithelium"     , "Fibroblasts" ,  "Bcell" ,"Tcells"  )
  }
  if("all" %in%  subset_target){
    subset_target <- c("NA"  ,"Plasma"  ,"Epithelium"     , "Fibroblasts" ,  "Bcell" ,"Tcells"  )
  }
  if("all" %in% fov){
    fov <- c(1:66)

  }
  #Select fov
  pat_cells <- meta_seu[meta_seu$fov %in% fov,]$cell

  chord <- chord[chord$id_source %in% pat_cells & chord$source_cell_type %in% subset_source & chord$target_cell_type %in% subset_target ,]

  chord <- as.data.frame(table(chord$source_cell_type, chord$target_cell_type))
  chord <- chord[chord$Freq > freq,]
  chord_1 <- chord %>% spread(key = Var2, value = Freq)
  rownames(chord_1) <- chord_1$Var1
  chord_1 <- chord_1[,2:ncol(chord_1)]

  if(down == TRUE){
    return(chord)
  }else{
    return(chord_1)
  }


}
##CHORD PLOTS
chord_plot <- function(meta_seu ,mat, cell_type,color_chord){
  cells <- unique(meta_seu$refined)
  names(color_chord) <- cells

  col_mat <- matrix("#F0F0F0", nrow = nrow(mat), ncol = ncol(mat))
  colnames(col_mat) <-colnames(mat)
  rownames(col_mat) <- rownames(mat)

  for(z in rownames(mat)){
    col_mat[z,] <- color_chord[z]

  }

  if (!("all" %in% cell_type)) {

    cols <- setdiff(colnames(col_mat), c(cell_type))

    # Get all row names except "hola23" and "hola333"
    rows <- setdiff(rownames(col_mat), c(cell_type))

    col_mat[rows, cols] <- "#F0F0F0"


  }

  chordDiagram(mat,
               annotationTrack = "grid", preAllocateTracks = 1, grid.col = color_chord, col = col_mat)


  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
    circos.axis(h = "top", labels.cex = 0.5, major.tick.length = 0.2, sector.index = sector.name, track.index = 2)
  }, bg.border = NA)




}

##CHORD PLOTS prot
chord_plot_prot <- function(meta_seu ,mat, cell_type,color_chord){
  cells <- unique(meta_seu$subset)
  names(color_chord) <- cells

  col_mat <- matrix("#F0F0F0", nrow = nrow(mat), ncol = ncol(mat))
  colnames(col_mat) <-colnames(mat)
  rownames(col_mat) <- rownames(mat)

  for(z in rownames(mat)){
    col_mat[z,] <- color_chord[z]

  }

  if (!("all" %in% cell_type)) {

    cols <- setdiff(colnames(col_mat), c(cell_type))

    # Get all row names except "hola23" and "hola333"
    rows <- setdiff(rownames(col_mat), c(cell_type))

    col_mat[rows, cols] <- "#F0F0F0"


  }

  chordDiagram(mat,
               annotationTrack = "grid", preAllocateTracks = 1, grid.col = color_chord, col = col_mat)


  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
    circos.axis(h = "top", labels.cex = 0.5, major.tick.length = 0.2, sector.index = sector.name, track.index = 2)
  }, bg.border = NA)




}
