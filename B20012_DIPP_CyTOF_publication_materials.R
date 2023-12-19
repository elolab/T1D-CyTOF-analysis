
plot_updated_cell_clusters_pdf <- function(filename,
                                       set,
                                       color_clusters,
                                       labels) {
  
  # TSNE with coloring by clustering, with updated cluster names
  dev.flush()
  pdf(file=filename)
  cluster_ids <- sort(unique(set$updated_clustering.cell))
  plot(set$dr$tSNE1, set$dr$tSNE2, bty="n", pch=19, cex=0.1, col=color_clusters[set$dr$updated_clustering.cell], axes=FALSE, xlab="", ylab="")
  if (labels) {
    text(set$updated_agg$tSNE1, set$updated_agg$tSNE2, paste(cluster_ids, set$clusternames[cluster_ids], sep=". "))
    legend("topright", legend=paste(cluster_ids, set$clusternames[cluster_ids], sep=". "), pch=15, col=color_clusters, ncol=2, bty="n")
  }
  while (!is.null(dev.list())) dev.off()
}


plot_subcell_clusters <- function(filename,
                                    set,
                                    color_clusters) {
  
  dev.flush()
  pdf(file=filename)
  plot(set$dr$tSNE1, set$dr$tSNE2, bty="n", pch=19, cex=0.1, col=color_clusters[set$clustering.cell[set$tsne_inds]], axes=FALSE, xlab="", ylab="")
  while (!is.null(dev.list())) dev.off()
}


marker_tSNE_plot <- function(filename,
                             cell_clustername, # for example B_Cells
                             set,
                             color_clusters,
                             batch,
                             panel,
                             dr,
                             tsne_expr,
                             marker) {
  
  dev.flush()

  png(file=filename, width = 2160*3, height = 2160, units = "px", pointsize = 12)
  
  layout(mat = matrix(c(1,1,1,1,2,6,3,7,4,8,5,9), nrow = 2, ncol = 6),
         heights = c(1, 1), # Heights of the rows
         widths = c(1, 1,1,1,1)) # Widths of the columns
  
  par(mar = c(15, 15, 15, 15), mgp = c(0, 1000, -3)) #, mgp = c(12, 7, 5), c(3, 1, 0)

  cell_cluster_index <- match(cell_clustername, full_set$clusternames)
  color_scheme <- rep("#20202030", max(set$dr$updated_clustering.cell))
  color_scheme[cell_cluster_index] <- color_clusters[cell_cluster_index] #"#FF2020"

  # WORKING HERE
  samples <- rownames(set$anno[!is.na(set$anno["Timepoint"]), ])
  cluster_ids <- sort(unique(set$updated_clustering.cell))
  dr_subsample <- dr$sampleids.cell %in% samples
  plot(set$dr$tSNE1[dr_subsample], set$dr$tSNE2[dr_subsample], bty="n", pch=19, cex=0.1, col=color_scheme[set$dr$updated_clustering.cell[dr_subsample]], xlab="", ylab="", main=cell_clustername, cex.main = 8, axes=FALSE)
  mtext("t-SNE component 1", side=1, line=5, cex=8)
  mtext("t-SNE component 2", side=2, line=5, cex=8)
  axis(1, lwd=8, lwd.tick=0, lab=F)
  axis(2, lwd=8, lwd.tick=0, lab=F)

  cname <- rownames(panel[panel$cname == marker,])
  all_colors <- colorRampPalette(c("darkgrey", "blue"))(64)[as.numeric(cut(tsne_expr[, cname], breaks=64))]

  titles <- c("Case", "Control")
  
  for (i in 1:length(titles)) {
    casectrl <- titles[[i]]
    
    for (timepoint in c(1,2,3,4)) {
      
      print (paste0("Rendering ", casectrl, ", timepoint ", timepoint))
    
      samples <- rownames(set$anno[(set$anno["CaseCtrl"] == casectrl) & (set$anno["Batch"] == batch) & (set$anno["Timepoint"] == timepoint), ])
      colors <- all_colors[dr$sampleids.cell %in% samples]
      plot(dr$tSNE1[dr$sampleids.cell %in% samples], dr$tSNE2[dr$sampleids.cell %in% samples], pch=19, bty="n", cex=0.4, col=colors, xlab="", ylab="", axes=FALSE, main=paste0(marker, " - ", titles[i], " (timepoint ", timepoint,")"), cex.main = 8)
      adjustcolor(colors,alpha.f=0.2)
      mtext("t-SNE component 1", side=1, line=5, cex=8)
      mtext("t-SNE component 2", side=2, line=5, cex=8)
      axis(1, lwd=8, lwd.tick=0, lab=F)
      axis(2, lwd=8, lwd.tick=0, lab=F)
      
    }        
  }

  while (!is.null(dev.list())) dev.off()

}




marker_tSNE_plot_sub <- function(filename,
                             cell_type,
                             subcluster_index, 
                             set,
                             color_clusters,
                             batch,
                             panel,
                             dr,
                             tsne_expr,
                             marker) {
  
  dev.flush()
  
  png(file=filename, width = 2160*3, height = 2160, units = "px", pointsize = 12)
  
  layout(mat = matrix(c(1,1,1,1,2,6,3,7,4,8,5,9), nrow = 2, ncol = 6),
         heights = c(1, 1), # Heights of the rows
         widths = c(1, 1,1,1,1)) # Widths of the columns
  
  par(mar = c(15, 15, 15, 15), mgp = c(0, 1000, -3)) #, mgp = c(12, 7, 5), c(3, 1, 0)
  
  color_scheme <- rep("#20202030", max(set$clustering.cell))
  color_scheme[subcluster_index] <- color_clusters[subcluster_index] #"#FF2020"
  
  samples <- rownames(set$anno[!is.na(set$anno["Timepoint"]), ])
  cluster_ids <- sort(unique(set$clustering.cell))
  dr_subsample <- dr$sampleids.cell %in% samples
  
  plot(set$dr$tSNE1[dr_subsample], set$dr$tSNE2[dr_subsample], bty="n", pch=19, cex=0.1, col=color_scheme[set$clustering.cell[set$tsne_inds[dr_subsample]]], xlab="", ylab="", axes=FALSE, main=paste(cell_type," Subcluster ", subcluster_index), cex.main = 8,)
  mtext("t-SNE component 1", side=1, line=5, cex=8)
  mtext("t-SNE component 2", side=2, line=5, cex=8)
  axis(1, lwd=8, lwd.tick=0, lab=F)
  axis(2, lwd=8, lwd.tick=0, lab=F)

  cname <- rownames(panel[panel$cname == marker,])
  all_colors <- colorRampPalette(c("darkgrey", "blue"))(64)[as.numeric(cut(tsne_expr[, cname], breaks=64))]
  
  titles <- c("Case", "Control")
  #   cname <- rownames(panel[panel$cname == marker,])
  
  for (i in 1:length(titles)) {
    casectrl <- titles[[i]]
    
    for (timepoint in c(1,2,3,4)) {
      
      samples <- rownames(set$anno[(set$anno["CaseCtrl"] == casectrl) & (set$anno["Batch"] == batch) & (set$anno["Timepoint"] == timepoint), ])
      colors <- all_colors[dr$sampleids.cell %in% samples]
      #colors <- colorRampPalette(c("lightgrey", "blue"))(64)[as.numeric(cut(tsne_expr[dr$sampleids.cell %in% samples, cname], breaks=64))]
      plot(dr$tSNE1[dr$sampleids.cell %in% samples], dr$tSNE2[dr$sampleids.cell %in% samples], pch=19, bty="n", cex=0.4, col=colors, xlab="", ylab="", axes=FALSE, main=paste0(marker, " - ", titles[i], " (timepoint ", timepoint,")"), cex.main = 8)
      adjustcolor(colors,alpha.f=0.2)
      mtext("t-SNE component 1", side=1, line=5, cex=8)
      mtext("t-SNE component 2", side=2, line=5, cex=8)
      axis(1, lwd=8, lwd.tick=0, lab=F)
      axis(2, lwd=8, lwd.tick=0, lab=F)
      
    }        
  }

  while (!is.null(dev.list())) dev.off()
  
}


OLS_plot_case_vs_control_per_type <- function(filename,
                                     table_,
                                     type_,
                                     y_axis_label,
                                     x_axis_label,
                                     ...) {
  
  if (FALSE) { 
    filename <- paste("out/OLS-Case-vs-Control_CellTypeProportions", ".pdf", sep="")
    table_ <- proportion_table
    types_ <- names(marker_intensity_tables)
    y_axis_label <- "Proportion"
    y_axis_label <- "Age (years)"  
  }
  
  # type: cell_cluster_names, markers
  
  # DEBUG
  #table_ <- table_[1:8,]
  #types_ <- c(types_[1])
  
  #png(file=filename, width = 4320, height = 2160, units = "px", pointsize = 12)
  #par(mfrow = determine_image_tiling(length(types_)), mar = c(16, 16, 16, 16), mgp = c(12, 7, 5))
  
  by_variable = "Age"
  draw_lines = FALSE
  
  pdf(file=filename, width=4, height=4)
  #for (type_ in types_) {
  
  group <- "Case"
  x1 <- table_[table_$CaseCtrl==group,by_variable]
  y1 <- table_[table_$CaseCtrl==group, type_]
  
  group <- "Control"
  x2 <- table_[table_$CaseCtrl==group,by_variable]
  y2 <- table_[table_$CaseCtrl==group, type_]
  #points(x2,y2,col="blue",pch=19, bty="n",cex=2)
  
  # main="Case (red) vs Control (blue)"   
  plot(table_[,by_variable], table_[, type_], col=ifelse(table_[, "CaseCtrl"]=="Case", "red", "blue"), pch=19, bty="l", xlab=x_axis_label, ylab=y_axis_label, cex.axis = 1, cex.lab = 1, ...)
  #plot(c(x1,x2), c(y1,y2), col=c(rep("red", length(y1)), rep("blue", length(y2)) ), pch=19, bty="n", xlab="Age", ylab=y_axis_label, cex=4, cex.lab=4, cex.axis=4)
  
  mtext(type_, cex=1, font = 2)
  
  model1 = lm(y1~x1)
  
  newx <- seq(min(x1), max(x1), length.out=length(x1))
  preds <- predict(model1, newdata = data.frame(x1=newx), interval = 'confidence')
  polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = "#10106010", border = NA)
  abline(model1, col="red", lwd=3)
  lines(newx, preds[ ,3], lty = 'dashed', col = 'red')
  lines(newx, preds[ ,2], lty = 'dashed', col = 'red')
  
  model2 = lm(y2~x2)
  
  newx <- seq(min(x2), max(x2), length.out=length(x2))
  preds <- predict(model2, newdata = data.frame(x2=newx), interval = 'confidence')
  polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = "#60101010", border = NA)
  abline(model2, col="blue", lwd=3)
  lines(newx, preds[ ,3], lty = 'dashed', col = 'blue')
  lines(newx, preds[ ,2], lty = 'dashed', col = 'blue')
  
  if (draw_lines) {
    
    linecolors <-  list("Case" = "red", "Control" = "blue")  
    
    ids = unique(paste(table_$Pair, table_$CaseCtrl))
    for(i in 1:length(ids)) {
      id <- ids[i]
      sel <- which(paste(table_$Pair, table_$CaseCtrl)==id)
      chunk <- table_[sel, ]
      chunk_order <- order(chunk[, by_variable])
      ordered_chunk <- chunk[chunk_order,]
      lines(ordered_chunk[, by_variable], ordered_chunk[, type_], col=linecolors[[ordered_chunk[1, "CaseCtrl"]]])
    }
  }
  
  
  #}
  
  
  dev.off()
}



B200012_DIPP_CyTOF_workflow_figures <- function(CyTOF_workflow_cwd,
                                            sample_metadata_tsv,
                                            markers_tsv) {
  

  color_clusters <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72", 
                      "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3", 
                      "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D", 
                      "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999", 
                      "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000", 
                      "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00",
                      "#AC012C", "#AB2052", "#698670", "#D8AFAA", "#28AE22", 
                      "#519BD6", "#BF7F80", "#6DA412", "#A7892A", "#578A23", 
                      "#63A08C", "#52AF8A", "#D5A141", "#2DD3E7", "#86768D")
  
  
  anno <- read.csv(sample_metadata_tsv, stringsAsFactors=FALSE, sep="\t")
  rownames(anno) <- anno[["Sample"]]
  
  panel <- read.table(markers_tsv, sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE)

  markers_panel <- panel[(panel["type"] == "functional") | (panel["type"] == "lineage/functional") | (panel["type"] == "lineage"), ]
  markers_cname <- markers_panel[["cname"]]
  markers <- rownames(markers_panel)


  lineage_marker_panel <- panel[(panel["type"] == "lineage/functional") | (panel["type"] == "lineage"), ]
  #lineage_markers_cname <- lineage_marker_panel[["cname"]]
  #lineage_markers <- rownames(lineage_marker_panel)
  
    
  functional_markers_panel <- panel[(panel["type"] == "functional") | (panel["type"] == "lineage/functional"), ]
  functional_markers_cname <- functional_markers_panel[["cname"]]
  #functional_markers <- rownames(functional_markers_panel)
  
  
  # PROCESS SET AS A WHOLE
  
  proportion_table_rdata_filename = paste0("RData/proportion_table", ".RData")
  load(proportion_table_rdata_filename)
  combined_proportion_table <- with(
    list(
      left = proportion_table,
      right = anno),
    {
      shared_names <- intersect(rownames(left),rownames(right))
      cbind(left[shared_names,], right[shared_names,])
    })
  
  
  clusternames <- colnames(proportion_table)[!(colnames(proportion_table) %in% c(colnames(anno), "removed"))]
  
  intensity_table_rdata_filename <- paste0("RData/marker_intensity_tables", ".RData")
  load(intensity_table_rdata_filename)

  dir.create("out/")
  dir.create("out/tables/")
  
#  for (j in seq_along(marker_intensity_tables)) {
#    celltype <- names(marker_intensity_tables)[[j]]
    
    celltype <- "NK_cells"
    marker <- "CD161"
    celltype_index <- match(celltype, names(marker_intensity_tables))
    
    
    OLS_plot_case_vs_control(paste("out/OLS-Case-vs-Control_", celltype, ".pdf", sep=""),
                             marker_intensity_tables[[celltype_index]], # cluster names
                             marker,
                             "Marker intensity",
                             "Age (years)") 
    
    
#  }
  
  # CELL TYPE PROPORTIONS
  
  OLS_plot_case_vs_control(paste("out/OLS-Case-vs-Control_CellTypeProportions", ".pdf", sep=""),
                           combined_proportion_table,
                           names(marker_intensity_tables), # TODO: get cluster names from proportion table
                           "Proportion",
                           "Age (years)") 
  
  #rm(proportion_table)
  
  # PROCESS BATCHES ()
  for (batch in unique(anno[["Batch"]])) {

    intensity_table_rdata_filename <- paste0("RData/marker_intensity_tables_", batch, ".RData")
    load(intensity_table_rdata_filename)
    
    proportion_table_rdata_filename = paste0("RData/proportion_table-", batch, ".RData")
    load(proportion_table_rdata_filename)

    combined_proportion_table <- with(
      list(
        left = proportion_table,
        right = anno),
      {
        shared_names <- intersect(rownames(left),rownames(right))
        cbind(left[shared_names,], right[shared_names,])
      })

    #for (j in seq_along(marker_intensity_tables)) {
    #celltype <- names(marker_intensity_tables)[[j]]

    celltype <- "NK_cells"
    marker <- "CD161"
    celltype_index <- match(celltype, names(marker_intensity_tables))

    OLS_plot_case_vs_control_per_type(paste("out/OLS-Case-vs-Control_", batch, "_", celltype, ".pdf", sep=""),
                                            marker_intensity_tables[[celltype_index]],
                                             marker,
                                             "Marker intensity",
                                             "Age (years)",
                                              ylim=c(1,3)) 
      
    #}

    # PROCESS CELL TYPE PROPORTIONS

    celltype <- "CD4_T_cells" 

    OLS_plot_case_vs_control_per_type(paste("out/OLS-Case-vs-Control_CellTypeProportions-", batch, "-", celltype, ".pdf", sep=""),
                             combined_proportion_table,
                             celltype, # TODO: get cluster names from proportion table
                             "Proportion",
                             "Age (years)",
                             ylim=c(0,0.7))

    
    celltype <- "NK_cells" 
    
    OLS_plot_case_vs_control_per_type(paste("out/OLS-Case-vs-Control_CellTypeProportions-", batch, "-", celltype, ".pdf", sep=""),
                                      combined_proportion_table,
                                      celltype, # TODO: get cluster names from proportion table
                                      "Proportion",
                                      "Age (years)",
                                      ylim=c(0,0.2))
    

    celltype <- "CD8_T_cells" 
    
    OLS_plot_case_vs_control_per_type(paste("out/OLS-Case-vs-Control_CellTypeProportions-", batch, "-", celltype, ".pdf", sep=""),
                                      combined_proportion_table,
                                      celltype, # TODO: get cluster names from proportion table
                                      "Proportion",
                                      "Age (years)",
                                      ylim=c(0,0.4))
    
    rm(marker_intensity_tables)
    rm(proportion_table)
    
  } # loop though batches
  

  
  # DRAW update clustering
  
  load("RData/updated_full_set-re-t-SNE.RData")
    
  plot_updated_cell_clusters_pdf("out/tSNE-updated-clustering_labeled.pdf",
                             full_set,
                             color_clusters,
                             TRUE)
  
  plot_updated_cell_clusters_pdf("out/tSNE-updated-clustering_unlabeled.pdf",
                             full_set,
                             color_clusters,
                             FALSE)

  hmap <- function() {
    # Colors for the heatmap
    color <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)
    #color <- colorRampPalette(c("blue", "red"))(100)

    expr_heat_filename <- "RData/expr_heat.RData"

    if (file.exists(expr_heat_filename)) {
      load(expr_heat_filename)
    } else {
      
      rng <- colQuantiles(full_set$expr, probs = c(0.01, 0.99))
      expr01 <- t((t(full_set$expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
      rm (rng)
      expr01[expr01 < 0] <- 0
      expr01[expr01 > 1] <- 1

      # Calculate the median expression
      expr_median <- data.frame(full_set$expr, cell_clustering = full_set$updated_clustering.cell) %>%
        group_by(cell_clustering) %>% 
        summarize_all(funs(median))
      
      expr01_median <- data.frame(expr01, cell_clustering = full_set$updated_clustering.cell) %>%
        group_by(cell_clustering) %>% 
        summarize_all(funs(median))
      
      # This clustering is based on the markers that were used for the main clustering
      d <- dist(expr_median[, colnames(full_set$expr)], method = "euclidean")
      cluster_rows <- hclust(d, method = "average")
      
      expr_heat <- as.matrix(expr01_median[, colnames(expr01)])
      rownames(expr_heat) <- expr01_median$cell_clustering
      
      save(expr_heat, file=expr_heat_filename)      
      
    }

    labels_row <- with(list(x = setNames(full_set$updated_clustering_names.cell,as.character(full_set$clustering.cell))), x[as.character(sort(as.numeric(unique(names(x)))))])
    
    pdf(file="out/HMAP-updated-clustering-all-samples.pdf", width = 5.5, height = 4)
    
    pheatmap(expr_heat[-2, rownames(lineage_marker_panel)], color = color, main="",
             cluster_cols = TRUE, cluster_rows = TRUE, 
             labels_col =lineage_marker_panel$cname, #labels_row = labels_row, 
             display_numbers = FALSE, number_color = "black", 
             border_color = "white") 
    
        
    dev.off()
    
  } # hmap_function 
  
  hmap()

  # PLOT UPDATED CELL CLUSTERS, WITHOUT THE CELLS marked to be removed

  plot_updated_cell_clusters_wo_cells_marked_to_be_removed_pdf <- function(filename,
                                                                           set,
                                                                           color_clusters,
                                                                           labels) {
    
    # TSNE with coloring by clustering, with updated cluster names
    dev.flush()
    pdf(file=filename)
    cluster_ids <- sort(unique(set$updated_clustering.cell))
    
    id_of_removed_cell <- match("removed", set$clusternames)
    bool_retained_cells <- (set$dr$updated_clustering.cell != id_of_removed_cell)
    
    plot(set$dr$tSNE1[bool_retained_cells], set$dr$tSNE2[bool_retained_cells], bty="n", pch=19, cex=0.05, col=color_clusters[set$dr$updated_clustering.cell[bool_retained_cells]], axes=FALSE, xlab="", ylab="")
    while (!is.null(dev.list())) dev.off()
  }

  plot_updated_cell_clusters_wo_cells_marked_to_be_removed_pdf("out/tSNE-updated-clustering_unlabeled_without_cells_to_be_removed.pdf",
                                                               full_set,
                                                               color_clusters,
                                                               FALSE)
  
  
  
  # DRAW tSNE about timepoints
  batch <- "MULTIPLE"
  cell_type <- "NK_cells"
  marker <- "CD161"
      
  marker_tSNE_plot(paste("out/functional-markers-case-vs-ctrl-fullset-tSNE-", batch,"-",cell_type,"-", marker, ".png", sep=""),
                   cell_type,
                   full_set,
                   color_clusters,
                   batch,
                   full_set$panel,
                   full_set$dr,
                   full_set$expr[full_set$tsne_inds, ],
                   marker)
    

  
  # PROCESS SUBCELLS
  
  for (subtype in c("CD4_T_cells", "CD8_T_cells")) {
    
    subset_file <- paste0("RData/updated-celltype-set-",subtype,".RData")
    
    load(subset_file)

    plot_updated_cell_clusters_pdf(paste0("out/", subtype, "/tSNE-subcell-clustering-", subtype,"-labeled",".pdf"),
                          celltype_set,
                          color_clusters,
                          TRUE)
    
    plot_updated_cell_clusters_pdf(paste0("out/", subtype,"/tSNE-subcell-clustering-", subtype,"-unlabeled",".pdf"),
                               celltype_set,
                               color_clusters,
                               FALSE)
    
  } # subtype
  
  
  for (subtype in c("CD4_T_cells", "CD8_T_cells")) {

    subtype_lineage_marker_panel <- panel[(panel[subtype] == "lineage/functional") | (panel[subtype] == "lineage"), ]
    subtype_lineage_markers_cname <- subtype_lineage_marker_panel[["cname"]]
    subtype_lineage_markers <- rownames(subtype_lineage_marker_panel)


    generate_expr_heat <- function(this_set) {
      
        rng <- colQuantiles(this_set$expr, probs = c(0.01, 0.99))
        expr01 <- t((t(this_set$expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
        rm (rng)
        expr01[expr01 < 0] <- 0
        expr01[expr01 > 1] <- 1
        
        # Calculate the median expression
        expr_median <- data.frame(this_set$expr, cell_clustering = this_set$updated_clustering.cell) %>%
          group_by(cell_clustering) %>% 
          summarize_all(funs(median))
        
        expr01_median <- data.frame(expr01, cell_clustering = this_set$updated_clustering.cell) %>%
          group_by(cell_clustering) %>% 
          summarize_all(funs(median))
        
        # This clustering is based on the markers that were used for the main clustering
        d <- dist(expr_median[, colnames(this_set$expr)], method = "euclidean")
        cluster_rows <- hclust(d, method = "average")
        
        expr_heat <- as.matrix(expr01_median[, colnames(expr01)])
        rownames(expr_heat) <- expr01_median$cell_clustering

        return (expr_heat)

    }

    expr_heat_filename <- paste0("RData/expr_heat_",subtype,".RData")
    
    if (file.exists(expr_heat_filename)) {
      load(expr_heat_filename)
    } else {
      subset_file <- paste0("RData/updated-celltype-set-",subtype,".RData")
      load(subset_file)
      expr_heat <- generate_expr_heat(celltype_set)
      save(expr_heat, file=expr_heat_filename)      
    }  
    
    hmap <- function(expr_heat, subtype, marker_panel) {
      color <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)
      #labels_row <- with(list(x = setNames(this_set$updated_clustering_names.cell,as.character(this_set$clustering.cell))), x[as.character(sort(as.numeric(unique(names(x)))))])

      pdf(file=paste0("out/",subtype,"/HMAP-updated-clustering-",subtype,".pdf"), width = 5.5, height = 4)

      pheatmap(expr_heat[, rownames(marker_panel)], color = color, main="",
               cluster_cols = TRUE, cluster_rows = TRUE,
               labels_col = marker_panel$cname, #labels_row = labels_row,
               display_numbers = FALSE, number_color = "black",
               #fontsize = 60, fontsize_number = 40,
               #annotation_row = annotation_row, annotation_colors = annotation_colors[-2],
               #annotation_legend = annotation_legend,
               border_color = "white")

      dev.off()
      
    } # hmap_function 
    
    hmap(expr_heat, subtype, subtype_lineage_marker_panel)
    
        
  } # subtype
  
  
  # GENERATE SUBCELLTYPE METACLUSTERING
  for (subtype in c("CD4_T_cells", "CD8_T_cells")) {
    
    subtype_lineage_marker_panel <- panel[(panel[subtype] == "lineage/functional") | (panel[subtype] == "lineage"), ]
    subtype_lineage_markers_cname <- subtype_lineage_marker_panel[["cname"]]
    subtype_lineage_markers <- rownames(subtype_lineage_marker_panel)
    
    generate_expr_heat_metaclustering <- function(this_set) {
      
      rng <- colQuantiles(this_set$expr, probs = c(0.01, 0.99))
      expr01 <- t((t(this_set$expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
      rm (rng)
      expr01[expr01 < 0] <- 0
      expr01[expr01 > 1] <- 1
      
      # Calculate the median expression
      expr_median <- data.frame(this_set$expr, cell_clustering = this_set$clustering.cell) %>%
        group_by(cell_clustering) %>% 
        summarize_all(funs(median))
      
      expr01_median <- data.frame(expr01, cell_clustering = this_set$clustering.cell) %>%
        group_by(cell_clustering) %>% 
        summarize_all(funs(median))
      
      # This clustering is based on the markers that were used for the main clustering
      d <- dist(expr_median[, colnames(this_set$expr)], method = "euclidean")
      cluster_rows <- hclust(d, method = "average")
      
      expr_heat <- as.matrix(expr01_median[, colnames(expr01)])
      rownames(expr_heat) <- expr01_median$cell_clustering
      
      return (expr_heat)
      
    }

    expr_heat_filename <- paste0("RData/expr_metaclustering_heat_",subtype,".RData")
    
    if (file.exists(expr_heat_filename)) {
      load(expr_heat_filename)
    } else {
      subset_file <- paste0("RData/celltype-set-",subtype,".RData")
      load(subset_file)
      expr_heat <- generate_expr_heat_metaclustering(celltype_set)
      save(expr_heat, file=expr_heat_filename)      
    }  
    
    hmap <- function(expr_heat, subtype, marker_panel) {
      color <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)
      #labels_row <- with(list(x = setNames(this_set$updated_clustering_names.cell,as.character(this_set$clustering.cell))), x[as.character(sort(as.numeric(unique(names(x)))))])
      
      pdf(file=paste0("out/",subtype,"/HMAP-meta-clustering-",subtype,".pdf"), width = 5.5, height = 4)
      
      pheatmap(expr_heat[, rownames(marker_panel)], color = color, main="",
               cluster_cols = TRUE, cluster_rows = TRUE,
               labels_col = marker_panel$cname, #labels_row = labels_row,
               display_numbers = FALSE, number_color = "black",
               #fontsize = 60, fontsize_number = 40,
               #annotation_row = annotation_row, annotation_colors = annotation_colors[-2],
               #annotation_legend = annotation_legend,
               border_color = "white")
      
      dev.off()
      
    } # hmap_function 
    
    hmap(expr_heat, subtype, subtype_lineage_marker_panel)
    
    
        
  } # subtype  


  # PLOT without cells in clusters marked to be deleted
  for (subtype in c("CD4_T_cells", "CD8_T_cells")) {
    
    subset_file <- paste0("RData/updated-celltype-set-",subtype,".RData")
    
    load(subset_file)

    plot_updated_cell_clusters_wo_cells_marked_to_be_removed_pdf(paste0("out/", subtype,"/tSNE-subcell-clustering-", subtype,"-unlabeled",".pdf"),
                                   celltype_set,
                                   color_clusters,
                                   FALSE)
    
  } # subtype
  

  plot_updated_cell_clusters_wo_cells_marked_to_be_removed_pdf <- function(filename,
                                                                           set,
                                                                           color_clusters,
                                                                           labels) {
    
    # TSNE with coloring by clustering, with updated cluster names
    dev.flush()
    pdf(file=filename)
    cluster_ids <- sort(unique(set$updated_clustering.cell))
    
    id_of_removed_cell <- match("CD45RAplusCCR7plus", set$clusternames)
    bool_retained_cells <- (set$dr$updated_clustering.cell != id_of_removed_cell)
    
    plot(set$dr$tSNE1[bool_retained_cells], set$dr$tSNE2[bool_retained_cells], bty="n", pch=19, cex=0.05, col=color_clusters[set$dr$updated_clustering.cell[bool_retained_cells]], axes=FALSE, xlab="", ylab="")
    # if (labels) {
    #   text(set$updated_agg$tSNE1, set$updated_agg$tSNE2, paste(cluster_ids, set$clusternames[cluster_ids], sep=". "))
    #   legend("topright", legend=paste(cluster_ids, set$clusternames[cluster_ids], sep=". "), pch=15, col=color_clusters, ncol=2, bty="n")
    # }
    while (!is.null(dev.list())) dev.off()
  }
  
    
  
  
  # OLS plots BATCH
  
  for (subtype in c("CD4_T_cells", "CD8_T_cells")) {
    
    subtype_functional_marker_panel <- panel[(panel[subtype] == "lineage/functional") | (panel[subtype] == "functional"), ]
    subtype_functional_markers_cname <- subtype_functional_marker_panel[["cname"]]
    subtype_functional_markers <- rownames(subtype_functional_marker_panel)

    for (batch in c("IAA", "GADA", "MULTIPLE")) {

      intensity_table_rdata_filename <- paste0("RData/marker_intensity_tables_", subtype, "_" , batch, ".RData")
      load(intensity_table_rdata_filename)
      
      proportion_table_rdata_filename = paste0("RData/proportion_table-", subtype, "_" , batch, ".RData")
      load(proportion_table_rdata_filename)
      
      combined_proportion_table <- with(
        list(
          left = proportion_table,
          right = anno),
        {
          shared_names <- intersect(rownames(left),rownames(right))
          cbind(left[shared_names,], right[shared_names,])
        })
      
#      for (j in seq_along(marker_intensity_tables)) {
#        celltype <- names(marker_intensity_tables)[[j]]

        if (subtype == "CD4_T_cells") {

          celltype <- "CCR4plus"
          marker <- "CD161"
          celltype_index <- match(celltype, names(marker_intensity_tables))

          OLS_plot_case_vs_control_per_type(paste0("out/",subtype,"/OLS-Case-vs-Control_", batch, "_", celltype, ".pdf"),
                                   marker_intensity_tables[[celltype_index]],
                                   marker,
                                   "Marker intensity",
                                   "Age (years)",
                                   ylim=c(0.3,2.5))
          
          # celltype <- "CD45RAplusCCR7plus"
          # marker <- "CD161"
          # celltype_index <- match(celltype, names(marker_intensity_tables))
          # 
          # OLS_plot_case_vs_control_per_type(paste0("out/",subtype,"/OLS-Case-vs-Control_", batch, "_", celltype, ".pdf"),
          #                          marker_intensity_tables[[celltype_index]],
          #                          marker,
          #                          "Marker intensity",
          #                          "Age (years)",
          #                          ylim=c(0,0.25)) 
          # 
          # 
          # celltype <- "CD57plus"
          # marker <- "TIGIT"
          # celltype_index <- match(celltype, names(marker_intensity_tables))
          # 
          # OLS_plot_case_vs_control_per_type(paste0("out/",subtype,"/OLS-Case-vs-Control_", batch, "_", celltype, ".pdf"),
          #                                   marker_intensity_tables[[celltype_index]],
          #                                   marker,
          #                                   "Marker intensity",
          #                                   "Age (years)",
          #                                   ylim=c(0,1.7)) 
          # 
          # celltype <- "CD25plusCD127minus"
          # marker <- "CD39"
          # celltype_index <- match(celltype, names(marker_intensity_tables))
          # 
          # OLS_plot_case_vs_control_per_type(paste0("out/",subtype,"/OLS-Case-vs-Control_", batch, "_", celltype, ".pdf"),
          #                                   marker_intensity_tables[[celltype_index]],
          #                                   marker,
          #                                   "Marker intensity",
          #                                   "Age (years)",
          #                                   ylim=c(0,3.6)) 
          # 
          # celltype <- "HLA_DRplusICOSplus"
          # marker <- "CD39"
          # celltype_index <- match(celltype, names(marker_intensity_tables))
          # 
          # OLS_plot_case_vs_control_per_type(paste0("out/",subtype,"/OLS-Case-vs-Control_", batch, "_", celltype, ".pdf"),
          #                                   marker_intensity_tables[[celltype_index]],
          #                                   marker,
          #                                   "Marker intensity",
          #                                   "Age (years)",
          #                                   ylim=c(0,2)) 
          
          
        } else {
          
          celltype <- "CD27plusPD_1plus"
          marker <- "CD39"
          celltype_index <- match(celltype, names(marker_intensity_tables))
          
          OLS_plot_case_vs_control_per_type(paste0("out/",subtype,"/OLS-Case-vs-Control_", batch, "_", celltype, ".pdf"),
                                            marker_intensity_tables[[celltype_index]],
                                            marker,
                                            "Marker intensity",
                                            "Age (years)",
                                            ylim=c(0,25)) 
          
          
          
        }
        
#      }
      
      # PROCESS CELL TYPE PROPORTIONS
      # OLS_plot_case_vs_control(paste0("out/",subtype,"/OLS-Case-vs-Control_CellTypeProportions-", batch, ".pdf"),
      #                          combined_proportion_table,
      #                          names(marker_intensity_tables), # TODO: get cluster names from proportion table
      #                          "Proportion",
      #                          "Age (years)") 
      
      rm(marker_intensity_tables)
      rm(proportion_table)
    
    } # batch
    
  } # subset (CD4 vs CD8)

  # 

  load("RData/updated-celltype-set-CD4_T_cells.RData")    
  
  #including CD25+CD127– with the memory Treg phenotype
  
  batch <- "IAA"
  marker <- "CD39"
  subtype <- "CD4_T_cells"
  subcluster_name <- "CD25plusCD127minus" #8
  subcluster_index <- match(subcluster_name, celltype_set$clusternames)
  
  samples_case <- rownames(celltype_set$anno[(celltype_set$anno["CaseCtrl"] == "Case") & (celltype_set$anno["Batch"] == batch), ])
  samples_ctrl <- rownames(celltype_set$anno[(celltype_set$anno["CaseCtrl"] == "Control") & (celltype_set$anno["Batch"] == batch), ])
  
  marker_tSNE_plot_sub(paste0("out/",subtype,"/functional-markers-case-vs-ctrl-subset-",subtype,"-tSNE-", batch,"-",subcluster_index,"-", marker, ".png"),
                       subtype,
                       subcluster_index,
                       celltype_set,
                       color_clusters,
                       batch,
                       celltype_set$panel,
                       celltype_set$dr,
                       celltype_set$expr[celltype_set$tsne_inds, ],
                       marker)


  load("RData/updated_full_set-re-t-SNE.RData")
  
  tSNE_color_by_sample_separate_images <- function(filename,
                                                   dr,
                                                   agg,
                                                   table) {
    
    # TSNE with coloring by samples, each sample separately
    print ("TSNE with coloring by per sample")
    while (!is.null(dev.list())) dev.off()
    
    #png(file=filename, width = 6000, height = 6000, units = "px", pointsize = 12)
    #par(mfrow = determine_image_tiling(max(dr$sampleindices.cell)))
  
    bc_samples <- rownames(anno)[grep("_BC_", rownames(anno))]
    
    for (sampleName in bc_samples) {
      png(file=paste0("out/tSNE-sample-distribution-in-metaclusters-",sampleName,".png"))
      plot(dr$tSNE1, dr$tSNE2, bty="n", pch=19, cex=0.5, col="#20202010", axes=FALSE, xlab="", ylab="", main=sampleName, cex.main=1)
      points(dr$tSNE1[dr$sampleids.cell == sampleName], dr$tSNE2[dr$sampleids.cell == sampleName], bty="n", pch=19, cex=0.5, col="#2020FF90", axes=FALSE, xlab="", ylab="")
      text(agg$tSNE1, agg$tSNE2, agg$markerTags, cex=4)
      while (!is.null(dev.list())) dev.off()
    }
  }
 
  

  tSNE_color_by_sample_separate_images("out/tSNE-sample-distribution-in-metaclusters.png",
                                       full_set$dr,
                                       full_set$agg,
                                       full_set$table)
  
  
  
  
  

}
