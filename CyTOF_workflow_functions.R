library(tsne)
library(Rtsne)

redo_tSNE <- function(curr_set,
                      clustering_markers_,
                      cells_per_sample_) {
  
  seed_number = 12345
  
  dups_ <- which(!duplicated(curr_set$expr[, clustering_markers_])) ## Find and skip duplicates
  inds_ <- split(1:length(curr_set$sampleids.cell), curr_set$sampleids.cell) ## Data subsampling: create indices by sample
  set.seed(seed_number)
  tsne_inds_ <- lapply(names(inds_), function(i) {
    intersect(sample(inds_[[i]], min(length(inds_[[i]]), cells_per_sample_), replace = FALSE), dups_)
  })
  tsne_inds_ <- unlist(tsne_inds_)
  rm(inds_)
  rm(dups_)
  
  set.seed(seed_number)
  tsne_out_ <- Rtsne(curr_set$expr[tsne_inds_, clustering_markers_], check_duplicates = FALSE, pca = FALSE)
  dr_ <- data.frame(tSNE1 = tsne_out_$Y[, 1], tSNE2 = tsne_out_$Y[, 2], curr_set$expr[tsne_inds_, clustering_markers_])

  dr_["updated_clustering.cell"] <- curr_set$updated_clustering.cell[tsne_inds_]
  dr_["updated_clustering_names.cell"] <- curr_set$updated_clustering_names.cell[tsne_inds_]
  
  updated_agg <- aggregate(dr_,
                           by = list(factor(curr_set$updated_clustering.cell[tsne_inds_], levels = 1:max(curr_set$updated_clustering.cell[tsne_inds_]))),
                           FUN = mean)

  dr_$sampleids.cell <- curr_set$sampleids.cell[tsne_inds_]
  dr_$sampleindices.cell <- curr_set$sampleindices.cell[tsne_inds_]
  
  curr_set$tsne_inds <- tsne_inds_
  curr_set$dr <- dr_
  curr_set$updated_agg <- updated_agg
    
  return (curr_set)

}

make_marker_intensity_tables <- function(expr,
                                         markers,
                                         updated_clustering_names.cell,
                                         table,
                                         sampleindices.cell,
                                         cell_cluster_names,
                                         anno,
                                         panel) {
  
  
  
  # Write per sample cell type marker intensity table
  marker_intensity_tables = list()
  for (i in 1:length(cell_cluster_names)) { # for each cell types
    marker_intensity_tables[[i]] <- NULL
    cl_name <- cell_cluster_names[i]
    print (cl_name)
    CLUSTER <- (updated_clustering_names.cell == cl_name)
    tbl <- NULL
    for (j in table[["row_index"]]) { # loop through samples
      SAMPLE <- (sampleindices.cell == j)
      cm <- colMeans(expr[SAMPLE & CLUSTER, markers, drop=FALSE])
      tbl <- rbind(tbl, cm)
    }
    # Rename metals with marker names
    colnames(tbl) <- panel[colnames(tbl), "cname"]
    
    tbl <- as.data.frame(cbind(tbl, anno))
    rownames(tbl) <- table[["sample"]]

    marker_intensity_tables[[i]] <- tbl
  }
  names(marker_intensity_tables) <- cell_cluster_names
  return (marker_intensity_tables)
}


make_proportion_table <- function(updated_clustering_names.cell,
                                  table,
                                  sampleindices.cell,
                                  cell_cluster_names,
                                  anno) {
  
  
  # Write per sample cell type proportion table
  proportion_table = NULL
  for (cl_name in cell_cluster_names) { # for each cell types
    CLUSTER <- (updated_clustering_names.cell == cl_name)
    sample_proportions = c()
    for (j in table[["row_index"]]) { # loop through samples
      SAMPLE <- (sampleindices.cell == j)
      sample_proportions <- c(sample_proportions, (sum(SAMPLE & CLUSTER, na.rm=TRUE) / sum(SAMPLE, na.rm=TRUE)))
    }
    proportion_table <- cbind(proportion_table, sample_proportions)
  }
  colnames(proportion_table) <- cell_cluster_names
  rownames(proportion_table) <- table[["sample"]]
  
  #proportion_table <- cbind(proportion_table, anno)
  
  # Drop samples without a proper timepoint (mostly BC samples).
  return (proportion_table)
}


