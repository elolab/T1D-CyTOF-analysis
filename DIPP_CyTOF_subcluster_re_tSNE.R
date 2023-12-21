  for (subtype in c("CD4_T", "CD8_T")) {
    
    subset_file <- paste0("RData/celltype-set-",subtype,"_cells.RData")
    
    if (file.exists(subset_file)) {
  
      load(subset_file)
      
    } else {
  
      cl_name <- paste0(subtype, "_cells")
      
      load(paste0("celltype-set-",subtype,"_cells.RData"))
      
      cluster_merging_df <- read.table(paste0("../","cluster_merging_", cl_name, ".tsv"), sep="\t", header=TRUE, stringsAsFactors=FALSE)
      celltype_set <- update_clustering(cluster_merging_df, celltype_set)
  
      selection <- (panel[cl_name] == "lineage" | panel[cl_name] == "lineage/functional")
      celltype_clustering_markers <- rownames(panel)[selection]
      
      celltype_set <- redo_tSNE(celltype_set,
                                celltype_clustering_markers,
                                4000);
      
      save(celltype_set, file = paste0("celltype-set-",subtype,"_cells-re-t.RData"))
  
    }    
