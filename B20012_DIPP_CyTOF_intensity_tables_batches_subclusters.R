
## BEGIN PIPELINE EXECUTION

B200012_DIPP_CyTOF_workflow_intensity_tables_batches_subclusters <- function(CyTOF_workflow_cwd,
                                                                 sample_metadata_tsv,
                                                                 subclusters) {
  
 
  dir.create("out/")
  dir.create("out/tables/")
  
  anno <- read.csv(sample_metadata_tsv, stringsAsFactors=FALSE, sep="\t")
  rownames(anno) <- anno[["Sample"]]
  
  load("RData/panel.RData")
  
  for (celltype in subclusters) {
  
  functional_markers_panel <- panel[(panel[celltype] == "functional") | (panel[celltype] == "lineage/functional") , ]
  functional_markers_cname <- functional_markers_panel[["cname"]]
  markers <- rownames(functional_markers_panel)
  
  
  load(paste0("RData/updated-celltype-set-",celltype,".RData"))

  for (batch in unique(anno[["Batch"]])) {
    
    batch_table <- celltype_set$table
    batch_table["row_index"] <- 1:nrow(batch_table)
    batch_table <- batch_table[celltype_set$anno[["Batch"]] == batch, ]
    batch_anno <- celltype_set$anno[celltype_set$anno[["Batch"]] == batch, ]
    
    intensity_table_rdata_filename <- paste0("RData/marker_intensity_tables_", celltype,"_", batch, ".RData")
    
    if (file.exists(intensity_table_rdata_filename)) {
      load(intensity_table_rdata_filename)
    } else {
      
      marker_intensity_tables <- make_marker_intensity_tables(celltype_set$expr,
                                                              markers,
                                                              celltype_set$updated_clustering_names.cell,
                                                              batch_table, 
                                                              celltype_set$sampleindices.cell,
                                                              celltype_set$clusternames,
                                                              batch_anno,
                                                              panel)
      
      save(marker_intensity_tables, file=intensity_table_rdata_filename)
    }
    
    proportion_table_rdata_filename = paste0("RData/proportion_table-", celltype, "_", batch, ".RData")
    if (file.exists(proportion_table_rdata_filename)) {
      load(proportion_table_rdata_filename)
    } else {
      
      proportion_table <- make_proportion_table(celltype_set$updated_clustering_names.cell,
                                                batch_table, #full_set$table,
                                                celltype_set$sampleindices.cell,
                                                celltype_set$clusternames,
                                                batch_anno) #full_set$anno)
      
      save(proportion_table, file = proportion_table_rdata_filename)
    }
    
    # WRITE TSV-tables
    write.csv(proportion_table, paste0("out/tables/proportion-table-", celltype, "-", batch, ".csv"))
  } # loop through batches
  
  } # loop through subclusters 
  
}
