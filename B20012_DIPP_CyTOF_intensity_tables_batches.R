
## BEGIN PIPELINE EXECUTION

B200012_DIPP_CyTOF_workflow_intensity_tables_batches <- function(CyTOF_workflow_cwd,
                                                                 sample_metadata_tsv) {
  
  dir.create("out/")
  dir.create("out/tables/")
  
  anno <- read.csv(sample_metadata_tsv, stringsAsFactors=FALSE, sep="\t")
  rownames(anno) <- anno[["Sample"]]
  
  load("RData/panel.RData")
  
  #  functional_markers_panel <- panel[(panel["type"] == "functional") | (panel["type"] == "lineage/functional"), ]
  #  functional_markers_cname <- functional_markers_panel[["cname"]]
  
  markers_panel <- panel[(panel["type"] == "functional") | (panel["type"] == "lineage/functional") | (panel["type"] == "lineage"), ]
  markers_cname <- markers_panel[["cname"]]
  markers <- rownames(markers_panel)
  
  load("RData/updated_full_set.RData")
  
  
  for (batch in unique(anno[["Batch"]])) {
    
    batch_table <- full_set$table
    batch_table["row_index"] <- 1:nrow(batch_table)
    batch_table <- batch_table[full_set$anno[["Batch"]] == batch, ]
    batch_anno <- full_set$anno[full_set$anno[["Batch"]] == batch, ]
    
    intensity_table_rdata_filename <- paste0("RData/marker_intensity_tables_", batch, ".RData")
    
    if (file.exists(intensity_table_rdata_filename)) {
      load(intensity_table_rdata_filename)
    } else {
      
      marker_intensity_tables <- make_marker_intensity_tables(full_set$expr,
                                                              markers,
                                                              full_set$updated_clustering_names.cell,
                                                              batch_table, # full_set$table,
                                                              full_set$sampleindices.cell,
                                                              full_set$clusternames,
                                                              batch_anno,
                                                              panel)
      
      save(marker_intensity_tables, file=intensity_table_rdata_filename)
    }
    
    proportion_table_rdata_filename = paste0("RData/proportion_table-", batch, ".RData")
    if (file.exists(proportion_table_rdata_filename)) {
      load(proportion_table_rdata_filename)
    } else {
      
      proportion_table <- make_proportion_table(full_set$updated_clustering_names.cell,
                                                batch_table, #full_set$table,
                                                full_set$sampleindices.cell,
                                                full_set$clusternames,
                                                batch_anno) #full_set$anno)
      
      save(proportion_table, file = proportion_table_rdata_filename)
    }
    
    # WRITE TSV-tables
    write.csv(proportion_table, paste0("out/tables/proportion-table-", batch, ".csv"))
  } # loop through batches
  
  
  
}
