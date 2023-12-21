
## BEGIN PIPELINE EXECUTION

DIPP_CyTOF_workflow_intensity_tables <- function(CyTOF_workflow_cwd,
                           sample_metadata_tsv) {

  setwd(CyTOF_workflow_cwd)
  
  anno <- read.csv(sample_metadata_tsv, stringsAsFactors=FALSE, sep="\t")
  rownames(anno) <- anno[["Sample"]]
  
  load("RData/panel.RData")
  
  markers_panel <- panel[(panel["type"] == "functional") | (panel["type"] == "lineage/functional") | (panel["type"] == "lineage"), ]
  markers_cname <- markers_panel[["cname"]]
  markers <- rownames(markers_panel)
  
  load("RData/updated_full_set.RData")
 
  batch_table <- full_set$table
  batch_table["row_index"] <- 1:nrow(batch_table)

  intensity_table_rdata_filename <- paste0("RData/marker_intensity_tables", ".RData")
  
  if (file.exists(intensity_table_rdata_filename)) {
    load(intensity_table_rdata_filename)
  } else {
    marker_intensity_tables <- make_marker_intensity_tables(full_set$expr,
                                                            markers,
                                                            full_set$updated_clustering_names.cell,
                                                            batch_table, # full_set$table,
                                                            full_set$sampleindices.cell,
                                                            full_set$clusternames,
                                                            full_set$anno,
                                                            panel)
    
    save(marker_intensity_tables, file=intensity_table_rdata_filename)
  }

  proportion_table_rdata_filename = paste0("RData/proportion_table",".RData")
  if (file.exists(proportion_table_rdata_filename)) {
    load(proportion_table_rdata_filename)
  } else {
    
    proportion_table <- make_proportion_table(full_set$updated_clustering_names.cell,
                                              batch_table, #full_set$table,
                                              full_set$sampleindices.cell,
                                              full_set$clusternames,
                                              full_set$anno)
    
    save(proportion_table, file = proportion_table_rdata_filename)
  }
  
  # WRITE TSV-tables
  write.csv(proportion_table, paste0("out/tables/proportion-table", ".csv"))
  

    
}
