
B200012_DIPP_CyTOF_workflow_re_tSNE <- function(CyTOF_workflow_cwd,
                                            markers_tsv) {
  
  setwd(CyTOF_workflow_cwd)
 
  panel <- read.table(markers_tsv, sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE)

  if (file.exists("RData/updated_full_set-re-t-SNE.RData")) {
    
    print ("RData/updated_full_set-re-t-SNE.RData already exists, skipping.")
    
  } else {
    
    load("RData/updated_full_set.RData")
    
    clustering_markers <- rownames(panel)[panel$type=="lineage" | panel$type=="lineage/functional"]
    
    full_set <- redo_tSNE(full_set,
                          clustering_markers,
                          3500);
    
    save(full_set, file = "RData/updated_full_set-re-t-SNE.RData")
  } 
  

}