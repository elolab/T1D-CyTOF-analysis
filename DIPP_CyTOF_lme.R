library(foreach)
library(writexl)

lmesummary_to_df <- function(lme_results) {
  # LME for each row
  lme.out <- foreach(i=names(lme_results), .combine=rbind, .multicombine=TRUE) %do% {
    lme_result <- lme_results[[i]]
    s <- summary(lme_result)
    coef <- s$coefficients[-1,1]
    names(coef) <- paste("coef.",names(coef),sep="")
    p <- s$coefficients[-1,5]
    names(p) <- paste("p.",names(p),sep="")
    c(coef, p)
  }
  fdr <- apply(lme.out[,grep("p.",colnames(lme.out))], 2, function(x) p.adjust(x, method="BH"))
  colnames(fdr) <- gsub("p\\.", "fdr\\.", colnames(fdr))
  lme.out <- data.frame(lme.out, fdr, stringsAsFactors=FALSE)
  rownames(lme.out) <- names(lme_results)
  return (lme.out)
}


DIPP_CyTOF_workflow_lme_significant_results <- function(CyTOF_workflow_cwd) {
  
  setwd(paste0(CyTOF_workflow_cwd, "/out/tables/lme/" ))
  
  files <- list.files("./", pattern="celltype_")
  variables <- c("CaseCtrlControl","Age.CaseCtrlControl","Age")
  
  # Gather all results
  results.all <- foreach(i=1:length(variables)) %do% {
    results.variable <- foreach(j=1:length(files), .combine=rbind) %do% {
      # Read results
      results <- read.delim(files[j],sep="\t")
      # Parse details from filename
      names <- unlist(strsplit(files[j],"-"))
      if(length(names)==3) {
        cells <- gsub(".tsv","",names[3])
        subset <- names[2]
      }
      if(length(names)==2) {
        cells <- gsub(".tsv","",names[2])
        subset <- "NONE"
      }
      # Generate matrix
      results <- cbind(subset,cells,marker=results[,1],results[,c(paste(c("coef","p","fdr"),variables[i],sep="."))])
      colnames(results) <- c("Subset","Cell type","Marker","Coefficient","P-value","FDR")
      # Pass matrix for merging
      results
    }
    results.variable
  }
  names(results.all) <- variables
  
  # Sort and filter
  for(i in 1:length(results.all)) {
    results.all[[i]] <- results.all[[i]][order(results.all[[i]]$`P-value`, decreasing=FALSE),]
    results.all[[i]] <- results.all[[i]][which(results.all[[i]]$`P-value`<=0.01),]
  }
  
  # Store results
  write_xlsx(results.all, path="DIPP_CyTOF.xlsx")
}

DIPP_CyTOF_workflow_lme_generate_venns_of_methods <- function(CyTOF_workflow_cwd) {
  
  library(readxl)
  library(ggvenn)
  library(UpSetR)
  
  # Folder names
  folders <- c("lme-hla-id", "lme-hla-pair", "lme-hla-pair-Nid", "lme-id", "lme-pair", "lme-pair-Nid")
  
  # Gather and plot summaries for group
  results <- foreach(f=folders) %do% {
    df <- data.frame(read_xlsx(paste(f,"/","DIPP_CyTOF.xlsx",sep=""), sheet=1))
    apply(df[,1:3], 1, function(x) paste(x, collapse=", "))
  }
  names(results) <- folders
  
  pdf("LME_overlaps_group.pdf", width=5, height=5)
  ggvenn(results[1:3], stroke_size=0.5, fill_color=rep("lightgrey", 3), set_name_size=6, text_size=3, show_elements=FALSE, 
         show_percentage=TRUE)
  ggvenn(results[4:6], stroke_size=0.5, fill_color=rep("lightgrey", 3), set_name_size=6, text_size=3, show_elements=FALSE, 
         show_percentage=TRUE)
  upset(fromList(results), nsets=length(results), sets=names(results), keep.order=FALSE, matrix.color="darkblue")
  dev.off()
  
  # Gather and plot summaries for group*age
  results <- foreach(f=folders) %do% {
    df <- data.frame(read_xlsx(paste(f,"/","DIPP_CyTOF.xlsx",sep=""), sheet=2))
    apply(df[,1:3], 1, function(x) paste(x, collapse=", "))
  }
  names(results) <- folders
  
  pdf("LME_overlaps_group_age.pdf", width=5, height=5)
  ggvenn(results[1:3], stroke_size=0.5, fill_color=rep("lightgrey", 3), set_name_size=6, text_size=3, show_elements=FALSE, 
         show_percentage=TRUE)
  ggvenn(results[4:6], stroke_size=0.5, fill_color=rep("lightgrey", 3), set_name_size=6, text_size=3, show_elements=FALSE, show_percentage=TRUE)
  upset(fromList(results), nsets=length(results), sets=names(results), keep.order=FALSE, matrix.color="darkblue")
  dev.off()
}


craft_lme_row <- function(f, pr_int_table, prefix_list) {

  if (FALSE) {  
    pr_int_table <- marker_intensity_tables[[j]]
    prefix_list <- c(clustername=clustername, marker=marker)
    
    pr_int_table <- combined_proportion_table
    prefix_list <- c()
  }
  
  suppressWarnings(suppressMessages(fit1 <- lmer(f, data = pr_int_table, na.action = na.omit, REML=FALSE)))

  out.coef <- coef(summary(fit1))[ ,"Estimate"]; names(out.coef) <- paste("coef",names(out.coef),sep=".")
  out.p <-coef(summary(fit1))[ ,"Pr(>|t|)"]; names(out.p) <- paste("p",names(out.p),sep=".")
  
  row <- c(prefix_list, out.coef, out.p)
  return (row)
}

add_fdr <- function(lme.table) {
  fdr <- apply(lme.table[,grep("^p\\.",colnames(lme.table))], 2, function(x) p.adjust(x, method="BH"))
  colnames(fdr) <- gsub("^p\\.", "fdr\\.", colnames(fdr))
  lme.table <- data.frame(lme.table, fdr, stringsAsFactors=FALSE)
  return (lme.table)
}


DIPP_CyTOF_workflow_lme <- function(CyTOF_workflow_cwd,
                                        sample_metadata_tsv,
                                        markers_tsv,
                                        sub_cluster,
                                        formula_str,
                                        formula_str_per_batch,
                                        whole_set_flag,
                                        batched_set_flag) {

  setwd(CyTOF_workflow_cwd)
  
  dir.create(paste0("out/",sub_cluster,"/tables/"))
  dir.create(paste0("out/",sub_cluster,"/tables/lme/"))
  
  fileConn<-file(paste0("out/",sub_cluster,"/tables/lme/formula.txt"))
  writeLines(c(formula_str, formula_str_per_batch), fileConn)
  close(fileConn)
  
  anno <- read.csv(sample_metadata_tsv, stringsAsFactors=FALSE, sep="\t")
  rownames(anno) <- anno[["Sample"]]
  
  #load("RData/panel.RData")
  panel <- read.table(markers_tsv, sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE)
  
  #coefs <- c("Age" , "SexM","CaseCtrlControl", "Age:CaseCtrlControl", "HLADR3_4", "HLADR4")
  
  #functional_markers_panel <- panel[(panel["type"] == "functional") | (panel["type"] == "lineage/functional"), ]
  #functional_markers_cname <- functional_markers_panel[["cname"]]

  if (sub_cluster != "") {  
    column <- sub_cluster
  } else {
    column <- "type" 
  }
  
  markers_panel <- panel[(panel[column] == "functional") | (panel[column] == "lineage/functional"), ]
  markers_cname <- markers_panel[["cname"]]
  markers <- rownames(markers_panel)
  
  # load("RData/updated_full_set.RData")

#  formula_str <- " ~ Age + Sex + Age*CaseCtrl + HLA + (1|Batch) + (1|Pair/Id)"
#  formula_str_per_batch <- " ~ Age + Sex + Age*CaseCtrl + HLA + (1|Pair/Id)"
  
#  formula_str <- " ~ Age + Sex + Age*CaseCtrl + (1|HLA) + (1|Batch) + (1|Pair/Id)"
#  formula_str_per_batch <- " ~ Age + Sex + Age*CaseCtrl + (1|HLA) + (1|Pair/Id)"

#  formula_str <- " ~ Age + Sex + Age*CaseCtrl + (1|HLA) + (1|Batch) + (1|Id)"
#  formula_str_per_batch <- " ~ Age + Sex + Age*CaseCtrl + (1|HLA) + (1|Id)"
  
  #formula_str <- " ~ Age + Sex + Age*CaseCtrl + (1|Pair:Id)"h
  #formula_str <- " ~ Age + Sex + Age*CaseCtrl + (1|HLA) + (1|Pair:Id)"
  
  # WHOLE SET
  
  if (whole_set_flag) {
  
    if (sub_cluster != "") {  
      proportion_table_rdata_filename = paste0("RData/proportion_table-", sub_cluster ,".RData")
      intensity_table_rdata_filename <- paste0("RData/marker_intensity_tables_", sub_cluster, ".RData")
    } else {
      proportion_table_rdata_filename = paste0("RData/proportion_table", ".RData")
      intensity_table_rdata_filename <- paste0("RData/marker_intensity_tables", ".RData")
    }
    
    load(proportion_table_rdata_filename)
    
    clusternames <- colnames(proportion_table)[!(colnames(proportion_table) %in% c(colnames(anno), "removed"))]
    
    load(intensity_table_rdata_filename)
    combined_proportion_table <- with(
      list(
        left = proportion_table,
        right = anno),
      {
        shared_names <- intersect(rownames(left),rownames(right))
        cbind(left[shared_names,], right[shared_names,])
      })
    
    # CELL PROPORTIONS LME 
    
    #lme.table <- matrix(nrow = 0, ncol=(12))
    lme.table <- NULL
    
    for (i in 1:length(clusternames)) {
      celltype <- clusternames[[i]]
      f <- formula(paste0(celltype, formula_str))

      row <- craft_lme_row(f, combined_proportion_table, c())
      #message(length(row))
      if (is.null(lme.table)) {
        lme.table <- row
      } else {
        lme.table <- rbind(lme.table, row)
      }
    }
    
    lme.table <- data.frame(lme.table)
    rownames(lme.table) <- clusternames
    lme.table <- add_fdr(lme.table)
    
    write.table(lme.table, file=paste0("out/",sub_cluster,"/tables/lme/LME_cell_type_proportions", ".tsv"), sep="\t", col.names = NA)
    
    
    # lme_results_cell_type_proportions <- list()
    # for (i in 1:length(clusternames)) {
    #   celltype <- clusternames[[i]]
    #   f <- formula(paste0(celltype, formula_str))
    #   lme_results_cell_type_proportions[[i]] <- lmerTest::lmer(f, data = combined_proportion_table, na.action = na.omit, REML=FALSE)
    # }
    # names(lme_results_cell_type_proportions) <- clusternames
    # 
    # df <- lmesummary_to_df(lme_results_cell_type_proportions)
    # write.table(df, file=paste0("out/",sub_cluster,"/tables/lme/LME_cell_type_proportions", ".tsv"), sep="\t", col.names = NA)
    
    # INTENSITY TABLE LME

    #lme.table <- matrix(nrow = 0, ncol=14)
    lme.table <- NULL
    
    for (j in 1:length(clusternames)) {
      clustername <- clusternames[[j]]
      
      marker_intensity_table <- marker_intensity_tables[[j]]
      marker_intensity_table <- marker_intensity_table[, markers_cname]
      
      combined_intensity_table <- with(
        list(
          left = marker_intensity_table,
          right = anno),
        {
          shared_names <- intersect(rownames(left),rownames(right))
          cbind(left[shared_names,], right[shared_names,])
        })
      
      
      for (i in 1:length(markers_cname)) {
        marker <- markers_cname[i]
        f <- formula(paste0(marker, formula_str))

        row <- craft_lme_row(f, combined_intensity_table, c(clustername=clustername, marker=marker))
        
        if (is.null(lme.table)) {
          lme.table <- row
        } else {
          lme.table <- rbind(lme.table, row)
        }
        
      }
    }
    
    lme.table <- data.frame(lme.table)
    lme.table <- add_fdr(lme.table)
    write.table(lme.table, file=paste0("out/",sub_cluster,"/tables/lme/celltype_LME_marker_intensities", ".tsv"), sep="\t", row.names = FALSE, col.names = TRUE)
    
  } # whole set flag
    
  # BATCH

  if (batched_set_flag) {
    
    for (batch in unique(anno[["Batch"]])) {
  
      if (sub_cluster != "") {  
        proportion_table_rdata_filename = paste0("RData/proportion_table-", sub_cluster,"_", batch, ".RData")
        intensity_table_rdata_filename <- paste0("RData/marker_intensity_tables_", sub_cluster, "_", batch, ".RData")
      } else {
        proportion_table_rdata_filename = paste0("RData/proportion_table-", batch, ".RData")
        intensity_table_rdata_filename <- paste0("RData/marker_intensity_tables_", batch, ".RData")
      }
  
      load(proportion_table_rdata_filename)
  
      clusternames <- colnames(proportion_table)[!(colnames(proportion_table) %in% c(colnames(anno), "removed"))]
      
      load(intensity_table_rdata_filename)
  
      combined_proportion_table <- with(
        list(
          left = proportion_table,
          right = anno),
        {
          shared_names <- intersect(rownames(left),rownames(right))
          cbind(left[shared_names,], right[shared_names,])
        })

      lme.table <- matrix(nrow = 0, ncol=12)
      for (i in 1:length(clusternames)) {
        celltype <- clusternames[[i]]
        f <- formula(paste0(celltype, formula_str_per_batch))
        
        row <- craft_lme_row(f, combined_proportion_table, c())
        lme.table <- rbind(lme.table, row)
      }
      
      lme.table <- data.frame(lme.table)
      rownames(lme.table) <- clusternames
      lme.table <- add_fdr(lme.table)
      
      if (sub_cluster != "") {
        LME_proportion_tsv_filename <- paste0("out/",sub_cluster,"/tables/lme/LME_cell_type_proportions-", batch,".tsv")
      } else {
        LME_proportion_tsv_filename <- paste0("out/tables/lme/LME_cell_type_proportions-", batch,".tsv")
      }

      write.table(lme.table, file=LME_proportion_tsv_filename, sep="\t", col.names = NA)
  
      # INTENSITY TABLE LME
      
      lme.table <- matrix(nrow = 0, ncol=(14))
      
      for (j in 1:length(clusternames)) {
        clustername <- clusternames[[j]]    
        for (i in 1:length(markers_cname)) {
          marker <- markers_cname[i]
          f <- formula(paste0(marker, formula_str_per_batch))
          
          row <- craft_lme_row(f, marker_intensity_tables[[j]], c(clustername=clustername, marker=marker))
          lme.table <- rbind(lme.table, row)
        }
        
      }
      
      lme.table <- data.frame(lme.table)
      lme.table <- add_fdr(lme.table)
      
      if (sub_cluster != "") {
        LME_intensity_tsv_filename <- paste0("out/",sub_cluster,"/tables/lme/celltype_LME_marker_intensities-", batch, ".tsv")
      } else {
        LME_intensity_tsv_filename <- paste0("out/tables/lme/celltype_LME_marker_intensities-",batch,"-", ".tsv")
      }

      write.table(lme.table, file=LME_intensity_tsv_filename, sep="\t", row.names = FALSE, col.names = TRUE)

    }
 
  } # batched set flag 
   
}
