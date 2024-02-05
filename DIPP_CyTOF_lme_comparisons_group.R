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


# craft_lme_row <- function(f, pr_int_table, prefix_list) {
#   
#   if (FALSE) {  
#     pr_int_table <- marker_intensity_tables[[j]]
#     prefix_list <- c(clustername=clustername, marker=marker)
#     
    pr_int_table <- combined_proportion_table
    prefix_list <- c()
#   }
#   
#   suppressWarnings(suppressMessages(fit1 <- lmer(f, data = pr_int_table, na.action = na.omit, REML=FALSE)))
#   
#   out.coef <- coef(summary(fit1))[ ,"Estimate"]; names(out.coef) <- paste("coef",names(out.coef),sep=".")
#   out.p <-coef(summary(fit1))[ ,"Pr(>|t|)"]; names(out.p) <- paste("p",names(out.p),sep=".")
#   
#   row <- c(prefix_list, out.coef, out.p)
#   return (row)
# }

do_lme_row <- function(f1, f2, pr_int_table, prefix_list) {

  if (FALSE) {  
    pr_int_table <- marker_intensity_tables[[j]]
    prefix_list <- c(clustername=clustername, marker=marker)
    pr_int_table <- combined_proportion_table
    prefix_list <- c()
  }
  
  suppressWarnings(suppressMessages(fit1 <- lmer(f1, data = pr_int_table, na.action = na.omit, REML=FALSE)))
  suppressWarnings(suppressMessages(fit2 <- lmer(f2, data = pr_int_table, na.action = na.omit, REML=FALSE)))
  
  aov <- suppressWarnings(suppressMessages(anova(fit1,fit2)))
  if (aov$`Pr(>Chisq)`[2]<0.05) {
    if(aov$AIC[2]<aov$AIC[1]) {
      better <- 2
    } else if (aov$AIC[1]<aov$AIC[2]) {
      better <- 1
    } else {
    better <- 0
    }
  } else {
    better <- 0
  }

  f1_out.coef <- coef(summary(fit1))[,"Estimate"]; names(f1_out.coef) <- paste("coef.f1",names(f1_out.coef),sep=".")
  f1_out.p <-coef(summary(fit1))[,"Pr(>|t|)"]; names(f1_out.p) <- paste("p.f1",names(f1_out.p),sep=".")
  
  f2_out.coef <- coef(summary(fit2))[,"Estimate"]; names(f2_out.coef) <- paste("coef.f2",names(f2_out.coef),sep=".")
  f2_out.p <-coef(summary(fit2))[,"Pr(>|t|)"]; names(f2_out.p) <- paste("p.f2",names(f2_out.p),sep=".")
  
  row <- c(prefix_list, better_model=better, AIC_f1=aov$AIC[1], AIC_f2=aov$AIC[2], Pr_Chisq=aov$`Pr(>Chisq)`[2], f1_out.coef, f1_out.p, f2_out.coef, f2_out.p)
  return (row)
}

add_fdr <- function(lme.comparison) {
  fdr <- apply(lme.comparison[,grep("p.",colnames(lme.comparison))], 2, function(x) p.adjust(x, method="BH"))
  colnames(fdr) <- gsub("p\\.", "fdr\\.", colnames(fdr))
  lme.comparison <- data.frame(lme.comparison, fdr, stringsAsFactors=FALSE)
  return (lme.comparison)
}

DIPP_CyTOF_workflow_lme_comparisons_group <- function(CyTOF_workflow_cwd,
                                                        sample_metadata_tsv,
                                                        markers_tsv,
                                                        formula_str_f1,
                                                        formula_str_f2) {
  
  setwd(CyTOF_workflow_cwd)
  
  dir.create("out/")
  dir.create("out/tables/")
  dir.create("out/tables/lmecomp/")
  
  fileConn<-file(paste0("out/tables/lmecomp/formula.txt"))
  writeLines(c(formula_str_f1, formula_str_f2), fileConn)

  close(fileConn)

  anno <- read.csv(sample_metadata_tsv, stringsAsFactors=FALSE, sep="\t")
  rownames(anno) <- anno[["Sample"]]
  
  load("RData/panel.RData")
  panel <- read.table(markers_tsv, sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE)
  
  markers_panel <- panel[(panel["type"] == "functional") | (panel["type"] == "lineage/functional"), ]
  markers_cname <- markers_panel[["cname"]]
  markers <- rownames(markers_panel)

  # WHOLE SET
  
  proportion_table_rdata_filename = paste0("RData/proportion_table", ".RData")
  load(proportion_table_rdata_filename)
  
  clusternames <- colnames(proportion_table)[!(colnames(proportion_table) %in% c(colnames(anno), "removed"))]
  
  intensity_table_rdata_filename <- paste0("RData/marker_intensity_tables", ".RData")
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
  
  lme.table <- NULL
  for (i in 1:length(clusternames)) {
    celltype <- clusternames[[i]]
    f1 <- formula(paste0(celltype, formula_str_f1))
    f2 <- formula(paste0(celltype, formula_str_f2))

    row <- do_lme_row(f1, f2, combined_proportion_table, c())
    
    if (is.null(lme.table)) {
      lme.table <- row
    } else {
      lme.table <- rbind(lme.table, row)
    }
    
  }

  lme.table <- data.frame(lme.table)
  rownames(lme.table) <- clusternames

  lme.table <- add_fdr(lme.table)
  
  write.table(lme.table, file=paste0("out/tables/lmecomp/LME_cell_type_proportions_comparison", ".tsv"), sep="\t", col.names = NA)
  
  # INTENSITY TABLE LME
  
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
      f1 <- formula(paste0(marker, formula_str_f1))
      f2 <- formula(paste0(marker, formula_str_f2))

      row <- do_lme_row(f1, f2, combined_intensity_table, c(clustername=clustername, marker=marker))
      
      if (is.null(lme.table)) {
        lme.table <- row
      } else {
        lme.table <- rbind(lme.table, row)
      }

    }

  }
  
  lme.table <- data.frame(lme.table)
  #rownames(lme.comparison) <- clusternames

  lme.table <- add_fdr(lme.table)

  write.table(lme.table, file=paste0("out/tables/lmecomp/LME_cell_type_marker_intensities_comparison",".tsv"), sep="\t", row.names = FALSE, col.names = TRUE)

}
