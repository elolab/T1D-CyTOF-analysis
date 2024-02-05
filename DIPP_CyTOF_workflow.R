.libPaths("Rlibs/R/library/4.3")

print ("THE START")  

source("DIPP_CyTOF_intensity_tables_batches.R")
source("DIPP_CyTOF_intensity_tables_batches_subclusters.R")
source("DIPP_CyTOF_intensity_tables.R")
source("DIPP_CyTOF_intensity_tables_subclusters.R")

source("DIPP_CyTOF_lme.R")
#source("DIPP_CyTOF_lme_comparisons.R")
source("DIPP_CyTOF_lme_comparisons_group.R")
source("DIPP_CyTOF_lme_models.R")

source("DIPP_CyTOF_re_tSNE.R")
source("DIPP_CyTOF_publication_materials.R")

# Prerequisite: CyTOF workflow based on tutorial: 
# https://www.bioconductor.org/help/course-materials/2017/BioC2017/Day2/Workshops/CyTOF/doc/cytofWorkflow_BioC2017workshop.html

project_folder = "./"

print ("CREATING PROPORTION AND INTENSITY TABLES")

DIPP_CyTOF_workflow_intensity_tables_batches(project_folder,
                                                     "cfg/sample_metadata.tsv")

DIPP_CyTOF_workflow_intensity_tables(project_folder,
                                                     "cfg/sample_metadata.tsv")


print ("RUNNING LME")

formula_str <- "~ Age + Sex + HLA + Group*Age + (1|Batch) + (1|Pair)"
formula_str_per_batch <- "" 
 
DIPP_CyTOF_workflow_lme(project_folder,
                                "cfg/sample_metadata.tsv",
                                "cfg/markers.tsv",
                                "",
                                formula_str,
                                formula_str_per_batch,
                                TRUE,
                                FALSE)
 
#cmp_formula_str <- " ~ Age + Sex + Age*CaseCtrl + (1|Batch) + (1|Pair) + HLA"
#cmp_formula_str_per_batch <- " ~ Age + Sex + Age*CaseCtrl + (1|Pair) + HLA"

#DIPP_CyTOF_workflow_lme_comparisons(project_folder,
#                                            "cfg/sample_metadata.tsv",
#                                            "cfg/markers.tsv",
#                                            cmp_formula_str,
#                                            paste0(cmp_formula_str, " + HLA"),
#                                            cmp_formula_str_per_batch,
#                                            paste0(cmp_formula_str_per_batch, " + HLA"))


print ("RUNNING LME - OUTPUTTING MODEL RDATA FILES")

DIPP_CyTOF_workflow_lme_models(project_folder,
                                      "cfg/sample_metadata.tsv",
                                      "cfg/markers.tsv")

print ("CREATING PROPORTION AND INTENSITY TABLES FOR SUBCLUSTERS")
DIPP_CyTOF_workflow_intensity_tables_batches_subclusters(project_folder,
                                                    "cfg/sample_metadata.tsv")

print ("CREATING PROPORTION AND INTENSITY TABLES FOR SUBCLUSTER TYPES")

DIPP_CyTOF_workflow_intensity_tables_batches_subclusters(project_folder,
                                                                  "cfg/sample_metadata.tsv",
                                                                  c("CD4_T_cells", "CD8_T_cells"))


cmp_formula_str1 <- "~ Age + Sex + HLA + CaseCtrl*Age + (1|Batch) + (1|Pair)"
cmp_formula_str2 <- "~ Age + Sex + HLA + Group*Age + (1|Batch) + (1|Pair)"

DIPP_CyTOF_workflow_lme_comparisons_group(project_folder,
                                           "cfg/sample_metadata.tsv",
                                           "cfg/markers.tsv",
                                           cmp_formula_str1,
                                           cmp_formula_str2)

print ("LME TESTING FOR SUBCLUSTER: CD4_T_cells")
DIPP_CyTOF_workflow_lme(project_folder,
                                "cfg/sample_metadata.tsv",
                                "cfg/markers.tsv",
                                "CD4_T_cells",
                                formula_str,
                                formula_str_per_batch,
                                FALSE,
                                TRUE)

print ("LME TESTING FOR SUBCLUSTER: CD8_T_cells")
DIPP_CyTOF_workflow_lme(project_folder,
                                "cfg/sample_metadata.tsv",
                                "cfg/markers.tsv",
                                "CD8_T_cells",
                                formula_str,
                                formula_str_per_batch,
                                FALSE,
                                TRUE)



print ("Re-sampling data for tSNE (change sampling depth)")
DIPP_CyTOF_workflow_re_tSNE(project_folder,
                                    "cfg/markers.tsv")
  
print ("RUNNING PUBLICATION FIGURES")
DIPP_CyTOF_workflow_figures(project_folder,
                                    "cfg/sample_metadata.tsv",
                                    "cfg/markers.tsv")

print ("THE END")  
