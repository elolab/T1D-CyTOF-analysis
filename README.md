# DIPP CyTOF

These codes present the dataset specific analysis steps after CyTOF workflow based on following tutorial: 
https://www.bioconductor.org/help/course-materials/2017/BioC2017/Day2/Workshops/CyTOF/doc/cytofWorkflow_BioC2017workshop.html

In the CyTOF workflow the set cells are clustered by their cell types defined by lineage marker profile. The subsets CD4_T and CD8_T cells where then subsequently divided into subclusters.

## DIPP CyTOF data-analysis

The outline of the data-analysis steps (presented in DIPP_CyTOF_workflow.R):

1. Calculate intensity and proportion tables (batch is IAA, GADA, MULTIPLE separately):

    function: DIPP_CyTOF_workflow_intensity_tables_batches
    function: DIPP_CyTOF_workflow_intensity_tables

2. Do LME comparisons between Cases and Controls

    formula_str <- " ~ Age + Sex + Age*CaseCtrl + (1|Batch) + (1|Pair) + HLA"

    Note: The significance of the sample subsets (IAA-first, GADA-first, ≥2 Aab first groups) was further assessed using the following formula.
    formula_str <- "~ Age + Sex + HLA + Age*Group + (1|Batch) + (1|Pair)"

    formula_str_per_batch <- " ~ Age + Sex + Age*CaseCtrl + (1|Pair) + HLA"
    
    function: DIPP_CyTOF_workflow_lme

4. Calculate intensity and proportion tables for celltypes

    function: DIPP_CyTOF_workflow_intensity_tables_batches_subclusters
    function: DIPP_CyTOF_workflow_intensity_tables_batches_subclusters

5. Do LME comparisons between Cases and Control for CD4_T and CD8_T

    function: DIPP_CyTOF_workflow_lme

6. Re-sample tSNE for obtaining publication quality tSNE plot

    function: DIPP_CyTOF_workflow_re_tSNE

7. Draw plots for the publication

    function DIPP_CyTOF_workflow_figures

## System Requirements
* The workflow requires R language version 4.x.
* The workflow has been tested on R 4.1 and 4.3.
* R installation guide: https://cran.r-project.org/doc/manuals/r-release/R-admin.html



   
