# DIPP CyTOF

These codes present the dataset specific analysis steps after CyTOF workflow based on following tutorial: 
https://www.bioconductor.org/help/course-materials/2017/BioC2017/Day2/Workshops/CyTOF/doc/cytofWorkflow_BioC2017workshop.html

In the CyTOF workflow the set cells are clustered by their cell types defined by lineage marker profile. The subsets CD4_T and CD8_T cells where then subsequently divided into subclusters.

## DIPP CyTOF data-analysis

The outline of the data-analysis steps (presented in B20012_DIPP_CyTOF_workflow.R):

1. Calculate intensity and proportion tables (batch is IAA, GADA, MULTIPLE separately):

    function: B200012_DIPP_CyTOF_workflow_intensity_tables_batches

    function: B200012_DIPP_CyTOF_workflow_intensity_tables

2. Do LME comparisons between Cases and Controls

    formula_str <- " ~ Age + Sex + Age*CaseCtrl + (1|Batch) + (1|Pair) + HLA"

    formula_str_per_batch <- " ~ Age + Sex + Age*CaseCtrl + (1|Pair) + HLA"
    
    function: B200012_DIPP_CyTOF_workflow_lme

3. Calculate intensity and proportion tables for celltypes

    function: B200012_DIPP_CyTOF_workflow_intensity_tables_batches_subclusters

4. Do LME comparisons between Cases and Control for CD4_T and CD8_T

    function: B200012_DIPP_CyTOF_workflow_lme

5. Re-sample tSNE for obtaining publication quality tSNE plot

    function: B200012_DIPP_CyTOF_workflow_re_tSNE

6. Draw plots for the publication

    function B200012_DIPP_CyTOF_workflow_figures
