# NEK2_Survival
Code for the survival analysis of patients with high and low NEK2 and CD274 expression [![DOI](https://zenodo.org/badge/610370369.svg)](https://zenodo.org/badge/latestdoi/610370369)

# CONTENTS OF Code

## Related packages
[survminer](https://github.com/kassambara/survminer) : The survminer R package provides functions for facilitating survival analysis and visualization.

[biomanager](https://github.com/Bioconductor/BiocManager) : The BiocManager package, as the modern successor package to BiocInstaller, allows users to install and manage packages from the Bioconductor project. 
Bioconductor focuses on the statistical analysis and comprehension of high-throughput genomic data.

[Biobase](https://github.com/Bioconductor/Biobase) : Biobase is an R/Bioconductor package that implements base functions for Bioconductor.

[GEOquery](https://github.com/seandavi/GEOquery) : The bridge between the NCBI Gene Expression Omnibus and Bioconductor

[Survival](https://github.com/therneau/survival) : Survival package for R

[tidyverse](https://github.com/tidyverse/tidyverse) : The tidyverse is a set of packages that work in harmony because they share common data representations and API design.

## PROJECT 1: UAMS Survival Data and NEK2 and CD274 (n1136) - (Line 47 - Line 917)
     SECTION 1: Load Packages
     SECTION 2:  Import UAMS data and create Expression Set for Survival Data
       Subsection 2A: Finding Event Free Survival (EFS) Cut-point for NEK2 and CD274 
         Sub-subsection 2AI:  Kaplan-Meier(EFS) for NEK2 (204641_at)
         Sub-subsection 2AII: Kaplan-Meier(EFS) for CD274 (223834_at)
       Subsection 2B: Finding Overall Survival (OS) Cut-point for NEK2 and CD274 
         Sub-subsection 2BI:  Kaplan-Meier(OS) for NEK2 (204641_at)
         Sub-subsection 2BII: Kaplan-Meier(OS) for CD274 (223834_at)
     SECTION 3 - Hazard Ratios for NEK2 and CD274 (EFS)
     SECTION 4 - Hazard Ratios for NEK2 and CD274 (OS)
     SECTION 5 - Import and prepare NEK2_CD274 SUBGROUP combinations for EFS analysis
       Subsection 5A: Hazard Ratio for NEK2Hi_CD274Hi_vs_NEK2Lo_CD274Hi (EFS)
         Subsection 5B: Hazard Ratio for NEK2Hi_CD274Lo_vs_NEK2Lo_CD274Hi (EFS)
         Subsection 5C: Hazard Ratio for NEK2Lo_CD274Lo_vs_NEK2Lo_CD274Hi (EFS)
     SECTION 6 - Import and prepare NEK2_CD274 SUBGROUP combinations for OS analysis
       Subsection 6A: Hazard Ratio for NEK2Hi_CD274Hi_vs_NEK2Lo_CD274Hi (OS)
       Subsection 6B: Hazard Ratio for NEK2Hi_CD274Lo_vs_NEK2Lo_CD274Hi (OS)
       Subsection 6C: Hazard Ratio for NEK2Lo_CD274Lo_vs_NEK2Lo_CD274Hi (OS)

## PROJECT 2: MMRF Survival Data and NEK2 and CD274 (n592) - (Line 921 - 1795)
     SECTION 1: Load Packages
     SECTION 2:  Import MMRF data and create Expression Set for Survival Data
       Subsection 2A: Finding Event Free Survival (EFS) Cut-point for NEK2 and CD274 
         Sub-subsection 2AI:  Kaplan-Meier(EFS) for NEK2 (ENSG00000117650)
         Sub-subsection 2AII: Kaplan-Meier(EFS) for CD274 (ENSG00000120217)
       Subsection 2B: Finding Overall Survival (OS) Cut-point for NEK2 and CD274 
          Sub-subsection 2BI:  Kaplan-Meier(OS) for NEK2 (ENSG00000117650)
          Sub-subsection 2BII: Kaplan-Meier(OS) for CD274 (ENSG00000120217)
     SECTION 3 - Hazard Ratios for NEK2 and CD274 (EFS)
     SECTION 4 - Hazard Ratios for NEK2 and CD274 (OS)
     SECTION 5 - Import and prepare NEK2_CD274 SUBGROUP combinations for EFS analysis
       Subsection 5A: Hazard Ratio for NEK2Hi_CD274Hi_vs_NEK2Lo_CD274Hi (EFS)
       Subsection 5B: Hazard Ratio for NEK2Hi_CD274Lo_vs_NEK2Lo_CD274Hi (EFS)
       Subsection 5C: Hazard Ratio for NEK2Lo_CD274Lo_vs_NEK2Lo_CD274Hi (EFS)
     SECTION 6 - Import and prepare NEK2_CD274 SUBGROUP combinations for OS analysis
       Subsection 6A: Hazard Ratio for NEK2Hi_CD274Hi_vs_NEK2Lo_CD274Hi (OS)
       Subsection 6B: Hazard Ratio for NEK2Hi_CD274Lo_vs_NEK2Lo_CD274Hi (OS)
       Subsection 6C: Hazard Ratio for NEK2Lo_CD274Lo_vs_NEK2Lo_CD274Hi (OS)  
