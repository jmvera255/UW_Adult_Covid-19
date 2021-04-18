# UW_Adult_Covid-19

This project repository contains R packages and scripts to generate results for the paper entitled
"The landscape of antibody binding in SARS-CoV-2 infection"

# Processing peptide arrays and muscle alignments

Please see https://dholk.primate.wisc.edu/project/dho/sequencing/Polypeptide%20Microarrays/public/COVID_19/begin.view for data
and instructions for preprocesing the peptide binding array data and muscle alignments.


# Finding Epitopes (For Figure 7)

Requirements:
- R or RStudio
- rmarkdown package (install.packages("rmarkdown")
- UW.Adult.Covid.19 - install.packages("devtools");devtools::install("UW.Adult.Covid.19")
- Matrix - install.packages("Matrix");
- data.table - install.packages("data.table")
- multtest - https://www.bioconductor.org/packages/release/bioc/html/multtest.html
- preprocessCore - http://bioconductor.org/packages/release/bioc/html/preprocessCore.html

See find_wu1_epitopes.Rmd for other requirements.


# Probe Heatmaps

<Insert Instructions here for heatmap>

