# UW_Adult_Covid-19

This project repository contains R packages and scripts to generate results for the paper entitled
"The landscape of antibody binding in SARS-CoV-2 infection"

# Processing peptide arrays and MUSCLE alignments

Please see https://dholk.primate.wisc.edu/project/dho/sequencing/Polypeptide%20Microarrays/public/COVID_19/begin.view for data
and instructions for preprocesing the peptide binding array data and MUSCLE (MUltiple Sequence Comparison by Log-Expectation) alignments for comparing each coronavirus to the Wu1-SARS2 sequence.  To run the find_wu1_epitopes.Rmd script, you don't need to be able to run pepMeld, you can download the "aggregated_data/df_stacked.tsv.gz" and "all_sequences_except_wi.tsv.gz" from the above website and run the find_wu1_epitopes.Rmd mentioned below.


# Identifying Epitopes (Figure 7)

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

