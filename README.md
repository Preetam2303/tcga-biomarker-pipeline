# tcga-biomarker-pipeline
R pipeline for TCGA transcriptomics biomarker discovery with ML &amp; SHAP
## Purpose of the Pipeline
This pipeline takes raw TCGA transcriptomics data from oral cancer samples, preprocesses and annotates gene expression values, performs differential expression and functional enrichment analyses, integrates multi-omics layers, and then trains machine learning models (e.g., Random Forest, SVM) to identify and validate immune-evasion biomarkers. Post hoc explainability with SHAP highlights the most predictive features for biological interpretation.

## Installation of Dependencies
In an R session, install the required packages:
```r
install.packages("Seurat")
install.packages("biomaRt")
install.packages("limma")
install.packages("clusterProfiler")
install.packages("randomForest")
install.packages("shapforxgboost")
# (Add any other CRAN/Bioconductor packages your pipeline uses)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"))
Rscript run_pipeline.R --input Tumor.csv --output results/
