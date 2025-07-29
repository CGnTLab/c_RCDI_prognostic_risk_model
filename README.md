# c-RCDI Prognostic Signature

## Research article details:

Title: Development of an integrative machine learning-based combined prognostic gene signature based on non-apoptotic regulated cell death pathways to improve patient survival in colorectal cancer

Authors: Avik Sengupta, Dr. Rahul Kumar

Affiliations: Department of Biotechnology, Indian Institute of Technology Hyderabad, Near NH-65, Vill. Kandi, Sangareddy, Telangana-502285, India

## Overview

We identified and developed an integrative machine learning based prognostic gene signature, derived from 13 combined non-apoptotic regulated cell death pathway genes, in colorectal cancer (CRC).

## Machine learning based gene selection of each of the 13 differentially expressed NARCD pathway genes

Gene selection was performed using 46 machine learning and generalized linear survival model combinations. The survival models used to generate model combinations are as follows: I) Univariate cox regression, II) LASSO regression, III) Ridge regression, IV) random forest survival (RSF), V) plsRcox, VI) CoxBoost, VII) StepCox, VIII) SuperPC, and IX) Gradient Boost Machine (GBM).

These model combinations were applied to each of the 13 DE-NARCD pathway genes to find the genes selected for optimal survival predictive power (C-index). Hence, risk scores for each pathway were calculated, derived from the best model combinations of each pathway.

## Identifying the best combined model of the NARCD pathways / development of c-RCDI

Logistic regression was applied to find the best 2-4 pathway combination by calculating the C-indices for each pathway(s) combination. The combined risk score was calculated further to get the combined risk score, which we named as the combined regulated cell death index (c-RCDI), formed from a signature of 43 genes.

## Predictive capability of the c-RCDI

- c-RCDI was used to distinguish patients into high- and low-risk groups, based on the optimal cutoff from the "survival" package in R.
- Kaplan-Meier (KM) and area under receiver operating characteristics (AUROC) analyses were performed for the predictive capability of the c-RCDI

## Hierarchical clustering of the patients based on the gene signature

Based on the gene expression of the 43 signature genes, we performed Hierarchical clustering of the patient samples to see the distinction based on the 43 genes. Also ,we extrapolated the clusters with the risk groups to compare and correlate them.

## 43 genes evaluation and validation of the 43 signature genes

We performed the following to evaluate and validate the 43 signature genes:
  
  - Multivariate Cox regression analysis of the signature genes
  - Protein-protein interaction (PPI) network analysis using "STRING database" and "Cytoscape"
  - Single-cell gene expression analysis using "Seurat" package in R
  - Gene dependency analysis using DepMap CRISPR/Cas9 knockout (KO) data
  - Gene expression correlation among each other
  - Weighted gene co-expression network analysis using "WGCNA" package in R to identify gene modules and module-module relationships

## c-RCDI as an independent prognostic clinical factor

To evaluate the prognostic potential of the c-RCDI as an independent clinical factor, we performed the following analysis:

  - Multivariate Cox regression analysis of c-RCDI and other clinical factors, age group, sex, and tumor stage
  - AUROC analysis of the clinical factors
  - Nomogram analysis

## Functional analysis of high- and low-risk groups

### Immune cell enrichment analysis

Using CIBERSORT algorithm from CIBSERSORTx website, we calculated the absolute immune enrichment scores of the 22 immune infiltrating cells. 

### Tumor immunotherapy response analysis

Using the TIDE algorithm from the TIDE website, we identified the responders vs non-responders patients among the risk groups, to further evaluate the differences in immunotherapy response.

### Drug sensitivity analysis

Using GDSC2 data and "oncoPredict()" tool in R, we calculated the drugs' IC50 values to determine their effectiveness against each patient of the high- and low-risk groups.

## Machine learning analysis

To evaluate our c-RCDI-based distinction of the patients, we further applied six classification machine learning models to the 43 genes' gene expression and c-RCDI to predict the classes.


Copyright | Computational Genomics & Transcriptomics (CG&T) Laboratory, Department of Biotechnology, Indian Institute of Technology Hyderabad, Telangana-502285, India
