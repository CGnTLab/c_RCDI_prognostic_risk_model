# c-RCDI Development

## Research article details:

Title: Development of an integrative machine learning-based combined prognostic gene signature based on non-apoptotic regulated cell death pathways to improve patient survival in colorectal cancer

Authors: Avik Sengupta, Dr. Rahul Kumar

Affiliations: Department of Biotechnology, Indian Institute of Technology Hyderabad, Near NH-65, Vill. Kandi, Sangareddy, Telangana-502285, India

## Overview

We identified and developed an integrative machine learning based prognostic gene signature, derived from 13 combined non-apoptotic regulated cell death pathway genes, in colorectal cancer (CRC).

## Machine learning based gene selection of each of the 13 differentially expressed NARCD pathway genes

Gene selection was performed using 46 machine learning and generalized linear survival model combinations. The survival models used to generate model combinations are as follows: I) Univariate cox regression, II) LASSO regression, III) Ridge regression, IV) random forest survival (RSF), V) plsRcox, VI) CoxBoost, VII) StepCox, VIII) SuperPC, and IX) Gradient Boost Machine (GBM).

These model combinations were applied to each of the 13 DE-NARCD pathway genes to find the genes selected for optimal survival predictive power (C-index). Hence, risk scores for each pathway were calculated, derived from the best model combinations of each pathway.

## Identifying the best combined model of the NARCD pathways + development of c-RCDI

Logistic regression was applied to find the best 2-4 pathway combination by calculating the C-indices for each pathway(s) combination. The combined risk score was calculated further to get the combined risk score, which we named as the combined regulated cell death index (c-RCDI), formed from a signature of 43 genes.

## Predictive capability of the c-RCDI

- c-RCDI was used to distinguish patients into high- and low-risk groups, based on the optimal cutoff from the "survival" package in R.
- Kaplan-Meier (KM) and area under receiver operating characteristics (AUROC) analyses were performed for the predictive capability of the c-RCDI

## Hierarchical clustering of the patients based on the gene signature
