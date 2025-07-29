setwd("/home2/New_Objective2_Biomarker_research/batch_corrected_rna_seq/new/")

rm(list=ls())
# install.packages("Rtsne")
library(Rtsne)
library(ggplot2)
library(ggpubr)

######### t-SNE implementation ################
# install.packages("plotly")
# install.packages("tsne")
library(tsne)
library(plotly)
gc()

library(mclust)   # For adjustedRandIndex
library(fossil)    # For rand.index
library(caret)     # For confusion matrix

######### CONCORDANCE ANALYSIS FUNCTION #########
calculate_concordance <- function(clustering_labels, risk_groups, dataset_name) {
  # Convert to factors for consistent handling
  clustering_factor <- as.factor(clustering_labels)
  risk_factor <- as.factor(risk_groups)
  
  # 1. Contingency Table
  contingency_table <- table(Risk_Group = risk_factor, Cluster = clustering_factor)
  cat("\nContingency Table for", dataset_name, ":\n")
  print(contingency_table)
  
  # 2. Adjusted Rand Index (ARI)
  ari_value <- adjustedRandIndex(clustering_factor, risk_factor)
  
  # 3. Standard Rand Index
  ri_value <- rand.index(as.numeric(clustering_factor), as.numeric(risk_factor))
  
  # 4. Overall Concordance Percentage
  n_total <- length(risk_factor)
  mapping1 <- sum(diag(contingency_table))  # Direct mapping
  mapping2 <- sum(diag(contingency_table[, c(2,1)]))  # Flipped mapping
  best_concordance <- max(mapping1, mapping2)
  concordance_percentage <- (best_concordance / n_total) * 100
  
  # 5. Cohen's Kappa
  kappa <- confusionMatrix(contingency_table)$overall['Kappa']
  
  # Results summary
  cat("\nConcordance Metrics for", dataset_name, ":\n")
  cat(sprintf("Adjusted Rand Index (ARI): %.4f\n", ari_value))
  cat(sprintf("Rand Index: %.4f\n", ri_value))
  cat(sprintf("Cohen's Kappa: %.4f\n", kappa))
  cat(sprintf("Concordance Percentage: %.2f%%\n", concordance_percentage))
  
  # Interpretation
  if(ari_value >= 0.70) {
    cat("Interpretation: EXCELLENT concordance (ARI >= 0.70)\n")
  } else if(ari_value >= 0.50) {
    cat("Interpretation: GOOD concordance (ARI >= 0.50)\n") 
  } else {
    cat("Interpretation: MODERATE/POOR concordance (ARI < 0.50)\n")
  }
  
  return(list(
    ari = ari_value,
    rand_index = ri_value,
    concordance_percentage = concordance_percentage,
    kappa = kappa,
    contingency_table = contingency_table
  ))
}

#####################################################
#### Unsupervised Clustering PCA and t-SNE plots ####
#####################################################
seed = 42
set.seed(seed)
library(ConsensusClusterPlus)
library(uwot)
library(cluster) # For silhouette
# signature genes
gs_ferr <- c("CDKN2A","CHMP6","TXNIP","CXCL2","TP63","BID","CDO1","PLIN4","SLC2A3","EGFR","ALOXE3","AURKA","GDF15","GPX2","RGS4")
gs_py <- c("CHMP6","TP63","GSDMC","IL18","GZMB")
gs_aut <- c("MAP1B","CLN3","ATP2B4","SUN2")
gs_net <- c("MYC","CDK4","LAMC2","ETS2","TCF7","NOX4","STAT1","YAP1","C1QC","C1QA","KNG1","SNAI1","IGHG1","CGAS","S100A11","CR1","ACTA2","LAMA2","CDK6","NFATC1","TRAF3IP2")
gs_RCDI <- unique(c(gs_ferr, gs_py, gs_aut, gs_net))

## input data creation (from z-scores) ##
library(data.table)
library(ConsensusClusterPlus)
library(Rtsne)
library(uwot)
library(ggplot2)
library(cluster)   # silhouette()

filenames <- c(
  "Train_RCDI_categorized.csv",
  "Test_RCDI_categorized.csv",
  "GSE39582_RCDI_categorized.csv",
  "GSE161158_RCDI_categorized.csv"
)

# Load files
file_list <- lapply(filenames, function(file) {
  data <- data.frame(fread(file, header = TRUE, sep = ","), check.names = FALSE, row.names = 1)
  return(data)
})
names(file_list) <- c("Train", "Test", "GSE39582", "GSE161158")

# Store clustering results
cluster_assignments_list <- list()

# Loop over datasets
lapply(names(file_list), function(name) {
  message("Running clustering for: ", name)
  
  df <- file_list[[name]]
  mat <- df[, colnames(df) %in% gs_RCDI]
  mat <- t(scale(t(mat), scale = TRUE, center = TRUE))
  mat <- t(mat)
  mat <- as.matrix(mat)
  
  # Consensus clustering
  results <- ConsensusClusterPlus(
    mat, maxK = 6, reps = 1000, pItem = 0.8, pFeature = 0.8,
    clusterAlg = "hc", distance = "pearson", seed = 42, plot = "png",
    title = paste0("clustering_", name)
  )
  calcICL(results, plot = "png", title = paste0("icl_", name), writeTable = TRUE)
  cluster_assignments <- results[[2]]$consensusClass
  cluster_assignments_list[[name]] <<- cluster_assignments
  
  #### UMAP ####
  umap_result <- umap(t(mat), n_neighbors = 10, min_dist = 0.01, metric = "euclidean")
  umap_data <- data.frame(UMAP1 = umap_result[,1], UMAP2 = umap_result[,2], Cluster = as.factor(cluster_assignments))
  umap_data$Sample <- colnames(mat)
  
  #### UMAP Plot (NO LABELS) ####
  umap_plot <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
    geom_point(size = 3) +
    xlab("UMAP1") + ylab("UMAP2") +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14, face = "bold"),
          legend.position = "right")
  
  #### Save plots ####
  tiff(paste0("PCA_", name, ".tif"), res = 600, width = 5, height = 5, units = "in", compression = "lzw")
  print(pca_plot)
  dev.off()
  
  tiff(paste0("TSNE_", name, ".tif"), res = 600, width = 5, height = 5, units = "in", compression = "lzw")
  print(tsne_plot)
  dev.off()
  
  tiff(paste0("UMAP_", name, ".tif"), res = 600, width = 5, height = 5, units = "in", compression = "lzw")
  print(umap_plot)
  dev.off()
  
  #### Save cluster coordinates ####
  fwrite(umap_data, paste0("umap_clusters_", name, ".csv"), sep = ",", row.names = FALSE)
})

## Risk class vs clusters -- alluvial plots

library(ggalluvial)
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)

lapply(names(merged_list), function(name) {
  m <- merged_list[[name]]
  
  features <- c("RiskClass", "Cluster")
  m <- m[, names(m) %in% features]
  m <- as.data.frame(m)
  
  df <- m %>%
    group_by(RiskClass, Cluster) %>%
    summarise(Freq = n(), .groups = "drop")
  
  tiff(paste0("alluvial_plot_risk_cluster_", name, ".tif"), 
       res = 600, height = 4, width = 4, units = "in", compression = "lzw")
  
  gg <- ggplot(df, aes(axis1 = RiskClass, axis2 = Cluster, y = Freq)) +
    geom_alluvium(aes(fill = RiskClass), curve_type = "sigmoid", alpha = 0.8, size = 0.4) +
    geom_stratum(fill = "white", color = "black", size = 0.4) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.5, family = "Helvetica") +
    scale_fill_manual(values = c("high" = "salmon", "low" = "skyblue")) +
    scale_x_discrete(limits = c("RiskClass", "Cluster"), expand = c(0.1, 0.05)) +
    theme_classic(base_family = "Helvetica") +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 10),
      axis.title = element_blank(),
      legend.position = "top",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10)
    ) +
    labs(title = paste("Alluvial Plot:", name), y = "Frequency", fill = "RiskClass")
  
  print(gg)
  dev.off()
})
