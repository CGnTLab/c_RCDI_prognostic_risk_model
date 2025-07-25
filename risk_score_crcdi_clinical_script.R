#!/usr/bin/Rscript
setwd("/home2/New_Objective2_Biomarker_research/batch_corrected_rna_seq/new/")
gc()
rm(list=ls())
set.seed(42)

library(dplyr)
library(data.table)
library(readxl)
library(BiocParallel)
library(parallel)
num_detected_cores <- detectCores(logical = FALSE) - 2 # logical = FALSE gives physical cores
register(MulticoreParam(workers = num_detected_cores))

#### Dataset preprocessing ####
merged <- data.frame(fread("combined_z_score_clinical_batch.csv", sep=','), row.names=1, check.names=F)

train  <- merged %>% dplyr::sample_frac(0.8) ##dim(train) #743 14558
test1   <- merged[!(rownames(merged) %in% rownames(train)), ]

input_cibersort_train <- train[, -c((ncol(train)-1):ncol(train))]
input_cibersort_train <- t(input_cibersort_train)
fwrite(input_cibersort_train, "input_cibersort_train.tsv", sep='\t', row.names = T)

input_cibersort_test <- test1[, -c((ncol(test1)-1):ncol(test1))]
input_cibersort_test <- t(input_cibersort_test)
fwrite(input_cibersort_test, "input_cibersort_test.tsv", sep='\t', row.names = T)

degs_emtab <- data.frame(fread("DEGs_emtab_0.5.txt", sep = "\t", header = T), row.names=1, check.names = F)
degs_emtab <- row.names(degs_emtab)
degs_tcga <- data.frame(fread("DEGs_TCGA_0.5.txt", sep = "\t", header = T), row.names=1, check.names = F)
degs_tcga <- row.names(degs_tcga)
common_degs <- intersect(degs_emtab, degs_tcga)
genes_to_select <- common_degs

test2 <- data.frame(fread("independent_ds_gse39582.csv", sep = ",", header = T), row.names=1, check.names=F)
test2 <- test2 %>%
  dplyr::rename(
    OS_MONTHS = os.months,
    OS_STATUS = os.event
  )
test3 <- data.frame(fread("independent_ds_gse161158.csv", sep=",", header=T), row.names = 1, check.names = F)
test3 <- test3 %>%
  dplyr::rename(
    OS_MONTHS = dfs_time,
    OS_STATUS = dfs_event
  )
test4 <- data.frame(fread("tcga_brca_dataset.csv", sep=",", header=T), row.names = 1, check.names = F)
test4 <- test4 %>%
  dplyr::rename(
    OS_MONTHS = OS.time,
    OS_STATUS = OS
  )
test5 <- data.frame(fread("tcga_luad_dataset.csv", sep=",", header=T), row.names = 1, check.names = F)
test5 <- test5 %>%
  dplyr::rename(
    OS_MONTHS = OS.time,
    OS_STATUS = OS
  ) 
test6 <- data.frame(fread("tcga_skcm_dataset.csv", sep=",", header=T), row.names = 1, check.names = F)
test6 <- test6 %>%
  dplyr::rename(
    OS_MONTHS = OS.time,
    OS_STATUS = OS
  ) 
test7 <- data.frame(fread("tcga_gbm_dataset.csv", sep=",", header=T), row.names = 1, check.names = F)
test7 <- test7 %>%
  dplyr::rename(
    OS_MONTHS = OS.time,
    OS_STATUS = OS
  ) 
test8 <- data.frame(fread("tcga_paad_dataset.csv", sep=",", header=T), row.names = 1, check.names = F)
test8 <- test8 %>%
  dplyr::rename(
    OS_MONTHS = OS.time,
    OS_STATUS = OS
  )

data_list <- list(train, test1, test2, test3, test4, test5, test6, test7, test8)

# Find common column names across all data frames
common_cols <- Reduce(intersect, lapply(data_list, colnames))

train <- train[, common_cols, drop = FALSE]
test1 <- test1[, common_cols, drop = FALSE]
test2 <- test2[, common_cols, drop = FALSE]
test3 <- test3[, common_cols, drop = FALSE]
test4 <- test4[, common_cols, drop = FALSE]
test5 <- test5[, common_cols, drop = FALSE]
test6 <- test6[, common_cols, drop = FALSE]
test7 <- test7[, common_cols, drop = FALSE]
test8 <- test8[, common_cols, drop = FALSE]

# subset for only DEGs for all datasets
train <- train[, names(train) %in% c(genes_to_select, "OS_MONTHS", "OS_STATUS")]
test1 <- test1[, names(test1) %in% c(genes_to_select, "OS_MONTHS", "OS_STATUS")]
test2 <- test2[, names(test2) %in% c(genes_to_select, "OS_MONTHS", "OS_STATUS")]
test3 <- test3[, names(test3) %in% c(genes_to_select, "OS_MONTHS", "OS_STATUS")]
test4 <- test4[, names(test4) %in% c(genes_to_select, "OS_MONTHS", "OS_STATUS")]
test5 <- test5[, names(test5) %in% c(genes_to_select, "OS_MONTHS", "OS_STATUS")]
test6 <- test6[, names(test6) %in% c(genes_to_select, "OS_MONTHS", "OS_STATUS")]
test7 <- test7[, names(test7) %in% c(genes_to_select, "OS_MONTHS", "OS_STATUS")]
test8 <- test8[, names(test8) %in% c(genes_to_select, "OS_MONTHS", "OS_STATUS")]


## gene sets of signature genes ##

gs_ferr <- c("CDKN2A","CHMP6","TXNIP","CXCL2","TP63","BID","CDO1","PLIN4","SLC2A3","EGFR","ALOXE3","AURKA","GDF15","GPX2","RGS4")
gs_py <- c("CHMP6","TP63","GSDMC","IL18","GZMB")
gs_aut <- c("MAP1B","CLN3","ATP2B4","SUN2")
gs_net <- c("MYC","CDK4","LAMC2","ETS2","TCF7","NOX4","STAT1","YAP1","C1QC","C1QA","KNG1","SNAI1","IGHG1","CGAS","S100A11","CR1","ACTA2","LAMA2","CDK6","NFATC1","TRAF3IP2")
gs_RCDI <- unique(c(gs_ferr, gs_py, gs_aut, gs_net))

train <- train[, names(train) %in% c(gs_RCDI, "OS_MONTHS", "OS_STATUS")]
test1 <- test1[, names(test1) %in% c(gs_RCDI, "OS_MONTHS", "OS_STATUS")]
test2 <- test2[, names(test2) %in% c(gs_RCDI, "OS_MONTHS", "OS_STATUS")]
test3 <- test3[, names(test3) %in% c(gs_RCDI, "OS_MONTHS", "OS_STATUS")]
test4 <- test4[, names(test4) %in% c(gs_RCDI, "OS_MONTHS", "OS_STATUS")]
test5 <- test5[, names(test5) %in% c(gs_RCDI, "OS_MONTHS", "OS_STATUS")]
test6 <- test6[, names(test6) %in% c(gs_RCDI, "OS_MONTHS", "OS_STATUS")]
test7 <- test7[, names(test7) %in% c(gs_RCDI, "OS_MONTHS", "OS_STATUS")]
test8 <- test8[, names(test8) %in% c(gs_RCDI, "OS_MONTHS", "OS_STATUS")]

library("survival")
library("survminer")

if (!requireNamespace("survivalROC", quietly = TRUE)) install.packages("survivalROC")
if (!requireNamespace("PRROC", quietly = TRUE)) install.packages("PRROC")

# Load required packages
library(survivalROC)   # For time-dependent ROC analysis
library(PRROC)   # For Precision-Recall (PR) AUC
library(timeROC)

covariates<- colnames(train)
covariates <- head(covariates, -2)
length(covariates)

## Make dataframe list--
df_list <- list(
  Train= train,
  Test= test1,
  GSE39582 = test2,
  GSE161158 = test3,
  TCGA_BRCA = test4,
  TCGA_LUAD = test5,
  TCGA_SKCM = test6,
  TCGA_GBM = test7,
  TCGA_PAAD = test8
)

df_list <- bplapply(df_list, function(x) {
  class(x$OS_STATUS) <- "numeric"
  class(x$OS_MONTHS) <- "numeric"
  x <- na.omit(x)
  x <- x[!(x$OS_MONTHS < 1),]
}, BPPARAM = MulticoreParam(num_detected_cores))


################################ RCDI Risk Score grouping - All datasets ###############################

filenames <- c(
  "train_RCDI_combined_risk_scores_2.csv",
  "test_RCDI_combined_risk_scores_2.csv",
  "gse39582_RCDI_combined_risk_scores_2.csv",
  "gse161158_RCDI_combined_risk_scores_2.csv",
  "tcga_brca_RCDI_combined_risk_scores_2.csv",
  "tcga_luad_RCDI_combined_risk_scores_2.csv",
  "tcga_skcm_RCDI_combined_risk_scores_2.csv",
  "tcga_gbm_RCDI_combined_risk_scores_2.csv",
  "tcga_paad_RCDI_combined_risk_scores_2.csv"
)
file_list <- lapply(filenames, function(file) {
  if (file.exists(file)) {
    data <- read.csv(file, row.names = 1, check.names = FALSE)  # Read CSV file
    return(data)
  } else {
    warning(paste("File not found:", file))
    return(NULL)
  }
})
names(file_list) <- c("Train", "Test", "GSE39582", "GSE161158", "TCGA_BRCA", "TCGA_LUAD", "TCGA_SKCM", "TCGA_GBM", "TCGA_PAAD")
print(file_list)

## Merged list with RS_scores and gene expression values ##
merged_list <- Map(function(file_data, df_data) {
  merged <- merge(file_data, df_data, by = "row.names", all = TRUE)
  rownames(merged) <- merged$Row.names
  merged <- merged[, -which(colnames(merged) == "Row.names")]
  merged <- merged[, !duplicated(colnames(merged))]
  return(merged)
}, file_list, df_list)
names(merged_list) <- names(file_list)
merged_list <- lapply(merged_list, function(df) {
  df <- df[, -which(colnames(df) %in% c("OS_STATUS.y", "OS_MONTHS.y"))]  # Remove unwanted columns
  names(df)[names(df) == "OS_STATUS.x"] <- "OS_STATUS"  # Rename OS_STATUS.x to OS_STATUS
  names(df)[names(df) == "OS_MONTHS.x"] <- "OS_MONTHS"  # Rename OS_MONTHS.x to OS_MONTHS
  return(df)
})
names(merged_list) <- c("Train", "Test", "GSE39582", "GSE161158", "TCGA_BRCA", "TCGA_LUAD", "TCGA_SKCM", "TCGA_GBM", "TCGA_PAAD")
lapply(merged_list, function(df) {
  head(df)  # View the first few rows of each merged dataset
})

## Risk grouping based on risk score
library(survival)
library(survminer)
library(survivalROC)
library(survMisc)
library(survAUC)

cat_list <- lapply(seq_along(merged_list), function(i) {
  x <- merged_list[[i]]
  dataset_name <- names(merged_list)[i]
  x$OS_STATUS <- as.numeric(x$OS_STATUS)
  x$OS_MONTHS <- as.numeric(x$OS_MONTHS)
  x <- na.omit(x)
  x <- x[!(x$OS_MONTHS < 1),]
  cols = c("RS_ferroptosis", "RS_autosis", "RS_pyroptosis", "RS_NETosis") # save the unwanted columns
  x <- x[, -which(colnames(x) %in% cols)] # Remove unwanted columns
  res.cut <- surv_cutpoint(x, time = "OS_MONTHS", event = "OS_STATUS",
                           variables = c("RCDI"))
  jpeg_filename <- paste0("RCDI_cutoff_", dataset_name, "_2.jpeg")
  jpeg(jpeg_filename, units = "in", width = 6, height = 5, res = 600)
  print(plot(res.cut, "RCDI", 
             palette = c("salmon", "skyblue"),  # Low=skyblue, High=salmon
             ggtheme = theme_minimal() + 
               theme(
                 axis.title = element_text(size = 24),  # Axis titles (X/Y labels)
                 axis.text = element_text(size = 14)    # Axis tick labels
               )
  ))
  dev.off()
  cutoff_value <- res.cut$cutpoint[[1]]
  res.cat <- surv_categorize(res.cut, variables = "RCDI")
  names(res.cat)[names(res.cat) == "RCDI"] <- "RiskClass"
  res.cat1 <- res.cat[, c("RiskClass"), drop=F]
  m <- merge(x,res.cat1, by=0)
  row.names(m) <- m[,1]
  m<- m[,-1]
  return(m)
})


names(cat_list) <- c("Train", "Test", "GSE39582", "GSE161158", "TCGA_BRCA", "TCGA_LUAD", "TCGA_SKCM", "TCGA_GBM", "TCGA_PAAD")
# save the full datasets
lapply(names(cat_list), function(dataset_name) {
  # Fetch the dataframe
  df <- cat_list[[dataset_name]]
  # Define the output file name based on the dataset name
  output_file <- paste0(dataset_name, "_RCDI_categorized_2.csv")
  # Save the dataframe to a CSV file
  fwrite(df, output_file, sep = ",", row.names = TRUE)
  # Print a message for confirmation
  print(paste("Saved dataset:", dataset_name, "to file:", output_file))
})


cat_list <- lapply(seq_along(merged_list), function(i) {
  x <- merged_list[[i]]
  dataset_name <- names(merged_list)[i]
  
  # Preprocessing
  x$OS_STATUS <- as.numeric(x$OS_STATUS)
  x$OS_MONTHS <- as.numeric(x$OS_MONTHS)
  x <- na.omit(x)
  x <- x[!(x$OS_MONTHS < 1),]
  
  # Drop unwanted risk score columns
  cols = c("RS_ferroptosis", "RS_autosis", "RS_pyroptosis", "RS_NETosis")
  x <- x[, !(colnames(x) %in% cols)]
  
  # -------------------------------
  # Use median to categorize RCDI
  # -------------------------------
  cutoff_value <- median(x$RCDI, na.rm = TRUE)
  x$RiskClass <- ifelse(x$RCDI > cutoff_value, "high", "low")
  
  # Save cutoff plot (optional, you can comment this out if not needed)
  jpeg_filename <- paste0("RCDI_cutoff_", dataset_name, "_median.jpeg")
  jpeg(jpeg_filename, units = "in", width = 6, height = 5, res = 600)
  hist(x$RCDI, breaks = 40, col = "gray80", main = paste0("RCDI - ", dataset_name),
       xlab = "RCDI", ylab = "Frequency")
  abline(v = cutoff_value, col = "red", lwd = 2, lty = 2)
  legend("topright", legend = paste("Median =", round(cutoff_value, 3)),
         col = "red", lty = 2, lwd = 2, cex = 1)
  dev.off()
  
  return(x)
})



# Rename datasets in the list
names(cat_list) <- c("Train", "Test", "GSE39582", "GSE161158", "TCGA_BRCA", "TCGA_LUAD", "TCGA_SKCM", "TCGA_GBM", "TCGA_PAAD")

# Save categorized datasets using median-based RiskClass assignment
lapply(names(cat_list), function(dataset_name) {
  df <- cat_list[[dataset_name]]
  output_file <- paste0(dataset_name, "_RCDI_categorized_median.csv")  # Use "_median" for clarity
  fwrite(df, output_file, sep = ",", row.names = TRUE)
  print(paste("Saved dataset:", dataset_name, "to file:", output_file))
})

############################## RCDI Risk groups plot - All ds ######################################
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
# devtools::install_github("zabore/condsurv")
library(condsurv)
library(ggplot2)
library("survival")
library("survminer")

# palette <- "jco"        # Journal of Clinical Oncology color palette
names(cat_list) <- c("Train", "Test", "GSE39582", "GSE161158", "TCGA_BRCA", "TCGA_LUAD", "TCGA_SKCM", "TCGA_GBM", "TCGA_PAAD")
# Run survival analysis and save each plot
km_plot_list <- bplapply(seq_along(cat_list), function(i) {
  data_i <- cat_list[[i]]
  dataset_name <- names(cat_list)[i]
  
  # Ensure RiskClass is treated as factor
  data_i <- data_i[, -which(colnames(data_i) %in% c("RCDI"))]
  
  # KM and Cox model
  kmcurve <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ RiskClass, data = data_i)
  coxph_model <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ RiskClass, data = data_i)
  hr <- summary(coxph_model)$coefficients[1, "exp(coef)"]
  hr <- 1 / hr  # Invert to express HR (High vs Low)
  
  # Output file path
  file_name <- paste0("KMplot_", dataset_name, "_RCDI_High_vs_Low_2.tiff")
  
  # Save TIFF
  tiff(file_name, res = 600, units = "in", height = 6, width = 6, compression = "lzw")
  
  g <- ggsurvplot(
    fit = kmcurve,
    data = data_i,
    risk.table = TRUE,
    conf.int = TRUE,
    pval = TRUE,
    palette = c("salmon", "skyblue"),
    xlab = "Time (Months)",
    ylab = "Survival Probability",
    legend.labs = c("High-Risk", "Low-Risk"),
    legend.title = "Risk Group",
    font.main = c(18, "bold"),
    font.x = c(16, "bold"),
    font.y = c(16, "bold"),
    font.tickslab = c(14, "plain"),
    font.legend = c(14, "plain"),
    risk.table.fontsize = 4.8,
    risk.table.height = 0.25,
    ggtheme = theme_classic(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        axis.ticks = element_line(color = "black"),
        legend.position = "bottom",
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold")
      )
  )
  
  # Annotate HR
  g$plot <- g$plot +
    ggplot2::annotate(
      "text",
      x = 0.75 * max(kmcurve$time),
      y = 0.97,
      label = paste0("HR: ", round(hr, 2)),
      size = 5.5,
      fontface = "bold",
      color = "black"
    )
  
  print(g)
  dev.off()
}, BPPARAM = MulticoreParam(num_detected_cores))


## median cutoff

library(survival)
library(survminer)
library(ggplot2)
library(BiocParallel)

# Optional: Set names if not already done
names(cat_list) <- c("Train", "Test", "GSE39582", "GSE161158", "TCGA_BRCA", "TCGA_LUAD", "TCGA_SKCM", "TCGA_GBM", "TCGA_PAAD")

# Parallel KM Plotting
km_plot_list <- bplapply(seq_along(cat_list), function(i) {
  data_i <- cat_list[[i]]
  dataset_name <- names(cat_list)[i]
  
  # Ensure RiskClass is a factor (low, high)
  data_i$RiskClass <- factor(data_i$RiskClass, levels = c("low", "high"))
  
  # Remove RCDI if present (not needed for KM)
  data_i <- data_i[, !(colnames(data_i) %in% "RCDI")]
  
  # KM & Cox model
  kmcurve <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ RiskClass, data = data_i)
  coxph_model <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ RiskClass, data = data_i)
  
  # Compute HR (High vs Low)
  hr <- summary(coxph_model)$coefficients[1, "exp(coef)"]
  
  # Plot file path
  file_name <- paste0("KMplot_", dataset_name, "_RCDI_High_vs_Low_median.tiff")
  
  # Save TIFF
  tiff(file_name, res = 600, units = "in", height = 6, width = 6, compression = "lzw")
  
  g <- ggsurvplot(
    fit = kmcurve,
    data = data_i,
    risk.table = TRUE,
    conf.int = TRUE,
    pval = TRUE,
    palette = c("skyblue", "salmon"),  # low = skyblue, high = salmon
    xlab = "Time (Months)",
    ylab = "Survival Probability",
    legend.labs = c("Low-Risk", "High-Risk"),
    legend.title = "Risk Group",
    font.main = c(18, "bold"),
    font.x = c(16, "bold"),
    font.y = c(16, "bold"),
    font.tickslab = c(14, "plain"),
    font.legend = c(14, "plain"),
    risk.table.fontsize = 4.8,
    risk.table.height = 0.25,
    ggtheme = theme_classic(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        axis.ticks = element_line(color = "black"),
        legend.position = "bottom",
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold")
      )
  )
  
  # Annotate HR
  g$plot <- g$plot +
    annotate("text",
             x = 0.75 * max(kmcurve$time),
             y = 0.97,
             label = paste0("HR: ", round(hr, 2)),
             size = 5.5,
             fontface = "bold",
             color = "black")
  
  print(g)
  dev.off()
}, BPPARAM = MulticoreParam(num_detected_cores))


############################ Boxplots for individual genes - all DS ####################################
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(data.table)
library(dplyr)

Train <- cat_list$Train
cols <- c("OS_MONTHS", "OS_STATUS", "RiskClass","RCDI")
gene_sig <- colnames(Train[, -which(colnames(Train) %in% cols)])

lapply(names(cat_list), function(ds_name) {
  df <- cat_list[[ds_name]]
  
  # Create a riskScore column from RiskClass for plotting
  df$RiskClass <- factor(df$RiskClass, levels = c("low", "high"))
  
  # Loop through each gene
  for (gene in gene_sig) {
    # Skip non-numeric gene columns
    if (!is.numeric(df[[gene]])) next
    
    message("Plotting: ", gene, " for ", ds_name)
    
    # Output filename
    file_name <- paste0("./Individual_genes_expr_high_vs_low_RCDI/",ds_name,"/boxplot_high_vs_low_risk_", ds_name, "_", gene, "_zscore.jpeg")
    
    # Save high-res JPEG
    jpeg(filename = file_name, units = "in", width = 7, height = 5, res = 600)
    
    # Create boxplot with ggpubr
    p <- ggboxplot(df, x = "RiskClass", y = gene, 
                   add = "jitter", color = "RiskClass", size = 0.8,
                   palette = "jco", add.params = list(size = 1.5, alpha = 0.5),
                   bxp.errorbar = TRUE, width = 0.3, legend = "none", order = c("low", "high")) +
      stat_compare_means(method = "t.test", label = "p.format", size = 5,
                         ref.group = "low",
                         symnum.args = list(
                           cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
                           symbols = c("****", "***", "**", "*", "ns")
                         )) +
      xlab('Risk Class') + 
      ylab("Z-Score") +
      ggtitle(paste0("Expression - ", gene, " (", ds_name, ")")) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    print(p)
    dev.off()
  }
})

# 
############################ KM plots for individual genes ####################################

Train <- cat_list$Train
cols <- c("OS_MONTHS", "OS_STATUS", "RiskClass","RCDI")
gene_sig <- colnames(Train[, -which(colnames(Train) %in% cols)])


library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
# devtools::install_github("zabore/condsurv")
library(condsurv)
library(ggplot2)
library("survival")
library("survminer")
lapply(names(cat_list), function(ds_name) {
  df <- cat_list[[ds_name]]
  df <- df[, -which(colnames(df) %in% c("RiskClass", "RCDI"))]
  for (i in 1:length(gene_sig)) {
  	gene <- gene_sig[i]
  	print(gene)
  	cutpoint <- surv_cutpoint(df, time = "OS_MONTHS", event = "OS_STATUS", variables = gene)
  	optimal_cut <- cutpoint$cutpoint[[1]]
  	df$group <- ifelse(df[,names(df) %in% gene, drop=F] >= optimal_cut, "High_expression", "Low_expression")
  	
  	jpeg(paste0("./KM_plots_high_vs_low_expr_RCDI/",ds_name,"/KM_plot_high_vs_low_", gene, "_", ds_name, ".jpeg"), 
  	     res=600, units="in", height=7, width=7)
  	kmcurve<- survfit(Surv(OS_MONTHS, OS_STATUS) ~ group, data=df)
  	coxph_model <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ group, data = df)
  	hr <- summary(coxph_model)$coefficients[1, "exp(coef)"]
  	hr<- 1/hr
  	# Kaplan-Meier Plot with Large Fonts
  	g <- ggsurvplot(
  	  kmcurve, 
  	  conf.int = TRUE,
  	  pval = TRUE,
  	  pval.size = 6,  # Larger p-value text
  	  risk.table = TRUE,
  	  risk.table.height = 0.2,  # Reduce risk table size
  	  risk.table.fontsize = 4,  # Smaller risk table text
  	  xlab = "Time (Months)",
  	  ylab = "Survival Probability",
  	  legend.title = "Expression",
  	  legend.labs = c("High", "Low"),
  	  title = paste0("KM Plot: ", gene)
  	    # Large base font for scientific publications
  	)
  	
  	# Adjust annotations
  	g$plot <- g$plot + 
  	  ggplot2::annotate("text", x = 0.75 * max(kmcurve$time), y = 0.95, 
  	                    label = paste0("HR: ", round(hr, 3)), size = 6, fontface = "bold") +
  	  theme(
  	    plot.title = element_text(size = 20, face = "bold"),
  	    axis.title = element_text(size = 18, face = "bold"),
  	    axis.text = element_text(size = 16),
  	    legend.text = element_text(size = 14),
  	    legend.key.size = unit(1.2, "cm")
  	  )
  	
  	print(g)
  	dev.off()
  }
})

################################## Forest plots ####################################

###### All RCDI geneset #######

library(survival)
library(survminer)
library(MASS)
# Do multivariate analysis on the geneset
data <- data.frame(fread("./Train_RCDI_categorized.csv", sep = ",", header = T), row.names = 1, check.names = F)
data$OS_MONTHS <- as.numeric(data$OS_MONTHS)
data$OS_STATUS <- as.numeric(data$OS_STATUS)
cols <- c("RCDI", "RiskClass")
data <- data[, -which(colnames(data) %in% cols)]
data <- na.omit(data)
data <- data[data$OS_MONTHS >= 1, ]
var <- setdiff(colnames(data), c("OS_MONTHS", "OS_STATUS"))

formula_string <- paste("Surv(OS_MONTHS, OS_STATUS) ~ ", paste(var, collapse = " + "))
formula <- as.formula(formula_string)
cox_model <- coxph(formula, data = data)

tiff("forest_plot_all_genes_RCDI_multi.tiff", res = 600, units = "in", height = 14, width = 9, compression = "lzw")
ggforest(cox_model,
         data = data,
         fontsize = 0.8,                 # Increase font size
         refLabel = "Reference",
         noDigits = 2) +                 # Add precision
  theme_void(base_size = 14) +       # Cleaner base theme
  theme(
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 13),
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    strip.text = element_text(face = "bold", size = 14)
  )
dev.off()

# Do univariate analysis on the geneset

library(survival)
library(data.table)
library(ggplot2)

# Load and clean the data
data <- data.frame(fread("./Train_RCDI_categorized.csv", sep = ",", header = T), row.names = 1, check.names = F)
data$OS_MONTHS <- as.numeric(data$OS_MONTHS)
data$OS_STATUS <- as.numeric(data$OS_STATUS)
cols <- c("RCDI", "RiskClass")
data <- data[, -which(colnames(data) %in% cols)]
data <- na.omit(data)
data <- data[data$OS_MONTHS >= 1, ]
var <- setdiff(colnames(data), c("OS_MONTHS", "OS_STATUS"))

# Run univariate Cox regressions
univ_results <- lapply(var, function(v) {
  f <- as.formula(paste("Surv(OS_MONTHS, OS_STATUS) ~", v))
  model <- coxph(f, data = data)
  summary_model <- summary(model)
  data.frame(
    Variable = v,
    HR = summary_model$coefficients[1, "exp(coef)"],
    CI_lower = summary_model$conf.int[1, "lower .95"],
    CI_upper = summary_model$conf.int[1, "upper .95"],
    p = summary_model$coefficients[1, "Pr(>|z|)"]
  )
})

# Combine results
univ_df <- do.call(rbind, univ_results)
univ_df$Significance <- ifelse(univ_df$p < 0.05, "*", "")

# Order by HR for plotting
univ_df <- univ_df[order(univ_df$HR, decreasing = TRUE), ]
univ_df$Variable <- factor(univ_df$Variable, levels = univ_df$Variable)

# Plot
tiff("forest_plot_all_genes_RCDI_univaraite.tiff", res = 600, units = "in", height = 12, width = 9, compression = "lzw")
ggplot(univ_df, aes(x = Variable, y = HR)) +
  geom_point(size = 3.5, color = "black") +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.4, color = "black") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray30", size = 0.8) +
  coord_flip() +
  geom_text(
    aes(label = sprintf("p = %.3g%s", p, Significance), y = max(CI_upper) + 0.6),
    hjust = 0, size = 4.2, color = "black"
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Univariate Cox Regression – RCDI Gene Set",
    x = NULL,
    y = "Hazard Ratio (95% CI)"
  ) +
  theme(
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 14, face = "bold"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  expand_limits(y = max(univ_df$CI_upper) + 1)  # Add space on x-axis for p-values
dev.off()

################ Clinical factors ###################
library(dplyr)
columns <- c("Sample","Sex", "Age_group", "Tumor_stage")
clinical_emtab <- read_xlsx("../../E-MTAB-12862/data/Supplementary_Table_01.xlsx", sheet=1)
clinical_emtab <- clinical_emtab[,c("RNA Tumor Sample Barcode", "Sex", "Age group","Tumour Stage")]
clinical_emtab <- na.omit(clinical_emtab)
names(clinical_emtab) <- columns
clinical_emtab <- as.data.frame(clinical_emtab)

tcga_coad <- read.table("../../UCSC_Xena_Browser/TCGA-COAD.clinical.tsv", sep='\t', check.names=F, header=T)
tcga_read <- read.table("../../UCSC_Xena_Browser/TCGA-READ.clinical.tsv", sep='\t', check.names=F, header=T)
clinical_cbio <- rbind(tcga_coad, tcga_read)
clinical_cbio <- clinical_cbio[, c("sample", "gender.demographic","age_at_index.demographic","ajcc_pathologic_stage.diagnoses")]
clinical_cbio <- na.omit(clinical_cbio)
names(clinical_cbio) <- columns

clinical_cbio <- clinical_cbio %>%
  mutate(Age_group = case_when(
    Age_group <= 65 ~ "<= 65",
    Age_group >= 66 & Age_group <= 79 ~ "66-79",
    Age_group >= 80 ~ ">= 80"
  ))
names(clinical_cbio) <- columns
clinical_emtab <- clinical_emtab %>%
  mutate(Age_group = case_when(
    Age_group == "≤ 65"~ "<= 65" ,
    Age_group == "≥ 80" ~ ">= 80",
    Age_group == "66-79" ~ "66-79",
  ))
clinical_cbio <- clinical_cbio %>%
  mutate(Sex = case_when(
   Sex == "male" ~ "Male",
   Sex == "female" ~ "Female"
    ))

clinical_combined <- rbind(clinical_emtab, clinical_cbio)

clinical_gse <- data.frame(fread("../../data_GSE39582_array/counts/clinical_data.csv", sep=',', header = T), check.names=F)
cols <- c("Sample","Sex","Years_at_diagnosis","TNM.stage")
clinical_gse <- clinical_gse[,names(clinical_gse) %in% cols]
clinical_gse <- na.omit(clinical_gse)
names(clinical_gse) <- columns
clinical_gse <- as.data.frame(clinical_gse)
clinical_gse <- clinical_gse %>%
  mutate(Age_group = case_when(
    Age_group <= 65 ~ "<= 65",
    Age_group >= 66 & Age_group <= 79 ~ "66-79",
    Age_group >= 80 ~ ">= 80"
  ))


clinical_combined <- rbind(clinical_combined, clinical_gse)

row.names(clinical_combined) <- clinical_combined[,1]
clinical_combined <- clinical_combined[,-1]
clinical_combined <- clinical_combined %>%
  mutate(Tumor_stage = case_when(
    Tumor_stage %in% c("Stage I", "Stage IA", "Stage II", "Stage IIA", "Stage IIB", "Stage IIC") ~ "Stage I - II",
    Tumor_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC", "Stage IV", "Stage IVA", "Stage IVB") ~ "Stage III - IV",
    Tumor_stage %in% c("1", "2") ~ "Stage I - II",
    Tumor_stage %in% c("3", "4") ~ "Stage III - IV",
    TRUE ~ Tumor_stage  # <- default: keep as-is if no match
  ))
clinical_combined <- clinical_combined %>%
  filter(across(everything(), ~ !(. %in% c("", "N/A", "0"))))
clinical_combined <- na.omit(clinical_combined)
write.table(clinical_combined, "clinical_combined_cbio_emtab__gse_43_genes.csv", sep=",", row.names=T, col.names=NA)

clinical_combined <- data.frame(fread("clinical_combined_cbio_emtab__gse_43_genes.csv", sep=',', header = T), row.names=1, check.names=F)
file_paths <- c("Train_RCDI_categorized.csv", "Test_RCDI_categorized.csv", "GSE39582_RCDI_categorized.csv")
dataset_names <- c("Train", "Test", "GSE39582")

auc_results <- data.frame(
  Dataset = character(),
  Time_months = numeric(),
  Nomogram_AUROC = numeric(),
  Sex_AUROC = numeric(),
  Age_group_AUROC = numeric(),
  Tumor_stage_AUROC = numeric(),
  RiskClass_AUROC = numeric(),
  stringsAsFactors = FALSE
)
lapply(seq_along(file_paths), function(i) {
  dataset <- dataset_names[i]
  file <- file_paths[i]
  
  message("Processing dataset: ", dataset)
  
  # Load expression data
  train_class <- data.frame(fread(file, sep = ',', header = TRUE), row.names = 1, check.names = FALSE)
  
  # Merge with clinical
  merged_data <- merge(train_class, clinical_combined, by = 0)
  row.names(merged_data) <- merged_data[, 1]
  merged_data <- merged_data[, -1]
  fwrite(merged_data, paste0(dataset, "_clinical_expr_combined_rcdi.csv"), sep = ",", row.names = TRUE)
  
  # Clean & format
  df <- merged_data
  df$OS_MONTHS <- as.numeric(df$OS_MONTHS)
  df$OS_STATUS <- as.numeric(df$OS_STATUS)
  df$RiskClass <- relevel(factor(df$RiskClass), ref = "low")
  df$Sex <- factor(df$Sex)
  df$Age_group <- factor(df$Age_group)
  df$Tumor_stage <- factor(df$Tumor_stage)
  
  # Cox formula and model
  cox_formula <- as.formula("Surv(OS_MONTHS, OS_STATUS) ~ `RCDI` + Sex + Age_group + Tumor_stage")
  cox_clinical <- coxph(cox_formula, data = df)
  summary_result <- summary(cox_clinical)
  fwrite(summary_result$coefficients, paste0("cox_model_coefficients_", dataset, ".csv"), sep = ",", row.names = T)
  
  # Forest plot
  svg(paste0("forest_plot_clinical_RCDI_", dataset, ".svg"), width = 10, height = 8)
  gg<- ggforest(cox_clinical,
           data = df,
           fontsize = 1,
           refLabel = "Reference",
           noDigits = 2,
           main = paste("Multivariate Cox: ", dataset))
  print(gg)
  dev.off()
  
  # Nomogram plot
  library(regplot)
  p <- regplot(cox_clinical,
          observation = df[1, ],
          failtime = c(12, 36, 60),
          prfail = TRUE,
          droplines = TRUE,
          title = paste("Nomogram: ", dataset),
          points = TRUE,
          rank = "range")
  class(p)
  dev.copy2pdf(file=paste0("nomogram_RCDI_43_genes_",dataset,".pdf"), width=10, height=7)
  
  # Nomogram risk score
  clinical_factors <- c("RCDI", "Sex", "Age_group", "Tumor_stage")
  df_data <- df[, clinical_factors]
  df$nomogram_score <- predict(cox_clinical, newdata = df_data, type = "risk")
  
  # Time-dependent ROC
  time_points <- c(12, 36, 60)
  roc_curve <- timeROC(T = df$OS_MONTHS,
                       delta = df$OS_STATUS,
                       marker = df$nomogram_score,
                       cause = 1,
                       times = time_points,
                       iid = TRUE)
  # Store AUC results
  auc_nomogram <- roc_curve$AUC

  # Save ROC plot
  svg(paste0("nomogram_ROC_AUC_RCDI_", dataset, ".svg"), height = 6, width = 6)
  plot(roc_curve, time = 12, col = "red", main = paste("ROC Curve -", dataset))
  plot(roc_curve, time = 36, col = "blue", add = TRUE)
  plot(roc_curve, time = 60, col = "green", add = TRUE)
  legend("bottomright",
         legend = c(paste("1-Year AUC:", round(roc_curve$AUC[1], 3)),
                    paste("3-Year AUC:", round(roc_curve$AUC[2], 3)),
                    paste("5-Year AUC:", round(roc_curve$AUC[3], 3))),
         col = c("red", "blue", "green"), lwd = 2)
  dev.off()
  
  ## Marginal ROC comparisons at 60 months
  time_points <- c(12, 36, 60)
  
  # Prepare the data (as in your original code)
  roc_data <- df %>%
    mutate(
      Age_group = case_when(
        Age_group == "<= 65" ~ 1,
        Age_group == "66-79" ~ 2,
        Age_group == ">= 80" ~ 3
      ),
      Tumor_stage = case_when(
        Tumor_stage == "Stage I - II" ~ 1,
        Tumor_stage == "Stage III - IV" ~ 2
      ),
      class = case_when(
        RiskClass == "high" ~ 2,
        RiskClass == "low" ~ 1
      ),
      Sex = case_when(
        Sex == "Female" ~ 1,
        Sex == "Male" ~ 2
      )
    )
  
  # Loop over each time point and generate plots
  for (t in seq_along(time_points)) {
    tp <- time_points[t]
    ROC.data.sex <- timeROC(
      T = roc_data$OS_MONTHS,
      delta = roc_data$OS_STATUS,
      marker = roc_data$Sex,
      other_markers = as.matrix(roc_data[, c("Age_group", "Tumor_stage", "class")]),
      cause = 1,
      weighting = "marginal",
      times = c(tp),
      iid = TRUE
    )
    ROC.data.age <- timeROC(
      T = roc_data$OS_MONTHS,
      delta = roc_data$OS_STATUS,
      marker = roc_data$Age_group,
      other_markers = as.matrix(roc_data[, c("Sex", "Tumor_stage", "class")]),
      cause = 1,
      weighting = "marginal",
      times = c(tp),
      iid = TRUE
    )
    ROC.data.stage <- timeROC(
      T = roc_data$OS_MONTHS,
      delta = roc_data$OS_STATUS,
      marker = roc_data$Tumor_stage,
      other_markers = as.matrix(roc_data[, c("Age_group", "Sex", "class")]),
      cause = 1,
      weighting = "marginal",
      times = c(tp),
      iid = TRUE
    )
    ROC.data.class <- timeROC(
      T = roc_data$OS_MONTHS,
      delta = roc_data$OS_STATUS,
      marker = roc_data$class,
      other_markers = as.matrix(roc_data[, c("Age_group", "Sex", "Tumor_stage")]),
      cause = 1,
      weighting = "marginal",
      times = c(tp),
      iid = TRUE
    )
    
    # Extract AUCs (use [2] if you have two time points, but here only one, so use [1])
    auc_sex <- ROC.data.sex$AUC[2]
    auc_age <- ROC.data.age$AUC[2]
    auc_stage <- ROC.data.stage$AUC[2]
    auc_class <- ROC.data.class$AUC[2]
    
    # Store AUC results
    auc_results <<- rbind(auc_results, data.frame(
      Dataset = dataset,
      Time_months = tp,
      Nomogram_AUROC = auc_nomogram[t],
      Sex_AUROC = auc_sex,
      Age_group_AUROC = auc_age,
      Tumor_stage_AUROC = auc_stage,
      RiskClass_AUROC = auc_class,
      stringsAsFactors = FALSE
    ))
    
    # Save ROC comparison plot for this time point
    jpeg(paste0("ROC_curves_marginal_factors_", dataset, "_", tp, "_months.jpeg"), width = 1800, height = 1800, res=300)
    plot(ROC.data.class, time = tp, col = "red", lwd = 2, title = FALSE)
    plot(ROC.data.sex, time = tp, col = "blue", add = TRUE, lwd = 2, lty = 4)
    plot(ROC.data.age, time = tp, col = "violet", add = TRUE, lwd = 2, lty = 3)
    plot(ROC.data.stage, time = tp, col = "green", add = TRUE, lwd = 2, lty = 3)
    abline(v = 0.5, col = "gray", lty = 2)
    legend("bottomright",
           legend = c(
             paste("Risk Class (AUC =", round(auc_class, 2), ")"),
             paste("Sex (AUC =", round(auc_sex, 2), ")"),
             paste("Age Group (AUC =", round(auc_age, 2), ")"),
             paste("Tumor Stage (AUC =", round(auc_stage, 2), ")")
           ),
           col = c("red", "blue", "violet", "green"),
           lty = c(1, 4, 3, 3), lwd = 2)
    title(paste("ROC Curves for Clinical Factors at", tp, "Months -", dataset))
    dev.off()
  }
})

# Save AUC results
write.csv(auc_results, "auroc_timepoints_summary.csv", row.names = FALSE)


### For GSE161158 dataset ######
columns <- c("Sample", "Age_group", "Tumor_stage")
clinical_gse161158 <- read.table("../../array_GSE161158/clinical_data.csv", sep=",", header=T, stringsAsFactors=F)
clinical_gse161158 <- clinical_gse161158[, c("Sample", "age", "Stage")]
names(clinical_gse161158) <- columns
clinical_gse161158 <- as.data.frame(clinical_gse161158)
clinical_gse161158 <- clinical_gse161158 %>%
  mutate(Age_group = case_when(
    Age_group <= 65 ~ "<= 65",
    Age_group >= 66 & Age_group <= 79 ~ "66-79",
    Age_group >= 80 ~ ">= 80"
  ))
clinical_gse161158 <- clinical_gse161158 %>%
  mutate(Tumor_stage = case_when(
    Tumor_stage %in% c("stage 1") ~ "Stage I - II",
    Tumor_stage %in% c("stage 2") ~ "Stage I - II",
    Tumor_stage %in% c("stage 3") ~ "Stage III - IV",
    Tumor_stage %in% c("stage 4") ~ "Stage III - IV",
    TRUE ~ Tumor_stage  # <- default: keep as-is if no match
  ))
clinical_gse161158 <- clinical_gse161158 %>%
  filter(across(everything(), ~ !(. %in% c("", "N/A", "0"))))
clinical_gse161158 <- na.omit(clinical_gse161158)

file_paths <- c("GSE161158_RCDI_categorized.csv")
dataset_names <- c("GSE161158")

auc_results <- data.frame(
  Dataset = character(),
  Time_months = numeric(),
  Nomogram_AUROC = numeric(),
  Age_group_AUROC = numeric(),
  Tumor_stage_AUROC = numeric(),
  RiskClass_AUROC = numeric(),
  stringsAsFactors = FALSE
)
lapply(seq_along(file_paths), function(i) {
  dataset <- dataset_names[i]
  file <- file_paths[i]
  
  message("Processing dataset: ", dataset)
  
  # Load expression data
  train_class <- data.frame(fread(file, sep = ',', header = TRUE), row.names = 1, check.names = FALSE)
  
  # Merge with clinical
  merged_data <- merge(train_class, clinical_gse161158, by.x = 0, by.y=1)
  row.names(merged_data) <- merged_data[, 1]
  merged_data <- merged_data[, -1]
  fwrite(merged_data, paste0(dataset, "_clinical_expr_combined_rcdi.csv"), sep = ",", row.names = TRUE)
  
  # Clean & format
  df <- merged_data
  df$OS_MONTHS <- as.numeric(df$OS_MONTHS)
  df$OS_STATUS <- as.numeric(df$OS_STATUS)
  df$RiskClass <- relevel(factor(df$RiskClass), ref = "low")
  # df$Sex <- factor(df$Sex)
  df$Age_group <- factor(df$Age_group)
  df$Tumor_stage <- factor(df$Tumor_stage)
  
  # Cox formula and model
  cox_formula <- as.formula("Surv(OS_MONTHS, OS_STATUS) ~ `RCDI` + Age_group + Tumor_stage")
  cox_clinical <- coxph(cox_formula, data = df)
  summary_result <- summary(cox_clinical)
  fwrite(summary_result$coefficients, paste0("cox_model_coefficients_", dataset, ".csv"), sep = ",", row.names = T)
  
  # Forest plot
  svg(paste0("forest_plot_clinical_RCDI_", dataset, ".svg"), width = 10, height = 8)
  gg<- ggforest(cox_clinical,
                data = df,
                fontsize = 1,
                refLabel = "Reference",
                noDigits = 2,
                main = paste("Multivariate Cox: ", dataset))
  print(gg)
  dev.off()
  
  # Nomogram plot
  library(regplot)
  p <- regplot(cox_clinical,
               observation = df[1, ],
               failtime = c(12, 36, 60),
               prfail = TRUE,
               droplines = TRUE,
               title = paste("Nomogram: ", dataset),
               points = TRUE,
               rank = "range")
  class(p)
  dev.copy2pdf(file=paste0("nomogram_RCDI_43_genes_",dataset,".pdf"), width=10, height=7)
  
  # Nomogram risk score
  clinical_factors <- c("RCDI", "Age_group", "Tumor_stage")
  df_data <- df[, clinical_factors]
  df$nomogram_score <- predict(cox_clinical, newdata = df_data, type = "risk")
  
  # Time-dependent ROC
  time_points <- c(12, 36, 60)
  roc_curve <- timeROC(T = df$OS_MONTHS,
                       delta = df$OS_STATUS,
                       marker = df$nomogram_score,
                       cause = 1,
                       times = time_points,
                       iid = TRUE)
  # Store AUC results
  auc_nomogram <- roc_curve$AUC
  
  # Save ROC plot
  svg(paste0("nomogram_ROC_AUC_RCDI_", dataset, ".svg"), height = 6, width = 6)
  plot(roc_curve, time = 12, col = "red", main = paste("ROC Curve -", dataset))
  plot(roc_curve, time = 36, col = "blue", add = TRUE)
  plot(roc_curve, time = 60, col = "green", add = TRUE)
  legend("bottomright",
         legend = c(paste("1-Year AUC:", round(roc_curve$AUC[1], 3)),
                    paste("3-Year AUC:", round(roc_curve$AUC[2], 3)),
                    paste("5-Year AUC:", round(roc_curve$AUC[3], 3))),
         col = c("red", "blue", "green"), lwd = 2)
  dev.off()
  
  ## Marginal ROC comparisons at 60 months
  time_points <- c(12, 36, 60)
  
  # Prepare the data (as in your original code)
  roc_data <- df %>%
    mutate(
      Age_group = case_when(
        Age_group == "<= 65" ~ 1,
        Age_group == "66-79" ~ 2,
        Age_group == ">= 80" ~ 3
      ),
      Tumor_stage = case_when(
        Tumor_stage == "Stage I - II" ~ 1,
        Tumor_stage == "Stage III - IV" ~ 2
      ),
      class = case_when(
        RiskClass == "high" ~ 2,
        RiskClass == "low" ~ 1
      )
    )
  
  # Loop over each time point and generate plots
  for (t in seq_along(time_points)) {
    tp <- time_points[t]
    ROC.data.age <- timeROC(
      T = roc_data$OS_MONTHS,
      delta = roc_data$OS_STATUS,
      marker = roc_data$Age_group,
      other_markers = as.matrix(roc_data[, c("Tumor_stage", "class")]),
      cause = 1,
      weighting = "marginal",
      times = c(tp),
      iid = TRUE
    )
    ROC.data.stage <- timeROC(
      T = roc_data$OS_MONTHS,
      delta = roc_data$OS_STATUS,
      marker = roc_data$Tumor_stage,
      other_markers = as.matrix(roc_data[, c("Age_group", "class")]),
      cause = 1,
      weighting = "marginal",
      times = c(tp),
      iid = TRUE
    )
    ROC.data.class <- timeROC(
      T = roc_data$OS_MONTHS,
      delta = roc_data$OS_STATUS,
      marker = roc_data$class,
      other_markers = as.matrix(roc_data[, c("Age_group", "Tumor_stage")]),
      cause = 1,
      weighting = "marginal",
      times = c(tp),
      iid = TRUE
    )
    
    # Extract AUCs (use [2] if you have two time points, but here only one, so use [1])
    auc_age <- ROC.data.age$AUC[2]
    auc_stage <- ROC.data.stage$AUC[2]
    auc_class <- ROC.data.class$AUC[2]
    
    # Store AUC results
    auc_results <<- rbind(auc_results, data.frame(
      Dataset = dataset,
      Time_months = tp,
      Nomogram_AUROC = auc_nomogram[t],
      Age_group_AUROC = auc_age,
      Tumor_stage_AUROC = auc_stage,
      RiskClass_AUROC = auc_class,
      stringsAsFactors = FALSE
    ))
    
    # Save ROC comparison plot for this time point
    jpeg(paste0("ROC_curves_marginal_factors_", dataset, "_", tp, "_months.jpeg"), width = 1800, height = 1800, res=300)
    plot(ROC.data.class, time = tp, col = "red", lwd = 2, title = FALSE)
    # plot(ROC.data.sex, time = tp, col = "blue", add = TRUE, lwd = 2, lty = 4)
    plot(ROC.data.age, time = tp, col = "violet", add = TRUE, lwd = 2, lty = 3)
    plot(ROC.data.stage, time = tp, col = "green", add = TRUE, lwd = 2, lty = 3)
    abline(v = 0.5, col = "gray", lty = 2)
    legend("bottomright",
           legend = c(
             paste("Risk Class (AUC =", round(auc_class, 2), ")"),
             # paste("Sex (AUC =", round(auc_sex, 2), ")"),
             paste("Age Group (AUC =", round(auc_age, 2), ")"),
             paste("Tumor Stage (AUC =", round(auc_stage, 2), ")")
           ),
           col = c("red", "violet", "green"),
           lty = c(1, 4, 3, 3), lwd = 2)
    title(paste("ROC Curves for Clinical Factors at", tp, "Months -", dataset))
    dev.off()
  }
})

write.csv(auc_results, "auroc_timepoints_summary_gse161158.csv", row.names = FALSE)


##### gene lists ##############
setwd("/home2/New_Objective2_Biomarker_research/gene_lists")
ferroptosis <- read.table("./ferroptosis_gene_list.txt", sep=",", header=F, stringsAsFactors=F)
ferroptosis <- unique(ferroptosis$V1) # 259
entosis <- read.table("./entosis_core_genes.txt", sep=",", header=F, stringsAsFactors=F)
entosis <- unique(entosis$V1) # 55
NETosis <- read.table("./NETosis_related_genes.txt", sep=",", header=F, stringsAsFactors=F)
NETosis <- unique(NETosis$V1) # 247
cuproptosis <- read.table("./cuproptosis_gene_list.txt", sep=",", header=F, stringsAsFactors=F)
cuproptosis <- unique(cuproptosis$V1) #64
ICD <- read.table("./ICD_gene_list.txt", sep=",", header=F, stringsAsFactors=F)
ICD <- unique(ICD$V1) #534
LCD <- read.table("./lcd_gene_list.txt", sep=",", header=F, stringsAsFactors=F)
LCD <- unique(LCD$V1) # 218
oxp <- read.table("./oxeiptosis_related_genes.txt", sep=",", header=F, stringsAsFactors=F)
oxp <- unique(oxp$V1) # 88
parthanatos <- read.table("./parthanatos_related_genes.txt", sep=",", header=F, stringsAsFactors=F)
parthanatos <- unique(parthanatos$V1) # 132
pyroptosis <- read.table("./pyroptosis_gene_list.txt", sep=",", header=F, stringsAsFactors=F)
pyroptosis <- unique(pyroptosis$V1) # 57
autosis <- read.table("./autosis_gene_list.txt", sep=",", header=F, stringsAsFactors=F)
autosis <- unique(autosis$V1) # 65
necroptosis <- read.table("./necroptosis_gene_list.txt", sep=",", header=F, stringsAsFactors=F)
necroptosis <- unique(necroptosis$V1) # 209
anoikis <- read.table("./anoikis_gene_list.txt", sep=",", header=F, stringsAsFactors=F)
anoikis <- unique(anoikis$V1) # 34
autophagy <- read.table("./autophagy_gene_list.txt", sep=",", header=F, stringsAsFactors=F)
autophagy <- unique(autophagy$V1) #380  

## common DEGs ##
TCGA_DEGs <- read.table("../batch_corrected_rna_seq/new/DEGs_TCGA_0.5.txt", sep="\t", header=T, stringsAsFactors=F)
EMTAB_DEGs <- read.table("../batch_corrected_rna_seq/new/DEGs_emtab_0.5.txt", sep="\t", header=T, stringsAsFactors=F)
common_DEGs <- intersect(TCGA_DEGs$X, EMTAB_DEGs$X)
common_DEGs <- unique(common_DEGs)

# Create a named list of RCD gene sets
rcd_gene_sets <- list(
  Ferroptosis = ferroptosis,
  Entosis = entosis,
  NETosis = NETosis,
  Cuproptosis = cuproptosis,
  ICD = ICD,
  LCD = LCD,
  Oxeiptosis = oxp,
  Parthanatos = parthanatos,
  Pyroptosis = pyroptosis,
  Autosis = autosis,
  Necroptosis = necroptosis,
  Anoikis = anoikis,
  Autophagy = autophagy
)

# Count how many common DEGs fall into each pathway
deg_counts_in_pathways <- sapply(rcd_gene_sets, function(gene_set) {
  length(intersect(common_DEGs, gene_set))
})

# Display result as a sorted table
deg_counts_in_pathways <- sort(deg_counts_in_pathways, decreasing = TRUE)
print(deg_counts_in_pathways)
