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




### This section is the steps for dataset loading ###
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
