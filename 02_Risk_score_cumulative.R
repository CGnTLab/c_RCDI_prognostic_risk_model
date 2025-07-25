#!/usr/bin/Rscript

#################### Basic settings #######################

setwd("/home2/New_Objective2_Biomarker_research/batch_corrected_rna_seq/new/")

options(repos = c(CRAN = "https://cran.rstudio.com/"))
gc()
rm(list=ls())
seed = 42
set.seed(seed)

ntree= 1000
nodesize=5
nodedepth=5

library(dplyr)
library(data.table)
library(randomForestSRC)
# install.packages("randomForestExplainer")
library(randomForestExplainer)

# library(readxl)


library(BiocParallel)
library(parallel)
num_detected_cores <- detectCores(logical = FALSE) - 2 # logical = FALSE gives physical cores
register(MulticoreParam(workers = num_detected_cores))


################ 1. For Parthanatos ##################


# RSF #

print("---------------------------------------")
print("Analysis started for Parthanatos: RSF")
print("---------------------------------------")
## Load datasets
merged <- data.frame(fread("combined_z_score_clinical_batch.csv", sep=','), row.names=1, check.names=F)

train  <- merged %>% dplyr::sample_frac(0.8) ##dim(train) #743 14558
test1   <- merged[!(rownames(merged) %in% rownames(train)), ]

degs_emtab <- data.frame(fread("DEGs_emtab_0.5.txt", sep = "\t", header = T), row.names=1, check.names = F)
degs_emtab <- row.names(degs_emtab)
degs_tcga <- data.frame(fread("DEGs_TCGA_0.5.txt", sep = "\t", header = T), row.names=1, check.names = F)
degs_tcga <- row.names(degs_tcga)
common_degs <- intersect(degs_emtab, degs_tcga)
genes_to_select <- common_degs

# train <- data.frame(fread("emtab_z_score_clinical.csv", sep=',', header=T), row.names=1, check.names=F)
# # # # train <- train[ (train$OS_MONTHS <= 140.38356164) ,]
# test1 <- data.frame(fread("tcga_z_score_clinical.csv", sep = ",", header = T), row.names=1, check.names=F)
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

gene_set <- read.table("../../gene_lists/parthanatos_related_genes.txt", sep = ",", header = FALSE, check.names = FALSE)$V1
gene_set <- unique(gene_set)

train <- train[, names(train) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test1 <- test1[, names(test1) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test2 <- test2[, names(test2) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test3 <- test3[, names(test3) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test4 <- test4[, names(test4) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test5 <- test5[, names(test5) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test6 <- test6[, names(test6) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test7 <- test7[, names(test7) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test8 <- test8[, names(test8) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]

library("survival")
library("survminer")

if (!requireNamespace("survivalROC", quietly = TRUE)) install.packages("survivalROC")
if (!requireNamespace("PRROC", quietly = TRUE)) install.packages("PRROC")

# Load required packages
library(survivalROC)   # For time-dependent ROC analysis
library(PRROC)   # For Precision-Recall (PR) AUC

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


# select the train dataset
train <- df_list$Train ## explicitly select the train dataset

# mention time points for survivalROC and PRROC
time_points <- c(12, 36, 60)

# Gene selection
genes <- covariates
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]
library(randomForestSRC)

rf_model_final <- rfsrc(Surv(OS_MONTHS, OS_STATUS) ~ ., data = train_data,
                        importance = T, ntree = ntree, # Use full ntree for final model
                        proximity = T, forest = T,
                        seed = seed, nodesize = nodesize,
                        bootstrap = "by.root", ncores= num_detected_cores, nodedepth = nodedepth) # Use topvars for final model
vs.rf_model <- var.select(object = rf_model_final, conservative = "high")
topvars <- vs.rf_model$topvars
length(topvars)

print(paste0("Length of selected genes: ",length(topvars)))
print(paste0("Genes selected: ",paste(topvars, collapse = ",")))

# --- Evaluate Final Model on Test Datasets (test1, test2) ---
RS_list_final <- lapply(df_list, function(x) { # Use df_list to evaluate on train, test1, test2
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")] # Use 'genes' or 'topvars' consistently
  predict(rf_model_final, newdata = df_data)$predicted
})

# Create a list of dataframes instead of combining
results_list <- lapply(names(df_list), function(dataset_name) {
  df <- df_list[[dataset_name]]
  rs <- RS_list_final[[dataset_name]]
  data.frame(
    Sample_Name = rownames(df),
    OS_MONTHS = df$OS_MONTHS,
    OS_STATUS = df$OS_STATUS,
    RS_parthanatos = rs,
    stringsAsFactors = FALSE
  )
})
names(results_list) <- names(df_list)  # Name the list elements

# Save each dataframe as a CSV
lapply(names(results_list), function(dataset_name) {
  fwrite(results_list[[dataset_name]], paste0(tolower(dataset_name), ".csv"), sep = ",")
})

# Keep results_list for further use
combined_df <- results_list  # Store as a named list instead of a single dataframe

################## 2. For Lysosomal Cell Death ######################

rm(list = setdiff(ls(), c("combined_df")))
seed = 42
set.seed(seed)

ntree= 1000
nodesize=5
nodedepth=5

library(dplyr)
library(data.table)
# library(readxl)


library(BiocParallel)
library(parallel)
num_detected_cores <- detectCores(logical = FALSE) - 2 # logical = FALSE gives physical cores
register(MulticoreParam(workers = num_detected_cores))


# RSF #

print("----------------------------------------------")
print("Analysis started for Lysosomal Cell Death: LASSO + RSF")
print("----------------------------------------------")

merged <- data.frame(fread("combined_z_score_clinical_batch.csv", sep=','), row.names=1, check.names=F)

train  <- merged %>% dplyr::sample_frac(0.8) ##dim(train) #743 14558
test1   <- merged[!(rownames(merged) %in% rownames(train)), ]

degs_emtab <- data.frame(fread("DEGs_emtab_0.5.txt", sep = "\t", header = T), row.names=1, check.names = F)
degs_emtab <- row.names(degs_emtab)
degs_tcga <- data.frame(fread("DEGs_TCGA_0.5.txt", sep = "\t", header = T), row.names=1, check.names = F)
degs_tcga <- row.names(degs_tcga)
common_degs <- intersect(degs_emtab, degs_tcga)
genes_to_select <- common_degs

# train <- data.frame(fread("emtab_z_score_clinical.csv", sep=',', header=T), row.names=1, check.names=F)
# # # # train <- train[ (train$OS_MONTHS <= 140.38356164) ,]
# test1 <- data.frame(fread("tcga_z_score_clinical.csv", sep = ",", header = T), row.names=1, check.names=F)
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

gene_set <- read.table("../../gene_lists/lcd_gene_list.txt", sep = ",", header = FALSE, check.names = FALSE)$V1
gene_set <- unique(gene_set)

train <- train[, names(train) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test1 <- test1[, names(test1) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test2 <- test2[, names(test2) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test3 <- test3[, names(test3) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test4 <- test4[, names(test4) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test5 <- test5[, names(test5) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test6 <- test6[, names(test6) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test7 <- test7[, names(test7) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test8 <- test8[, names(test8) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]

library("survival")
library("survminer")

if (!requireNamespace("survivalROC", quietly = TRUE)) install.packages("survivalROC")
if (!requireNamespace("PRROC", quietly = TRUE)) install.packages("PRROC")

# Load required packages
library(survivalROC)   # For time-dependent ROC analysis
library(PRROC)   # For Precision-Recall (PR) AUC

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
  # x <- x[ (x$OS_MONTHS <= 140.38356164) ,]
}, BPPARAM = MulticoreParam(num_detected_cores))


# select the train dataset
train <- df_list$Train ## explicitly select the train dataset

# mention time points for survivalROC and PRROC
time_points <- c(12, 36, 60)


## Gene selection

library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
# install.packages("devtools")
# devtools::install_github("zabore/condsurv")
library(condsurv)
library("glmnet")
genes <- covariates
data <-train[, colnames(train) %in% c(genes, "OS_STATUS", "OS_MONTHS")]
set.seed(seed)

cv.lasso <- cv.glmnet(x = as.matrix(data[, -c((ncol(data) - 1):ncol(data))]),  # Exclude OS_Months and OS_Status
                      y = Surv(data$OS_MONTHS, data$OS_STATUS),
                      family = "cox",
                      nfolds = 10,
                      alpha=1)

# Extract optimal lambda
best_lambda_lasso <- cv.lasso$lambda.min  # Or lambda.1se for more conservative model
# best_lambda_lasso <- cv.lasso$lambda.1se

# Fit the Lasso model with the optimal lambda
lasso_fit <- glmnet(x = as.matrix(data[, -c((ncol(data) - 1):ncol(data))]),  # Exclude OS_Months and OS_Status
                    y = Surv(data$OS_MONTHS, data$OS_STATUS),
                    family = "cox",
                    lambda = best_lambda_lasso,
                    alpha=1)

cf_lasso<- coef(lasso_fit)
cf_lasso_df<- as.data.frame(as.matrix(cf_lasso))
cf_lasso_df_non_zero<- cf_lasso_df[!(cf_lasso_df$s0==0), , drop=F]
dim(cf_lasso_df_non_zero) #10 1
genes <- row.names(cf_lasso_df_non_zero)
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]
library(randomForestSRC)

rf_model_final <- rfsrc(Surv(OS_MONTHS, OS_STATUS) ~ ., data = train_data,
                        importance = T, ntree = ntree, # Use full ntree for final model
                        proximity = T, forest = T,
                        seed = seed, nodesize = nodesize,
                        bootstrap = "by.root", ncores= num_detected_cores, nodedepth = nodedepth) # Use topvars for final model
vs.rf_model <- var.select(object = rf_model_final, conservative = "high")
topvars <- vs.rf_model$topvars
length(topvars)

print(paste0("Length of selected genes: ",length(topvars)))
print(paste0("Genes selected: ",paste(topvars, collapse = ",")))

# --- Evaluate Final Model on Test Datasets (test1, test2) ---
RS_list_final <- lapply(df_list, function(x) { # Use df_list to evaluate on train, test1, test2
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")] # Use 'genes' or 'topvars' consistently
  predict(rf_model_final, newdata = df_data)$predicted
})

# Add RS_LCD to each dataframe in results_list
results_list <- lapply(names(df_list), function(dataset_name) {
  df <- combined_df[[dataset_name]]
  rs <- RS_list_final[[dataset_name]]
  df$RS_LCD <- rs
  return(df)
})
names(results_list) <- names(df_list)

# Save each dataframe as a CSV
lapply(names(results_list), function(dataset_name) {
  fwrite(results_list[[dataset_name]], paste0(tolower(dataset_name), ".csv"), sep = ",")
})

# Update combined_df to remain a list
combined_df <- results_list

########################### 3. For ICD ##############################

rm(list = setdiff(ls(), c("combined_df")))
seed = 42
set.seed(seed)

ntree= 1000
nodesize=5
nodedepth=5

library(dplyr)
library(data.table)
# library(readxl)


library(BiocParallel)
library(parallel)
num_detected_cores <- detectCores(logical = FALSE) - 2 # logical = FALSE gives physical cores
register(MulticoreParam(workers = num_detected_cores))


# LASSO + RSF #
print("----------------------------------------------")
print("Analysis started for ICD: LASSO + RSF")
print("----------------------------------------------")

merged <- data.frame(fread("combined_z_score_clinical_batch.csv", sep=','), row.names=1, check.names=F)

train  <- merged %>% dplyr::sample_frac(0.8) ##dim(train) #743 14558
test1   <- merged[!(rownames(merged) %in% rownames(train)), ]

degs_emtab <- data.frame(fread("DEGs_emtab_0.5.txt", sep = "\t", header = T), row.names=1, check.names = F)
degs_emtab <- row.names(degs_emtab)
degs_tcga <- data.frame(fread("DEGs_TCGA_0.5.txt", sep = "\t", header = T), row.names=1, check.names = F)
degs_tcga <- row.names(degs_tcga)
common_degs <- intersect(degs_emtab, degs_tcga)
genes_to_select <- common_degs

# train <- data.frame(fread("emtab_z_score_clinical.csv", sep=',', header=T), row.names=1, check.names=F)
# # # # train <- train[ (train$OS_MONTHS <= 140.38356164) ,]
# test1 <- data.frame(fread("tcga_z_score_clinical.csv", sep = ",", header = T), row.names=1, check.names=F)
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

gene_set <- read.table("../../gene_lists/ICD_gene_list.txt", sep = ",", header = FALSE, check.names = FALSE)$V1
gene_set <- unique(gene_set)

train <- train[, names(train) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test1 <- test1[, names(test1) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test2 <- test2[, names(test2) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test3 <- test3[, names(test3) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test4 <- test4[, names(test4) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test5 <- test5[, names(test5) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test6 <- test6[, names(test6) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test7 <- test7[, names(test7) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test8 <- test8[, names(test8) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]

library("survival")
library("survminer")

if (!requireNamespace("survivalROC", quietly = TRUE)) install.packages("survivalROC")
if (!requireNamespace("PRROC", quietly = TRUE)) install.packages("PRROC")

# Load required packages
library(survivalROC)   # For time-dependent ROC analysis
library(PRROC)   # For Precision-Recall (PR) AUC

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
  # x <- x[ (x$OS_MONTHS <= 140.38356164) ,]
}, BPPARAM = MulticoreParam(num_detected_cores))


# select the train dataset
train <- df_list$Train ## explicitly select the train dataset

# mention time points for survivalROC and PRROC
time_points <- c(12, 36, 60)

## Gene Selection ##
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
# install.packages("devtools")
# devtools::install_github("zabore/condsurv")
library(condsurv)
library("glmnet")
genes <- covariates
data <-train[, colnames(train) %in% c(genes, "OS_STATUS", "OS_MONTHS")]
set.seed(seed)
cv.lasso <- cv.glmnet(x = as.matrix(data[, -c((ncol(data) - 1):ncol(data))]),  # Exclude OS_Months and OS_Status
                      y = Surv(data$OS_MONTHS, data$OS_STATUS),
                      family = "cox",
                      nfolds = 10,
                      alpha=1)

# Extract optimal lambda
best_lambda_lasso <- cv.lasso$lambda.min  # Or lambda.1se for more conservative model
# best_lambda_lasso <- cv.lasso$lambda.1se

# Fit the Lasso model with the optimal lambda
lasso_fit <- glmnet(x = as.matrix(data[, -c((ncol(data) - 1):ncol(data))]),  # Exclude OS_Months and OS_Status
                    y = Surv(data$OS_MONTHS, data$OS_STATUS),
                    family = "cox",
                    lambda = best_lambda_lasso,
                    alpha=1)

cf_lasso<- coef(lasso_fit)
cf_lasso_df<- as.data.frame(as.matrix(cf_lasso))
cf_lasso_df_non_zero<- cf_lasso_df[!(cf_lasso_df$s0==0), , drop=F]
dim(cf_lasso_df_non_zero) #10 1
genes <- row.names(cf_lasso_df_non_zero)
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

rf_model_final <- rfsrc(Surv(OS_MONTHS, OS_STATUS) ~ ., data = train_data,
                        importance = T, ntree = ntree, # Use full ntree for final model
                        proximity = T, forest = T,
                        seed = seed, nodesize = nodesize,
                        bootstrap = "by.root", ncores= num_detected_cores, nodedepth = nodedepth) # Use topvars for final model
vs.rf_model <- var.select(object = rf_model_final, conservative = "high")
topvars <- vs.rf_model$topvars
length(topvars)

print(paste0("Length of selected genes: ",length(topvars)))
print(paste0("Genes selected: ",paste(topvars, collapse = ",")))

# --- Evaluate Final Model on Test Datasets (test1, test2) ---
RS_list_final <- lapply(df_list, function(x) { # Use df_list to evaluate on train, test1, test2
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")] # Use 'genes' or 'topvars' consistently
  predict(rf_model_final, newdata = df_data)$predicted
})

results_list <- lapply(names(df_list), function(dataset_name) {
  df <- combined_df[[dataset_name]]
  rs <- RS_list_final[[dataset_name]]
  df$RS_ICD <- rs
  return(df)
})
names(results_list) <- names(df_list)

# Save each dataframe as a CSV
lapply(names(results_list), function(dataset_name) {
  fwrite(results_list[[dataset_name]], paste0(tolower(dataset_name), ".csv"), sep = ",")
})

# Update combined_df to remain a list
combined_df <- results_list

########################### 4. For Ferroptosis #############################

rm(list = setdiff(ls(), c("combined_df")))
seed = 42
set.seed(seed)

ntree= 1000
nodesize=5
nodedepth=5

library(dplyr)
library(data.table)
# library(readxl)


library(BiocParallel)
library(parallel)
num_detected_cores <- detectCores(logical = FALSE) - 2 # logical = FALSE gives physical cores
register(MulticoreParam(workers = num_detected_cores))


# LASSO + RSF #

print("----------------------------------------------")
print("Analysis started for Ferroptosis: LASSO + RSF")
print("----------------------------------------------")

merged <- data.frame(fread("combined_z_score_clinical_batch.csv", sep=','), row.names=1, check.names=F)

train  <- merged %>% dplyr::sample_frac(0.8) ##dim(train) #743 14558
test1   <- merged[!(rownames(merged) %in% rownames(train)), ]

degs_emtab <- data.frame(fread("DEGs_emtab_0.5.txt", sep = "\t", header = T), row.names=1, check.names = F)
degs_emtab <- row.names(degs_emtab)
degs_tcga <- data.frame(fread("DEGs_TCGA_0.5.txt", sep = "\t", header = T), row.names=1, check.names = F)
degs_tcga <- row.names(degs_tcga)
common_degs <- intersect(degs_emtab, degs_tcga)
genes_to_select <- common_degs

# train <- data.frame(fread("emtab_z_score_clinical.csv", sep=',', header=T), row.names=1, check.names=F)
# # # # train <- train[ (train$OS_MONTHS <= 140.38356164) ,]
# test1 <- data.frame(fread("tcga_z_score_clinical.csv", sep = ",", header = T), row.names=1, check.names=F)
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

gene_set <- read.table("../../gene_lists/ferroptosis_gene_list.txt", sep = ",", header = FALSE, check.names = FALSE)$V1
gene_set <- unique(gene_set)

train <- train[, names(train) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test1 <- test1[, names(test1) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test2 <- test2[, names(test2) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test3 <- test3[, names(test3) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test4 <- test4[, names(test4) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test5 <- test5[, names(test5) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test6 <- test6[, names(test6) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test7 <- test7[, names(test7) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test8 <- test8[, names(test8) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]

library("survival")
library("survminer")

if (!requireNamespace("survivalROC", quietly = TRUE)) install.packages("survivalROC")
if (!requireNamespace("PRROC", quietly = TRUE)) install.packages("PRROC")

# Load required packages
library(survivalROC)   # For time-dependent ROC analysis
library(PRROC)   # For Precision-Recall (PR) AUC

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
  # x <- x[ (x$OS_MONTHS <= 140.38356164) ,]
}, BPPARAM = MulticoreParam(num_detected_cores))


# select the train dataset
train <- df_list$Train ## explicitly select the train dataset

# mention time points for survivalROC and PRROC
time_points <- c(12, 36, 60)


## Gene Selection ##
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
# install.packages("devtools")
# devtools::install_github("zabore/condsurv")
library(condsurv)
library("glmnet")
genes <- covariates
data <-train[, colnames(train) %in% c(genes, "OS_STATUS", "OS_MONTHS")]
set.seed(seed)

cv.lasso <- cv.glmnet(x = as.matrix(data[, -c((ncol(data) - 1):ncol(data))]),  # Exclude OS_Months and OS_Status
                      y = Surv(data$OS_MONTHS, data$OS_STATUS),
                      family = "cox",
                      nfolds = 10,
                      alpha=1)

# Extract optimal lambda
best_lambda_lasso <- cv.lasso$lambda.min  # Or lambda.1se for more conservative model
# best_lambda_lasso <- cv.lasso$lambda.1se

# Fit the Lasso model with the optimal lambda
lasso_fit <- glmnet(x = as.matrix(data[, -c((ncol(data) - 1):ncol(data))]),  # Exclude OS_Months and OS_Status
                    y = Surv(data$OS_MONTHS, data$OS_STATUS),
                    family = "cox",
                    lambda = best_lambda_lasso,
                    alpha=1)

cf_lasso<- coef(lasso_fit)
cf_lasso_df<- as.data.frame(as.matrix(cf_lasso))
cf_lasso_df_non_zero<- cf_lasso_df[!(cf_lasso_df$s0==0), , drop=F]
dim(cf_lasso_df_non_zero) #10 1
genes <- row.names(cf_lasso_df_non_zero)
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

rf_model_final <- rfsrc(Surv(OS_MONTHS, OS_STATUS) ~ ., data = train_data,
                        importance = T, ntree = ntree, # Use full ntree for final model
                        proximity = T, forest = T,
                        seed = seed, nodesize = nodesize,
                        bootstrap = "by.root", ncores= num_detected_cores, nodedepth = nodedepth) # Use topvars for final model
vs.rf_model <- var.select(object = rf_model_final, conservative = "high")
topvars <- vs.rf_model$topvars
length(topvars)

print(paste0("Length of selected genes: ",length(topvars)))
print(paste0("Genes selected: ",paste(topvars, collapse = ",")))

# --- Evaluate Final Model on Test Datasets (test1, test2) ---
RS_list_final <- lapply(df_list, function(x) { # Use df_list to evaluate on train, test1, test2
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")] # Use 'genes' or 'topvars' consistently
  predict(rf_model_final, newdata = df_data)$predicted
})

results_list <- lapply(names(df_list), function(dataset_name) {
  df <- combined_df[[dataset_name]]
  rs <- RS_list_final[[dataset_name]]
  df$RS_ferroptosis <- rs
  return(df)
})
names(results_list) <- names(df_list)

# Save each dataframe as a CSV
lapply(names(results_list), function(dataset_name) {
  fwrite(results_list[[dataset_name]], paste0(tolower(dataset_name), ".csv"), sep = ",")
})

# Update combined_df to remain a list
combined_df <- results_list

########################## 5. For Necroptosis ##############################

rm(list = setdiff(ls(), c("combined_df")))
seed = 42
set.seed(seed)

ntree= 1000
nodesize=5
nodedepth=5

library(dplyr)
library(data.table)
# library(readxl)


library(BiocParallel)
library(parallel)
num_detected_cores <- detectCores(logical = FALSE) - 2 # logical = FALSE gives physical cores
register(MulticoreParam(workers = num_detected_cores))


# Uni + LASSO + RSF #
print("----------------------------------------------")
print("Analysis started for Necroptosis: LASSO + RSF")
print("----------------------------------------------")
merged <- data.frame(fread("combined_z_score_clinical_batch.csv", sep=','), row.names=1, check.names=F)

train  <- merged %>% dplyr::sample_frac(0.8) ##dim(train) #743 14558
test1   <- merged[!(rownames(merged) %in% rownames(train)), ]

degs_emtab <- data.frame(fread("DEGs_emtab_0.5.txt", sep = "\t", header = T), row.names=1, check.names = F)
degs_emtab <- row.names(degs_emtab)
degs_tcga <- data.frame(fread("DEGs_TCGA_0.5.txt", sep = "\t", header = T), row.names=1, check.names = F)
degs_tcga <- row.names(degs_tcga)
common_degs <- intersect(degs_emtab, degs_tcga)
genes_to_select <- common_degs

# train <- data.frame(fread("emtab_z_score_clinical.csv", sep=',', header=T), row.names=1, check.names=F)
# # # # train <- train[ (train$OS_MONTHS <= 140.38356164) ,]
# test1 <- data.frame(fread("tcga_z_score_clinical.csv", sep = ",", header = T), row.names=1, check.names=F)
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

gene_set <- read.table("../../gene_lists/necroptosis_gene_list.txt", sep = ",", header = FALSE, check.names = FALSE)$V1
gene_set <- unique(gene_set)

train <- train[, names(train) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test1 <- test1[, names(test1) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test2 <- test2[, names(test2) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test3 <- test3[, names(test3) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test4 <- test4[, names(test4) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test5 <- test5[, names(test5) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test6 <- test6[, names(test6) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test7 <- test7[, names(test7) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test8 <- test8[, names(test8) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]

library("survival")
library("survminer")

if (!requireNamespace("survivalROC", quietly = TRUE)) install.packages("survivalROC")
if (!requireNamespace("PRROC", quietly = TRUE)) install.packages("PRROC")

# Load required packages
library(survivalROC)   # For time-dependent ROC analysis
library(PRROC)   # For Precision-Recall (PR) AUC

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
  # x <- x[ (x$OS_MONTHS <= 140.38356164) ,]
}, BPPARAM = MulticoreParam(num_detected_cores))


# select the train dataset
train <- df_list$Train ## explicitly select the train dataset

# mention time points for survivalROC and PRROC
time_points <- c(12, 36, 60)

## Gene Selection ##
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
# install.packages("devtools")
# devtools::install_github("zabore/condsurv")
library(condsurv)
library("glmnet")
genes <- covariates
data <-train[, colnames(train) %in% c(genes, "OS_STATUS", "OS_MONTHS")]
set.seed(seed)

cv.lasso <- cv.glmnet(x = as.matrix(data[, -c((ncol(data) - 1):ncol(data))]),  # Exclude OS_Months and OS_Status
                      y = Surv(data$OS_MONTHS, data$OS_STATUS),
                      family = "cox",
                      nfolds = 10,
                      alpha=1)

# Extract optimal lambda
best_lambda_lasso <- cv.lasso$lambda.min  # Or lambda.1se for more conservative model
# best_lambda_lasso <- cv.lasso$lambda.1se

# Fit the Lasso model with the optimal lambda
lasso_fit <- glmnet(x = as.matrix(data[, -c((ncol(data) - 1):ncol(data))]),  # Exclude OS_Months and OS_Status
                    y = Surv(data$OS_MONTHS, data$OS_STATUS),
                    family = "cox",
                    lambda = best_lambda_lasso,
                    alpha=1)

cf_lasso<- coef(lasso_fit)
cf_lasso_df<- as.data.frame(as.matrix(cf_lasso))
cf_lasso_df_non_zero<- cf_lasso_df[!(cf_lasso_df$s0==0), , drop=F]
dim(cf_lasso_df_non_zero) #10 1
genes <- row.names(cf_lasso_df_non_zero)
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

rf_model_final <- rfsrc(Surv(OS_MONTHS, OS_STATUS) ~ ., data = train_data,
                        importance = T, ntree = ntree, # Use full ntree for final model
                        proximity = T, forest = T,
                        seed = seed, nodesize = nodesize,
                        bootstrap = "by.root", ncores= num_detected_cores, nodedepth = nodedepth) # Use topvars for final model
vs.rf_model <- var.select(object = rf_model_final, conservative = "high")
topvars <- vs.rf_model$topvars
length(topvars)

print(paste0("Length of selected genes: ",length(topvars)))
print(paste0("Genes selected: ",paste(topvars, collapse = ",")))

# --- Evaluate Final Model on Test Datasets (test1, test2) ---
RS_list_final <- lapply(df_list, function(x) { # Use df_list to evaluate on train, test1, test2
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")] # Use 'genes' or 'topvars' consistently
  predict(rf_model_final, newdata = df_data)$predicted
})

results_list <- lapply(names(df_list), function(dataset_name) {
  df <- combined_df[[dataset_name]]
  rs <- RS_list_final[[dataset_name]]
  df$RS_necroptosis <- rs
  return(df)
})
names(results_list) <- names(df_list)

# Save each dataframe as a CSV
lapply(names(results_list), function(dataset_name) {
  fwrite(results_list[[dataset_name]], paste0(tolower(dataset_name), ".csv"), sep = ",")
})

# Update combined_df to remain a list
combined_df <- results_list

########################## 6. For Oxeiptosis ##############################

rm(list = setdiff(ls(), c("combined_df")))
seed = 42
set.seed(seed)

ntree= 1000
nodesize=5
nodedepth=5

library(dplyr)
library(data.table)
# library(readxl)


library(BiocParallel)
library(parallel)
num_detected_cores <- detectCores(logical = FALSE) - 2 # logical = FALSE gives physical cores
register(MulticoreParam(workers = num_detected_cores))


# LASSO + RSF #

print("----------------------------------------------")
print("Analysis started for Oxeiptosis: RSF")
print("----------------------------------------------")

merged <- data.frame(fread("combined_z_score_clinical_batch.csv", sep=','), row.names=1, check.names=F)

train  <- merged %>% dplyr::sample_frac(0.8) ##dim(train) #743 14558
test1   <- merged[!(rownames(merged) %in% rownames(train)), ]

degs_emtab <- data.frame(fread("DEGs_emtab_0.5.txt", sep = "\t", header = T), row.names=1, check.names = F)
degs_emtab <- row.names(degs_emtab)
degs_tcga <- data.frame(fread("DEGs_TCGA_0.5.txt", sep = "\t", header = T), row.names=1, check.names = F)
degs_tcga <- row.names(degs_tcga)
common_degs <- intersect(degs_emtab, degs_tcga)
genes_to_select <- common_degs

# train <- data.frame(fread("emtab_z_score_clinical.csv", sep=',', header=T), row.names=1, check.names=F)
# # # # train <- train[ (train$OS_MONTHS <= 140.38356164) ,]
# test1 <- data.frame(fread("tcga_z_score_clinical.csv", sep = ",", header = T), row.names=1, check.names=F)
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

gene_set <- read.table("../../gene_lists/oxeiptosis_related_genes.txt", sep = ",", header = FALSE, check.names = FALSE)$V1
gene_set <- unique(gene_set)

train <- train[, names(train) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test1 <- test1[, names(test1) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test2 <- test2[, names(test2) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test3 <- test3[, names(test3) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test4 <- test4[, names(test4) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test5 <- test5[, names(test5) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test6 <- test6[, names(test6) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test7 <- test7[, names(test7) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test8 <- test8[, names(test8) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]

library("survival")
library("survminer")

if (!requireNamespace("survivalROC", quietly = TRUE)) install.packages("survivalROC")
if (!requireNamespace("PRROC", quietly = TRUE)) install.packages("PRROC")

# Load required packages
library(survivalROC)   # For time-dependent ROC analysis
library(PRROC)   # For Precision-Recall (PR) AUC

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
  # x <- x[ (x$OS_MONTHS <= 140.38356164) ,]
}, BPPARAM = MulticoreParam(num_detected_cores))


# select the train dataset
train <- df_list$Train ## explicitly select the train dataset

# mention time points for survivalROC and PRROC
time_points <- c(12, 36, 60)

## Gene Selection ##

genes <- covariates
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

rf_model_final <- rfsrc(Surv(OS_MONTHS, OS_STATUS) ~ ., data = train_data,
                        importance = T, ntree = ntree, # Use full ntree for final model
                        proximity = T, forest = T,
                        seed = seed, nodesize = nodesize,
                        bootstrap = "by.root", ncores= num_detected_cores, nodedepth = nodedepth) # Use topvars for final model
vs.rf_model <- var.select(object = rf_model_final, conservative = "high")
topvars <- vs.rf_model$topvars
length(topvars)

print(paste0("Length of selected genes: ",length(topvars)))
print(paste0("Genes selected: ",paste(topvars, collapse = ",")))

# --- Evaluate Final Model on Test Datasets (test1, test2) ---
RS_list_final <- lapply(df_list, function(x) { # Use df_list to evaluate on train, test1, test2
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")] # Use 'genes' or 'topvars' consistently
  predict(rf_model_final, newdata = df_data)$predicted
})

results_list <- lapply(names(df_list), function(dataset_name) {
  df <- combined_df[[dataset_name]]
  rs <- RS_list_final[[dataset_name]]
  df$RS_oxeiptosis <- rs
  return(df)
})
names(results_list) <- names(df_list)

# Save each dataframe as a CSV
lapply(names(results_list), function(dataset_name) {
  fwrite(results_list[[dataset_name]], paste0(tolower(dataset_name), ".csv"), sep = ",")
})

# Update combined_df to remain a list
combined_df <- results_list


########################## 7. For NETosis ##############################

rm(list = setdiff(ls(), c("combined_df")))
seed = 42
set.seed(seed)

ntree= 1000
nodesize=5
nodedepth=5

library(dplyr)
library(data.table)
# library(readxl)


library(BiocParallel)
library(parallel)
num_detected_cores <- detectCores(logical = FALSE) - 2 # logical = FALSE gives physical cores
register(MulticoreParam(workers = num_detected_cores))


# LASSO + RSF #

print("----------------------------------------------")
print("Analysis started for NETosis: LASSO + RSF")
print("----------------------------------------------")

merged <- data.frame(fread("combined_z_score_clinical_batch.csv", sep=','), row.names=1, check.names=F)

train  <- merged %>% dplyr::sample_frac(0.8) ##dim(train) #743 14558
test1   <- merged[!(rownames(merged) %in% rownames(train)), ]

degs_emtab <- data.frame(fread("DEGs_emtab_0.5.txt", sep = "\t", header = T), row.names=1, check.names = F)
degs_emtab <- row.names(degs_emtab)
degs_tcga <- data.frame(fread("DEGs_TCGA_0.5.txt", sep = "\t", header = T), row.names=1, check.names = F)
degs_tcga <- row.names(degs_tcga)
common_degs <- intersect(degs_emtab, degs_tcga)
genes_to_select <- common_degs

# train <- data.frame(fread("emtab_z_score_clinical.csv", sep=',', header=T), row.names=1, check.names=F)
# # # # train <- train[ (train$OS_MONTHS <= 140.38356164) ,]
# test1 <- data.frame(fread("tcga_z_score_clinical.csv", sep = ",", header = T), row.names=1, check.names=F)
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

gene_set <- read.table("../../gene_lists/NETosis_related_genes.txt", sep = ",", header = FALSE, check.names = FALSE)$V1
gene_set <- unique(gene_set)

train <- train[, names(train) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test1 <- test1[, names(test1) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test2 <- test2[, names(test2) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test3 <- test3[, names(test3) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test4 <- test4[, names(test4) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test5 <- test5[, names(test5) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test6 <- test6[, names(test6) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test7 <- test7[, names(test7) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test8 <- test8[, names(test8) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]

library("survival")
library("survminer")

if (!requireNamespace("survivalROC", quietly = TRUE)) install.packages("survivalROC")
if (!requireNamespace("PRROC", quietly = TRUE)) install.packages("PRROC")

# Load required packages
library(survivalROC)   # For time-dependent ROC analysis
library(PRROC)   # For Precision-Recall (PR) AUC

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
  # x <- x[ (x$OS_MONTHS <= 140.38356164) ,]
}, BPPARAM = MulticoreParam(num_detected_cores))


# select the train dataset
train <- df_list$Train ## explicitly select the train dataset

# mention time points for survivalROC and PRROC
time_points <- c(12, 36, 60)


## Gene Selection ##
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
# install.packages("devtools")
# devtools::install_github("zabore/condsurv")
library(condsurv)
library("glmnet")
genes <- covariates
data <-train[, colnames(train) %in% c(genes, "OS_STATUS", "OS_MONTHS")]
set.seed(seed)

cv.lasso <- cv.glmnet(x = as.matrix(data[, -c((ncol(data) - 1):ncol(data))]),  # Exclude OS_Months and OS_Status
                      y = Surv(data$OS_MONTHS, data$OS_STATUS),
                      family = "cox",
                      nfolds = 10,
                      alpha=1)

# Extract optimal lambda
best_lambda_lasso <- cv.lasso$lambda.min  # Or lambda.1se for more conservative model
# best_lambda_lasso <- cv.lasso$lambda.1se

# Fit the Lasso model with the optimal lambda
lasso_fit <- glmnet(x = as.matrix(data[, -c((ncol(data) - 1):ncol(data))]),  # Exclude OS_Months and OS_Status
                    y = Surv(data$OS_MONTHS, data$OS_STATUS),
                    family = "cox",
                    lambda = best_lambda_lasso,
                    alpha=1)

cf_lasso<- coef(lasso_fit)
cf_lasso_df<- as.data.frame(as.matrix(cf_lasso))
cf_lasso_df_non_zero<- cf_lasso_df[!(cf_lasso_df$s0==0), , drop=F]
dim(cf_lasso_df_non_zero) #10 1
genes <- row.names(cf_lasso_df_non_zero)
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

rf_model_final <- rfsrc(Surv(OS_MONTHS, OS_STATUS) ~ ., data = train_data,
                        importance = T, ntree = ntree, # Use full ntree for final model
                        proximity = T, forest = T,
                        seed = seed, nodesize = nodesize,
                        bootstrap = "by.root", ncores= num_detected_cores, nodedepth = nodedepth) # Use topvars for final model
vs.rf_model <- var.select(object = rf_model_final, conservative = "high")
topvars <- vs.rf_model$topvars
length(topvars)

print(paste0("Length of selected genes: ",length(topvars)))
print(paste0("Genes selected: ",paste(topvars, collapse = ",")))

# --- Evaluate Final Model on Test Datasets (test1, test2) ---
RS_list_final <- lapply(df_list, function(x) { # Use df_list to evaluate on train, test1, test2
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")] # Use 'genes' or 'topvars' consistently
  predict(rf_model_final, newdata = df_data)$predicted
})

results_list <- lapply(names(df_list), function(dataset_name) {
  df <- combined_df[[dataset_name]]
  rs <- RS_list_final[[dataset_name]]
  df$RS_NETosis <- rs
  return(df)
})
names(results_list) <- names(df_list)

# Save each dataframe as a CSV
lapply(names(results_list), function(dataset_name) {
  fwrite(results_list[[dataset_name]], paste0(tolower(dataset_name), ".csv"), sep = ",")
})

# Update combined_df to remain a list
combined_df <- results_list

########################## 8. For Entosis ##############################

rm(list = setdiff(ls(), c("combined_df")))
seed = 42
set.seed(seed)

ntree= 1000
nodesize=5
nodedepth=5

library(dplyr)
library(data.table)
# library(readxl)


library(BiocParallel)
library(parallel)
num_detected_cores <- detectCores(logical = FALSE) - 2 # logical = FALSE gives physical cores
register(MulticoreParam(workers = num_detected_cores))


# RSF #
print("----------------------------------------------")
print("Analysis started for Entosis: RSF")
print("----------------------------------------------")

merged <- data.frame(fread("combined_z_score_clinical_batch.csv", sep=','), row.names=1, check.names=F)

train  <- merged %>% dplyr::sample_frac(0.8) ##dim(train) #743 14558
test1   <- merged[!(rownames(merged) %in% rownames(train)), ]

degs_emtab <- data.frame(fread("DEGs_emtab_0.5.txt", sep = "\t", header = T), row.names=1, check.names = F)
degs_emtab <- row.names(degs_emtab)
degs_tcga <- data.frame(fread("DEGs_TCGA_0.5.txt", sep = "\t", header = T), row.names=1, check.names = F)
degs_tcga <- row.names(degs_tcga)
common_degs <- intersect(degs_emtab, degs_tcga)
genes_to_select <- common_degs

# train <- data.frame(fread("emtab_z_score_clinical.csv", sep=',', header=T), row.names=1, check.names=F)
# # # # train <- train[ (train$OS_MONTHS <= 140.38356164) ,]
# test1 <- data.frame(fread("tcga_z_score_clinical.csv", sep = ",", header = T), row.names=1, check.names=F)
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

gene_set <- read.table("../../gene_lists/entosis_core_genes.txt", sep = ",", header = FALSE, check.names = FALSE)$V1
gene_set <- unique(gene_set)

train <- train[, names(train) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test1 <- test1[, names(test1) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test2 <- test2[, names(test2) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test3 <- test3[, names(test3) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test4 <- test4[, names(test4) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test5 <- test5[, names(test5) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test6 <- test6[, names(test6) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test7 <- test7[, names(test7) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test8 <- test8[, names(test8) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]

library("survival")
library("survminer")

if (!requireNamespace("survivalROC", quietly = TRUE)) install.packages("survivalROC")
if (!requireNamespace("PRROC", quietly = TRUE)) install.packages("PRROC")

# Load required packages
library(survivalROC)   # For time-dependent ROC analysis
library(PRROC)   # For Precision-Recall (PR) AUC

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
  # x <- x[ (x$OS_MONTHS <= 140.38356164) ,]
}, BPPARAM = MulticoreParam(num_detected_cores))


# select the train dataset
train <- df_list$Train ## explicitly select the train dataset

# mention time points for survivalROC and PRROC
time_points <- c(12, 36, 60)

## Gene Selection ##

genes <- covariates
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

rf_model_final <- rfsrc(Surv(OS_MONTHS, OS_STATUS) ~ ., data = train_data,
                        importance = T, ntree = ntree, # Use full ntree for final model
                        proximity = T, forest = T,
                        seed = seed, nodesize = nodesize,
                        bootstrap = "by.root", ncores= num_detected_cores, nodedepth = nodedepth) # Use topvars for final model
vs.rf_model <- var.select(object = rf_model_final, conservative = "high")
topvars <- vs.rf_model$topvars
print(paste0("Length of selected genes: ",length(topvars)))
print(paste0("Genes selected: ",paste(topvars, collapse = ",")))


# --- Evaluate Final Model on Test Datasets (test1, test2) ---
RS_list_final <- lapply(df_list, function(x) { # Use df_list to evaluate on train, test1, test2
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")] # Use 'genes' or 'topvars' consistently
  predict(rf_model_final, newdata = df_data)$predicted
})

results_list <- lapply(names(df_list), function(dataset_name) {
  df <- combined_df[[dataset_name]]
  rs <- RS_list_final[[dataset_name]]
  df$RS_entosis <- rs
  return(df)
})
names(results_list) <- names(df_list)

# Save each dataframe as a CSV
lapply(names(results_list), function(dataset_name) {
  fwrite(results_list[[dataset_name]], paste0(tolower(dataset_name), ".csv"), sep = ",")
})

# Update combined_df to remain a list
combined_df <- results_list


########################## 9. For Cuproptosis ##############################

rm(list = setdiff(ls(), c("combined_df")))
seed = 42
set.seed(seed)

ntree= 1000
nodesize=5
nodedepth=5

library(dplyr)
library(data.table)
# library(readxl)


library(BiocParallel)
library(parallel)
num_detected_cores <- detectCores(logical = FALSE) - 2 # logical = FALSE gives physical cores
register(MulticoreParam(workers = num_detected_cores))


# RSF #
print("----------------------------------------------")
print("Analysis started for Cuproptosis: RSF")
print("----------------------------------------------")

merged <- data.frame(fread("combined_z_score_clinical_batch.csv", sep=','), row.names=1, check.names=F)

train  <- merged %>% dplyr::sample_frac(0.8) ##dim(train) #743 14558
test1   <- merged[!(rownames(merged) %in% rownames(train)), ]

degs_emtab <- data.frame(fread("DEGs_emtab_0.5.txt", sep = "\t", header = T), row.names=1, check.names = F)
degs_emtab <- row.names(degs_emtab)
degs_tcga <- data.frame(fread("DEGs_TCGA_0.5.txt", sep = "\t", header = T), row.names=1, check.names = F)
degs_tcga <- row.names(degs_tcga)
common_degs <- intersect(degs_emtab, degs_tcga)
genes_to_select <- common_degs

# train <- data.frame(fread("emtab_z_score_clinical.csv", sep=',', header=T), row.names=1, check.names=F)
# # # # train <- train[ (train$OS_MONTHS <= 140.38356164) ,]
# test1 <- data.frame(fread("tcga_z_score_clinical.csv", sep = ",", header = T), row.names=1, check.names=F)
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

gene_set <- read.table("../../gene_lists/cuproptosis_gene_list.txt", sep = ",", header = FALSE, check.names = FALSE)$V1
gene_set <- unique(gene_set)

train <- train[, names(train) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test1 <- test1[, names(test1) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test2 <- test2[, names(test2) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test3 <- test3[, names(test3) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test4 <- test4[, names(test4) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test5 <- test5[, names(test5) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test6 <- test6[, names(test6) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test7 <- test7[, names(test7) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test8 <- test8[, names(test8) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]

library("survival")
library("survminer")

if (!requireNamespace("survivalROC", quietly = TRUE)) install.packages("survivalROC")
if (!requireNamespace("PRROC", quietly = TRUE)) install.packages("PRROC")

# Load required packages
library(survivalROC)   # For time-dependent ROC analysis
library(PRROC)   # For Precision-Recall (PR) AUC

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
  # x <- x[ (x$OS_MONTHS <= 140.38356164) ,]
}, BPPARAM = MulticoreParam(num_detected_cores))


# select the train dataset
train <- df_list$Train ## explicitly select the train dataset

# mention time points for survivalROC and PRROC
time_points <- c(12, 36, 60)

## Gene Selection ##

genes <- covariates
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

rf_model_final <- rfsrc(Surv(OS_MONTHS, OS_STATUS) ~ ., data = train_data,
                        importance = T, ntree = ntree, # Use full ntree for final model
                        proximity = T, forest = T,
                        seed = seed, nodesize = nodesize,
                        bootstrap = "by.root", ncores= num_detected_cores, nodedepth = nodedepth) # Use topvars for final model
vs.rf_model <- var.select(object = rf_model_final, conservative = "high")
topvars <- vs.rf_model$topvars
length(topvars)
print(paste0("Length of selected genes: ",length(topvars)))
print(paste0("Genes selected: ",paste(topvars, collapse = ",")))

# --- Evaluate Final Model on Test Datasets (test1, test2) ---
RS_list_final <- lapply(df_list, function(x) { # Use df_list to evaluate on train, test1, test2
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")] # Use 'genes' or 'topvars' consistently
  predict(rf_model_final, newdata = df_data)$predicted
})

results_list <- lapply(names(df_list), function(dataset_name) {
  df <- combined_df[[dataset_name]]
  rs <- RS_list_final[[dataset_name]]
  df$RS_cuproptosis <- rs
  return(df)
})
names(results_list) <- names(df_list)

# Save each dataframe as a CSV
lapply(names(results_list), function(dataset_name) {
  fwrite(results_list[[dataset_name]], paste0(tolower(dataset_name), ".csv"), sep = ",")
})

# Update combined_df to remain a list
combined_df <- results_list


########################## 10. For Pyroptosis ##############################

rm(list = setdiff(ls(), c("combined_df")))
seed = 42
set.seed(seed)

ntree= 1000
nodesize=5
nodedepth=5

library(dplyr)
library(data.table)
# library(readxl)


library(BiocParallel)
library(parallel)
num_detected_cores <- detectCores(logical = FALSE) - 2 # logical = FALSE gives physical cores
register(MulticoreParam(workers = num_detected_cores))


# RSF #
print("----------------------------------------------")
print("Analysis started for Pyroptosis: RSF")
print("----------------------------------------------")

merged <- data.frame(fread("combined_z_score_clinical_batch.csv", sep=','), row.names=1, check.names=F)

train  <- merged %>% dplyr::sample_frac(0.8) ##dim(train) #743 14558
test1   <- merged[!(rownames(merged) %in% rownames(train)), ]

degs_emtab <- data.frame(fread("DEGs_emtab_0.5.txt", sep = "\t", header = T), row.names=1, check.names = F)
degs_emtab <- row.names(degs_emtab)
degs_tcga <- data.frame(fread("DEGs_TCGA_0.5.txt", sep = "\t", header = T), row.names=1, check.names = F)
degs_tcga <- row.names(degs_tcga)
common_degs <- intersect(degs_emtab, degs_tcga)
genes_to_select <- common_degs

# train <- data.frame(fread("emtab_z_score_clinical.csv", sep=',', header=T), row.names=1, check.names=F)
# # # # train <- train[ (train$OS_MONTHS <= 140.38356164) ,]
# test1 <- data.frame(fread("tcga_z_score_clinical.csv", sep = ",", header = T), row.names=1, check.names=F)
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

gene_set <- read.table("../../gene_lists/pyroptosis_gene_list.txt", sep = ",", header = FALSE, check.names = FALSE)$V1
gene_set <- unique(gene_set)

train <- train[, names(train) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test1 <- test1[, names(test1) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test2 <- test2[, names(test2) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test3 <- test3[, names(test3) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test4 <- test4[, names(test4) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test5 <- test5[, names(test5) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test6 <- test6[, names(test6) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test7 <- test7[, names(test7) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test8 <- test8[, names(test8) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]

library("survival")
library("survminer")

if (!requireNamespace("survivalROC", quietly = TRUE)) install.packages("survivalROC")
if (!requireNamespace("PRROC", quietly = TRUE)) install.packages("PRROC")

# Load required packages
library(survivalROC)   # For time-dependent ROC analysis
library(PRROC)   # For Precision-Recall (PR) AUC

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
  # x <- x[ (x$OS_MONTHS <= 140.38356164) ,]
}, BPPARAM = MulticoreParam(num_detected_cores))


# select the train dataset
train <- df_list$Train ## explicitly select the train dataset

# mention time points for survivalROC and PRROC
time_points <- c(12, 36, 60)

## Gene Selection ##

genes <- covariates
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

rf_model_final <- rfsrc(Surv(OS_MONTHS, OS_STATUS) ~ ., data = train_data,
                        importance = T, ntree = ntree, # Use full ntree for final model
                        proximity = T, forest = T,
                        seed = seed, nodesize = nodesize,
                        bootstrap = "by.root", ncores= num_detected_cores, nodedepth = nodedepth) # Use topvars for final model
vs.rf_model <- var.select(object = rf_model_final, conservative = "high")
topvars <- vs.rf_model$topvars
length(topvars)

print(paste0("Length of selected genes: ",length(topvars)))
print(paste0("Genes selected: ",paste(topvars, collapse = ",")))

# --- Evaluate Final Model on Test Datasets (test1, test2) ---
RS_list_final <- lapply(df_list, function(x) { # Use df_list to evaluate on train, test1, test2
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")] # Use 'genes' or 'topvars' consistently
  predict(rf_model_final, newdata = df_data)$predicted
})

results_list <- lapply(names(df_list), function(dataset_name) {
  df <- combined_df[[dataset_name]]
  rs <- RS_list_final[[dataset_name]]
  df$RS_pyroptosis <- rs
  return(df)
})
names(results_list) <- names(df_list)

# Save each dataframe as a CSV
lapply(names(results_list), function(dataset_name) {
  fwrite(results_list[[dataset_name]], paste0(tolower(dataset_name), ".csv"), sep = ",")
})

# Update combined_df to remain a list
combined_df <- results_list


########################## 11. For Autophagy ##############################

rm(list = setdiff(ls(), c("combined_df")))
seed = 42
set.seed(seed)

ntree= 1000
nodesize=5
nodedepth=5

library(dplyr)
library(data.table)
# library(readxl)


library(BiocParallel)
library(parallel)
num_detected_cores <- detectCores(logical = FALSE) - 2 # logical = FALSE gives physical cores
register(MulticoreParam(workers = num_detected_cores))


# LASSO + RSF #

print("----------------------------------------------")
print("Analysis started for Autophagy: LASSO + RSF")
print("----------------------------------------------")

merged <- data.frame(fread("combined_z_score_clinical_batch.csv", sep=','), row.names=1, check.names=F)

train  <- merged %>% dplyr::sample_frac(0.8) ##dim(train) #743 14558
test1   <- merged[!(rownames(merged) %in% rownames(train)), ]

degs_emtab <- data.frame(fread("DEGs_emtab_0.5.txt", sep = "\t", header = T), row.names=1, check.names = F)
degs_emtab <- row.names(degs_emtab)
degs_tcga <- data.frame(fread("DEGs_TCGA_0.5.txt", sep = "\t", header = T), row.names=1, check.names = F)
degs_tcga <- row.names(degs_tcga)
common_degs <- intersect(degs_emtab, degs_tcga)
genes_to_select <- common_degs

# train <- data.frame(fread("emtab_z_score_clinical.csv", sep=',', header=T), row.names=1, check.names=F)
# # # # train <- train[ (train$OS_MONTHS <= 140.38356164) ,]
# test1 <- data.frame(fread("tcga_z_score_clinical.csv", sep = ",", header = T), row.names=1, check.names=F)
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

gene_set <- read.table("../../gene_lists/autophagy_gene_list.txt", sep = ",", header = FALSE, check.names = FALSE)$V1
gene_set <- unique(gene_set)

train <- train[, names(train) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test1 <- test1[, names(test1) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test2 <- test2[, names(test2) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test3 <- test3[, names(test3) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test4 <- test4[, names(test4) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test5 <- test5[, names(test5) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test6 <- test6[, names(test6) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test7 <- test7[, names(test7) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test8 <- test8[, names(test8) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]

library("survival")
library("survminer")

if (!requireNamespace("survivalROC", quietly = TRUE)) install.packages("survivalROC")
if (!requireNamespace("PRROC", quietly = TRUE)) install.packages("PRROC")

# Load required packages
library(survivalROC)   # For time-dependent ROC analysis
library(PRROC)   # For Precision-Recall (PR) AUC

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
  # x <- x[ (x$OS_MONTHS <= 140.38356164) ,]
}, BPPARAM = MulticoreParam(num_detected_cores))


# select the train dataset
train <- df_list$Train ## explicitly select the train dataset

# mention time points for survivalROC and PRROC
time_points <- c(12, 36, 60)

## Gene Selection ##
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
# install.packages("devtools")
# devtools::install_github("zabore/condsurv")
library(condsurv)
library("glmnet")
genes <- covariates
data <-train[, colnames(train) %in% c(genes, "OS_STATUS", "OS_MONTHS")]
set.seed(seed)

cv.lasso <- cv.glmnet(x = as.matrix(data[, -c((ncol(data) - 1):ncol(data))]),  # Exclude OS_Months and OS_Status
                      y = Surv(data$OS_MONTHS, data$OS_STATUS),
                      family = "cox",
                      nfolds = 10,
                      alpha=1)

# Extract optimal lambda
best_lambda_lasso <- cv.lasso$lambda.min  # Or lambda.1se for more conservative model
# best_lambda_lasso <- cv.lasso$lambda.1se

# Fit the Lasso model with the optimal lambda
lasso_fit <- glmnet(x = as.matrix(data[, -c((ncol(data) - 1):ncol(data))]),  # Exclude OS_Months and OS_Status
                    y = Surv(data$OS_MONTHS, data$OS_STATUS),
                    family = "cox",
                    lambda = best_lambda_lasso,
                    alpha=1)

cf_lasso<- coef(lasso_fit)
cf_lasso_df<- as.data.frame(as.matrix(cf_lasso))
cf_lasso_df_non_zero<- cf_lasso_df[!(cf_lasso_df$s0==0), , drop=F]
dim(cf_lasso_df_non_zero) #10 1
genes <- row.names(cf_lasso_df_non_zero)
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

rf_model_final <- rfsrc(Surv(OS_MONTHS, OS_STATUS) ~ ., data = train_data,
                        importance = T, ntree = ntree, # Use full ntree for final model
                        proximity = T, forest = T,
                        seed = seed, nodesize = nodesize,
                        bootstrap = "by.root", ncores= num_detected_cores, nodedepth = nodedepth) # Use topvars for final model
vs.rf_model <- var.select(object = rf_model_final, conservative = "high")
topvars <- vs.rf_model$topvars
length(topvars)

print(paste0("Length of selected genes: ",length(topvars)))
print(paste0("Genes selected: ",paste(topvars, collapse = ",")))

# --- Evaluate Final Model on Test Datasets (test1, test2) ---
RS_list_final <- lapply(df_list, function(x) { # Use df_list to evaluate on train, test1, test2
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")] # Use 'genes' or 'topvars' consistently
  predict(rf_model_final, newdata = df_data)$predicted
})

results_list <- lapply(names(df_list), function(dataset_name) {
  df <- combined_df[[dataset_name]]
  rs <- RS_list_final[[dataset_name]]
  df$RS_autophagy <- rs
  return(df)
})
names(results_list) <- names(df_list)

# Save each dataframe as a CSV
lapply(names(results_list), function(dataset_name) {
  fwrite(results_list[[dataset_name]], paste0(tolower(dataset_name), ".csv"), sep = ",")
})

# Update combined_df to remain a list
combined_df <- results_list


########################## 12. For Anoikis ##############################

rm(list = setdiff(ls(), c("combined_df")))
seed = 42
set.seed(seed)

ntree= 1000
nodesize=5
nodedepth=5

library(dplyr)
library(data.table)
# library(readxl)


library(BiocParallel)
library(parallel)
num_detected_cores <- detectCores(logical = FALSE) - 2 # logical = FALSE gives physical cores
register(MulticoreParam(workers = num_detected_cores))


# LASSO + RSF #

print("----------------------------------------------")
print("Analysis started for Anoikis: RSF")
print("----------------------------------------------")

merged <- data.frame(fread("combined_z_score_clinical_batch.csv", sep=','), row.names=1, check.names=F)

train  <- merged %>% dplyr::sample_frac(0.8) ##dim(train) #743 14558
test1   <- merged[!(rownames(merged) %in% rownames(train)), ]

degs_emtab <- data.frame(fread("DEGs_emtab_0.5.txt", sep = "\t", header = T), row.names=1, check.names = F)
degs_emtab <- row.names(degs_emtab)
degs_tcga <- data.frame(fread("DEGs_TCGA_0.5.txt", sep = "\t", header = T), row.names=1, check.names = F)
degs_tcga <- row.names(degs_tcga)
common_degs <- intersect(degs_emtab, degs_tcga)
genes_to_select <- common_degs

# train <- data.frame(fread("emtab_z_score_clinical.csv", sep=',', header=T), row.names=1, check.names=F)
# # # # train <- train[ (train$OS_MONTHS <= 140.38356164) ,]
# test1 <- data.frame(fread("tcga_z_score_clinical.csv", sep = ",", header = T), row.names=1, check.names=F)
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

gene_set <- read.table("../../gene_lists/anoikis_gene_list.txt", sep = ",", header = FALSE, check.names = FALSE)$V1
gene_set <- unique(gene_set)

train <- train[, names(train) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test1 <- test1[, names(test1) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test2 <- test2[, names(test2) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test3 <- test3[, names(test3) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test4 <- test4[, names(test4) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test5 <- test5[, names(test5) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test6 <- test6[, names(test6) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test7 <- test7[, names(test7) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test8 <- test8[, names(test8) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]

library("survival")
library("survminer")

if (!requireNamespace("survivalROC", quietly = TRUE)) install.packages("survivalROC")
if (!requireNamespace("PRROC", quietly = TRUE)) install.packages("PRROC")

# Load required packages
library(survivalROC)   # For time-dependent ROC analysis
library(PRROC)   # For Precision-Recall (PR) AUC

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
  # x <- x[ (x$OS_MONTHS <= 140.38356164) ,]
}, BPPARAM = MulticoreParam(num_detected_cores))


# select the train dataset
train <- df_list$Train ## explicitly select the train dataset

# mention time points for survivalROC and PRROC
time_points <- c(12, 36, 60)

## Gene Selection ##

genes <- covariates
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

rf_model_final <- rfsrc(Surv(OS_MONTHS, OS_STATUS) ~ ., data = train_data,
                        importance = T, ntree = ntree, # Use full ntree for final model
                        proximity = T, forest = T,
                        seed = seed, nodesize = nodesize,
                        bootstrap = "by.root", ncores= num_detected_cores, nodedepth = nodedepth) # Use topvars for final model
vs.rf_model <- var.select(object = rf_model_final, conservative = "high")
topvars <- vs.rf_model$topvars
length(topvars)

print(paste0("Length of selected genes: ",length(topvars)))
print(paste0("Genes selected: ",paste(topvars, collapse = ",")))

# --- Evaluate Final Model on Test Datasets (test1, test2) ---
RS_list_final <- lapply(df_list, function(x) { # Use df_list to evaluate on train, test1, test2
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")] # Use 'genes' or 'topvars' consistently
  predict(rf_model_final, newdata = df_data)$predicted
})

results_list <- lapply(names(df_list), function(dataset_name) {
  df <- combined_df[[dataset_name]]
  rs <- RS_list_final[[dataset_name]]
  df$RS_anoikis <- rs
  return(df)
})
names(results_list) <- names(df_list)

# Save each dataframe as a CSV
lapply(names(results_list), function(dataset_name) {
  fwrite(results_list[[dataset_name]], paste0(tolower(dataset_name), ".csv"), sep = ",")
})

# Update combined_df to remain a list
combined_df <- results_list

########################## 13. For Autosis ##############################

rm(list = setdiff(ls(), c("combined_df")))
seed = 42
set.seed(seed)

ntree= 1000
nodesize=5
nodedepth=5

library(dplyr)
library(data.table)
# library(readxl)


library(BiocParallel)
library(parallel)
num_detected_cores <- detectCores(logical = FALSE) - 2 # logical = FALSE gives physical cores
register(MulticoreParam(workers = num_detected_cores))


# LASSO + RSF #

print("----------------------------------------------")
print("Analysis started for Autosis: LASSO + RSF")
print("----------------------------------------------")

merged <- data.frame(fread("combined_z_score_clinical_batch.csv", sep=','), row.names=1, check.names=F)

train  <- merged %>% dplyr::sample_frac(0.8) ##dim(train) #743 14558
test1   <- merged[!(rownames(merged) %in% rownames(train)), ]

degs_emtab <- data.frame(fread("DEGs_emtab_0.5.txt", sep = "\t", header = T), row.names=1, check.names = F)
degs_emtab <- row.names(degs_emtab)
degs_tcga <- data.frame(fread("DEGs_TCGA_0.5.txt", sep = "\t", header = T), row.names=1, check.names = F)
degs_tcga <- row.names(degs_tcga)
common_degs <- intersect(degs_emtab, degs_tcga)
genes_to_select <- common_degs

# train <- data.frame(fread("emtab_z_score_clinical.csv", sep=',', header=T), row.names=1, check.names=F)
# # # # train <- train[ (train$OS_MONTHS <= 140.38356164) ,]
# test1 <- data.frame(fread("tcga_z_score_clinical.csv", sep = ",", header = T), row.names=1, check.names=F)
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

gene_set <- read.table("../../gene_lists/autosis_gene_list.txt", sep = ",", header = FALSE, check.names = FALSE)$V1
gene_set <- unique(gene_set)

train <- train[, names(train) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test1 <- test1[, names(test1) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test2 <- test2[, names(test2) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test3 <- test3[, names(test3) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test4 <- test4[, names(test4) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test5 <- test5[, names(test5) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test6 <- test6[, names(test6) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test7 <- test7[, names(test7) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test8 <- test8[, names(test8) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]

library("survival")
library("survminer")

if (!requireNamespace("survivalROC", quietly = TRUE)) install.packages("survivalROC")
if (!requireNamespace("PRROC", quietly = TRUE)) install.packages("PRROC")

# Load required packages
library(survivalROC)   # For time-dependent ROC analysis
library(PRROC)   # For Precision-Recall (PR) AUC

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
  # x <- x[ (x$OS_MONTHS <= 140.38356164) ,]
}, BPPARAM = MulticoreParam(num_detected_cores))


# select the train dataset
train <- df_list$Train ## explicitly select the train dataset

# mention time points for survivalROC and PRROC
time_points <- c(12, 36, 60)

## Gene Selection ##
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
# install.packages("devtools")
# devtools::install_github("zabore/condsurv")
library(condsurv)
library("glmnet")
genes <- covariates
data <-train[, colnames(train) %in% c(genes, "OS_STATUS", "OS_MONTHS")]
set.seed(seed)

cv.lasso <- cv.glmnet(x = as.matrix(data[, -c((ncol(data) - 1):ncol(data))]),  # Exclude OS_Months and OS_Status
                      y = Surv(data$OS_MONTHS, data$OS_STATUS),
                      family = "cox",
                      nfolds = 10,
                      alpha=1)

# Extract optimal lambda
best_lambda_lasso <- cv.lasso$lambda.min  # Or lambda.1se for more conservative model
# best_lambda_lasso <- cv.lasso$lambda.1se

# Fit the Lasso model with the optimal lambda
lasso_fit <- glmnet(x = as.matrix(data[, -c((ncol(data) - 1):ncol(data))]),  # Exclude OS_Months and OS_Status
                    y = Surv(data$OS_MONTHS, data$OS_STATUS),
                    family = "cox",
                    lambda = best_lambda_lasso,
                    alpha=1)

cf_lasso<- coef(lasso_fit)
cf_lasso_df<- as.data.frame(as.matrix(cf_lasso))
cf_lasso_df_non_zero<- cf_lasso_df[!(cf_lasso_df$s0==0), , drop=F]
dim(cf_lasso_df_non_zero) #10 1
genes <- row.names(cf_lasso_df_non_zero)
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

rf_model_final <- rfsrc(Surv(OS_MONTHS, OS_STATUS) ~ ., data = train_data,
                        importance = T, ntree = ntree, # Use full ntree for final model
                        proximity = T, forest = T,
                        seed = seed, nodesize = nodesize,
                        bootstrap = "by.root", ncores= num_detected_cores, nodedepth = nodedepth) # Use topvars for final model
vs.rf_model <- var.select(object = rf_model_final, conservative = "high")
topvars <- vs.rf_model$topvars
length(topvars)

print(paste0("Length of selected genes: ",length(topvars)))
print(paste0("Genes selected: ",paste(topvars, collapse = ",")))

# --- Evaluate Final Model on Test Datasets (test1, test2) ---
RS_list_final <- lapply(df_list, function(x) { # Use df_list to evaluate on train, test1, test2
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")] # Use 'genes' or 'topvars' consistently
  predict(rf_model_final, newdata = df_data)$predicted
})

results_list <- lapply(names(df_list), function(dataset_name) {
  df <- combined_df[[dataset_name]]
  rs <- RS_list_final[[dataset_name]]
  df$RS_autosis <- rs
  return(df)
})
names(results_list) <- names(df_list)

# Save each dataframe as a CSV
lapply(names(results_list), function(dataset_name) {
  fwrite(results_list[[dataset_name]], paste0(tolower(dataset_name), ".csv"), sep = ",")
})

# Update combined_df to remain a list
combined_df <- results_list


################################# SAVE RESULTS #######################################
print("Saving final files")
# fwrite(combined_df, "All_RS_models.csv", sep=',', row.names=F)


########################## Combined Risk Score Analysis ##############################
library(survcomp)

set.seed(seed)
print("----------------------------------------------")
print("Analysis started for Combined Risk Score Models")
print("----------------------------------------------")

# Load the saved datasets
datasets <- c("train", "test", "gse39582", "gse161158", "tcga_brca", "tcga_luad", "tcga_skcm", "tcga_gbm", "tcga_paad")
data_list <- lapply(datasets, function(ds) {
  df <- data.frame(fread(paste0(ds, ".csv"), sep = ","), row.names=1, check.names = F)
  df$OS_STATUS <- as.numeric(df$OS_STATUS)
  df$OS_MONTHS <- as.numeric(df$OS_MONTHS)
  return(df)
})
names(data_list) <- datasets

# Define all risk score columns
risk_scores <- c("RS_parthanatos", "RS_LCD", "RS_ICD", "RS_ferroptosis",
                 "RS_necroptosis", "RS_oxeiptosis", "RS_NETosis", "RS_entosis",
                 "RS_cuproptosis", "RS_pyroptosis", "RS_autophagy", "RS_anoikis", "RS_autosis")

# risk_scores <- c("RS_LCD", "RS_ferroptosis", "RCD_ICD",
#                  "RS_cuproptosis", "RS_pyroptosis", "RS_necroptosis", "RS_NETosis")

# Generate combinations of 2, 3, and 4 risk scores
combinations <- unlist(lapply(2:4, function(k) {
  combn(risk_scores, k, simplify = FALSE)
}), recursive = FALSE)

# Function to calculate C-index for a logistic regression model
calculate_cindex_bootstrap <- function(data, formula_str, n_boot = 1000) {
  set.seed(seed)  # For reproducibility
  fit <- glm(as.formula(formula_str), data = data, family = binomial())
  risk_score <- predict(fit, type = "response")
  
  # Base c-index
  base_cindex <- concordance.index(
    x = risk_score,
    surv.time = data$OS_MONTHS,
    surv.event = data$OS_STATUS,
    method = "noether"
  )$c.index
  
  # Bootstrapping
  boot_cindices <- replicate(n_boot, {
    idx <- sample(seq_len(nrow(data)), replace = TRUE)
    boot_data <- data[idx, ]
    boot_fit <- glm(as.formula(formula_str), data = boot_data, family = binomial())
    boot_score <- predict(boot_fit, type = "response")
    concordance.index(
      x = boot_score,
      surv.time = boot_data$OS_MONTHS,
      surv.event = boot_data$OS_STATUS,
      method = "noether"
    )$c.index
  })
  
  ci_lower <- quantile(boot_cindices, 0.025, na.rm = TRUE)
  ci_upper <- quantile(boot_cindices, 0.975, na.rm = TRUE)
  
  return(list(cindex = base_cindex, lower = ci_lower, upper = ci_upper))
}

# Evaluate each combination
results <- lapply(combinations, function(combo) {
  formula_str <- paste("OS_STATUS ~", paste(combo, collapse = " + "))
  
  cindex_info <- lapply(names(data_list), function(ds) {
    data <- data_list[[ds]]
    if (all(combo %in% colnames(data))) {
      calculate_cindex_bootstrap(data, formula_str)
    } else {
      list(cindex = NA, lower = NA, upper = NA)
    }
  })
  names(cindex_info) <- names(data_list)
  
  data.frame(
    Combination = paste(combo, collapse = ","),
    Train_Cindex = cindex_info$train$cindex,
    Train_CI_Lower = cindex_info$train$lower,
    Train_CI_Upper = cindex_info$train$upper,
    
    Test_Cindex = cindex_info$test$cindex,
    Test_CI_Lower = cindex_info$test$lower,
    Test_CI_Upper = cindex_info$test$upper,
    
    gse39582_Cindex = cindex_info$gse39582$cindex,
    gse39582_CI_Lower = cindex_info$gse39582$lower,
    gse39582_CI_Upper = cindex_info$gse39582$upper,
    
    gse161158_Cindex = cindex_info$gse161158$cindex,
    gse161158_CI_Lower = cindex_info$gse161158$lower,
    gse161158_CI_Upper = cindex_info$gse161158$upper,
    
    tcga_brca_Cindex = cindex_info$tcga_brca$cindex,
    tcga_brca_CI_Lower = cindex_info$tcga_brca$lower,
    tcga_brca_CI_Upper = cindex_info$tcga_brca$upper,
    
    tcga_luad_Cindex = cindex_info$tcga_luad$cindex,
    tcga_luad_CI_Lower = cindex_info$tcga_luad$lower,
    tcga_luad_CI_Upper = cindex_info$tcga_luad$upper,
    
    tcga_skcm_Cindex = cindex_info$tcga_skcm$cindex,
    tcga_skcm_CI_Lower = cindex_info$tcga_skcm$lower,
    tcga_skcm_CI_Upper = cindex_info$tcga_skcm$upper,
    
    tcga_gbm_Cindex = cindex_info$tcga_gbm$cindex,
    tcga_gbm_CI_Lower = cindex_info$tcga_gbm$lower,
    tcga_gbm_CI_Upper = cindex_info$tcga_gbm$upper,
    
    tcga_paad_Cindex = cindex_info$tcga_paad$cindex,
    tcga_paad_CI_Lower = cindex_info$tcga_paad$lower,
    tcga_paad_CI_Upper = cindex_info$tcga_paad$upper,
    
    
    Avg_Cindex = mean(sapply(cindex_info, function(x) x$cindex), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})

# Combine results into a data frame
results_df_cindex <- do.call(rbind, results)

# Sort by average C-index (descending)
results_df_cindex <- results_df_cindex[order(-results_df_cindex$Avg_Cindex), ]

# Save results
fwrite(results_df_cindex, "combined_risk_score_cindices_2.csv", sep = ",")

# Print the top combinations
print("Top 5 Combined Risk Score Models by Average C-index:")
print(head(results_df_cindex, 5))
print("Calculated and saved the best combination of C-indices")

#### Plot RCDI C-indices for the top3 combinations for all the 4 datasets ####

print("----------------------")
print("Plotting C-indices for top 3 combinations .....")
print("----------------------")

library(ggplot2)
library(dplyr)
library(tidyr)

# Extract top 3 combinations
top3 <- results_df_cindex %>%
  arrange(desc(Avg_Cindex)) %>%
  slice(1:3)

top3_long <- top3 %>%
  select(Combination,
         Train_Cindex, Train_CI_Lower, Train_CI_Upper,
         Test_Cindex, Test_CI_Lower, Test_CI_Upper,
         Ind1_Cindex, Ind1_CI_Lower, Ind1_CI_Upper,
         Ind2_Cindex, Ind2_CI_Lower, Ind2_CI_Upper) %>%
  pivot_longer(
    cols = -Combination,
    names_to = c("Dataset", ".value"),
    names_pattern = "(Train|Test|Ind1|Ind2)_(Cindex|CI_Lower|CI_Upper)"
  ) %>%
  rename(Cindex = Cindex, Lower = CI_Lower, Upper = CI_Upper)


# Set dataset factor levels in desired order
top3_long$Dataset <- factor(top3_long$Dataset, levels = c("Train", "Test", "Ind1", "Ind2"))

# Calculate mean C-index per combination for reference line & annotation
mean_lines <- top3_long %>%
  group_by(Combination) %>%
  summarise(mean_cindex = mean(Cindex, na.rm = TRUE))

# Merge for annotation
top3_long <- left_join(top3_long, mean_lines, by = "Combination")

# Change x axis labels to be more readable
top3_long$Combination <- gsub("RS_", "", top3_long$Combination)
top3_long$Combination <- gsub(",", " + ", top3_long$Combination)
top3_long$Combination <- gsub(" ", "\n", top3_long$Combination)

# capitalize the first letter of each word
top3_long$Combination <- tools::toTitleCase(top3_long$Combination)

# Reorder the x-axis by descending mean C-index values
top3_long$Combination <- reorder(top3_long$Combination, -top3_long$mean_cindex)

# Create the bar plot
p <- ggplot(top3_long, aes(x = Combination, y = Cindex, fill = Dataset)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), 
                width = 0.2, position = position_dodge(width = 0.8)) +
  geom_text(aes(label = round(Cindex, 3)), 
            position = position_dodge(width = 0.8), 
            vjust = -0.6, size = 3.5, family = "Helvetica") +
  stat_summary(aes(group = Combination), fun = mean, geom = "crossbar",
               width = 0.5, color = "black", fatten = 2,
               position = position_nudge(x = 0),
               show.legend = FALSE) +
  stat_summary(
    aes(group = Combination, label = after_stat(paste("Mean:", round(y, 3)))), 
    fun = mean, geom = "text", vjust = -2, size = 4, 
    position = position_nudge(x = 0), family = "Helvetica"
  ) +
  scale_fill_manual(values = c("Train" = "#1b9e77", 
                               "Test" = "#d95f02", 
                               "Ind1" = "#7570b3", 
                               "Ind2" = "#e7298a")) +
  labs(title = "Top 3 RCDI Models Across Datasets",
       y = "C-index",
       x = "Combination",
       fill = "Dataset") +
  theme_bw(base_family = "Helvetica") +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 10, face = "bold", family = "Helvetica"),
    axis.text.y = element_text(size = 10, hjust = 1, face = "bold", family = "Helvetica"),
    axis.title = element_text(size = 12, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.5, 0.95))
print(p)

# Save the plot
tiff("Top3_RCDI_C_index.tiff", units = "in", width = 7, height = 5, res = 600, compression = "lzw")
print(p)
dev.off()

print("----------------------")
print("Plotted C-indices for top 3 combinations successfully")
print("----------------------")

# Step 2: Compute RCDI for the best combination
# best_combination <- results_df_cindex$Combination[1]  # Extract the best combination
best_combination <- results_df_cindex$Combination[563]  # Extract the Fe,NET,Py,autosis combination
best_combo_scores <- unlist(strsplit(best_combination, ","))  # Split into individual pathways
best_combo_scores <- as.vector(best_combo_scores)

# Logistic regression formula for the best combination
formula_str <- paste("OS_STATUS ~", paste(best_combo_scores, collapse = " + "))

# Function to calculate RCDI for all datasets
calculate_RCDI <- function(data, formula_str) {
  # Ensure all required risk scores exist in the dataset
  required_cols <- c(best_combo_scores, "OS_MONTHS", "OS_STATUS")
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(paste("Dataset is missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  # Fit logistic regression
  fit <- glm(as.formula(formula_str), data = data, family = binomial())
  # Use predict() to calculate RCDI (Combined Risk Score)
  data$RCDI <- predict(fit, type = "response")  # Logistic regression predicted probabilities
  # Return the dataset with RCDI
  return(data)
}

# Calculate RCDI for each dataset
RCDI_results <- lapply(names(data_list), function(ds) {
  data <- data_list[[ds]]
  result <- calculate_RCDI(data, formula_str)
  # Select relevant columns for saving
  result <- result[, c("OS_MONTHS", "OS_STATUS", best_combo_scores, "RCDI"), drop = FALSE]
  return(result)
})

names(RCDI_results) <- names(data_list)

# Save RCDI results for each dataset
lapply(names(RCDI_results), function(ds) {
  result <- RCDI_results[[ds]]
  output_file <- paste0(ds, "_RCDI_combined_risk_scores_2.csv")
  fwrite(result, output_file, sep = ",", row.names = T)
  print(paste("Saved RCDI results for", ds, "to", output_file))
})

################## Barplot of C-indices of RCDI and F, P, Ne, NE RCD pathways ######################

print("---------------------------------------------------")
print("Creating Barplot of C-indices of RCDI and F, Py, NET, Autosis RCD pathways")
print("---------------------------------------------------")

library(ggplot2)
library(readxl)

# Step 1: Read the Excel file
file_path <- "algorithms_chosen.xlsx"
data <- read_excel(file_path)

# Step 2: Filter the relevant pathways and the best CCDI pathway
selected_pathways <- c("Ferroptosis", "Pyroptosis", "Autosis", "NETosis")
filtered_data <- data[data$Gene_list %in% selected_pathways, ]
filtered_data <- filtered_data[, c("Gene_list","C-index Train", "Test", "Ind (GSE39582)", "Ind (GSE161158)")]
colnames(filtered_data) <- c("Gene_list","Train", "Test", "GSE39582", "GSE161158")

# Add a row for the best CCDI (from the top combination)
best_rcdi <- top3[1, ]
best_rcdi <- best_rcdi[, c("Combination","Train_Cindex", "Test_Cindex", "Ind1_Cindex", "Ind2_Cindex")]
colnames(best_rcdi) <- c("Gene_list","Train", "Test", "GSE39582", "GSE161158")

# Combine the filtered pathways and the best CCDI into one data frame
plot_data <- rbind(filtered_data, best_rcdi)

# Step 3: Create a new data frame for the barplot
plot_data <- plot_data[, c("Gene_list", "Train", "Test", "GSE39582", "GSE161158")]
colnames(plot_data) <- c("Pathway", "Training", "Test", "GSE39582", "GSE161158")
plot_data$Pathway[5] <- "c-RCDI"

# Step 4: Melt the data for ggplot
library(reshape2)
melted_plot_data <- melt(plot_data, id.vars = "Pathway", variable.name = "Dataset", value.name = "C_Index")
melted_plot_data$Pathway <- factor(melted_plot_data$Pathway, 
                                   levels = c("c-RCDI", "Ferroptosis", "Autosis", "NETosis", "Pyroptosis"))

# Reorder the x-axis by descending mean C-index values
melted_plot_data$Pathway <- reorder(melted_plot_data$Pathway, -melted_plot_data$C_Index)

# Step 5: Create the barplot using ggplot2
barplot <- ggplot(melted_plot_data, aes(x = Pathway, y = C_Index, fill = Dataset)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black") +  # Use geom_col for barplot
  geom_errorbar(aes(ymin = C_Index - 0.02, ymax = C_Index + 0.02),  # Add error bars (example: +/- 0.02 as placeholder)
                width = 0.2, position = position_dodge(width = 0.8)) +
  geom_text(aes(label = round(C_Index, 3)), 
            position = position_dodge(width = 0.8), 
            vjust = -2, size = 3.5, family = "Helvetica") +
  stat_summary(aes(group = Pathway), fun = mean, geom = "crossbar",
               width = 0.5, color = "black", fatten = 2,
               position = position_nudge(x = 0),
               show.legend = FALSE) +
  stat_summary(
    aes(group = Pathway, label = after_stat(paste("Mean:", round(y, 3)))), 
    fun = mean, geom = "text", vjust = -2, size = 4, face = "bold",
    position = position_nudge(x = 0), family = "Helvetica"
  ) +
  scale_fill_manual(values = c("Training" = "#1b9e77", 
                               "Test" = "#d95f02", 
                               "GSE39582" = "#7570b3", 
                               "GSE161158" = "#e7298a")) +
  labs(
    title = NULL,
    x = "Signature",
    y = "C-Index",
    fill = "Dataset"
  ) +
  theme_bw(base_family = "Helvetica") +  # Use theme_bw for a clean look
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = "bold", family = "Helvetica", colour = "black"),
    axis.text.y = element_text(size = 10, hjust = 0.5, face = "bold", family = "Helvetica", color = "black"),
    axis.title = element_text(size = 12, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.5, 0.95))  # Set y-axis limits for better visualization

# Step 4: Print the barplot
print(barplot)

# Save the plot as a TIF file
ggsave("cindex_comparison_signatures.tiff", device = "tiff",barplot, units = "in", width = 8, height = 8, dpi = 600, bg = "white", compression = "lzw")


print("Barplot of C-indices for RCDI and selected pathways created successfully")

################## Calculate AUC and PR-AUC for time-dependent ROC at 12, 36, and 60 months ######################

print("---------------------------------------------------")
print("Calculating AUC and PR-AUC for time-dependent ROC at 12, 36 and 60 months")
print("---------------------------------------------------")

RCDI_results <- list(
  Train = read.csv("train_RCDI_combined_risk_scores.csv", check.names = FALSE, row.names = 1),
  Test  = read.csv("test_RCDI_combined_risk_scores.csv", check.names = FALSE, row.names = 1),
  Ind1  = read.csv("ind1_RCDI_combined_risk_scores.csv", check.names = FALSE, row.names = 1),
  Ind2  = read.csv("ind2_RCDI_combined_risk_scores.csv", check.names = FALSE, row.names = 1)
)

# To verify the list
names(RCDI_results)

library(timeROC)
library(PRROC)
time_points <- c(12, 36, 60)
performance_results <- lapply(names(RCDI_results), function(ds) {
  data <- RCDI_results[[ds]]
  # Calculate time-dependent ROC
  roc <- timeROC(
    T = data$OS_MONTHS,
    delta = data$OS_STATUS,
    marker = data$RCDI,
    cause = 1,
    times = time_points,
    iid = TRUE
  )
  # Extract AUC values for specified time points
  auc_values <- roc$AUC
  # Calculate PR-AUC for specific time points
  pr_auc_values <- sapply(time_points, function(tp) {
    # Define binary outcome for PR-AUC (event within time point vs. not)
    binary_outcome <- ifelse(data$OS_MONTHS <= tp & data$OS_STATUS == 1, 1, 0)
    prroc_obj <- pr.curve(scores.class0 = data$RCDI, weights.class0 = binary_outcome, curve = FALSE)
    return(prroc_obj$auc.integral)
  })
  # Combine results into a data frame
  data.frame(
    Dataset = ds,
    Time_12_AUC = auc_values[1],
    Time_36_AUC = auc_values[2],
    Time_60_AUC = auc_values[3],
    Time_12_PR_AUC = pr_auc_values[1],
    Time_36_PR_AUC = pr_auc_values[2],
    Time_60_PR_AUC = pr_auc_values[3]
  )
})

# Combine AUC and PR-AUC results into a single data frame
performance_results_df <- do.call(rbind, performance_results)

# Save performance results to a CSV
fwrite(performance_results_df, "time_dependent_performance_RCDI_results.csv", sep = ",")

# Print AUC and PR-AUC results
print("AUC and PR-AUC Results for Time-Dependent ROC of RCDI:")
print(performance_results_df)

print("ROC-AUC calculated and saved successfully")


########################## Scientific AUC Plots for RCDI ##############################

print("---------------------------------------------------")
print("Creating ROC-AUC plots")
print("---------------------------------------------------")

if (!require("patchwork")) install.packages("patchwork")
library(timeROC)
library(ggplot2)
library(cowplot)
library(patchwork)

# Set time points and base color palette
time_points <- c(12, 36, 60)
base_colors <- c("1-year" = "#1b9e77", "3-years" = "#d95f02", "5-years" = "#7570b3")

# Container for plots
roc_plots <- list()

# Loop through each dataset in the RCDI results
for (ds in names(RCDI_results)) {
  data <- RCDI_results[[ds]]
  
  # Compute ROC
  roc_obj <- timeROC(
    T = data$OS_MONTHS,
    delta = data$OS_STATUS,
    marker = data$RCDI,
    cause = 1,
    times = time_points,
    iid = TRUE
  )
  
  # Extract ROC coordinates
  roc_data <- data.frame(
    FPR = c(roc_obj$FP[, "t=12"], roc_obj$FP[, "t=36"], roc_obj$FP[, "t=60"]),
    TPR = c(roc_obj$TP[, "t=12"], roc_obj$TP[, "t=36"], roc_obj$TP[, "t=60"]),
    Time = factor(rep(c("1-year", "3-years", "5-years"), each = nrow(roc_obj$FP)))
  )
  
  # Format AUC into legend
  auc_vals <- round(roc_obj$AUC, 2)
  legend_labels <- c(
    paste0("1-year (AUC = ", auc_vals[1], ")"),
    paste0("3-years (AUC = ", auc_vals[2], ")"),
    paste0("5-years (AUC = ", auc_vals[3], ")")
  )
  names(base_colors) <- legend_labels
  roc_data$Time <- factor(roc_data$Time,
                          levels = c("1-year", "3-years", "5-years"),
                          labels = legend_labels)
  
  # Plot
  p <- ggplot(roc_data, aes(x = FPR, y = TPR, color = Time)) +
    geom_line(size = 1.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray70") +
    scale_color_manual(values = base_colors) +
    labs(
      # title = paste0("Time-dependent ROC (", toupper(ds), ")"),
      title = NULL,
      x = "1 - Specificity",
      y = "Sensitivity",
      color = "Time points"
    ) +
    theme_classic(base_family = "Helvetica") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right",
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 14, color= "black"),
      axis.title = element_text(size = 16, face = "bold", color= "black"),
    )
  
  # Save individual plot
  # ggsave(paste0("TimeROC_RCDI_", ds, ".pdf"), p, width = 7, height = 6)
  ggsave(paste0("TimeROC_RCDI_", ds, ".tif"), p, dpi = 600, width = 7, height = 5, units = "in", compression = "lzw")
  
  roc_plots[[ds]] <- p
}

#  Combine all into a single panel #
panel_plot <- (roc_plots$train + roc_plots$test) / (roc_plots$ind1 + roc_plots$ind2) +
  plot_annotation(title = "Time-dependent ROC Curves (RCDI Model)",
                  theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)))

# Save full-panel figure
ggsave("RCDI_TimeROC_4Panel.pdf", panel_plot, width = 14, height = 12)
ggsave("RCDI_TimeROC_4Panel.png", panel_plot, dpi = 600, width = 14, height = 12)

# Show
print(panel_plot)

print("Successfully created ROC-AUC plots")
########################## END OF AUC PLOTS ##############################

############################################### END #############################################
