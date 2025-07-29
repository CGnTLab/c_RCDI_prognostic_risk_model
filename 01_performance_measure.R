#!/usr/bin/Rscript

setwd("/home2/New_Objective2_Biomarker_research/batch_corrected_rna_seq/")

#### ====================================================== ###
################# Identification of DEGs ######################
###=========================================================###

rm(list= ls())
#rm(list=ls())
set.seed(42)

library(dplyr)
library(data.table)
library(readxl)
# BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

# EMTAB dataset DEGs -------

### Reading the dataframes - counts values ###
library(data.table)
emtab <- data.frame(fread("/home2/New_Objective2_Biomarker_research/E-MTAB-12862/data/CRC.SW.mRNA.symbol.count.txt", sep='\t'), row.names=1, check.names=F)
dim(emtab)
emtab_tumor <- emtab[, grepl("\\.T.*$", colnames(emtab))]
dim(emtab_tumor)
emtab_normal <- emtab[, grepl("\\.N.*$", colnames(emtab))]
dim(emtab_normal)
emtab <- cbind(emtab_tumor, emtab_normal)


###### DEG analysis using limma ## --------
#BiocManager::install("DESeq2")
# library(DESeq2)

# Load required libraries
library(edgeR)
library(limma)
library(data.table)

# Define group labels for each sample (e.g., 1063 tumor and 120 normal samples)
group <- c(rep("Tumor", 1063), rep("Normal", 120))
if(length(group) != ncol(emtab)) stop("The number of samples in 'data' does not match the length of 'group'.")

# Create a DGEList object
dge <- DGEList(counts = emtab)
# Filter out lowly expressed genes (optional; adjust cutoff as needed)
keep <- filterByExpr(dge, group = group)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Normalize library sizes
dge <- calcNormFactors(dge)

# Create design matrix
# First, set the group factor with "Normal" as the reference level
group <- factor(group, levels = c("Normal", "Tumor"))
design <- model.matrix(~ group)

# Apply voom transformation to account for mean-variance relationship
v <- voom(dge, design, plot = TRUE)

# Fit the linear model using limma
fit <- lmFit(v, design)

# Compute moderated t-statistics, p-values, and log-odds of differential expression
fit <- eBayes(fit)

# Extract differential expression results for the tumor effect
results_emtab <- topTable(fit, coef = "groupTumor", number = Inf, adjust.method = "BH", sort.by = "P")
print(head(results_emtab))

# Save the results to a CSV file
write.csv(results_emtab, file = "DEGs_Tumor_vs_Normal_EMTAB.csv", row.names = TRUE)
results_emtab <- read.table("DEGs_Tumor_vs_Normal_EMTAB.csv", sep=',', header=T, check.names=F, row.names=1)
degs_emtab <- results_emtab[(abs(results_emtab$logFC) > 0.5) & (results_emtab$adj.P.Val) < 0.05,] ## 4005 6
write.table(degs_emtab, "DEGs_emtab_0.5.txt", sep="\t", row.names = T, col.names = NA)

top10_genes <- rownames(degs_emtab[order(degs_emtab$adj.P.Val), ][1:10, ])

# Create Volcano Plot using EnhancedVolcano
library(EnhancedVolcano)
tiff("volcano_plot_emtab.tif", width = 8, height = 6, units = "in", bg = "white", res=600, compression = "lzw")
EnhancedVolcano(results_emtab,
                lab = rownames(results_emtab),
                x = 'logFC',
                y = 'adj.P.Val',
                title = NULL,
                subtitle = NULL,
                xlab = bquote(~Log[2]~ 'Fold Change'),
                ylab = bquote(~-Log[10]~ 'Adjusted P-value'),
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 1,
                labSize = 4.0,
                colAlpha = 0.4,
                legendLabels = c('ns', 'log2 FC', 'p.adj', 'significant'),
                # legendLabels = c('NS', 'Significant'),
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 3.0,
                drawConnectors = F,
                selectLab = top10_genes,  # Highlight top 10 genes
                axisLabSize = 16,
                caption = NULL
                )
dev.off()

# TCGA dataset DEGs -------
cbio_tumor <- data.frame(fread("/home2/New_Objective2_Biomarker_research/UCSC_Xena_Browser/TCGA_CRC_cancer_count_batch.csv", sep=',', header=T), check.names=F, row.names = 1)
cbio_normal <- data.frame(fread("/home2/New_Objective2_Biomarker_research/UCSC_Xena_Browser/TCGA_CRC_normal_count_batch.csv", sep=',', header=T), check.names=F, row.names = 1)
cbio <- merge(cbio_tumor, cbio_normal, by=0)
row.names(cbio) <- cbio[, 1]
cbio <- cbio[,-1]
# dim(cbio)

group <- c(rep("Tumor", 618), rep("Normal", 51))
if(length(group) != ncol(cbio)) stop("The number of samples in 'data' does not match the length of 'group'.")

# Create a DGEList object
dge <- DGEList(counts = cbio)
# Filter out lowly expressed genes (optional; adjust cutoff as needed)
keep <- filterByExpr(dge, group = group)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Normalize library sizes
dge <- calcNormFactors(dge)

# Create design matrix
# First, set the group factor with "Normal" as the reference level
group <- factor(group, levels = c("Normal", "Tumor"))
design <- model.matrix(~ group)
# The coefficient 'groupTumor' will reflect the log2 fold change of Tumor vs Normal

# Apply voom transformation to account for mean-variance relationship
v <- voom(dge, design, plot = TRUE)

# Fit the linear model using limma
fit <- lmFit(v, design)

# Compute moderated t-statistics, p-values, and log-odds of differential expression
fit <- eBayes(fit)

# Extract differential expression results for the tumor effect
results_tcga <- topTable(fit, coef = "groupTumor", number = Inf, adjust.method = "BH", sort.by = "P")
print(head(results_tcga))

# Save the results to a CSV file
write.csv(results_tcga, file = "DEGs_Tumor_vs_Normal_TCGA.csv", row.names = TRUE)
results_tcga <- read.table("DEGs_Tumor_vs_Normal_TCGA.csv", row.names=1, sep=',', check.names = F, header = T)
degs_cbio <- results_tcga[(abs(results_tcga$logFC) > 0.5) & (results_tcga$adj.P.Val) < 0.05,] # 7061  6
write.table(degs_cbio, "DEGs_TCGA_0.5.txt", sep="\t", row.names = T, col.names=NA)

# Create Volcano Plot using EnhancedVolcano
tiff("volcano_plot_tcga.tif", width = 8, height = 8, units = "in", bg = "white", res=600, compression = "lzw")
EnhancedVolcano(results_tcga,
                lab = rownames(results_tcga),
                x = 'logFC',
                y = 'adj.P.Val',
                title = 'Volcano Plot of Differential Expression : TCGA',
                subtitle = 'Cutoffs: |log2FC| > 0.5 and adj.P.Val < 0.05',
                xlab = bquote(~Log[2]~ 'Fold Change'),
                ylab = bquote(~-Log[10]~ 'Adjusted P-value'),
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 0.5,
                labSize = 6.0,
                colAlpha = 0.6,
                legendLabels = c('NS', 'Log2 FC', 'Adj. P-value', 'Adj. P-value & Log2 FC'),
                legendPosition = 'top',
                legendLabSize = 15,
                legendIconSize = 6.0,
                drawConnectors = F
                # widthConnectors = 0.5,
                # colConnectors = "grey50"
                )
dev.off()

#### ====================================================== ###
################# Dataset loading and preprocesses ######################
###=========================================================###

options(repos = c(CRAN = "https://cran.rstudio.com/"))
gc()
rm(list=ls())
seed = 42
set.seed(seed)


start_time <- Sys.time()

library(dplyr)
library(data.table)
# library(readxl)


library(BiocParallel)
library(parallel)
# num_detected_cores <- detectCores(logical = FALSE) - 2 # logical = FALSE gives physical cores
num_detected_cores = 7 
register(MulticoreParam(workers = num_detected_cores))

# merged<- data.frame(fread("combined_z_score_clinical.csv", sep=','), row.names=1, check.names=F)
merged <- data.frame(fread("combined_z_score_clinical_batch.csv", sep=','), row.names=1, check.names=F)
# merged_full <- data.frame(fread("combined_z_score_clinical_batch.csv", sep=','), row.names=1, check.names=F)
# print(paste0("No. of samples in original dataset (E-MTAB-12862 + TCGA) :", nrow(merged_full)))
# merged <- merged_full[ (merged_full$OS_MONTHS <= 140.38356164) ,]
# print(paste0("No. of samples in dataset (E-MTAB-12862 + TCGA) filtered by OS_MONTHS :", nrow(merged)))

# merged<- data.frame(fread("emtab_z_score_clinical.csv", sep=','), row.names=1, check.names=F)
# train   <- merged
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
test4 <- data.frame(fread("independent_ds_gse14333.csv", sep=",", header=T), row.names = 1, check.names = F)
test4 <- test4 %>%
  dplyr::rename(
    OS_MONTHS = DFS_Time,
    OS_STATUS = DFS_Cens
  )

data_list <- list(train, test1, test2, test3, test4)

# Find common column names across all data frames
common_cols <- Reduce(intersect, lapply(data_list, colnames))

train <- train[, common_cols, drop = FALSE]
test1 <- test1[, common_cols, drop = FALSE]
test2 <- test2[, common_cols, drop = FALSE]
test3 <- test3[, common_cols, drop = FALSE]
test4 <- test4[, common_cols, drop = FALSE]


# subset for only DEGs for all datasets
train <- train[, names(train) %in% c(genes_to_select, "OS_MONTHS", "OS_STATUS")]
test1 <- test1[, names(test1) %in% c(genes_to_select, "OS_MONTHS", "OS_STATUS")]
test2 <- test2[, names(test2) %in% c(genes_to_select, "OS_MONTHS", "OS_STATUS")]
test3 <- test3[, names(test3) %in% c(genes_to_select, "OS_MONTHS", "OS_STATUS")]
test4 <- test4[, names(test4) %in% c(genes_to_select, "OS_MONTHS", "OS_STATUS")]


gene_set <- read.table("../gene_lists/parthanatos_related_genes.txt", sep = ",", header = FALSE, check.names = FALSE)$V1
gene_set <- unique(gene_set)


train <- train[, names(train) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test1 <- test1[, names(test1) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test2 <- test2[, names(test2) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test3 <- test3[, names(test3) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]
test4 <- test4[, names(test4) %in% c(gene_set, "OS_MONTHS", "OS_STATUS")]


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
  Ind1 = test2,
  Ind2 = test3
)

df_list <- bplapply(df_list, function(x) {
  class(x$OS_STATUS) <- "numeric"
  class(x$OS_MONTHS) <- "numeric"
  x <- na.omit(x)
  x <- x[!(x$OS_MONTHS < 1),]
  # x <- x[ (x$OS_MONTHS <= 140.38356164) ,]
}, BPPARAM = MulticoreParam(num_detected_cores))

# df_list <- lapply(df_list, function(x) {
#   # Identify feature columns (excluding OS_MONTHS and OS_STATUS)
#   feature_cols <- setdiff(colnames(x), c("OS_MONTHS", "OS_STATUS"))
#   
#   # Apply min-max scaling to each feature column
#   for (col in feature_cols) {
#     min_val <- min(x[[col]], na.rm = TRUE)
#     max_val <- max(x[[col]], na.rm = TRUE)
#     
#     # Check if min and max are the same (to avoid division by zero)
#     if (max_val > min_val) {
#       x[[col]] <- (x[[col]] - min_val) / (max_val - min_val)
#     } else {
#       # If min and max are the same (e.g., all values are the same), set all scaled values to 0.5 (or any constant)
#       x[[col]] <- 0.5 # or you could set to 0, or keep original values, depending on requirement
#       warning(paste("Column", col, "has constant values. Min-max scaling resulted in constant values set to 0.5."))
#     }
#   }
#   return(x)
# })

#scaling of the datasets again
# df_list <- bplapply(df_list, function(x) {
#   survival_data <- x[, c("OS_MONTHS", "OS_STATUS"), drop=FALSE]
#   scaled_data <- scale(x[, !(names(x) %in% c("OS_MONTHS", "OS_STATUS"))], center = TRUE, scale = TRUE)
#   x <- data.frame(scaled_data, survival_data)
#   return(x)
# }, BPPARAM = MulticoreParam(num_detected_cores))

# select the train dataset
train <- df_list$Train ## explicitly select the train dataset

# mention time points for survivalROC and PRROC
time_points <- c(12, 36, 60)

nodesize = 5
nodedepth = 5
ntree = 1000

###### Univariate Analysis  --
univ_formulas <- sapply(covariates,
                        function(x) {
                          formula_string <- paste0("Surv(OS_MONTHS, OS_STATUS) ~ ", x)
                          as.formula(formula_string)
                        })
univ_models <- bplapply( univ_formulas, function(x){coxph(x, data = train)}, BPPARAM = MulticoreParam(num_detected_cores))
univ_results <- bplapply(univ_models,
                       function(x){ 
                         s <- summary(x)
                         p.value<-signif(s$wald["pvalue"], digits=2)
                         wald.test<-signif(s$wald["test"], digits=2)
                         beta<-signif(s$coef[1], digits=2);#coeficient beta
                         HR <-signif(s$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(s$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(s$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       }, BPPARAM = MulticoreParam(num_detected_cores))
res <- t(as.data.frame(univ_results, check.names = FALSE))
res<- as.data.frame(res)
res$p.adjust <- p.adjust(res$p.value, method = "BH")
class(res$p.value) <- "numeric"
res_mod<- res[res$p.value < 0.05,] #pick univariate results with significant p values
res_mod <- res_mod[order(res_mod$p.adjust),]
dim(res_mod) #40  4
row.names(res_mod)

results <- data.frame() # Initialize an empty results dataframe

######################### 1.1 Univariate + Multivariate  ################################
# gc()
# genes <- row.names(res_mod)
# train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]
# formula_string <- paste("Surv(OS_MONTHS, OS_STATUS) ~ ", paste(genes, collapse = " + "))
# formula <- as.formula(formula_string)
# res.cox <- coxph(formula, data = train_data)
# # step.model <- stepAIC(res.cox, direction = "both", trace = FALSE)
# summary(res.cox)
# res.cox.coeff <- summary(res.cox)$coefficients
# res.cox.coeff.sig <- res.cox.coeff[res.cox.coeff[,ncol(res.cox.coeff)] < 0.05,]
# dim(res.cox.coeff.sig)
# res.cox.coeff.sig
# 
# genes <- row.names(res.cox.coeff.sig)
# train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]
# formula_string <- paste("Surv(OS_MONTHS, OS_STATUS) ~ ", paste(genes, collapse = " + "))
# formula <- as.formula(formula_string)
# res.cox <- coxph(formula, data = train_data)
# # step.model <- stepAIC(res.cox, direction = "both", trace = FALSE)
# summary(res.cox)
# res.cox.coeff <- summary(res.cox)$coefficients
# res.cox.coeff.sig <- res.cox.coeff[res.cox.coeff[,ncol(res.cox.coeff)] < 0.05,]
# dim(res.cox.coeff.sig)
# res.cox.coeff.sig
# 
# genes <- row.names(res.cox.coeff)
# RS_list <- bplapply(df_list, function(x) {
#   df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
#   predict(res.cox, type = "risk", newdata = df_data)
# }, BPPARAM = MulticoreParam(num_detected_cores))
# 
# 
# cc_list <- bplapply(names(df_list), function(x) {
#   RS <- RS_list[[x]]
#   df <- df_list[[x]]
#   cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
#   # Calculate Time-dependent ROC-AUC at specified time points
#   auc_values <- sapply(time_points, function(tp) {
#     roc_obj <- survivalROC(Stime = df$OS_MONTHS,
#                            status = df$OS_STATUS,
#                            marker = RS,
#                            predict.time = tp,
#                            method = "KM") # or "NNE", method for handling censoring
#     return(roc_obj$AUC)
#   })
#   names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
#   # Calculate Time-dependent PR-AUC at specified time points
#   pr_auc_values <- sapply(time_points, function(tp) {
#     binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
#     pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
#     return(pr_obj$auc.integral)
#   })
#   names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
#   
#   data.frame(
#     Model_combination = "Univariate + Multivariate",
#     Dataset_Name = x,
#     C_index = as.numeric(cc),
#     AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
#     AUC_36months = as.numeric(auc_values["AUC_36months"]),
#     AUC_60months = as.numeric(auc_values["AUC_60months"]),
#     PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
#     PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
#     PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
#     gene_set = paste(genes, collapse = ","),
#     No_of_genes_selected = length(genes)
#   )
# }, BPPARAM = MulticoreParam(num_detected_cores))
# 
# results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe
# 
# print(results[, 1:3])
#Concordance Index (C-index) for RFC model on training data:  0.5922707 


######################### 1.2 Univariate + StepCox  ################################
genes <- row.names(res_mod)
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]
# formula_string <- paste("Surv(OS_MONTHS, OS_STATUS) ~ ", paste(genes, collapse = " + "))
# formula <- as.formula(formula_string)
# res.cox <- coxph(formula, data = train_data)
step.model <- step(coxph(Surv(OS_MONTHS,OS_STATUS)~.,train_data), direction = "both", trace = FALSE)
summary(step.model)
res.cox.coeff <- summary(step.model)$coefficients

genes <- row.names(res.cox.coeff)

RS_list <- bplapply(df_list, function(x) {
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
  predict(step.model, type = "risk", newdata = df_data)
}, BPPARAM = MulticoreParam(num_detected_cores))


cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]                                                                                
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "Univariate + StepCox",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(genes, collapse = ","),
    No_of_genes_selected = length(genes)
  )
}, BPPARAM = MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])

######################### 1.3 Univariate + LASSO  ################################
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
# install.packages("devtools")
# devtools::install_github("zabore/condsurv")
library(condsurv)
library("glmnet")
genes <- row.names(res_mod)
data <-train[, colnames(train) %in% c(genes, "OS_STATUS", "OS_MONTHS")]

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

selected_genes <- row.names(cf_lasso_df_non_zero)
# rm(RS_list)
RS_list <- bplapply(df_list, function(x) {
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
  # newx <- as.matrix(df_data[, -c((ncol(df_data) - 1), ncol(df_data))])
  predict(cv.lasso, type = "link", newx = as.matrix(df_data[, -c((ncol(df_data)-1):ncol(df_data))]), s = cv.lasso$lambda.min)
}, BPPARAM = MulticoreParam(num_detected_cores))


# rm(cc_list)
cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "Univariate + LASSO",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(selected_genes, collapse = ","),
    No_of_genes_selected = length(selected_genes)
  )
}, BPPARAM = MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])

######################### 1.4 Univariate + RSF  ################################
gc()
library(randomForestSRC)
genes <- row.names(res_mod)
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

rf_model_final <- rfsrc(Surv(OS_MONTHS, OS_STATUS) ~ ., data = train_data,
                        importance = T, ntree = ntree, # Use full ntree for final model
                        proximity = T, forest = T,
                        seed = seed, nodesize = nodesize,
                        splitrule = "logrank", ncores= num_detected_cores,
                        nodedepth=nodedepth) 
# Use topvars for final model
vs.rf_model <- var.select(object = rf_model_final, conservative = "high", refit=T)
topvars <- vs.rf_model$topvars
length(topvars)

# --- Evaluate Final Model on Test Datasets (test1, test2) ---
RS_list_final <- lapply(df_list, function(x) { # Use df_list to evaluate on train, test1, test2
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")] # Use 'genes' or 'topvars' consistently
  predict(rf_model_final, newdata = df_data)$predicted
})

cc_list_final <- lapply(names(df_list), function(x) {
  RS <- RS_list_final[[x]]
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "Univariate + RSF",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(topvars	, collapse = ","),
    No_of_genes_selected = length(topvars)
  )
})

results <- bind_rows(results, bind_rows(cc_list_final)) # Append to results dataframe

print(results[, 1:3])


########################## 1.5 Univariate + SVM  ################################
#gc()

#library(survival)         # Survival analysis
#library(survivalsvm)      # Survival SVM
#library(caret)            # Feature selection
#library(e1071)            # SVM model
#library(dplyr)            # Data manipulation

#genes <- row.names(res_mod)
#train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

#gc()
#set.seed(seed)
#fit = survivalsvm(Surv(OS_MONTHS, OS_STATUS) ~., data= train_data, gamma.mu = 1)

#RS_list <- bplapply(df_list, function(x) {
#  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
#  as.numeric(predict(fit, df_data)$predicted)
#}, BPPARAM = MulticoreParam(num_detected_cores))

#cc_list <- bplapply(names(df_list), function(x) {
#  RS <- RS_list[[x]]                                                                                
#  df <- df_list[[x]]
#  cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
#  auc_values <- sapply(time_points, function(tp) {
#    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
#                           status = df$OS_STATUS,
#                           marker = RS,
#                           predict.time = tp,
#                           method = "KM") # or "NNE", method for handling censoring
#    return(roc_obj$AUC)
#  })
#  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
#  # Calculate Time-dependent PR-AUC at specified time points
#  pr_auc_values <- sapply(time_points, function(tp) {
#    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
#    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
#    return(pr_obj$auc.integral)
#  })
#  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
#  
#  data.frame(
#    Model_combination = "Univariate + SVM",
#    Dataset_Name = x,
#    C_index = as.numeric(cc),
#    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
#    AUC_36months = as.numeric(auc_values["AUC_36months"]),
#    AUC_60months = as.numeric(auc_values["AUC_60months"]),
#    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
#    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
#    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
#    gene_set = paste(genes, collapse = ","),
#    No_of_genes_selected = length(genes)
#  )
#}, BPPARAM = MulticoreParam(num_detected_cores))

#results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

#print(results[, 1:3])

######################### 1.6 Univariate + GBM  ################################
gc()
library(gbm)        # For Gradient Boosting Machine
genes <- row.names(res_mod)
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

gbm_model <- gbm(Surv(OS_MONTHS, OS_STATUS) ~ ., data = train_data, distribution = "coxph",
                 n.trees = 10000, interaction.depth = 5, shrinkage = 0.001,
                 cv.folds = 10, n.cores = num_detected_cores, verbose = FALSE, n.minobsinnode = 10)
best <- which.min(gbm_model$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS_MONTHS,OS_STATUS)~.,data = train_data,distribution = 'coxph',
           n.trees = best,
           interaction.depth = 5,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = num_detected_cores)

RS_list <- bplapply(df_list, function(x) {
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
  predict(fit, df_data, n.trees=best, type="link")
}, BPPARAM= MulticoreParam(num_detected_cores))

cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]                                                                                
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS,OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "Univariate + GBM",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(genes, collapse = ","),
    No_of_genes_selected = length(genes)
  )
}, BPPARAM= MulticoreParam(num_detected_cores))
cc_list

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])

######################### 1.7 Uni + XGB ##################################
#gc()
## install.packages("mlr3")
#library(mlr3)        # For machine learning tasks
## remotes::install_github("mlr-org/mlr3extralearners")
#library(mlr3extralearners)  # For survival learners, such as xgboost Cox
#library(xgboost)
## remotes::install_github("mlr-org/mlr3proba")
#library(mlr3proba)

#genes <- row.names(res_mod)
#train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

## Define the learner
#learner = mlr3::lrn("surv.xgboost.cox")
#print(learner)

#train_task <- as_task_surv(train_data, time="OS_MONTHS", event="OS_STATUS")

## install.packages("ml3tuning")
#library(mlr3tuning)
#library(paradox)

## Define a search space for hyperparameters
#search_space <- ps(
#  nrounds = p_int(50, 500),
#  eta = p_dbl(0.01, 0.3),
#  max_depth = p_int(2, 10),
#  min_child_weight = p_dbl(1, 10),
#  subsample = p_dbl(0.5, 1),
#  colsample_bytree = p_dbl(0.5, 1),
#  lambda = p_dbl(0, 1),
#  alpha = p_dbl(0, 1)
#)

## Define the tuning instance
#instance <- mlr3tuning::TuningInstanceBatchSingleCrit$new(
#  task = train_task,
#  learner = learner,
#  resampling = rsmp("cv", folds = 5),  # Cross-validation
#  measure = msr("surv.cindex"),       # Optimization metric
#  search_space = search_space,
#  terminator = trm("evals", n_evals = 50)  # Stop after 100 evaluations
#)

## Perform random search or grid search
#set.seed(seed)  # Ensure reproducibility
#tuner <- tnr("random_search")  # Or "grid_search"
#tuner$optimize(instance)

## Use the best hyperparameters
#learner$param_set$values <- instance$result_learner_param_vals

## Create train and test indices
#set.seed(seed)  # Ensure reproducibility

## Train the learner on the training data
#learner$train(train_task, seq_len(nrow(train_data)))

## Print model details
#print(learner$model)

## Get feature importance
#print(learner$importance())

## Assuming `df_list` is a list of datasets you want to apply the model to:
#RS_list <- bplapply(df_list, function(x) {
#  # Prepare the data (ensure to include relevant columns)
#  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
#  predictions = learner$predict_newdata(df_data)
#  # Predict survival using the trained xgboost model
#  as.numeric(predictions$lp)
#}, BPPARAM = MulticoreParam(num_detected_cores))

## Calculate C-index, AUC, and PR-AUC for each dataset
#cc_list <- bplapply(names(df_list), function(x) {
#  RS <- RS_list[[x]]  # Get the predicted risk scores
#  df <- df_list[[x]]
#  # Calculate Concordance Index (C-index)
#  cc <- learner$predict_newdata(df)$score()
#  # Calculate AUC for specific time points
#  auc_values <- sapply(time_points, function(tp) {
#    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
#                           status = df$OS_STATUS,
#                           marker = RS,
#                           predict.time = tp,
#                           method = "KM") # or "NNE", method for handling censoring
#    return(roc_obj$AUC)
#  })
#  names(auc_values) <- paste0("AUC_", time_points, "months")  # Naming AUC values
#  
#  # Calculate Time-dependent PR-AUC at specified time points
#  pr_auc_values <- sapply(time_points, function(tp) {
#    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
#    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
#    return(pr_obj$auc.integral)
#  })
#  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months")  # Naming PR-AUC values
#  
#  # Create a data frame with results
#  data.frame(
#    Model_combination = "Univariate + XGBoost",
#    Dataset_Name = x,
#    C_index = as.numeric(cc),
#    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
#    AUC_36months = as.numeric(auc_values["AUC_36months"]),
#    AUC_60months = as.numeric(auc_values["AUC_60months"]),
#    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
#    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
#    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
#    gene_set = paste(genes, collapse = ","),
#    No_of_genes_selected = length(genes)
#  )
#}, BPPARAM = MulticoreParam(num_detected_cores))

## Combine the results into a single dataframe
#results <- bind_rows(results, bind_rows(cc_list))

## Print results
#print(results[, 1:3])


######################### 1.8 Uni + CoxBoost ################################
gc()
# library(devtools)
# devtools::install_github("binderh/CoxBoost")
# install.packages("snowfall")
library(snowfall)
library(CoxBoost)
set.seed(seed)
genes <- row.names(res_mod)
train_data <- train[, names(train) %in% c(genes, "OS_STATUS", "OS_MONTHS")]

# CoxBoost in snowfall computing
sfInit(parallel = TRUE, cpus = num_detected_cores)
sfLibrary(CoxBoost) # Make CoxBoost library available on workers
sfExport("train_data") # Export train_data to worker environments
sfExport("genes") # Export genes to worker environments
sfExport("seed")
pen <- optimCoxBoostPenalty(train_data[,'OS_MONTHS'],train_data[,'OS_STATUS'], as.matrix(train_data[, !(names(train_data) %in% c("OS_MONTHS","OS_STATUS"))]),
                            trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(train_data[,'OS_MONTHS'],train_data[,'OS_STATUS'], as.matrix(train_data[, !(names(train_data) %in% c("OS_MONTHS","OS_STATUS"))]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty, parallel = T)
fit <- CoxBoost(train_data[,'OS_MONTHS'],train_data[,'OS_STATUS'], as.matrix(train_data[, !(names(train_data) %in% c("OS_MONTHS","OS_STATUS"))]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)
sfStop() # stop snowfall processing

## Genes selected
nonzero_coef <- coef(fit, type = "nonzero")
genes_selected <- names(nonzero_coef[nonzero_coef != 0])


RS_list <- bplapply(df_list, function(x) {
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
  as.numeric(predict(fit, type = "lp", newdata = df_data[, !(names(df_data) %in% c("OS_MONTHS","OS_STATUS"))], newtime= df_data$OS_MONTHS, newstatus=df_data$OS_STATUS))
}, BPPARAM= MulticoreParam(num_detected_cores))

cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]                                                                                
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "Univariate + CoxBoost",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(genes_selected, collapse = ","),
    No_of_genes_selected = length(genes_selected)
  )
}, BPPARAM= MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])


######################### 1.9 Uni + SuperPC  ################################

# install.packages("superpc")
library(superpc)
gc()
genes <- row.names(res_mod)
train_data <-train[, colnames(train) %in% c(genes, "OS_STATUS", "OS_MONTHS")]
data <- list(x=t(train_data[,-c((ncol(train_data)-1):ncol(train_data))]),y=train_data$OS_MONTHS,censoring.status=train_data$OS_STATUS,featurenames=colnames(train_data)[-c((ncol(train_data)-1):ncol(train_data))])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                     n.fold = 10,
                     n.components=3,
                     min.features=2,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)

RS_list <- bplapply(df_list, function(w) {
  df_data <- w[, c(genes, "OS_MONTHS", "OS_STATUS")]
  df_data <- list(x=t(df_data[,-c((ncol(df_data)-1):ncol(df_data))]),y=df_data$OS_MONTHS,censoring.status=df_data$OS_STATUS,featurenames=colnames(df_data)[-c((ncol(df_data)-1):ncol(df_data))])
  ff <- superpc.predict(fit,data,df_data,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  as.numeric(ff$v.pred)
}, BPPARAM= MulticoreParam(num_detected_cores))

cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]                                                                                
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "Univariate + SuperPC",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(genes, collapse = ","),
    No_of_genes_selected = length(genes)
  )
}, BPPARAM= MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])


######################### 1.10 Uni + Ridge  ################################
gc()
genes <- row.names(res_mod)
data <-train[, colnames(train) %in% c(genes, "OS_STATUS", "OS_MONTHS")]

cv.lasso <- cv.glmnet(x = as.matrix(data[, -c((ncol(data) - 1):ncol(data))]),  # Exclude OS_Months and OS_Status
                      y = Surv(data$OS_MONTHS, data$OS_STATUS),
                      family = "cox",
                      nfolds = 10,
                      alpha=0)

# Extract optimal lambda
best_lambda_lasso <- cv.lasso$lambda.min  # Or lambda.1se for more conservative model
# best_lambda_lasso <- cv.lasso$lambda.1se

# Fit the Lasso model with the optimal lambda
lasso_fit <- glmnet(x = as.matrix(data[, -c((ncol(data) - 1):ncol(data))]),  # Exclude OS_Months and OS_Status
                    y = Surv(data$OS_MONTHS, data$OS_STATUS),
                    family = "cox",
                    lambda = best_lambda_lasso,
                    alpha=0)

cf_lasso<- coef(lasso_fit)
cf_lasso_df<- as.data.frame(as.matrix(cf_lasso))
cf_lasso_df_non_zero<- cf_lasso_df[!(cf_lasso_df$s0==0), , drop=F]
dim(cf_lasso_df_non_zero) #10 1

RS_list <- bplapply(df_list, function(x) {
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
  predict(cv.lasso, type = "link", newx = as.matrix(df_data[, -c((ncol(df_data)-1):ncol(df_data))]), s = cv.lasso$lambda.min)
}, BPPARAM= MulticoreParam(num_detected_cores))

cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "Univariate + Ridge",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(genes, collapse = ","),
    No_of_genes_selected = length(genes)
  )
}, BPPARAM= MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])


######################### 1.10 Uni + plsRcox  ################################
# install.packages("devtools")

# devtools::install_github("fbertran/plsRcox")
library(plsRcox)

genes <- row.names(res_mod)

train_data <- train[, names(train) %in% c(genes, "OS_STATUS", "OS_MONTHS")]

set.seed(seed)
cv.plsRcox.res=cv.plsRcox(list(x=train_data[,genes],time=train_data$OS_MONTHS,status=train_data$OS_STATUS),nt=5,verbose = FALSE)
fit <- plsRcox(train_data[,genes],time=train_data$OS_MONTHS,event=train_data$OS_STATUS,nt=as.numeric(cv.plsRcox.res[5]))

RS_list <- bplapply(df_list, function(x) {
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
  as.numeric(predict(fit, newdata= df_data[, !(names(df_data) %in% c("OS_MONTHS","OS_STATUS"))], type="lp"))
}, BPPARAM= MulticoreParam(num_detected_cores))

cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]                                                                                
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS,OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "Univariate + plsRcox",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(genes, collapse = ","),
    No_of_genes_selected = length(genes)
  )
}, BPPARAM= MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])



########################  2.1 Multivariate  ################################

# genes <- covariates
# train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]
# formula_string <- paste("Surv(OS_MONTHS, OS_STATUS) ~ ", paste(genes, collapse = " + "))
# formula <- as.formula(formula_string)
# res.cox <- coxph(formula, data = train_data)
# # step.model <- stepAIC(res.cox, direction = "both", trace = FALSE)
# summary(res.cox)
# res.cox.coeff <- summary(res.cox)$coefficients
# res.cox.coeff.sig <- res.cox.coeff[res.cox.coeff[,ncol(res.cox.coeff)] < 0.05,]
# dim(res.cox.coeff.sig)
# res.cox.coeff.sig
# 
# genes <- row.names(res.cox.coeff.sig)
# train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]
# formula_string <- paste("Surv(OS_MONTHS, OS_STATUS) ~ ", paste(genes, collapse = " + "))
# formula <- as.formula(formula_string)
# res.cox <- coxph(formula, data = train_data)
#  # step.model <- stepAIC(res.cox, direction = "both", trace = FALSE)
# summary(res.cox)
# res.cox.coeff <- summary(res.cox)$coefficients
# res.cox.coeff.sig <- res.cox.coeff[res.cox.coeff[,ncol(res.cox.coeff)] < 0.05,]
# dim(res.cox.coeff.sig)
# res.cox.coeff.sig
# 
# genes <- row.names(res.cox.coeff)
# RS_list <- bplapply(df_list, function(x) {
#   df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
#   predict(res.cox, type = "risk", newdata = df_data)
# }, BPPARAM= MulticoreParam(num_detected_cores))
# cc_list <- bplapply(names(df_list), function(x) {
#   RS <- RS_list[[x]]
#   df <- df_list[[x]]
#   cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
#   auc_values <- sapply(time_points, function(tp) {
#     roc_obj <- survivalROC(Stime = df$OS_MONTHS,
#                            status = df$OS_STATUS,
#                            marker = RS,
#                            predict.time = tp,
#                            method = "KM") # or "NNE", method for handling censoring
#     return(roc_obj$AUC)
#   })
#   names(auc_values) <- paste0("AUC_", time_points, "months")  #Naming AUC values
#    #Calculate Time-dependent PR-AUC at specified time points
#   pr_auc_values <- sapply(time_points, function(tp) {
#     binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0)  #Event within time point vs. not
#     pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
#     return(pr_obj$auc.integral)
#   })
#   names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months")  #Naming PR-AUC values
# 
#   data.frame(
#     Model_combination = "Multivariate",
#     Dataset_Name = x,
#     C_index = as.numeric(cc),
#     AUC_12months = as.numeric(auc_values["AUC_12months"]),  #Extract AUC values
#     AUC_36months = as.numeric(auc_values["AUC_36months"]),
#     AUC_60months = as.numeric(auc_values["AUC_60months"]),
#     PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]),  #Extract PR-AUC values
#     PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
#     PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
#     gene_set = paste(genes, collapse = ","),
#     No_of_genes_selected = length(genes)
#   )
# }, BPPARAM= MulticoreParam(num_detected_cores))
# results <- bind_rows(results, bind_rows(cc_list))  #Append to results dataframe
# print(results[, 1:3])
# Concordance Index (C-index) for RFC model on training data:  0.5922707 


######################## 2.2 StepCox  ################################
genes <- covariates
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]
formula_string <- paste("Surv(OS_MONTHS, OS_STATUS) ~ ", paste(genes, collapse = " + "))
formula <- as.formula(formula_string)
res.cox <- coxph(formula, data = train_data)
step.model <- step(coxph(Surv(OS_MONTHS,OS_STATUS)~.,train_data), direction = "both", trace = FALSE)
summary(step.model)
res.cox.coeff <- summary(step.model)$coefficients

genes <- row.names(res.cox.coeff)

RS_list <- bplapply(df_list, function(x) {
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
  predict(step.model, type = "risk", newdata = df_data)
}, BPPARAM= MulticoreParam(num_detected_cores))

cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM")  #or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months")  #Naming AUC values
   #Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0)  #Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months")  #Naming PR-AUC values

  data.frame(
    Model_combination = "StepCox",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]),  #Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]),  #Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(genes, collapse = ","),
    No_of_genes_selected = length(genes)
  )
}, BPPARAM= MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list))  #Append to results dataframe

print(results[, 1:3])

######################### 2.3 LASSO  ################################
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

selected_genes <- row.names(cf_lasso_df_non_zero)
# rm(RS_list)
RS_list <- bplapply(df_list, function(x) {
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
  # newx <- as.matrix(df_data[, -c((ncol(df_data) - 1), ncol(df_data))])
  predict(cv.lasso, type = "link", newx = as.matrix(df_data[, -c((ncol(df_data)-1):ncol(df_data))]), s = cv.lasso$lambda.min)
}, BPPARAM= MulticoreParam(num_detected_cores))
# rm(cc_list)
cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "LASSO",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(selected_genes, collapse = ","),
    No_of_genes_selected = length(selected_genes)
  )
}, BPPARAM= MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])

######################### 2.4 RSF  ################################
gc()
library(randomForestSRC)
genes <- covariates
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

rf_model_final <- rfsrc(Surv(OS_MONTHS, OS_STATUS) ~ ., data = train_data,
                        importance = T, ntree = ntree, # Use full ntree for final model
                        proximity = T, forest = T,
                        seed = seed, nodesize = nodesize,
                        splitrule = "logrank",
                        ncores= num_detected_cores, nodedepth = nodedepth) 
# Use topvars for final model
vs.rf_model <- var.select(object = rf_model_final, conservative = "high")
topvars <- vs.rf_model$topvars
length(topvars)

# --- Evaluate Final Model on Test Datasets (test1, test2) ---
RS_list_final <- lapply(df_list, function(x) { # Use df_list to evaluate on train, test1, test2
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")] # Use 'genes' or 'topvars' consistently
  predict(rf_model_final, newdata = df_data)$predicted
})

cc_list_final <- lapply(names(df_list), function(x) {
  RS <- RS_list_final[[x]]
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "RSF",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(topvars, collapse = ","),
    No_of_genes_selected = length(topvars)
  )
})

results <- bind_rows(results, bind_rows(cc_list_final)) # Append to results dataframe

print(results[, 1:3])


########################## 2.5 CoxBoost  ################################
gc()
#CoxBoost
library(CoxBoost)
library(snowfall)
set.seed(seed)

sfInit(parallel = TRUE, cpus = num_detected_cores)

genes <- row.names(res_mod)
train_data <- train[, names(train) %in% c(genes, "OS_STATUS", "OS_MONTHS")]
sfLibrary(CoxBoost) # Make CoxBoost library available on workers
sfExport("train_data") # Export train_data to worker environments
sfExport("genes") # Export genes to worker environments
pen <- optimCoxBoostPenalty(train_data[,'OS_MONTHS'],train_data[,'OS_STATUS'], as.matrix(train_data[, !(names(train_data) %in% c("OS_MONTHS","OS_STATUS"))]),
                            trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(train_data[,'OS_MONTHS'],train_data[,'OS_STATUS'], as.matrix(train_data[, !(names(train_data) %in% c("OS_MONTHS","OS_STATUS"))]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty, parallel = T)
fit <- CoxBoost(train_data[,'OS_MONTHS'],train_data[,'OS_STATUS'], as.matrix(train_data[, !(names(train_data) %in% c("OS_MONTHS","OS_STATUS"))]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)
sfStop()

# gene selection
nonzero_coef <- coef(fit, type = "nonzero")
genes_selected <- names(nonzero_coef[nonzero_coef != 0])

RS_list <- bplapply(df_list, function(x) {
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
  as.numeric(predict(fit, type = "lp", newdata = df_data[, !(names(df_data) %in% c("OS_MONTHS","OS_STATUS"))], newtime= df_data$OS_MONTHS, newstatus=df_data$OS_STATUS))
}, BPPARAM= MulticoreParam(num_detected_cores))

cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]                                                                                
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "CoxBoost",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(genes_selected, collapse = ","),
    No_of_genes_selected = length(genes_selected)
  )
}, BPPARAM= MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])

########################## 2.6 SVM  ################################
#gc()
#library(survival)         # Survival analysis
#library(survivalsvm)      # Survival SVM
#library(caret)            # Feature selection
#library(e1071)            # SVM model
#library(dplyr)            # Data manipulation

#genes <- covariates
#train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

#set.seed(seed)
#fit = survivalsvm(Surv(OS_MONTHS, OS_STATUS) ~., data= train_data, gamma.mu = 1)

#RS_list <- bplapply(df_list, function(x) {
#  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
#  as.numeric(predict(fit, df_data)$predicted)
#}, BPPARAM= MulticoreParam(num_detected_cores))

#cc_list <- bplapply(names(df_list), function(x) {
#  RS <- RS_list[[x]]                                                                                
#  df <- df_list[[x]]
#  cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
#  auc_values <- sapply(time_points, function(tp) {
#    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
#                           status = df$OS_STATUS,
#                           marker = RS,
#                           predict.time = tp,
#                           method = "KM") # or "NNE", method for handling censoring
#    return(roc_obj$AUC)
#  })
#  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
#  # Calculate Time-dependent PR-AUC at specified time points
#  pr_auc_values <- sapply(time_points, function(tp) {
#    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
#    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
#    return(pr_obj$auc.integral)
#  })
#  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
#  
#  data.frame(
#    Model_combination = "SVM",
#    Dataset_Name = x,
#    C_index = as.numeric(cc),
#    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
#    AUC_36months = as.numeric(auc_values["AUC_36months"]),
#    AUC_60months = as.numeric(auc_values["AUC_60months"]),
#    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
#    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
#    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
#    gene_set = paste(genes, collapse = ","),
#    No_of_genes_selected = length(genes)
#  )
#}, BPPARAM= MulticoreParam(num_detected_cores))
#           shrinkage = 0.0001,
#           cv.folds = 10,n.cores = 8)

#RS_list <- bplapply(df_list, function(x) {
#  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
#  predict(fit, df_data, n.trees=best, type="link")
#}, BPPARAM= MulticoreParam(num_detected_cores))

#cc_list <- bplapply(names(df_list), function(x) {
#  RS <- RS_list[[x]]                                                                                
#  df <- df_list[[x]]
#  cc <- summary(coxph(Surv(OS_MONTHS,OS_STATUS) ~ RS, data = df))$concordance[1]
#  auc_values <- sapply(time_points, function(tp) {
#    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
#                           status = df$OS_STATUS,
#                           marker = RS,
#                           predict.time = tp,
#                           method = "KM") # or "NNE", method for handling censoring
#    return(roc_obj$AUC)
#  })
#  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
#  # Calculate Time-dependent PR-AUC at specified time points
#  pr_auc_values <- sapply(time_points, function(tp) {
#    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
#    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
#    return(pr_obj$auc.integral)
#  })
#  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
#  
#  data.frame(
#    Model_combination = "GBM",
#    Dataset_Name = x,
#    C_index = as.numeric(cc),
#    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
#    AUC_36months = as.numeric(auc_values["AUC_36months"]),
#    AUC_60months = as.numeric(auc_values["AUC_60months"]),
#    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
#    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
#    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
#    gene_set = paste(genes, collapse = ","),
#    No_of_genes_selected = length(genes)
#  )
#}, BPPARAM= MulticoreParam(num_detected_cores))

#results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

#print(results[, 1:3])


######################### 3.1 Uni + LASSO + StepCox  ################################
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
# install.packages("devtools")
# devtools::install_github("zabore/condsurv")
library(condsurv)
library("glmnet")
genes <- row.names(res_mod)
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

# genes <- row.names(cf_lasso_df_non_zero)
# rm(RS_list)
genes <- row.names(cf_lasso_df_non_zero)
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]
# formula_string <- paste("Surv(OS_MONTHS, OS_STATUS) ~ ", paste(genes, collapse = " + "))
# formula <- as.formula(formula_string)
# res.cox <- coxph(formula, data = train_data)
step.model <- step(coxph(Surv(OS_MONTHS,OS_STATUS)~.,train_data), direction = "both", trace = FALSE)
summary(step.model)
res.cox.coeff <- summary(step.model)$coefficients

genes <- row.names(res.cox.coeff)

RS_list <- bplapply(df_list, function(x) {
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
  predict(step.model, type = "risk", newdata = df_data)
}, BPPARAM= MulticoreParam(num_detected_cores))

cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]                                                                                
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "Univariate + LASSO + StepCox",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(genes, collapse = ","),
    No_of_genes_selected = length(genes)
  )
}, BPPARAM= MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])



########################## 3.2 Uni + LASSO + SVM  ################################
# 
# 
# gc()
# genes <- row.names(cf_lasso_df_non_zero)
# train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]
# 
# set.seed(seed)
# fit = survivalsvm(Surv(OS_MONTHS, OS_STATUS) ~., data= train_data, gamma.mu = 1)
# 
# RS_list <- bplapply(df_list, function(x) {
#   df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
#   as.numeric(predict(fit, df_data)$predicted)
# }, BPPARAM= MulticoreParam(num_detected_cores))
# 
# cc_list <- bplapply(names(df_list), function(x) {
#   RS <- RS_list[[x]]                                                                                
#   df <- df_list[[x]]
#   cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
#   auc_values <- sapply(time_points, function(tp) {
#     roc_obj <- survivalROC(Stime = df$OS_MONTHS,
#                            status = df$OS_STATUS,
#                            marker = RS,
#                            predict.time = tp,
#                            method = "KM") # or "NNE", method for handling censoring
#     return(roc_obj$AUC)
#   })
#   names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
#   # Calculate Time-dependent PR-AUC at specified time points
#   pr_auc_values <- sapply(time_points, function(tp) {
#     binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
#     pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
#     return(pr_obj$auc.integral)
#   })
#   names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
#   
#   data.frame(
#     Model_combination = "Univariate + LASSO + SVM",
#     Dataset_Name = x,
#     C_index = as.numeric(cc),
#     AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
#     AUC_36months = as.numeric(auc_values["AUC_36months"]),
#     AUC_60months = as.numeric(auc_values["AUC_60months"]),
#     PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
#     PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
#     PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
#     gene_set = paste(genes, collapse = ","),
#     No_of_genes_selected = length(genes)
#   )
# }, BPPARAM= MulticoreParam(num_detected_cores))
# 
# results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe
# 
# print(results[, 1:3])



######################### 3.3 Uni + LASSO + GBM  ################################

genes <- row.names(cf_lasso_df_non_zero)
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]
gbm_model <- gbm(Surv(OS_MONTHS, OS_STATUS) ~ ., data = train_data, distribution = "coxph",
                 n.trees = 10000, interaction.depth = 5, shrinkage = 0.001,
                 cv.folds = 10, n.cores = num_detected_cores, verbose = FALSE)
best <- which.min(gbm_model$cv.error)
fit <- gbm(formula = Surv(OS_MONTHS,OS_STATUS)~.,data = train_data,distribution = 'coxph',
           n.trees = best,
           interaction.depth = 5,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = num_detected_cores)

RS_list <- bplapply(df_list, function(x) {
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
  predict(fit, df_data, n.trees=best, type="link")
}, BPPARAM= MulticoreParam(num_detected_cores))

cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]                                                                                
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS,OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "Univariate + LASSO + GBM",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(genes, collapse = ","),
    No_of_genes_selected = length(genes)
  )
}, BPPARAM= MulticoreParam(num_detected_cores))
cc_list

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])


######################### 3.4 Uni + LASSO + CoxBoost  ################################

gc()
# CoxBoost
set.seed(seed)
library(CoxBoost)
genes <- row.names(cf_lasso_df_non_zero)
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

sfInit(parallel = TRUE, cpus = num_detected_cores)

train_data <- train[, names(train) %in% c(genes, "OS_STATUS", "OS_MONTHS")]
sfLibrary(CoxBoost) # Make CoxBoost library available on workers
sfExport("train_data") # Export train_data to worker environments
sfExport("genes") # Export genes to worker environments

pen <- optimCoxBoostPenalty(train_data[,'OS_MONTHS'],train_data[,'OS_STATUS'], as.matrix(train_data[, !(names(train_data) %in% c("OS_MONTHS","OS_STATUS"))]),
                            trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(train_data[,'OS_MONTHS'],train_data[,'OS_STATUS'], as.matrix(train_data[, !(names(train_data) %in% c("OS_MONTHS","OS_STATUS"))]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
fit <- CoxBoost(train_data[,'OS_MONTHS'],train_data[,'OS_STATUS'], as.matrix(train_data[, !(names(train_data) %in% c("OS_MONTHS","OS_STATUS"))]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)
sfStop()

nonzero_coef <- coef(fit, type = "nonzero")
genes_selected <- names(nonzero_coef[nonzero_coef != 0])

RS_list <- bplapply(df_list, function(x) {
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
  as.numeric(predict(fit, type = "lp", newdata = df_data[, !(names(df_data) %in% c("OS_MONTHS","OS_STATUS"))], newtime= df_data$OS_MONTHS, newstatus=df_data$OS_STATUS))
}, BPPARAM= MulticoreParam(num_detected_cores))

cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]                                                                                
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "Univariate + LASSO + CoxBoost",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(genes_selected, collapse = ","),
    No_of_genes_selected = length(genes_selected)
  )
}, BPPARAM= MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])



######################### 3.6 Uni + LASSO + RSF  ################################

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

# --- Evaluate Final Model on Test Datasets (test1, test2) ---
RS_list_final <- lapply(df_list, function(x) { # Use df_list to evaluate on train, test1, test2
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")] # Use 'genes' or 'topvars' consistently
  predict(rf_model_final, newdata = df_data)$predicted
})

cc_list_final <- lapply(names(df_list), function(x) {
  RS <- RS_list_final[[x]]
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "Uni + LASSO + RSF",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(topvars, collapse = ","),
    No_of_genes_selected = length(topvars)
  )
})

results <- bind_rows(results, bind_rows(cc_list_final)) # Append to results dataframe

print(results[, 1:3])

######################### 3.7 Uni + LASSO + plsRcox  ################################

genes <- row.names(cf_lasso_df_non_zero)
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

set.seed(seed)
cv.plsRcox.res=cv.plsRcox(list(x=train_data[,genes],time=train_data$OS_MONTHS,status=train_data$OS_STATUS),nt=5,verbose = FALSE)
fit <- plsRcox(train_data[,genes],time=train_data$OS_MONTHS,event=train_data$OS_STATUS,nt=as.numeric(cv.plsRcox.res[5]))

RS_list <- bplapply(df_list, function(x) {
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
  as.numeric(predict(fit, newdata= df_data[, !(names(df_data) %in% c("OS_MONTHS","OS_STATUS"))], type="lp"))
}, BPPARAM= MulticoreParam(num_detected_cores))

cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]                                                                                
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS,OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "Univariate + LASSO + plsRcox",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(genes, collapse = ","),
    No_of_genes_selected = length(genes)
  )
}, BPPARAM= MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])

######################### 3.8 Uni + LASSO + Ridge  ################################

genes <- row.names(cf_lasso_df_non_zero)
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

cv.lasso <- cv.glmnet(x = as.matrix(train_data[, -c((ncol(train_data) - 1):ncol(train_data))]),  # Exclude OS_Months and OS_Status
                      y = Surv(train_data$OS_MONTHS, train_data$OS_STATUS),
                      family = "cox",
                      nfolds = 10,
                      alpha=0)

# Extract optimal lambda
best_lambda_lasso <- cv.lasso$lambda.min  # Or lambda.1se for more conservative model
# best_lambda_lasso <- cv.lasso$lambda.1se

# Fit the Lasso model with the optimal lambda
lasso_fit <- glmnet(x = as.matrix(train_data[, -c((ncol(train_data) - 1):ncol(train_data))]),  # Exclude OS_Months and OS_Status
                    y = Surv(train_data$OS_MONTHS, train_data$OS_STATUS),
                    family = "cox",
                    lambda = best_lambda_lasso,
                    alpha=0)

cf_lasso<- coef(lasso_fit)
cf_lasso_df<- as.data.frame(as.matrix(cf_lasso))
cf_lasso_df_non_zero<- cf_lasso_df[!(cf_lasso_df$s0==0), , drop=F]
dim(cf_lasso_df_non_zero) #10 1

RS_list <- bplapply(df_list, function(x) {
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
  predict(cv.lasso, type = "link", newx = as.matrix(df_data[, -c((ncol(df_data)-1):ncol(df_data))]), s = cv.lasso$lambda.min)
}, BPPARAM= MulticoreParam(num_detected_cores))

cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "Univariate + LASSO + Ridge",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(genes, collapse = ","),
    No_of_genes_selected = length(genes)
  )
}, BPPARAM= MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])

######################### 3.9 Uni + LASSO + SuperPC  ################################

genes <- row.names(cf_lasso_df_non_zero)
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]
data <- list(x=t(train_data[,-c((ncol(train_data)-1):ncol(train_data))]),y=train_data$OS_MONTHS,censoring.status=train_data$OS_STATUS,featurenames=colnames(train_data)[-c((ncol(train_data)-1):ncol(train_data))])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                     n.fold = 10,
                     n.components=3,
                     min.features=2,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)

RS_list <- bplapply(df_list, function(w) {
  df_data <- w[, c(genes, "OS_MONTHS", "OS_STATUS")]
  df_data <- list(x=t(df_data[,-c((ncol(df_data)-1):ncol(df_data))]),y=df_data$OS_MONTHS,censoring.status=df_data$OS_STATUS,featurenames=colnames(df_data)[-c((ncol(df_data)-1):ncol(df_data))])
  ff <- superpc.predict(fit,data,df_data,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  as.numeric(ff$v.pred)
}, BPPARAM= MulticoreParam(num_detected_cores))

cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]                                                                                
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "Univariate + LASSO + SuperPC",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(genes, collapse = ","),
    No_of_genes_selected = length(genes)
  )
}, BPPARAM= MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])


########################## 3.10 Uni + LASSO + XGB ##################################
#gc()
## install.packages("mlr3")
#library(mlr3)        # For machine learning tasks
## remotes::install_github("mlr-org/mlr3extralearners")
#library(mlr3extralearners)  # For survival learners, such as xgboost Cox
#library(xgboost)
## remotes::install_github("mlr-org/mlr3proba")
#library(mlr3proba)

#genes <- row.names(cf_lasso_df_non_zero)
#train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

## Define the learner
#learner = mlr3::lrn("surv.xgboost.cox")
#print(learner)

#train_task <- as_task_surv(train_data, time="OS_MONTHS", event="OS_STATUS")

## install.packages("ml3tuning")
#library(mlr3tuning)
#library(paradox)

## Define a search space for hyperparameters
#search_space <- ps(
#  nrounds = p_int(50, 500),
#  eta = p_dbl(0.01, 0.3),
#  max_depth = p_int(2, 10),
#  min_child_weight = p_dbl(1, 10),
#  subsample = p_dbl(0.5, 1),
#  colsample_bytree = p_dbl(0.5, 1),
#  lambda = p_dbl(0, 1),
#  alpha = p_dbl(0, 1)
#)

## Define the tuning instance
#instance <- mlr3tuning::TuningInstanceBatchSingleCrit$new(
#  task = train_task,
#  learner = learner,
#  resampling = rsmp("cv", folds = 5),  # Cross-validation
#  measure = msr("surv.cindex"),       # Optimization metric
#  search_space = search_space,
#  terminator = trm("evals", n_evals = 100)  # Stop after 100 evaluations
#)

## Perform random search or grid search
#tuner <- tnr("random_search")  # Or "grid_search"
#tuner$optimize(instance)

## Use the best hyperparameters
#learner$param_set$values <- instance$result_learner_param_vals



## Create train and test indices
#set.seed(seed)  # Ensure reproducibility

## Train the learner on the training data
#learner$train(train_task, seq_len(nrow(train_data)))

## Print model details
#print(learner$model)

## Get feature importance
#print(learner$importance())

## Assuming `df_list` is a list of datasets you want to apply the model to:
#RS_list <- bplapply(df_list, function(x) {
#  # Prepare the data (ensure to include relevant columns)
#  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
#  predictions = learner$predict_newdata(df_data)
#  # Predict survival using the trained xgboost model
#  as.numeric(predictions$lp)
#}, BPPARAM = MulticoreParam(num_detected_cores))

## Calculate C-index, AUC, and PR-AUC for each dataset
#cc_list <- bplapply(names(df_list), function(x) {
#  RS <- RS_list[[x]]  # Get the predicted risk scores
#  df <- df_list[[x]]
#  # Calculate Concordance Index (C-index)
#  cc <- learner$predict_newdata(df)$score()
#  # Calculate AUC for specific time points
#  auc_values <- sapply(time_points, function(tp) {
#    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
#                           status = df$OS_STATUS,
#                           marker = RS,
#                           predict.time = tp,
#                           method = "KM") # or "NNE", method for handling censoring
#    return(roc_obj$AUC)
#  })
#  names(auc_values) <- paste0("AUC_", time_points, "months")  # Naming AUC values
#  
#  # Calculate Time-dependent PR-AUC at specified time points
#  pr_auc_values <- sapply(time_points, function(tp) {
#    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
#    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
#    return(pr_obj$auc.integral)
#  })
#  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months")  # Naming PR-AUC values
#  
#  # Create a data frame with results
#  data.frame(
#    Model_combination = "Univariate + XGBoost",
#    Dataset_Name = x,
#    C_index = as.numeric(cc),
#    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
#    AUC_36months = as.numeric(auc_values["AUC_36months"]),
#    AUC_60months = as.numeric(auc_values["AUC_60months"]),
#    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
#    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
#    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
#    gene_set = paste(genes, collapse = ","),
#    No_of_genes_selected = length(genes)
#  )
#}, BPPARAM = MulticoreParam(num_detected_cores))

## Combine the results into a single dataframe
#results <- bind_rows(results, bind_rows(cc_list))

## Print results
#print(results[, 1:3])



######################### 4.1 Uni + CoxBoost + LASSO  ################################

gc()
# library(devtools)
# devtools::install_github("binderh/CoxBoost")
# install.packages("snowfall")
library(snowfall)
library(CoxBoost)
set.seed(seed)
genes <- row.names(res_mod)
train_data <- train[, names(train) %in% c(genes, "OS_STATUS", "OS_MONTHS")]

# CoxBoost in snowfall computing
sfInit(parallel = TRUE, cpus = num_detected_cores)
sfLibrary(CoxBoost) # Make CoxBoost library available on workers
sfExport("train_data") # Export train_data to worker environments
sfExport("genes") # Export genes to worker environments
sfExport("seed")
pen <- optimCoxBoostPenalty(train_data[,'OS_MONTHS'],train_data[,'OS_STATUS'], as.matrix(train_data[, !(names(train_data) %in% c("OS_MONTHS","OS_STATUS"))]),
                            trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(train_data[,'OS_MONTHS'],train_data[,'OS_STATUS'], as.matrix(train_data[, !(names(train_data) %in% c("OS_MONTHS","OS_STATUS"))]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty, parallel = T)
fit <- CoxBoost(train_data[,'OS_MONTHS'],train_data[,'OS_STATUS'], as.matrix(train_data[, !(names(train_data) %in% c("OS_MONTHS","OS_STATUS"))]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)
sfStop() # stop snowfall processing

## Genes selected
nonzero_coef <- coef(fit, type = "nonzero")
genes <- names(nonzero_coef[nonzero_coef != 0])

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

# genes <- row.names(cf_lasso_df_non_zero)
# rm(RS_list)
RS_list <- bplapply(df_list, function(x) {
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
  # newx <- as.matrix(df_data[, -c((ncol(df_data) - 1), ncol(df_data))])
  predict(cv.lasso, type = "link", newx = as.matrix(df_data[, -c((ncol(df_data)-1):ncol(df_data))]), s = cv.lasso$lambda.min)
}, BPPARAM = MulticoreParam(num_detected_cores))


# rm(cc_list)
cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "Univariate + CoxBoost + LASSO",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(genes, collapse = ","),
    No_of_genes_selected = length(genes)
  )
}, BPPARAM = MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])

######################### 4.2 Uni + CoxBoost + StepCox  ################################

genes <- names(nonzero_coef[nonzero_coef != 0])
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

step.model <- step(coxph(Surv(OS_MONTHS,OS_STATUS)~.,train_data), direction = "both", trace = FALSE)
summary(step.model)
res.cox.coeff <- summary(step.model)$coefficients

genes <- row.names(res.cox.coeff)

RS_list <- bplapply(df_list, function(x) {
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
  predict(step.model, type = "risk", newdata = df_data)
}, BPPARAM= MulticoreParam(num_detected_cores))

cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]                                                                                
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "Univariate + CoxBoost + StepCox",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(genes, collapse = ","),
    No_of_genes_selected = length(genes)
  )
}, BPPARAM= MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])

######################### 4.3 Uni + CoxBoost + RSF  ################################

genes <- names(nonzero_coef[nonzero_coef != 0])
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

rf_model_final <- rfsrc(Surv(OS_MONTHS, OS_STATUS) ~ ., data = train_data,
                        importance = T, ntree = ntree, # Use full ntree for final model
                        proximity = T, forest = T,
                        seed = seed, nodesize = nodesize,
                        bootstrap = "by.root", ncores= num_detected_cores, nodedepth = nodedepth) # Use topvars for final model
vs.rf_model <- var.select(object = rf_model_final, conservative = "high")
topvars <- vs.rf_model$topvars
length(topvars)

# --- Evaluate Final Model on Test Datasets (test1, test2) ---
RS_list_final <- lapply(df_list, function(x) { # Use df_list to evaluate on train, test1, test2
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")] # Use 'genes' or 'topvars' consistently
  predict(rf_model_final, newdata = df_data)$predicted
})

cc_list_final <- lapply(names(df_list), function(x) {
  RS <- RS_list_final[[x]]
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "Uni + CoxBoost + RSF",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(topvars, collapse = ","),
    No_of_genes_selected = length(topvars)
  )
})

results <- bind_rows(results, bind_rows(cc_list_final)) # Append to results dataframe

print(results[, 1:3])

######################### 4.4 Uni + CoxBoost + SuperPC  ################################

genes <- names(nonzero_coef[nonzero_coef != 0])
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

data <- list(x=t(train_data[,-c((ncol(train_data)-1):ncol(train_data))]),y=train_data$OS_MONTHS,censoring.status=train_data$OS_STATUS,featurenames=colnames(train_data)[-c((ncol(train_data)-1):ncol(train_data))])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                     n.fold = 10,
                     n.components=3,
                     min.features=2,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)

RS_list <- bplapply(df_list, function(w) {
  df_data <- w[, c(genes, "OS_MONTHS", "OS_STATUS")]
  df_data <- list(x=t(df_data[,-c((ncol(df_data)-1):ncol(df_data))]),y=df_data$OS_MONTHS,censoring.status=df_data$OS_STATUS,featurenames=colnames(df_data)[-c((ncol(df_data)-1):ncol(df_data))])
  ff <- superpc.predict(fit,data,df_data,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  as.numeric(ff$v.pred)
}, BPPARAM= MulticoreParam(num_detected_cores))

cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]                                                                                
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "Univariate + CoxBoost + SuperPC",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(genes, collapse = ","),
    No_of_genes_selected = length(genes)
  )
}, BPPARAM= MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])

######################### 4.5 Uni + CoxBoost + GBM  ################################

library(gbm)
genes <- names(nonzero_coef[ nonzero_coef !=0])
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

gbm_model <- gbm(Surv(OS_MONTHS, OS_STATUS) ~ ., data = train_data, distribution = "coxph",
                 n.trees = 10000, interaction.depth = 5, shrinkage = 0.001,
                 cv.folds = 10, n.cores = num_detected_cores, verbose = FALSE)
best <- which.min(gbm_model$cv.error)
fit <- gbm(formula = Surv(OS_MONTHS,OS_STATUS)~.,data = train_data,distribution = 'coxph',
           n.trees = best,
           interaction.depth = 5,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = num_detected_cores)

RS_list <- bplapply(df_list, function(x) {
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
  predict(fit, df_data, n.trees=best, type="link")
}, BPPARAM= MulticoreParam(num_detected_cores))

cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]                                                                                
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS,OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "Univariate + CoxBoost + GBM",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(genes, collapse = ","),
    No_of_genes_selected = length(genes)
  )
}, BPPARAM= MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])


######################### 4.6 Uni + CoxBoost + plsRcox  ################################

# plsRcox
library(plsRcox)
genes <- names(nonzero_coef[nonzero_coef != 0])

train_data <- train[, names(train) %in% c(genes, "OS_STATUS", "OS_MONTHS")]

set.seed(seed)
cv.plsRcox.res=cv.plsRcox(list(x=train_data[,genes],time=train_data$OS_MONTHS,status=train_data$OS_STATUS),nt=5,verbose = FALSE)
fit <- plsRcox(train_data[,genes],time=train_data$OS_MONTHS,event=train_data$OS_STATUS,nt=as.numeric(cv.plsRcox.res[5]))

RS_list <- bplapply(df_list, function(x) {
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
  as.numeric(predict(fit, newdata= df_data[, !(names(df_data) %in% c("OS_MONTHS","OS_STATUS"))], type="lp"))
}, BPPARAM= MulticoreParam(num_detected_cores))

cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]                                                                                
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS,OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "Univariate + CoxBoost + plsRcox",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(genes, collapse = ","),
    No_of_genes_selected = length(genes)
  )
}, BPPARAM= MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])


######################### 4.7 Uni + CoxBoost + Ridge  ################################

gc()
genes <- names(nonzero_coef[nonzero_coef != 0])
data <-train[, colnames(train) %in% c(genes, "OS_STATUS", "OS_MONTHS")]

cv.lasso <- cv.glmnet(x = as.matrix(data[, -c((ncol(data) - 1):ncol(data))]),  # Exclude OS_Months and OS_Status
                      y = Surv(data$OS_MONTHS, data$OS_STATUS),
                      family = "cox",
                      nfolds = 10,
                      alpha=0)

# Extract optimal lambda
best_lambda_lasso <- cv.lasso$lambda.min  # Or lambda.1se for more conservative model
# best_lambda_lasso <- cv.lasso$lambda.1se

# Fit the Lasso model with the optimal lambda
lasso_fit <- glmnet(x = as.matrix(data[, -c((ncol(data) - 1):ncol(data))]),  # Exclude OS_Months and OS_Status
                    y = Surv(data$OS_MONTHS, data$OS_STATUS),
                    family = "cox",
                    lambda = best_lambda_lasso,
                    alpha=0)

cf_lasso<- coef(lasso_fit)
cf_lasso_df<- as.data.frame(as.matrix(cf_lasso))
cf_lasso_df_non_zero<- cf_lasso_df[!(cf_lasso_df$s0==0), , drop=F]
dim(cf_lasso_df_non_zero) #10 1

RS_list <- bplapply(df_list, function(x) {
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
  predict(cv.lasso, type = "link", newx = as.matrix(df_data[, -c((ncol(df_data)-1):ncol(df_data))]), s = cv.lasso$lambda.min)
}, BPPARAM= MulticoreParam(num_detected_cores))

cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "Univariate + CoxBoost + Ridge",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(genes, collapse = ","),
    No_of_genes_selected = length(genes)
  )
}, BPPARAM= MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])


######################### 4.8 Uni + CoxBoost +SVM  ################################
# gc()
# 
# genes <- names(nonzero_coef)
# train_data <- train[, names(train) %in% c(genes, "OS_STATUS", "OS_MONTHS")]
# set.seed(seed)
# fit = survivalsvm(Surv(OS_MONTHS, OS_STATUS) ~., data= train_data, gamma.mu = 1)
# 
# RS_list <- bplapply(df_list, function(x) {
#   df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
#   as.numeric(predict(fit, df_data)$predicted)
# }, BPPARAM= MulticoreParam(num_detected_cores))
# 
# cc_list <- bplapply(names(df_list), function(x) {
#   RS <- RS_list[[x]]
#   df <- df_list[[x]]
#   cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
#   auc_values <- sapply(time_points, function(tp) {
#     roc_obj <- survivalROC(Stime = df$OS_MONTHS,
#                            status = df$OS_STATUS,
#                            marker = RS,
#                            predict.time = tp,
#                            method = "KM") # or "NNE", method for handling censoring
#     return(roc_obj$AUC)
#   })
#   names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
#   # Calculate Time-dependent PR-AUC at specified time points
#   pr_auc_values <- sapply(time_points, function(tp) {
#     binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
#     pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
#     return(pr_obj$auc.integral)
#   })
#   names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
# 
#   data.frame(
#     Model_combination = "Univariate + CoxBoost + SVM",
#     Dataset_Name = x,
#     C_index = as.numeric(cc),
#     AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
#     AUC_36months = as.numeric(auc_values["AUC_36months"]),
#     AUC_60months = as.numeric(auc_values["AUC_60months"]),
#     PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
#     PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
#     PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
#     gene_set = paste(genes, collapse = ","),
#     No_of_genes_selected = length(genes)
#   )
# }, BPPARAM= MulticoreParam(num_detected_cores))
# 
# results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe
# 
# print(results[, 1:3])
# 


######################### 4.9 Uni + CoxBoost + XGB ##################################
#gc()
## install.packages("mlr3")
#library(mlr3)        # For machine learning tasks
## remotes::install_github("mlr-org/mlr3extralearners")
#library(mlr3extralearners)  # For survival learners, such as xgboost Cox
#library(xgboost)
## remotes::install_github("mlr-org/mlr3proba")
#library(mlr3proba)

#genes <- names(nonzero_coef[nonzero_coef != 0])
#train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

## Define the learner
#learner = mlr3::lrn("surv.xgboost.cox")
#print(learner)

#train_task <- as_task_surv(train_data, time="OS_MONTHS", event="OS_STATUS")

## install.packages("ml3tuning")
#library(mlr3tuning)
#library(paradox)

## Define a search space for hyperparameters
#search_space <- ps(
#  nrounds = p_int(50, 500),
#  eta = p_dbl(0.01, 0.3),
#  max_depth = p_int(2, 10),
#  min_child_weight = p_dbl(1, 10),
#  subsample = p_dbl(0.5, 1),
#  colsample_bytree = p_dbl(0.5, 1),
#  lambda = p_dbl(0, 1),
#  alpha = p_dbl(0, 1)
#)

## Define the tuning instance
#instance <- mlr3tuning::TuningInstanceBatchSingleCrit$new(
#  task = train_task,
#  learner = learner,
#  resampling = rsmp("cv", folds = 5),  # Cross-validation
#  measure = msr("surv.cindex"),       # Optimization metric
#  search_space = search_space,
#  terminator = trm("evals", n_evals = 100)  # Stop after 100 evaluations
#)

## Perform random search or grid search
#tuner <- tnr("random_search")  # Or "grid_search"
#tuner$optimize(instance)

## Use the best hyperparameters
#learner$param_set$values <- instance$result_learner_param_vals



## Create train and test indices
#set.seed(seed)  # Ensure reproducibility

## Train the learner on the training data
#learner$train(train_task, seq_len(nrow(train_data)))

## Print model details
#print(learner$model)

## Get feature importance
#print(learner$importance())

## Assuming `df_list` is a list of datasets you want to apply the model to:
#RS_list <- bplapply(df_list, function(x) {
#  # Prepare the data (ensure to include relevant columns)
#  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
#  predictions = learner$predict_newdata(df_data)
#  # Predict survival using the trained xgboost model
#  as.numeric(predictions$lp)
#}, BPPARAM = MulticoreParam(num_detected_cores))

## Calculate C-index, AUC, and PR-AUC for each dataset
#cc_list <- bplapply(names(df_list), function(x) {
#  RS <- RS_list[[x]]  # Get the predicted risk scores
#  df <- df_list[[x]]
#  # Calculate Concordance Index (C-index)
#  cc <- learner$predict_newdata(df)$score()
#  # Calculate AUC for specific time points
#  auc_values <- sapply(time_points, function(tp) {
#    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
#                           status = df$OS_STATUS,
#                           marker = RS,
#                           predict.time = tp,
#                           method = "KM") # or "NNE", method for handling censoring
#    return(roc_obj$AUC)
#  })
#  names(auc_values) <- paste0("AUC_", time_points, "months")  # Naming AUC values
#  
#  # Calculate Time-dependent PR-AUC at specified time points
#  pr_auc_values <- sapply(time_points, function(tp) {
#    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
#    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
#    return(pr_obj$auc.integral)
#  })
#  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months")  # Naming PR-AUC values
#  
#  # Create a data frame with results
#  data.frame(
#    Model_combination = "Univariate + XGBoost",
#    Dataset_Name = x,
#    C_index = as.numeric(cc),
#    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
#    AUC_36months = as.numeric(auc_values["AUC_36months"]),
#    AUC_60months = as.numeric(auc_values["AUC_60months"]),
#    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
#    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
#    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
#    gene_set = paste(genes, collapse = ","),
#    No_of_genes_selected = length(genes)
#  )
#}, BPPARAM = MulticoreParam(num_detected_cores))

## Combine the results into a single dataframe
#results <- bind_rows(results, bind_rows(cc_list))

## Print results
#print(results[, 1:3])



######################### 5.1 Uni + RSF + SVM  ################################
gc()
genes <- row.names(res_mod)
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

rf_model_final <- rfsrc(Surv(OS_MONTHS, OS_STATUS) ~ ., data = train_data,
                        importance = T, ntree = ntree, # Use full ntree for final model
                        proximity = T, forest = T,
                        seed = seed, nodesize = nodesize,
                        bootstrap = "by.root", ncores= num_detected_cores, nodedepth = nodedepth) # Use topvars for final model
vs.rf_model <- var.select(object = rf_model_final, conservative = "high")
topvars <- vs.rf_model$topvars
length(topvars)
# 
# genes <- topvars
# train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]
# gc()
# set.seed(seed)
# fit = survivalsvm(Surv(OS_MONTHS, OS_STATUS) ~., data= train_data, gamma.mu = 1)
# 
# RS_list <- bplapply(df_list, function(x) {
#   df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
#   as.numeric(predict(fit, df_data)$predicted)
# }, BPPARAM= MulticoreParam(num_detected_cores))
# 
# cc_list <- bplapply(names(df_list), function(x) {
#   RS <- RS_list[[x]]                                                                                
#   df <- df_list[[x]]
#   cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
#   auc_values <- sapply(time_points, function(tp) {
#     roc_obj <- survivalROC(Stime = df$OS_MONTHS,
#                            status = df$OS_STATUS,
#                            marker = RS,
#                            predict.time = tp,
#                            method = "KM") # or "NNE", method for handling censoring
#     return(roc_obj$AUC)
#   })
#   names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
#   # Calculate Time-dependent PR-AUC at specified time points
#   pr_auc_values <- sapply(time_points, function(tp) {
#     binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
#     pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
#     return(pr_obj$auc.integral)
#   })
#   names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
#   
#   data.frame(
#     Model_combination = "Univariate + RSF + SVM",
#     Dataset_Name = x,
#     C_index = as.numeric(cc),
#     AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
#     AUC_36months = as.numeric(auc_values["AUC_36months"]),
#     AUC_60months = as.numeric(auc_values["AUC_60months"]),
#     PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
#     PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
#     PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
#     gene_set = paste(genes, collapse = ","),
#     No_of_genes_selected = length(genes)
#   )
# }, BPPARAM= MulticoreParam(num_detected_cores))
# 
# results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe
# 
# print(results[, 1:3])


######################### 5.2 Uni + RSF + LASSO  ################################

genes <- topvars
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]
set.seed(seed)
cv.lasso <- cv.glmnet(x = as.matrix(train_data[, -c((ncol(train_data) - 1):ncol(train_data))]),  # Exclude OS_Months and OS_Status
                      y = Surv(train_data$OS_MONTHS, train_data$OS_STATUS),
                      family = "cox",
                      nfolds = 10,
                      alpha=1)

# Extract optimal lambda
best_lambda_lasso <- cv.lasso$lambda.min  # Or lambda.1se for more conservative model
# best_lambda_lasso <- cv.lasso$lambda.1se

# Fit the Lasso model with the optimal lambda
lasso_fit <- glmnet(x = as.matrix(train_data[, -c((ncol(train_data) - 1):ncol(train_data))]),  # Exclude OS_Months and OS_Status
                    y = Surv(train_data$OS_MONTHS, train_data$OS_STATUS),
                    family = "cox",
                    lambda = best_lambda_lasso,
                    alpha=1)

cf_lasso<- coef(lasso_fit)
cf_lasso_df<- as.data.frame(as.matrix(cf_lasso))
cf_lasso_df_non_zero<- cf_lasso_df[!(cf_lasso_df$s0==0), , drop=F]
dim(cf_lasso_df_non_zero) #10 1

selelected_genes <- row.names(cf_lasso_df_non_zero)
RS_list <- bplapply(df_list, function(x) {
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
  # newx <- as.matrix(df_data[, -c((ncol(df_data) - 1), ncol(df_data))])
  predict(cv.lasso, type = "link", newx = as.matrix(df_data[, -c((ncol(df_data)-1):ncol(df_data))]), s = cv.lasso$lambda.min)
}, BPPARAM = MulticoreParam(num_detected_cores))


# rm(cc_list)
cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "Univariate + RSF + LASSO",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(selelected_genes, collapse = ","),
    No_of_genes_selected = length(selelected_genes)
  )
}, BPPARAM = MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])

######################### 5.3 Uni + RSF + StepCox  ################################

genes <- topvars
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]
step.model <- step(coxph(Surv(OS_MONTHS,OS_STATUS)~.,train_data), direction = "both", trace = FALSE)
summary(step.model)
res.cox.coeff <- summary(step.model)$coefficients

genes <- row.names(res.cox.coeff)

RS_list <- bplapply(df_list, function(x) {
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
  predict(step.model, type = "risk", newdata = df_data)
}, BPPARAM= MulticoreParam(num_detected_cores))

cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]                                                                                
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "Univariate + RSF + StepCox",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(genes, collapse = ","),
    No_of_genes_selected = length(genes)
  )
}, BPPARAM= MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])

######################### 5.4 Uni + RSF + CoxBoost  ################################

genes <- topvars
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

sfInit(parallel = TRUE, cpus = num_detected_cores)
sfLibrary(CoxBoost) # Make CoxBoost library available on workers
sfExport("train_data") # Export train_data to worker environments
sfExport("genes") # Export genes to worker environments

pen <- optimCoxBoostPenalty(train_data[,'OS_MONTHS'],train_data[,'OS_STATUS'], as.matrix(train_data[, !(names(train_data) %in% c("OS_MONTHS","OS_STATUS"))]),
                            trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(train_data[,'OS_MONTHS'],train_data[,'OS_STATUS'], as.matrix(train_data[, !(names(train_data) %in% c("OS_MONTHS","OS_STATUS"))]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
fit <- CoxBoost(train_data[,'OS_MONTHS'],train_data[,'OS_STATUS'], as.matrix(train_data[, !(names(train_data) %in% c("OS_MONTHS","OS_STATUS"))]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)
sfStop()

nonzero_coef <- coef(fit, type = "nonzero")
genes_selected <- names(nonzero_coef[nonzero_coef != 0])

RS_list <- bplapply(df_list, function(x) {
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
  as.numeric(predict(fit, type = "lp", newdata = df_data[, !(names(df_data) %in% c("OS_MONTHS","OS_STATUS"))], newtime= df_data$OS_MONTHS, newstatus=df_data$OS_STATUS))
}, BPPARAM= MulticoreParam(num_detected_cores))

cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]                                                                                
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "Univariate + RSF + CoxBoost",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(genes_selected, collapse = ","),
    No_of_genes_selected = length(genes_selected)
  )
}, BPPARAM= MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])

######################### 5.5 Uni + RSF + GBM  ################################

genes <- topvars
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]
gbm_model <- gbm(Surv(OS_MONTHS, OS_STATUS) ~ ., data = train_data, distribution = "coxph",
                 n.trees = 10000, interaction.depth = 5, shrinkage = 0.001,
                 cv.folds = 10, n.cores = num_detected_cores, verbose = FALSE)
best <- which.min(gbm_model$cv.error)
fit <- gbm(formula = Surv(OS_MONTHS,OS_STATUS)~.,data = train_data,distribution = 'coxph',
           n.trees = best,
           interaction.depth = 5,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = num_detected_cores)

RS_list <- bplapply(df_list, function(x) {
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
  predict(fit, df_data, n.trees=best, type="link")
}, BPPARAM= MulticoreParam(num_detected_cores))

cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]                                                                                
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS,OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "Univariate + RSF + GBM",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(genes, collapse = ","),
    No_of_genes_selected = length(genes)
  )
}, BPPARAM= MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])


######################### 5.6 Uni + RSF + plsRcox  ################################

genes <- topvars
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

set.seed(seed)
cv.plsRcox.res=cv.plsRcox(list(x=train_data[,genes],time=train_data$OS_MONTHS,status=train_data$OS_STATUS),nt=5,verbose = FALSE)
fit <- plsRcox(train_data[,genes],time=train_data$OS_MONTHS,event=train_data$OS_STATUS,nt=as.numeric(cv.plsRcox.res[5]))

RS_list <- bplapply(df_list, function(x) {
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
  as.numeric(predict(fit, newdata= df_data[, !(names(df_data) %in% c("OS_MONTHS","OS_STATUS"))], type="lp"))
}, BPPARAM= MulticoreParam(num_detected_cores))

cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]                                                                                
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS,OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "Univariate + RSF + plsRcox",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(genes, collapse = ","),
    No_of_genes_selected = length(genes)
  )
}, BPPARAM= MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])


######################### 5.7 Uni + RSF + Ridge  ################################

genes <- topvars
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]
cv.lasso <- cv.glmnet(x = as.matrix(train_data[, -c((ncol(train_data) - 1):ncol(train_data))]),  # Exclude OS_Months and OS_Status
                      y = Surv(train_data$OS_MONTHS, train_data$OS_STATUS),
                      family = "cox",
                      nfolds = 10,
                      alpha=0)

# Extract optimal lambda
best_lambda_lasso <- cv.lasso$lambda.min  # Or lambda.1se for more conservative model
# best_lambda_lasso <- cv.lasso$lambda.1se

# Fit the Lasso model with the optimal lambda
lasso_fit <- glmnet(x = as.matrix(train_data[, -c((ncol(train_data) - 1):ncol(train_data))]),  # Exclude OS_Months and OS_Status
                    y = Surv(train_data$OS_MONTHS, train_data$OS_STATUS),
                    family = "cox",
                    lambda = best_lambda_lasso,
                    alpha=0)

cf_lasso<- coef(lasso_fit)
cf_lasso_df<- as.data.frame(as.matrix(cf_lasso))
cf_lasso_df_non_zero<- cf_lasso_df[!(cf_lasso_df$s0==0), , drop=F]
dim(cf_lasso_df_non_zero) #10 1

RS_list <- bplapply(df_list, function(x) {
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
  predict(cv.lasso, type = "link", newx = as.matrix(df_data[, -c((ncol(df_data)-1):ncol(df_data))]), s = cv.lasso$lambda.min)
}, BPPARAM= MulticoreParam(num_detected_cores))

cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "Univariate + RSF + Ridge",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(genes, collapse = ","),
    No_of_genes_selected = length(genes)
  )
}, BPPARAM= MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])





######################### 5.8 Uni + RSF + SuperPC  ################################

genes <- topvars
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]
data <- list(x=t(train_data[,-c((ncol(train_data)-1):ncol(train_data))]),y=train_data$OS_MONTHS,censoring.status=train_data$OS_STATUS,featurenames=colnames(train_data)[-c((ncol(train_data)-1):ncol(train_data))])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                     n.fold = 10,
                     n.components=3,
                     min.features=2,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)

RS_list <- bplapply(df_list, function(w) {
  df_data <- w[, c(genes, "OS_MONTHS", "OS_STATUS")]
  df_data <- list(x=t(df_data[,-c((ncol(df_data)-1):ncol(df_data))]),y=df_data$OS_MONTHS,censoring.status=df_data$OS_STATUS,featurenames=colnames(df_data)[-c((ncol(df_data)-1):ncol(df_data))])
  ff <- superpc.predict(fit,data,df_data,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  as.numeric(ff$v.pred)
}, BPPARAM= MulticoreParam(num_detected_cores))

cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]                                                                                
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "Univariate + RSF + SuperPC",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(genes, collapse = ","),
    No_of_genes_selected = length(genes)
  )
}, BPPARAM= MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])



########################## 5.9 Uni + RSF + XGB ##################################
#gc()
## install.packages("mlr3")
#library(mlr3)        # For machine learning tasks
## remotes::install_github("mlr-org/mlr3extralearners")
#library(mlr3extralearners)  # For survival learners, such as xgboost Cox
#library(xgboost)
## remotes::install_github("mlr-org/mlr3proba")
#library(mlr3proba)

#genes <- topvars
#train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

## Define the learner
#learner = mlr3::lrn("surv.xgboost.cox")
#print(learner)

#train_task <- as_task_surv(train_data, time="OS_MONTHS", event="OS_STATUS")

## install.packages("ml3tuning")
#library(mlr3tuning)
#library(paradox)

## Define a search space for hyperparameters
#search_space <- ps(
#  nrounds = p_int(50, 500),
#  eta = p_dbl(0.01, 0.3),
#  max_depth = p_int(2, 10),
#  min_child_weight = p_dbl(1, 10),
#  subsample = p_dbl(0.5, 1),
#  colsample_bytree = p_dbl(0.5, 1),
#  lambda = p_dbl(0, 1),
#  alpha = p_dbl(0, 1)
#)

## Define the tuning instance
#instance <- mlr3tuning::TuningInstanceBatchSingleCrit$new(
#  task = train_task,
#  learner = learner,
#  resampling = rsmp("cv", folds = 5),  # Cross-validation
#  measure = msr("surv.cindex"),       # Optimization metric
#  search_space = search_space,
#  terminator = trm("evals", n_evals = 100)  # Stop after 100 evaluations
#)

## Perform random search or grid search
#tuner <- tnr("random_search")  # Or "grid_search"
#tuner$optimize(instance)

## Use the best hyperparameters
#learner$param_set$values <- instance$result_learner_param_vals



## Create train and test indices
#set.seed(seed)  # Ensure reproducibility

## Train the learner on the training data
#learner$train(train_task, seq_len(nrow(train_data)))

## Print model details
#print(learner$model)

## Get feature importance
#print(learner$importance())

## Assuming `df_list` is a list of datasets you want to apply the model to:
#RS_list <- bplapply(df_list, function(x) {
#  # Prepare the data (ensure to include relevant columns)
#  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
#  predictions = learner$predict_newdata(df_data)
#  # Predict survival using the trained xgboost model
#  as.numeric(predictions$lp)
#}, BPPARAM = MulticoreParam(num_detected_cores))

## Calculate C-index, AUC, and PR-AUC for each dataset
#cc_list <- bplapply(names(df_list), function(x) {
#  RS <- RS_list[[x]]  # Get the predicted risk scores
#  df <- df_list[[x]]
#  # Calculate Concordance Index (C-index)
#  cc <- learner$predict_newdata(df)$score()
#  # Calculate AUC for specific time points
#  auc_values <- sapply(time_points, function(tp) {
#    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
#                           status = df$OS_STATUS,
#                           marker = RS,
#                           predict.time = tp,
#                           method = "KM") # or "NNE", method for handling censoring
#    return(roc_obj$AUC)
#  })
#  names(auc_values) <- paste0("AUC_", time_points, "months")  # Naming AUC values
#  
#  # Calculate Time-dependent PR-AUC at specified time points
#  pr_auc_values <- sapply(time_points, function(tp) {
#    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
#    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
#    return(pr_obj$auc.integral)
#  })
#  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months")  # Naming PR-AUC values
#  
#  # Create a data frame with results
#  data.frame(
#    Model_combination = "Univariate + XGBoost",
#    Dataset_Name = x,
#    C_index = as.numeric(cc),
#    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
#    AUC_36months = as.numeric(auc_values["AUC_36months"]),
#    AUC_60months = as.numeric(auc_values["AUC_60months"]),
#    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
#    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
#    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
#    gene_set = paste(genes, collapse = ","),
#    No_of_genes_selected = length(genes)
#  )
#}, BPPARAM = MulticoreParam(num_detected_cores))

## Combine the results into a single dataframe
#results <- bind_rows(results, bind_rows(cc_list))

## Print results
#print(results[, 1:3])




######################### 6.1 Uni + StepCox + GBM  ################################

#StepCox
genes <- row.names(res_mod)
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]
step.model <- step(coxph(Surv(OS_MONTHS,OS_STATUS)~.,train_data), direction = "both", trace = FALSE)
summary(step.model)
res.cox.coeff <- summary(step.model)$coefficients

gc()
# GBM
genes <- row.names(res.cox.coeff)
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

gbm_model <- gbm(Surv(OS_MONTHS, OS_STATUS) ~ ., data = train_data, distribution = "coxph",
                 n.trees = 10000, interaction.depth = 5, shrinkage = 0.001,
                 cv.folds = 10, n.cores = num_detected_cores, verbose = FALSE)
best <- which.min(gbm_model$cv.error)
fit <- gbm(formula = Surv(OS_MONTHS,OS_STATUS)~.,data = train_data,distribution = 'coxph',
           n.trees = best,
           interaction.depth = 5,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = num_detected_cores)

RS_list <- bplapply(df_list, function(x) {
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
  predict(fit, df_data, n.trees=best, type="link")
}, BPPARAM= MulticoreParam(num_detected_cores))

cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]                                                                                
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS,OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "Univariate + StepCox + GBM",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(genes, collapse = ","),
    No_of_genes_selected = length(genes)
  )
}, BPPARAM= MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])

######################### 6.2 Uni + StepCox + LASSO  ################################

genes <- row.names(res.cox.coeff)
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]
cv.lasso <- cv.glmnet(x = as.matrix(train_data[, -c((ncol(train_data) - 1):ncol(train_data))]),  # Exclude OS_Months and OS_Status
                      y = Surv(train_data$OS_MONTHS, train_data$OS_STATUS),
                      family = "cox",
                      nfolds = 10,
                      alpha=1)

# Extract optimal lambda
best_lambda_lasso <- cv.lasso$lambda.min  # Or lambda.1se for more conservative model
# best_lambda_lasso <- cv.lasso$lambda.1se

# Fit the Lasso model with the optimal lambda
lasso_fit <- glmnet(x = as.matrix(train_data[, -c((ncol(train_data) - 1):ncol(train_data))]),  # Exclude OS_Months and OS_Status
                    y = Surv(train_data$OS_MONTHS, train_data$OS_STATUS),
                    family = "cox",
                    lambda = best_lambda_lasso,
                    alpha=1)

cf_lasso<- coef(lasso_fit)
cf_lasso_df<- as.data.frame(as.matrix(cf_lasso))
cf_lasso_df_non_zero<- cf_lasso_df[!(cf_lasso_df$s0==0), , drop=F]
dim(cf_lasso_df_non_zero) #10 1
selected_genes <- row.names(cf_lasso_df_non_zero)
# rm(RS_list)
RS_list <- bplapply(df_list, function(x) {
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
  # newx <- as.matrix(df_data[, -c((ncol(df_data) - 1), ncol(df_data))])
  predict(cv.lasso, type = "link", newx = as.matrix(df_data[, -c((ncol(df_data)-1):ncol(df_data))]), s = cv.lasso$lambda.min)
}, BPPARAM = MulticoreParam(num_detected_cores))


# rm(cc_list)
cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "Univariate + StepCox + LASSO",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(selected_genes, collapse = ","),
    No_of_genes_selected = length(selected_genes)
  )
}, BPPARAM = MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])

######################### 6.3 Uni + StepCox + CoxBoost  ################################
genes <- row.names(res.cox.coeff)
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

sfInit(parallel = TRUE, cpus = num_detected_cores)
sfLibrary(CoxBoost) # Make CoxBoost library available on workers
sfExport("train_data") # Export train_data to worker environments
sfExport("genes") # Export genes to worker environments

pen <- optimCoxBoostPenalty(train_data[,'OS_MONTHS'],train_data[,'OS_STATUS'], as.matrix(train_data[, !(names(train_data) %in% c("OS_MONTHS","OS_STATUS"))]),
                            trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(train_data[,'OS_MONTHS'],train_data[,'OS_STATUS'], as.matrix(train_data[, !(names(train_data) %in% c("OS_MONTHS","OS_STATUS"))]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
fit <- CoxBoost(train_data[,'OS_MONTHS'],train_data[,'OS_STATUS'], as.matrix(train_data[, !(names(train_data) %in% c("OS_MONTHS","OS_STATUS"))]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)
sfStop()

nonzero_coef <- coef(fit, type = "nonzero")
genes_selected <- names(nonzero_coef[nonzero_coef != 0])

RS_list <- bplapply(df_list, function(x) {
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
  as.numeric(predict(fit, type = "lp", newdata = df_data[, !(names(df_data) %in% c("OS_MONTHS","OS_STATUS"))], newtime= df_data$OS_MONTHS, newstatus=df_data$OS_STATUS))
}, BPPARAM= MulticoreParam(num_detected_cores))

cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]                                                                                
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "Univariate + StepCox + CoxBoost",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(genes_selected, collapse = ","),
    No_of_genes_selected = length(genes_selected)
  )
}, BPPARAM= MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])

######################### 6.4 Uni + StepCox + SVM  ################################
# genes <- row.names(res.cox.coeff)
# train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]
# 
# set.seed(seed)
# fit = survivalsvm(Surv(OS_MONTHS, OS_STATUS) ~., data= train_data, gamma.mu = 1)
# 
# RS_list <- bplapply(df_list, function(x) {
#   df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
#   as.numeric(predict(fit, df_data)$predicted)
# }, BPPARAM= MulticoreParam(num_detected_cores))
# 
# cc_list <- bplapply(names(df_list), function(x) {
#   RS <- RS_list[[x]]                                                                                
#   df <- df_list[[x]]
#   cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
#   auc_values <- sapply(time_points, function(tp) {
#     roc_obj <- survivalROC(Stime = df$OS_MONTHS,
#                            status = df$OS_STATUS,
#                            marker = RS,
#                            predict.time = tp,
#                            method = "KM") # or "NNE", method for handling censoring
#     return(roc_obj$AUC)
#   })
#   names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
#   # Calculate Time-dependent PR-AUC at specified time points
#   pr_auc_values <- sapply(time_points, function(tp) {
#     binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
#     pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
#     return(pr_obj$auc.integral)
#   })
#   names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
#   
#   data.frame(
#     Model_combination = "Univariate + StepCox + SVM",
#     Dataset_Name = x,
#     C_index = as.numeric(cc),
#     AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
#     AUC_36months = as.numeric(auc_values["AUC_36months"]),
#     AUC_60months = as.numeric(auc_values["AUC_60months"]),
#     PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
#     PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
#     PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
#     gene_set = paste(genes, collapse = ","),
#     No_of_genes_selected = length(genes)
#   )
# }, BPPARAM= MulticoreParam(num_detected_cores))
# 
# results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe
# 
# print(results[, 1:3])

######################### 6.5 Uni + StepCox + RSF  ################################
genes <- row.names(res.cox.coeff)
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

rf_model_final <- rfsrc(Surv(OS_MONTHS, OS_STATUS) ~ ., data = train_data,
                        importance = T, ntree = ntree, # Use full ntree for final model
                        proximity = T, forest = T,
                        seed = seed, nodesize = nodesize,
                        bootstrap = "by.root", ncores= num_detected_cores, nodedepth = nodedepth) # Use topvars for final model
vs.rf_model <- var.select(object = rf_model_final, conservative = "high")
topvars <- vs.rf_model$topvars
length(topvars)

# --- Evaluate Final Model on Test Datasets (test1, test2) ---
RS_list_final <- lapply(df_list, function(x) { # Use df_list to evaluate on train, test1, test2
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")] # Use 'genes' or 'topvars' consistently
  predict(rf_model_final, newdata = df_data)$predicted
})

cc_list_final <- lapply(names(df_list), function(x) {
  RS <- RS_list_final[[x]]
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "Uni + StepCox + RSF",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(topvars, collapse = ","),
    No_of_genes_selected = length(topvars)
  )
})

results <- bind_rows(results, bind_rows(cc_list_final)) # Append to results dataframe

print(results[, 1:3])

######################### 6.6 Uni + StepCox + plsRcox  ################################

genes <- row.names(res.cox.coeff)
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

set.seed(seed)
cv.plsRcox.res=cv.plsRcox(list(x=train_data[,genes],time=train_data$OS_MONTHS,status=train_data$OS_STATUS),nt=5,verbose = FALSE)
fit <- plsRcox(train_data[,genes],time=train_data$OS_MONTHS,event=train_data$OS_STATUS,nt=as.numeric(cv.plsRcox.res[5]))

RS_list <- bplapply(df_list, function(x) {
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
  as.numeric(predict(fit, newdata= df_data[, !(names(df_data) %in% c("OS_MONTHS","OS_STATUS"))], type="lp"))
}, BPPARAM= MulticoreParam(num_detected_cores))

cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]                                                                                
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS,OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "Univariate + StepCox + plsRcox",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(genes, collapse = ","),
    No_of_genes_selected = length(genes)
  )
}, BPPARAM= MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])

######################### 6.7 Uni + StepCox + Ridge  ################################
genes <- row.names(res.cox.coeff)
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

cv.lasso <- cv.glmnet(x = as.matrix(train_data[, -c((ncol(train_data) - 1):ncol(train_data))]),  # Exclude OS_Months and OS_Status
                      y = Surv(train_data$OS_MONTHS, train_data$OS_STATUS),
                      family = "cox",
                      nfolds = 10,
                      alpha=0)

# Extract optimal lambda
best_lambda_lasso <- cv.lasso$lambda.min  # Or lambda.1se for more conservative model
# best_lambda_lasso <- cv.lasso$lambda.1se

# Fit the Lasso model with the optimal lambda
lasso_fit <- glmnet(x = as.matrix(train_data[, -c((ncol(train_data) - 1):ncol(train_data))]),  # Exclude OS_Months and OS_Status
                    y = Surv(train_data$OS_MONTHS, train_data$OS_STATUS),
                    family = "cox",
                    lambda = best_lambda_lasso,
                    alpha=0)

cf_lasso<- coef(lasso_fit)
cf_lasso_df<- as.data.frame(as.matrix(cf_lasso))
cf_lasso_df_non_zero<- cf_lasso_df[!(cf_lasso_df$s0==0), , drop=F]
dim(cf_lasso_df_non_zero) #10 1

RS_list <- bplapply(df_list, function(x) {
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
  predict(cv.lasso, type = "link", newx = as.matrix(df_data[, -c((ncol(df_data)-1):ncol(df_data))]), s = cv.lasso$lambda.min)
}, BPPARAM= MulticoreParam(num_detected_cores))

cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "Univariate + StepCox + Ridge",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(genes, collapse = ","),
    No_of_genes_selected = length(genes)
  )
}, BPPARAM= MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])

######################### 6.8 Uni + StepCox + SuperPC  ################################
genes <- row.names(res.cox.coeff)
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

data <- list(x=t(train_data[,-c((ncol(train_data)-1):ncol(train_data))]),y=train_data$OS_MONTHS,censoring.status=train_data$OS_STATUS,featurenames=colnames(train_data)[-c((ncol(train_data)-1):ncol(train_data))])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                     n.fold = 10,
                     n.components=3,
                     min.features=2,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)

RS_list <- bplapply(df_list, function(w) {
  df_data <- w[, c(genes, "OS_MONTHS", "OS_STATUS")]
  df_data <- list(x=t(df_data[,-c((ncol(df_data)-1):ncol(df_data))]),y=df_data$OS_MONTHS,censoring.status=df_data$OS_STATUS,featurenames=colnames(df_data)[-c((ncol(df_data)-1):ncol(df_data))])
  ff <- superpc.predict(fit,data,df_data,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  as.numeric(ff$v.pred)
}, BPPARAM= MulticoreParam(num_detected_cores))

cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]                                                                                
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "Univariate + StepCox + SuperPC",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(genes, collapse = ","),
    No_of_genes_selected = length(genes)
  )
}, BPPARAM= MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])


######################### 6.9 Uni + StepCox + XGB ##################################
#gc()
## install.packages("mlr3")
#library(mlr3)        # For machine learning tasks
## remotes::install_github("mlr-org/mlr3extralearners")
#library(mlr3extralearners)  # For survival learners, such as xgboost Cox
#library(xgboost)
## remotes::install_github("mlr-org/mlr3proba")
#library(mlr3proba)

#genes <- row.names(res.cox.coeff)
#train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

## Define the learner
#learner = mlr3::lrn("surv.xgboost.cox")
#print(learner)

#train_task <- as_task_surv(train_data, time="OS_MONTHS", event="OS_STATUS")

## install.packages("ml3tuning")
#library(mlr3tuning)
#library(paradox)

## Define a search space for hyperparameters
#search_space <- ps(
#  nrounds = p_int(50, 500),
#  eta = p_dbl(0.01, 0.3),
#  max_depth = p_int(2, 10),
#  min_child_weight = p_dbl(1, 10),
#  subsample = p_dbl(0.5, 1),
#  colsample_bytree = p_dbl(0.5, 1),
#  lambda = p_dbl(0, 1),
#  alpha = p_dbl(0, 1)
#)

## Define the tuning instance
#instance <- mlr3tuning::TuningInstanceBatchSingleCrit$new(
#  task = train_task,
#  learner = learner,
#  resampling = rsmp("cv", folds = 5),  # Cross-validation
#  measure = msr("surv.cindex"),       # Optimization metric
#  search_space = search_space,
#  terminator = trm("evals", n_evals = 100)  # Stop after 100 evaluations
#)

## Perform random search or grid search
#tuner <- tnr("random_search")  # Or "grid_search"
#tuner$optimize(instance)

## Use the best hyperparameters
#learner$param_set$values <- instance$result_learner_param_vals



## Create train and test indices
#set.seed(seed)  # Ensure reproducibility

## Train the learner on the training data
#learner$train(train_task, seq_len(nrow(train_data)))

## Print model details
#print(learner$model)

## Get feature importance
#print(learner$importance())

## Assuming `df_list` is a list of datasets you want to apply the model to:
#RS_list <- bplapply(df_list, function(x) {
#  # Prepare the data (ensure to include relevant columns)
#  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
#  predictions = learner$predict_newdata(df_data)
#  # Predict survival using the trained xgboost model
#  as.numeric(predictions$lp)
#}, BPPARAM = MulticoreParam(num_detected_cores))

## Calculate C-index, AUC, and PR-AUC for each dataset
#cc_list <- bplapply(names(df_list), function(x) {
#  RS <- RS_list[[x]]  # Get the predicted risk scores
#  df <- df_list[[x]]
#  # Calculate Concordance Index (C-index)
#  cc <- learner$predict_newdata(df)$score()
#  # Calculate AUC for specific time points
#  auc_values <- sapply(time_points, function(tp) {
#    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
#                           status = df$OS_STATUS,
#                           marker = RS,
#                           predict.time = tp,
#                           method = "KM") # or "NNE", method for handling censoring
#    return(roc_obj$AUC)
#  })
#  names(auc_values) <- paste0("AUC_", time_points, "months")  # Naming AUC values
#  
#  # Calculate Time-dependent PR-AUC at specified time points
#  pr_auc_values <- sapply(time_points, function(tp) {
#    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
#    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
#    return(pr_obj$auc.integral)
#  })
#  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months")  # Naming PR-AUC values
#  
#  # Create a data frame with results
#  data.frame(
#    Model_combination = "Univariate + XGBoost",
#    Dataset_Name = x,
#    C_index = as.numeric(cc),
#    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
#    AUC_36months = as.numeric(auc_values["AUC_36months"]),
#    AUC_60months = as.numeric(auc_values["AUC_60months"]),
#    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
#    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
#    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
#    gene_set = paste(genes, collapse = ","),
#    No_of_genes_selected = length(genes)
#  )
#}, BPPARAM = MulticoreParam(num_detected_cores))

## Combine the results into a single dataframe
#results <- bind_rows(results, bind_rows(cc_list))

## Print results
#print(results[, 1:3])


######################### 7.1 LASSO + StepCox  ################################


gc()
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

# genes <- row.names(cf_lasso_df_non_zero)
# rm(RS_list)
genes <- row.names(cf_lasso_df_non_zero)
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]
# formula_string <- paste("Surv(OS_MONTHS, OS_STATUS) ~ ", paste(genes, collapse = " + "))
# formula <- as.formula(formula_string)
# res.cox <- coxph(formula, data = train_data)
step.model <- step(coxph(Surv(OS_MONTHS,OS_STATUS)~.,train_data), direction = "both", trace = FALSE)
summary(step.model)
res.cox.coeff <- summary(step.model)$coefficients

genes <- row.names(res.cox.coeff)

RS_list <- bplapply(df_list, function(x) {
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
  predict(step.model, type = "risk", newdata = df_data)
}, BPPARAM= MulticoreParam(num_detected_cores))

cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]                                                                                
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "LASSO + StepCox",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(genes, collapse = ","),
    No_of_genes_selected = length(genes)
  )
}, BPPARAM= MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])



######################### 7.2 LASSO + SVM  ################################


# genes <- row.names(cf_lasso_df_non_zero)
# train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]
# 
# set.seed(seed)
# fit = survivalsvm(Surv(OS_MONTHS, OS_STATUS) ~., data= train_data, gamma.mu = 1)
# 
# RS_list <- bplapply(df_list, function(x) {
#   df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
#   as.numeric(predict(fit, df_data)$predicted)
# }, BPPARAM= MulticoreParam(num_detected_cores))
# 
# cc_list <- bplapply(names(df_list), function(x) {
#   RS <- RS_list[[x]]                                                                                
#   df <- df_list[[x]]
#   cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
#   auc_values <- sapply(time_points, function(tp) {
#     roc_obj <- survivalROC(Stime = df$OS_MONTHS,
#                            status = df$OS_STATUS,
#                            marker = RS,
#                            predict.time = tp,
#                            method = "KM") # or "NNE", method for handling censoring
#     return(roc_obj$AUC)
#   })
#   names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
#   # Calculate Time-dependent PR-AUC at specified time points
#   pr_auc_values <- sapply(time_points, function(tp) {
#     binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
#     pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
#     return(pr_obj$auc.integral)
#   })
#   names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
#   
#   data.frame(
#     Model_combination = "LASSO + SVM",
#     Dataset_Name = x,
#     C_index = as.numeric(cc),
#     AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
#     AUC_36months = as.numeric(auc_values["AUC_36months"]),
#     AUC_60months = as.numeric(auc_values["AUC_60months"]),
#     PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
#     PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
#     PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
#     gene_set = paste(genes, collapse = ","),
#     No_of_genes_selected = length(genes)
#   )
# }, BPPARAM= MulticoreParam(num_detected_cores))
# 
# results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe
# 
# print(results[, 1:3])



######################### 7.3 LASSO + GBM  ################################
gc()

genes <- row.names(cf_lasso_df_non_zero)
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]
gbm_model <- gbm(Surv(OS_MONTHS, OS_STATUS) ~ ., data = train_data, distribution = "coxph",
                 n.trees = 10000, interaction.depth = 5, shrinkage = 0.001,
                 cv.folds = 10, n.cores = 2, verbose = FALSE, n.minobsinnode = 10)
best <- which.min(gbm_model$cv.error)
fit <- gbm(formula = Surv(OS_MONTHS,OS_STATUS)~.,data = train_data,distribution = 'coxph',
           n.trees = best,
           interaction.depth = 5,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = num_detected_cores)

RS_list <- bplapply(df_list, function(x) {
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
  predict(fit, df_data, n.trees=best, type="link")
}, BPPARAM= MulticoreParam(num_detected_cores))

cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS,OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values

  data.frame(
    Model_combination = "LASSO + GBM",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(genes, collapse = ","),
    No_of_genes_selected = length(genes)
  )
}, BPPARAM= MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])


######################### 7.4 LASSO + CoxBoost  ################################


set.seed(seed)
library(CoxBoost)
genes <- row.names(cf_lasso_df_non_zero)
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

sfInit(parallel = TRUE, cpus = num_detected_cores)

train_data <- train[, names(train) %in% c(genes, "OS_STATUS", "OS_MONTHS")]
sfLibrary(CoxBoost) # Make CoxBoost library available on workers
sfExport("train_data") # Export train_data to worker environments
sfExport("genes") # Export genes to worker environments

pen <- optimCoxBoostPenalty(train_data[,'OS_MONTHS'],train_data[,'OS_STATUS'], as.matrix(train_data[, !(names(train_data) %in% c("OS_MONTHS","OS_STATUS"))]),
                            trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(train_data[,'OS_MONTHS'],train_data[,'OS_STATUS'], as.matrix(train_data[, !(names(train_data) %in% c("OS_MONTHS","OS_STATUS"))]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
fit <- CoxBoost(train_data[,'OS_MONTHS'],train_data[,'OS_STATUS'], as.matrix(train_data[, !(names(train_data) %in% c("OS_MONTHS","OS_STATUS"))]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)
sfStop()

nonzero_coef <- coef(fit, type = "nonzero")
genes_selected <- names(nonzero_coef[nonzero_coef != 0])

RS_list <- bplapply(df_list, function(x) {
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
  as.numeric(predict(fit, type = "lp", newdata = df_data[, !(names(df_data) %in% c("OS_MONTHS","OS_STATUS"))], newtime= df_data$OS_MONTHS, newstatus=df_data$OS_STATUS))
}, BPPARAM= MulticoreParam(num_detected_cores))

cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]                                                                                
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "LASSO + CoxBoost",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(genes_selected, collapse = ","),
    No_of_genes_selected = length(genes_selected)
  )
}, BPPARAM= MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])



######################### 7.5 LASSO + RSF  ################################

library(randomForestSRC)
genes <- row.names(cf_lasso_df_non_zero)
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

rf_model_final <- rfsrc(Surv(OS_MONTHS, OS_STATUS) ~ ., data = train_data,
                        importance = T, ntree = ntree, # Use full ntree for final model
                        proximity = T, forest = T,
                        seed = seed, nodesize = nodesize,
                        bootstrap = "by.root", ncores= num_detected_cores, nodedepth=nodedepth) # Use topvars for final model
vs.rf_model <- var.select(object = rf_model_final, conservative = "high")
topvars <- vs.rf_model$topvars
length(topvars)

# --- Evaluate Final Model on Test Datasets (test1, test2) ---
RS_list_final <- lapply(df_list, function(x) { # Use df_list to evaluate on train, test1, test2
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")] # Use 'genes' or 'topvars' consistently
  predict(rf_model_final, newdata = df_data)$predicted
})

cc_list_final <- lapply(names(df_list), function(x) {
  RS <- RS_list_final[[x]]
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "LASSO + RSF",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(topvars, collapse = ","),
    No_of_genes_selected = length(topvars)
  )
})

results <- bind_rows(results, bind_rows(cc_list_final)) # Append to results dataframe

print(results[, 1:3])

######################### 7.6 LASSO + plsRcox  ################################

genes <- row.names(cf_lasso_df_non_zero)
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

set.seed(seed)
cv.plsRcox.res=cv.plsRcox(list(x=train_data[,genes],time=train_data$OS_MONTHS,status=train_data$OS_STATUS),nt=5,verbose = FALSE)
fit <- plsRcox(train_data[,genes],time=train_data$OS_MONTHS,event=train_data$OS_STATUS,nt=as.numeric(cv.plsRcox.res[5]))

RS_list <- bplapply(df_list, function(x) {
  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
  as.numeric(predict(fit, newdata= df_data[, !(names(df_data) %in% c("OS_MONTHS","OS_STATUS"))], type="lp"))
}, BPPARAM= MulticoreParam(num_detected_cores))

cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]                                                                                
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS,OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "LASSO + plsRcox",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(genes, collapse = ","),
    No_of_genes_selected = length(genes)
  )
}, BPPARAM= MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])


######################### 7.7  LASSO + SuperPC  ################################

genes <- row.names(cf_lasso_df_non_zero)
train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]
data <- list(x=t(train_data[,-c((ncol(train_data)-1):ncol(train_data))]),y=train_data$OS_MONTHS,censoring.status=train_data$OS_STATUS,featurenames=colnames(train_data)[-c((ncol(train_data)-1):ncol(train_data))])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                     n.fold = 10,
                     n.components=3,
                     min.features=2,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)

RS_list <- bplapply(df_list, function(w) {
  df_data <- w[, c(genes, "OS_MONTHS", "OS_STATUS")]
  df_data <- list(x=t(df_data[,-c((ncol(df_data)-1):ncol(df_data))]),y=df_data$OS_MONTHS,censoring.status=df_data$OS_STATUS,featurenames=colnames(df_data)[-c((ncol(df_data)-1):ncol(df_data))])
  ff <- superpc.predict(fit,data,df_data,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  as.numeric(ff$v.pred)
}, BPPARAM= MulticoreParam(num_detected_cores))

cc_list <- bplapply(names(df_list), function(x) {
  RS <- RS_list[[x]]                                                                                
  df <- df_list[[x]]
  cc <- summary(coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = df))$concordance[1]
  auc_values <- sapply(time_points, function(tp) {
    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
                           status = df$OS_STATUS,
                           marker = RS,
                           predict.time = tp,
                           method = "KM") # or "NNE", method for handling censoring
    return(roc_obj$AUC)
  })
  names(auc_values) <- paste0("AUC_", time_points, "months") # Naming AUC values
  # Calculate Time-dependent PR-AUC at specified time points
  pr_auc_values <- sapply(time_points, function(tp) {
    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
    return(pr_obj$auc.integral)
  })
  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months") # Naming PR-AUC values
  
  data.frame(
    Model_combination = "LASSO + SuperPC",
    Dataset_Name = x,
    C_index = as.numeric(cc),
    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
    AUC_36months = as.numeric(auc_values["AUC_36months"]),
    AUC_60months = as.numeric(auc_values["AUC_60months"]),
    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
    gene_set = paste(genes, collapse = ","),
    No_of_genes_selected = length(genes)
  )
}, BPPARAM= MulticoreParam(num_detected_cores))

results <- bind_rows(results, bind_rows(cc_list)) # Append to results dataframe

print(results[, 1:3])


######################### 7.8 LASSO + XGB ##################################
#gc()
## install.packages("mlr3")
#library(mlr3)        # For machine learning tasks
## remotes::install_github("mlr-org/mlr3extralearners")
#library(mlr3extralearners)  # For survival learners, such as xgboost Cox
#library(xgboost)
## remotes::install_github("mlr-org/mlr3proba")
#library(mlr3proba)

#genes <- row.names(cf_lasso_df_non_zero)
#train_data <-train[, colnames(train) %in% c(genes, "OS_MONTHS", "OS_STATUS")]

## Define the learner
#learner = mlr3::lrn("surv.xgboost.cox")
#print(learner)

#train_task <- as_task_surv(train_data, time="OS_MONTHS", event="OS_STATUS")

## install.packages("ml3tuning")
#library(mlr3tuning)
#library(paradox)

## Define a search space for hyperparameters
#search_space <- ps(
#  nrounds = p_int(50, 500),
#  eta = p_dbl(0.01, 0.3),
#  max_depth = p_int(2, 10),
#  min_child_weight = p_dbl(1, 10),
#  subsample = p_dbl(0.5, 1),
#  colsample_bytree = p_dbl(0.5, 1),
#  lambda = p_dbl(0, 1),
#  alpha = p_dbl(0, 1)
#)

## Define the tuning instance
#instance <- mlr3tuning::TuningInstanceBatchSingleCrit$new(
#  task = train_task,
#  learner = learner,
#  resampling = rsmp("cv", folds = 5),  # Cross-validation
#  measure = msr("surv.cindex"),       # Optimization metric
#  search_space = search_space,
#  terminator = trm("evals", n_evals = 100)  # Stop after 100 evaluations
#)

## Perform random search or grid search
#tuner <- tnr("random_search")  # Or "grid_search"
#tuner$optimize(instance)

## Use the best hyperparameters
#learner$param_set$values <- instance$result_learner_param_vals

## Create train and test indices
#set.seed(seed)  # Ensure reproducibility

## Train the learner on the training data
#learner$train(train_task, seq_len(nrow(train_data)))

## Print model details
#print(learner$model)

## Get feature importance
#print(learner$importance())

## Assuming `df_list` is a list of datasets you want to apply the model to:
#RS_list <- bplapply(df_list, function(x) {
#  # Prepare the data (ensure to include relevant columns)
#  df_data <- x[, c(genes, "OS_MONTHS", "OS_STATUS")]
#  predictions = learner$predict_newdata(df_data)
#  # Predict survival using the trained xgboost model
#  as.numeric(predictions$lp)
#}, BPPARAM = MulticoreParam(num_detected_cores))

## Calculate C-index, AUC, and PR-AUC for each dataset
#cc_list <- bplapply(names(df_list), function(x) {
#  RS <- RS_list[[x]]  # Get the predicted risk scores
#  df <- df_list[[x]]
#  # Calculate Concordance Index (C-index)
#  cc <- learner$predict_newdata(df)$score()
#  # Calculate AUC for specific time points
#  auc_values <- sapply(time_points, function(tp) {
#    roc_obj <- survivalROC(Stime = df$OS_MONTHS,
#                           status = df$OS_STATUS,
#                           marker = RS,
#                           predict.time = tp,
#                           method = "KM") # or "NNE", method for handling censoring
#    return(roc_obj$AUC)
#  })
#  names(auc_values) <- paste0("AUC_", time_points, "months")  # Naming AUC values
#  
#  # Calculate Time-dependent PR-AUC at specified time points
#  pr_auc_values <- sapply(time_points, function(tp) {
#    binary_outcome <- ifelse(df$OS_MONTHS <= tp & df$OS_STATUS == 1, 1, 0) # Event within time point vs. not
#    pr_obj <- pr.curve(scores.class0 = RS, weights.class0 = binary_outcome, curve = FALSE)
#    return(pr_obj$auc.integral)
#  })
#  names(pr_auc_values) <- paste0("PR_AUC_", time_points, "months")  # Naming PR-AUC values
#  
#  # Create a data frame with results
#  data.frame(
#    Model_combination = "Univariate + XGBoost",
#    Dataset_Name = x,
#    C_index = as.numeric(cc),
#    AUC_12months = as.numeric(auc_values["AUC_12months"]), # Extract AUC values
#    AUC_36months = as.numeric(auc_values["AUC_36months"]),
#    AUC_60months = as.numeric(auc_values["AUC_60months"]),
#    PR_AUC_12months = as.numeric(pr_auc_values["PR_AUC_12months"]), # Extract PR-AUC values
#    PR_AUC_36months = as.numeric(pr_auc_values["PR_AUC_36months"]),
#    PR_AUC_60months = as.numeric(pr_auc_values["PR_AUC_60months"]),
#    gene_set = paste(genes, collapse = ","),
#    No_of_genes_selected = length(genes)
#  )
#}, BPPARAM = MulticoreParam(num_detected_cores))

## Combine the results into a single dataframe
#results <- bind_rows(results, bind_rows(cc_list))

## Print results
#print(results[, 1:3])


########################## SAVE RESULTS ################################


print(paste("Start time: ", start_time))
end_time <- Sys.time()
print(paste("End time: ", end_time))

fwrite(results, paste0("complete_gene_selection_results_0.5_parthanatos_alltrain_common_degs.csv"), sep = ",", row.names = T)

library(tidyr)
results2 <- results %>% dplyr::select(Model_combination, Dataset_Name, gene_set, No_of_genes_selected, C_index)

results2 <- results2 %>% pivot_wider(names_from = Dataset_Name, values_from = C_index)
res2 <- results %>%
  group_by(Model_combination) %>%
  summarise(Mean_C_index = mean(unlist(C_index)))
results2 <- merge(results2, res2, by=1)
results2 <- results2 %>% arrange(desc(Mean_C_index))
fwrite(results2, paste0("gene_selection_best_results_0.5_parthanatos_alltrain_common_degs.csv"), sep=',')
# }

############################################### END #############################################
