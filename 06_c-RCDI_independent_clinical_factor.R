#!/usr/bin/Rscript

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
