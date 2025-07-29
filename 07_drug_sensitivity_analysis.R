#!/usr/bin/Rscript

setwd("/home2/New_Objective2_Biomarker_research/Drug_Response")
rm(list=ls())

# BiocManager::install(TCGAbiolinks, force=T)
library(TCGAbiolinks)
# install.packages("oncoPredict")
library(oncoPredict)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(readr)

seed = 42
set.seed(seed)

trainingPtype = readRDS(file = "GDSC2_Res.rds")
trainingPtype <- exp(trainingPtype) #since we will transform it again in pcalctype function
trainingExprData=readRDS(file='GDSC2_Expr_short.rds')
trainingExprData <- t(scale(t(trainingExprData), center = T, scale = T)) #z-score transformation

testExprData=data.frame(fread("../batch_corrected_rna_seq/new/complete_z_scores_Train.csv", header=T, sep=","), row.names=1, check.names=F) %>%
  t() %>% as.matrix()

#Additional parameters. 
#_______________________________________________________
#batchCorrect options: "eb" for ComBat, "qn" for quantiles normalization, "standardize", or "none"
#"eb" is good to use when you use microarray training data to build models on microarray testing data.
#"standardize is good to use when you use microarray training data to build models on RNA-seq testing data (this is what Paul used in the 2017 IDWAS paper that used GDSC microarray to impute in TCGA RNA-Seq data, see methods section of that paper for rationale)
batchCorrect<-"standardize"

#Determine whether or not to power transform the phenotype data.
#Default is TRUE.
powerTransformPhenotype<-T

#Determine percentage of low varying genes to remove.
#Default is 0.2 (seemingly arbitrary).
removeLowVaryingGenes<-0.2

#Determine method to remove low varying genes.
#Options are 'homogenizeData' and 'rawData'
#homogenizeData is likely better if there is ComBat batch correction, raw data was used in the 2017 IDWAS paper that used GDSC microarray to impute in TCGA RNA-Seq data.
removeLowVaringGenesFrom<-"rawData"

#Determine the minimum number of training samples required to train on.
#Note: this shouldn't be an issue if you train using GDSC or CTRP because there are many samples in both training datasets.
#10, I believe, is arbitrary and testing could be done to get a better number.
minNumSamples=10

#Determine how you would like to deal with duplicate gene IDs.
#Sometimes based on how you clean the data, there shouldn't be any duplicates to deal with.
#Options are -1 for ask user, 1 for summarize by mean, and 2 for disregard duplicates
selection<- 1

#Determine if you'd like to print outputs.
#Default is TRUE.
printOutput=TRUE

#Indicate whether or not you'd like to use PCA for feature/gene reduction. Options are 'TRUE' and 'FALSE'.
#Note: If you indicate 'report_pca=TRUE' you need to also indicate 'pca=TRUE'
pcr=FALSE

#Indicate whether you want to output the principal components. Options are 'TRUE' and 'FALSE'.
report_pc=FALSE

#Indicate if you want correlation coefficients for biomarker discovery. These are the correlations between a given gene of interest across all samples vs. a given drug response across samples.
#These correlations can be ranked to obtain a ranked correlation to determine highly correlated drug-gene associations.
cc=FALSE

#Indicate whether or not you want to output the R^2 values for the data you train on from true and predicted values.
#These values represent the percentage in which the optimal model accounts for the variance in the training data.
#Options are 'TRUE' and 'FALSE'.
rsq=FALSE

#Indicate percent variability (of the training data) you'd like principal components to reflect if pcr=TRUE. Default is .80
percent=80

#Run the calcPhenotype() function using the parameters you specified above.
#__________________________________________________________________________________________________________________________________
# wd<-tempdir()
# savedir<-setwd(wd)

calcPhenotype(trainingExprData=trainingExprData,
              trainingPtype=trainingPtype,
              testExprData=testExprData,
              batchCorrect=batchCorrect,
              powerTransformPhenotype=powerTransformPhenotype,
              removeLowVaryingGenes=removeLowVaryingGenes,
              minNumSamples=minNumSamples,
              selection=selection,
              printOutput=printOutput,
              pcr=pcr,
              removeLowVaringGenesFrom=removeLowVaringGenesFrom,
              report_pc=report_pc,
              cc=cc,
              percent=percent,
              rsq=rsq)

res_pred <- data.frame(fread("./calcPhenotype_Output/DrugPredictions.csv", sep=',', header=T), row.names=1, check.names=F)
# res_pred <- (1e-3)*res_pred # Convert to mM scale
res_pred <- log(res_pred) #natural log for drug sensitivity analysis

# res_pred <- na.omit(res_pred)
a <- read.table("../batch_corrected_rna_seq/new/Train_RCDI_categorized.csv", sep=',', row.names=1, check.names=F, header=T)
a <- a %>% dplyr::select(RiskClass)

# we show examples of plotting using p value criteria and significance. One such example is shown below, logIC50 thershold is taken to be 1.5, p-val< 0.05.
# Single pdf for the plot
pdf("drug_plots_gdsc2/drugs_high_vs_low_gdsc2_psig_train_1.5_low>high.pdf")
remove_outliers_iqr <- function(df, column = "Drug") {
  df_clean <- df %>% 
    group_by(RiskClass) %>% 
    mutate(Q1 = quantile(.data[[column]], 0.25, na.rm = TRUE),
           Q3 = quantile(.data[[column]], 0.75, na.rm = TRUE),
           IQR = Q3 - Q1,
           lower = Q1 - 1.5 * IQR,
           upper = Q3 + 1.5 * IQR) %>%
    filter(.data[[column]] >= lower & .data[[column]] <= upper) %>%
    ungroup() %>%
    dplyr::select(-Q1, -Q3, -IQR, -lower, -upper)
  return(df_clean)
}
drugs <- names(res_pred)
for (i in 1:length(drugs)) {
  drug <- drugs[i]
  # print(drug)
  res_pred_sub <- res_pred[, i, drop=F]
  merged <- merge(res_pred_sub, a, by=0)
  row.names(merged) <- merged[,1]
  merged<- merged[,-1]
  names(merged)[1] <- "Drug"
  merged <- remove_outliers_iqr(merged)
  p <- ggboxplot(merged, x = "RiskClass", y = "Drug",
                 color = "RiskClass", palette =c("low"="skyblue","high"="salmon"),
                 add = c("jitter"), add.params=list(size=1, alpha=0.2), 
                 size=0.9, bxp.errorbar = FALSE, width=0.5)+
    geom_hline(yintercept = 1.5, linetype = "dotted", color = "black", size = 1) + 
    #coord_cartesian(ylim = c(0, 5))+
    # stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9) +
    labs(x="RiskClass", y= expression(bold(logIC[50])), title = paste0(drug))+
    theme(
      axis.title.x = element_text(face="bold", colour="Black", size=30),
      axis.title.y = element_text(face="bold", colour="Black", size=30),
      axis.text.x=element_text(size=20, face="bold"),
      panel.border = element_rect(color = "black",fill = NA,size = 1),
      panel.background = element_blank(),
      axis.text.y=element_text(angle=90, size=20),
      legend.position= "none",
      plot.title = element_text(size=30, face="bold", hjust = 0.5),
      legend.title = element_text(size = 20, face = "bold"),
      legend.text = element_text(size = 18)
    )
  #p <- p + coord_cartesian(ylim = y_range_1)	
  p<- p+ stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),
                            method = "t.test", ref.group = "low", hide.ns=TRUE, size=8)
  
  #extract median value
  merged_high <- merged[merged$RiskClass == "high",]
  median_value_high <- median(merged_high$Drug, na.rm =T)
  merged_low <- merged[merged$RiskClass == "low",]
  median_value_low <- median(merged_low$Drug, na.rm = T)
  # Extract p-value
  p_value <- compare_means(Drug ~ RiskClass, data = merged, method = "t.test")$p
  # Print plot only if p-value is significant
  if ((p_value < 0.05) & (median_value_low < 1.5) & (median_value_high < median_value_low)) {
    print(p)
  }
}
dev.off()
