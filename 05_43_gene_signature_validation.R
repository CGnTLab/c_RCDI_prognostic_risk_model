############ 1. Multivariate Cox regression analysis #############################

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


############ 2. Correlation plots #############################
setwd("/home2/New_Objective2_Biomarker_research/batch_corrected_rna_seq/new/")

library(data.table)
library(corrplot)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(reshape2)

### List of all the expression tables -----------

file_list <- list("Train_RCDI_categorized.csv",
                  "Test_RCDI_categorized.csv",
                  "GSE39582_RCDI_categorized.csv",
                  "GSE161158_RCDI_categorized.csv")
df_list <- lapply(file_list, function(x) {
  data.frame(fread(x, header = TRUE, sep = ","), check.names = FALSE, row.names=1)
})
names(df_list) <- c("Train", "Test", "GSE39582", "GSE161158")

# Assuming you have an expression matrix 'expr_matrix' and two gene sets 'gene_set1' and 'gene_set2'

### List of genes -----------

gs_ferr <- c("CDKN2A","CHMP6","TXNIP","CXCL2","TP63","BID","CDO1","PLIN4","SLC2A3","EGFR","ALOXE3","AURKA","GDF15","GPX2","RGS4")
gs_py <- c("CHMP6","TP63","GSDMC","IL18","GZMB")
gs_aut <- c("MAP1B","CLN3","ATP2B4","SUN2")
gs_net <- c("MYC","CDK4","LAMC2","ETS2","TCF7","NOX4","STAT1","YAP1","C1QC","C1QA","KNG1","SNAI1","IGHG1","CGAS","S100A11","CR1","ACTA2","LAMA2","CDK6","NFATC1","TRAF3IP2")
gs_RCDI <- unique(c(gs_ferr, gs_py, gs_aut, gs_net))

# signature genes
sig_genes <- gs_RCDI

# create directory

output_dir <- "Correlation_Plots"
if (!dir.exists(output_dir)) dir.create(output_dir)


plot_cross_correlation_ggplot <- function(data, x_genes, filename) {
  # Data filtering and preparation
  numeric_data <- data[, sapply(data, is.numeric)]
  x_data <- numeric_data[, intersect(colnames(numeric_data), x_genes), drop = FALSE]
  # y_data <- numeric_data[, intersect(colnames(numeric_data), y_genes), drop = FALSE]
  
  # Ensure consistent order
  x_data <- x_data[, x_genes[x_genes %in% colnames(x_data)], drop = FALSE]
  # y_data <- y_data[, y_genes[y_genes %in% colnames(y_data)], drop = FALSE]
  
  # Compute correlations and p-values
  results <- expand.grid(Y = colnames(x_data), X = colnames(x_data)) %>%
    rowwise() %>%
    mutate(
      cor = cor(x_data[[Y]], x_data[[X]], use = "pairwise.complete.obs"),
      pval = cor.test(x_data[[Y]], x_data[[X]], use = "pairwise.complete.obs")$p.value,
      sig = case_when(
        pval <= 0.001 ~ "***",
        pval <= 0.01  ~ "**",
        pval <= 0.05  ~ "*",
        TRUE ~ ""
      )
    ) %>%
    ungroup()
  
  # Start TIFF device
  tiff(filename, width = 8, height = 8, units = "in", res = 600, compression = "lzw")
  
  # Plot
  p <- ggplot(results, aes(x = X, y = Y, fill = cor)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sig), color = "black", size = 2.5) +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0, limit = c(-1, 1), name = "PCC"
    ) +
    labs(x = "c-RCDI signature genes", y = "c-RCDI signature genes") +
    theme_bw(base_size = 12, base_family = "Helvetica") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid = element_blank(),
      plot.margin = margin(10, 10, 10, 10),
      axis.title = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 10, colour = "black"),
      legend.title = element_text(size = 12, face = "bold", colour = "black"),
      legend.text = element_text(size = 10, colour = "black")
    )
  
  print(p)  # Important: render the plot in the tiff device
  dev.off()
}

for (name in names(df_list)) {
  plot_cross_correlation_ggplot(df_list[[name]], sig_genes, paste0(output_dir, "/", name, "_corrplot.tiff"))
}

############ 3. CRISPR/Cas9- knockout (KO) #############################

# gene dependency data extraction from depmap portal through R
crispr_data <- depmap::depmap_crispr()
crc_crispr <- crispr_data %>%
  filter(depmap_id %in% crc_lines &
           gene_name %in% genes)

gg <- ggplot(crc_crispr, aes(x = gene_name, y = dependency)) +
  geom_boxplot(color = "red", fill = NA, outlier.shape = NA, width = 0.6) +  # Boxplot with salmon border and no fill
  geom_jitter(width = 0.1, size = 0.1, alpha = 1, color = "black") +  # Jitter points in salmon
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5, alpha = 0.8) +
  labs(title = NULL,
       x = "Gene",
       y = "CERES Dependency Score",
       caption = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(color = "black", angle= 0),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, face = "bold", color = "black"),
        legend.position = "none",
        title = element_text(face = "bold", colour = "black", size = 14, vjust = 0.5, hjust = 0.5)) +
  coord_flip()

# Save as TIFF
tiff("CRISPR_gene_dependency_43_genes.tif", res = 600, width = 5, height = 9, units = "in", compression = "lzw")
print(gg)
dev.off()

############ 4. Single cell expression analysis #############################
