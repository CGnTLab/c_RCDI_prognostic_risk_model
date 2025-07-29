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
setwd("/data/avik/single_cell_crc/single_cell_analysis")
library(Seurat)
library(ggplot2)
library(SingleR)
library(dplyr)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)
library(data.table)
seed = 42
set.seed(seed)

sc.data <- data.frame(fread("GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt", sep="\t"), check.names=F, row.names=1, nThread=18)

metadata <- data.frame(fread("GSE132465_GEO_processed_CRC_10X_cell_annotation.txt", sep="\t"), check.names=F)
meta1 <- metadata[metadata$Class == "Tumor",]
row.names(meta1)<- meta1[,1]
meta1 <- meta1[,-1]
columns <- rownames(meta1)
sc.data <- sc.data[, colnames(sc.data) %in% columns]
sc.data <- sc.data[rowSums(sc.data) > 10, ]
srat <- CreateSeuratObject(counts = sc.data, project = "smc", meta.data=meta1) #adding metadata as meta.data
srat
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(srat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

##subset features as a part of QC
srat <- subset(srat, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & percent.mt < 15 & nCount_RNA < 40000)

#RNA study
srat[["RNA"]] <- split(srat[["RNA"]], f = srat$orig.ident)

srat <- SCTransform(srat, vars.to.regress = "percent.mt", verbose = FALSE)

##dimensionality reductions
#linear
srat <- RunPCA(srat, features = VariableFeatures(object = srat))

# Examine and visualize PCA results a few different ways
print(srat[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(srat, dims = 1:2, reduction = "pca")
DimPlot(srat, reduction = "pca") + NoLegend()

DimHeatmap(srat, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(srat, dims = 1:30, cells = 500, balanced = TRUE)

#determine the dimentionality
ElbowPlot(srat, ndims = 50)


##check for batch correction and optimize
srat <- RunUMAP(srat, dims = 1:50)

DimPlot(srat, reduction="umap", group.by="orig.ident", label=T)

srat <- FindNeighbors(srat, dims = 1:50)
srat <- FindClusters(srat, resolution = 1.2)
# Look at cluster IDs of the first 5 cells
head(Idents(srat), 5)

#run non-linear dimensionality reduction

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(srat, reduction = "umap")
saveRDS(srat, file = "srat_crc_19may25.rds")

srat <- PrepSCTFindMarkers(srat)
# install.packages("devtools")
devtools::install_github('immunogenomics/presto')
#Finding differentially expressed features (cluster biomarkers)
# find all markers of cluster 2
cluster2.markers <- FindMarkers(srat, ident.1 = 2)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(srat, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
srat.markers <- FindAllMarkers(srat, only.pos = F)

#Find the top 5 differential genes from each cluster
srat.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)

#Find the conserved genes from each cluster
DefaultAssay(srat) <- "RNA"

Idents(srat) <- "Cell_type"
library(ggplot2)

# Our gene list of 43 signature genes
gene_list <- c("CDKN2A",
"CHMP6",
"TXNIP",
"CXCL2",
"TP63",
"BID",
"CDO1",
"PLIN4",
"SLC2A3",
"EGFR",
"ALOXE3",
"AURKA",
"GDF15",
"GPX2",
"RGS4",
"GSDMC",
"IL18",
"GZMB",
"MAP1B",
"CLN3",
"ATP2B4",
"SUN2",
"MYC",
"CDK4",
"AMC2",
"ETS2",
"TCF7",
"NOX4",
"STAT1",
"YAP1",
"C1QC",
"C1QA",
"KNG1",
"SNAI1",
"IGHG1",
"CGAS",
"S100A11",
"CR1",
"ACTA2",
"LAMA2",
"CDK6",
"NFATC1",
"TRAF3IP2")

genes_present <- gene_list[gene_list %in% rownames(srat)]

for (gene in genes_present) {
  p <- FeaturePlot(srat, features = gene, label = TRUE, repel = TRUE) +
       ggtitle(gene) +
       theme(plot.title = element_text(hjust = 0.5, size = 14))
  
  tiff(paste0("featureplot_", gene,".tif"), res=600, units="in", width=7, height=6, compression="lzw")
  print(p)
  dev.off()
}

############ 5. WGCNA #############################
if (!("DESeq2" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("DESeq2", update = FALSE)
}

if (!("impute" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("impute")
}

if (!("ggforce" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  install.packages("ggforce")
}

if (!("ComplexHeatmap" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  install.packages("ComplexHeatmap")
}

if (!("WGCNA" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  install.packages("WGCNA")
}

if (!("flashClust" %in% installed.packages())) {
BiocManager::install('flashClust')
}

## load the libraries
library(dplyr)
library(WGCNA)
library(flashClust)
library(curl)
library(ggplot2)
library(readxl)
library(data.table)

enableWGCNAThreads()
allowWGCNAThreads()

data <- data.frame(fread("complete_z_scores_Train.csv", header = TRUE, sep = ","), row.names = 1, check.names= FALSE)

##### identify & remove outlier genes -------
gsg <-goodSamplesGenes(data)
summary(gsg)
gsg$allOK

if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(data)[!gsg$goodGenes], collapse = ", "))); #Identifies and prints outlier genes
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(data)[!gsg$goodSamples], collapse = ", "))); #Identifies and prints oulier samples
  data <- data[gsg$goodSamples == TRUE, gsg$goodGenes == TRUE] # Removes the offending genes and samples from the data
}

##### identify & remove outlier samples -------
sampleTree <- hclust(dist(data), method = "average") #Clustering samples based on distance 

#Setting the graphical parameters
par(cex = 0.6);
par(mar = c(0,4,2,0))

#Plotting the cluster dendogram
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# remove the outlier(s) using cutree function
#Setting the graphical parameters
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
#draw on line to show cutoff height
abline(h = 350, col = "red");

# retain the rest of the samples
cut.sampleTree <- cutreeStatic(sampleTree, cutHeight = 350, minSize = 10) #returns numeric vector
#Remove outlier
data <- data[cut.sampleTree==1, ]

#### Network Construction ---------

# Note: We should be maximizing the R^2 value and minimizing mean connectivity.

# select beta values
spt <- pickSoftThreshold(data)

# plot the R^2 vs soft thresholds
tiff("WGCNA_results/r2_vs_soft_thres.tif", res=600, units="in", compression="lzw", height=5, width=5)
par(mar=c(1,1,1,1))
plot(spt$fitIndices[,1],spt$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(spt$fitIndices[,1],spt$fitIndices[,2],col="red")
abline(h=0.80,col="red")
dev.off()

# plot mean connectivity as a function of soft thresholds
tiff("WGCNA_results/mean_connectivity.tif", res=600, units="in", compression="lzw", height=5, width=5)
par(mar=c(1,1,1,1))
plot(spt$fitIndices[,1], spt$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(spt$fitIndices[,1], spt$fitIndices[,5], labels= spt$fitIndices[,1],col="red")
dev.off()

# Calling the adjacency function

softPower <- 6 # from above plots
adjacency <- adjacency(data, power = softPower)

##### Module Construction ---------

# Topological Overlap Matrix for dissimilarity determination
gc()
TOM <- TOMsimilarity(adjacency) # To convert the adjacency matrix into a TOM similarity matrix
TOM.dissimilarity <- 1-TOM # To convert this matrix into a dissimilarity matrix you can subtract the TOM object from 1.

# Hierarchical clustering
geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average")  # clustering on the basis of dissimilarity

#plotting the dendrogram
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", 
     labels = FALSE, hang = 0.04)

## NOTE: To identify modules from this gene dendrogram, you can use the cutreeDynamic() function. This will allow you to set a 
## minimum cluster size. For genomic data like this it is more beneficial to set minimum module sizes relatively high as 
## you are working with high loads of data. The authors of WGCNA recommend to start at a minClusterSize = 30.

Modules <- cutreeDynamic(dendro = geneTree, distM = TOM.dissimilarity, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)

table(Modules) #returns a table of the counts of factor levels in an object. In this case how many genes are assigned to each created module. 

# Plot the module assignments under gene dendogram for visualization

ModuleColors <- labels2colors(Modules) #assigns each module number a color
table(ModuleColors) #returns the counts for each color (aka the number of genes within each module)

#plots the gene dendrogram with the module colors
plotDendroAndColors(geneTree, ModuleColors,"Module",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

## Module Eigengene identification

MElist <- moduleEigengenes(data, colors = ModuleColors) 
MEs <- MElist$eigengenes 
head(MEs)

#### Module Merging --------

ME.dissimilarity = 1-cor(MElist$eigengenes, use="complete") #Calculate eigengene dissimilarity
METree = hclust(as.dist(ME.dissimilarity), method = "average") #Clustering eigengenes 
tiff("WGCNA_results/cluster_dendogram_before_merge.tif", res=600, compression="lzw", height=5, width=5, units="in")
par(mar = c(0,4,2,0)) #seting margin sizes
par(cex = 0.6);#scaling the graphic
plot(METree)
abline(h=.25, col = "red") #a height of .25 corresponds to correlation of .75
dev.off()
merge <- mergeCloseModules(data, ModuleColors, cutHeight = .25)

# The merged module colors, assigning one color to each module
mergedColors = merge$colors
# Eigengenes of the new merged modules
mergedMEs = merge$newMEs


## We take your 43 gene list and check where do they lie in modules
geneList <- read.table("./43_gene_list.txt", header = F, sep = ",")
geneList <- as.character(geneList$V1)
gene_indices = match(geneList, colnames(data))
module_assignments = mergedColors[gene_indices]
# Create results table
geneList_results = data.frame(
  Gene = geneList,
  MergedModule = module_assignments
)
write.table(geneList_results, "WGCNA_results/module_genes_43.csv", sep=',', row.names=T, col.names=T)
table(module_assignments)

## Plotting both original and merged module dendograms
tiff("WGCNA_results/cluster_dendogram_before_and after_merge.tif", res=600, compression="lzw", height=6, width=9, units="in")
plotDendroAndColors(geneTree, cbind(ModuleColors, mergedColors), 
                    c("Original Module", "Merged Module"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors for original and merged modules")
dev.off()

###### External Trait Matching -------

# read clinical_data

df_clinical <- read.csv("Train_clinical_expr_combined_rcdi.csv", sep = ",", row.names = 1, check.names = FALSE)
df_clinical <- df_clinical %>%
  select(RCDI, RiskClass, Sex, Age_group, Tumor_stage) %>%
  dplyr::mutate(RiskClass = recode(RiskClass,
                            "low" = 1,
                            "high" = 2),
         Sex = recode(Sex,
                 "Male" = 2,
                 "Female" = 1),
    Age_group = recode(Age_group,
                       "<= 65" = 1,
                       "66-79" = 2,
                       ">= 80" = 3),
    Tumor_stage = recode(Tumor_stage,
                         "Stage I - II" = 1,
                         "Stage III - IV" = 2))

# MAtch the trait data to expression data by sample name
expr_data <- data[row.names(data) %in% row.names(df_clinical),]
expr_data <- expr_data[row.names(df_clinical),]

# Define numbers of genes and samples
nGenes = ncol(expr_data)
nSamples = nrow(expr_data)
mod_merged_MEs = mergedMEs[row.names(mergedMEs) %in% row.names(df_clinical),] # subset for the samples in the clinical data
mod_merged_MEs <- mod_merged_MEs[row.names(df_clinical),] # same order of samples as df_clinical
module.trait.correlation = cor(mod_merged_MEs, df_clinical, use = "p") #p for pearson correlation coefficient 
module.trait.Pvalue = corPvalueStudent(module.trait.correlation, nSamples) #calculate the p-value associated with the correlation

# Will display correlations and their p-values
textMatrix = paste(signif(module.trait.correlation, 2), "\n(",
                   signif(module.trait.Pvalue, 1), ")", sep = "");
dim(textMatrix) = dim(module.trait.correlation)
par(mar = c(6, 8.5, 3, 1))
# Display the correlation values within a heatmap plot
tiff("WGCNA_results/heatmap_RCDI_correlation_new.tif", res=600, compression="lzw", height=7, width=5, units="in")

labeledHeatmap(Matrix = module.trait.correlation,
               xLabels = names(df_clinical),
               yLabels = names(mod_merged_MEs),
               ySymbols = names(mod_merged_MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = T,
               cex.text = 0.4,
               cex.lab.y = 0.6,
               cex.lab.x = 0.6,
               zlim = c(-1,1),
               main = NULL)
dev.off()

####### Target Gene Identification -------
# Define variable weight containing the RCDI column of df_clinical
rcdi = as.data.frame(df_clinical$RCDI)
names(rcdi) = "RCDI"

modNames = substring(names(mod_merged_MEs), 3) #extract module names

#Calculate the module membership and the associated p-values
geneModuleMembership = as.data.frame(cor(expr_data, mod_merged_MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

#Calculate the gene significance and associated p-values
geneTraitSignificance = as.data.frame(cor(expr_data, rcdi, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(rcdi), sep="")
names(GSPvalue) = paste("p.GS.", names(rcdi), sep="")
head(GSPvalue)

hubGenes = names(which((abs(GSPvalue) > 0.2) & (abs(MMPvalue) > 0.8)))


## Scatter plot of gene significance vs module membership
par(mar=c(1,1,1,1))
module = "purple" #choose a module
column = match(module, modNames)
moduleGenes = mergedColors==module
verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),
                   abs(geneTraitSignificance[moduleGenes,1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

par(mar=c(1,1,1,1))
module = "brown" #choose a module
column = match(module, modNames)
moduleGenes = mergedColors==module
verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),
                   abs(geneTraitSignificance[moduleGenes,1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

# Network Visualization of Eigengenes
# Isolate rcdi from the clinical traits
rcdi = as.data.frame(df_clinical$RCDI)
names(rcdi) = "RCDI"
# Add the weight to existing module eigengenes
mod_MEs <- MEs[row.names(MEs) %in% row.names(df_clinical),]
mod_MEs <- mod_MEs[row.names(df_clinical),]
MET = orderMEs(cbind(mod_MEs, rcdi))
# Plot the relationships among the eigengenes and the trait
tiff("WGCNA_results/corr_RCDI_modules.tif", res=600, compression="lzw", height=6, width=9, units="in")
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(5,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)
dev.off()
# Plot the dendrogram
tiff("WGCNA_results/dendogram.tif", res=600, compression="lzw", height=6, width=9, units="in")
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
dev.off()

# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
tiff("WGCNA_results/corr_heatmap.tif", res=600, compression="lzw", height=6, width=9, units="in")
par(cex = 1.0)
par(cex = 1.0, mar = c(1,1,1,1))
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(5,5,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()
