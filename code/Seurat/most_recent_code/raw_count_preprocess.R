library(dplyr)
library(Seurat)
library(tidyverse)

Zheng_data <- read.table("original_data/raw_counts/GSE98638_HCC.TCell.S5063.count.txt",
                         header = TRUE, sep = "", dec = ".")

Guo_data <- read.table("original_data/raw_counts/GSE99254_NSCLC.TCell.S12346.count.txt",
                       header = TRUE, sep = "", dec = ".")

###### Zheng data preprocessing  ########

# remove duplicate gene symbols 
Zheng_1 <- Zheng_data
Zheng_1 <- Zheng_1[!duplicated(Zheng_1$symbol), ]
dim(Zheng_1)

# remove geneID column
Zheng_1[1] <- list(NULL)
Zheng_1[1:5,1:5]

# remove NA rownames from Zheng_1
which(is.na(Zheng_1[[1]])) #62
Zheng_1 <- Zheng_1[-c(62),]
anyNA(Zheng_1[[1]])

# set row names as gene symbols
rownames(Zheng_1) <- Zheng_1[[1]]
# col names already are the patient cell names

Zheng_seurat <- CreateSeuratObject(counts = Zheng_1, project = "Zheng", min.cells = 3, min.features = 200)
Zheng_seurat

# QC
VlnPlot(Zheng_seurat, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(Zheng_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
Zheng_seurat <- subset(Zheng_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & nCount_RNA < 2800000)

VlnPlot(Zheng_seurat, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(Zheng_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Normalization
Zheng_seurat <- NormalizeData(Zheng_seurat, normalization.method = "LogNormalize", scale.factor = 10000)

# Scaling
all.genes <- rownames(Zheng_seurat)
Zheng_seurat <- ScaleData(Zheng_seurat, features = all.genes)

###### Guo data preprocessing  ########

# remove duplicate gene symbols 
Guo_1 <- Guo_data
Guo_1 <- Guo_1[!duplicated(Guo_1$symbol), ]
dim(Guo_1)

# remove geneID column
Guo_1[1] <- list(NULL)
Guo_1[1:5,1:5]

# remove NA rownames from Guo_1
Guo_1 <- na.omit(Guo_1)
anyNA(Guo_1[[1]])

# set row names as gene symbols
rownames(Guo_1) <- Guo_1[[1]]
Guo_1[1:5,1:5]
# col names already are the patient cell names

Guo_seurat <- CreateSeuratObject(counts = Guo_1, project = "Guo", min.cells = 3, min.features = 200)
Guo_seurat

# QC
VlnPlot(Guo_seurat, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(Guo_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
Guo_seurat <- subset(Guo_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & nCount_RNA < 2000000)

VlnPlot(Guo_seurat, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(Guo_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Normalization
Guo_seurat <- NormalizeData(Guo_seurat, normalization.method = "LogNormalize", scale.factor = 10000)

# Scaling
all.genes.Guo <- rownames(Guo_seurat)
Guo_seurat <- ScaleData(Guo_seurat, features = all.genes.Guo)
Guo_seurat@assays$RNA@scale.data[1:5,1:5]


dim(Zheng_seurat)
dim(Guo_seurat)

# Get overlapping genes (for merging seurat objects later)
# prune matrices: only keep rows with the shared genes

overlap_gene <- intersect(rownames(Zheng_seurat@assays$RNA@scale.data),
                          rownames(Guo_seurat@assays$RNA@scale.data))

length(overlap_gene)

subset_genes_Zheng_seurat <- Zheng_seurat@assays$RNA@scale.data[overlap_gene, ]
pruned_Zheng_seurat <- CreateSeuratObject(subset_genes_Zheng_seurat) # Create a new Seurat object with just the genes of interest

subset_genes_Guo_seurat <- Guo_seurat@assays$RNA@scale.data[overlap_gene, ]
pruned_Guo_seurat <- CreateSeuratObject(subset_genes_Guo_seurat) # Create a new Seurat object with just the genes of interest

# check they have the same number of features (genes)
pruned_Zheng_seurat
pruned_Guo_seurat



