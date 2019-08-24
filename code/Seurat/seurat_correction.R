# This file takes in preprocessed data after the R script "preprocess.R"
# gets run, and performs Seurat correction and visualization on the combined 
# datasets Zheng et al and Guo et al.

library(Seurat)
library(cowplot)


######## PART 1 : Seurat correction and visualization #########

# Set up Zheng object
Zheng_seurat <- CreateSeuratObject(counts = Zheng_4)
Zheng_seurat@meta.data$batch <- "Batch1"

# add in metadata: cell id, patient, cluster and sample type 
#Zheng_seurat@meta.data$cell.id <- 

# Set up Guo object
Guo_seurat <- CreateSeuratObject(counts = Guo_4)
Guo_seurat@meta.data$batch <- "Batch2"

# merge two seurat objects
merged_uncorrected <- merge(x=Zheng_seurat, y=Guo_seurat)

# Plot clusters before correction
merged_uncorrected <- FindVariableFeatures(merged_uncorrected, selection.method = "vst", nfeatures = 2000)

# Run the standard workflow for visualization and clustering
merged_uncorrected <- ScaleData(merged_uncorrected, verbose = FALSE)

merged_uncorrected <- RunPCA(merged_uncorrected, npcs = 30, verbose = FALSE)

# t-SNE and Clustering
merged_uncorrected <- RunUMAP(merged_uncorrected, reduction = "pca", dims = 1:20)
merged_uncorrected <- FindNeighbors(merged_uncorrected, reduction = "pca", dims = 1:20)
merged_uncorrected <- FindClusters(merged_uncorrected, resolution = 0.5)

# Visualization
# group.by: color scheme
p1 <- DimPlot(merged_uncorrected, reduction = "umap", group.by = "batch")
p2 <- DimPlot(merged_uncorrected, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

# split.by: name of a metadata column to split plot by (side by side view)
DimPlot(merged_uncorrected, reduction = "umap", split.by = "batch", label = TRUE)


######## Batch Correction: Integrate two datasets ########## 

immune.anchors <- FindIntegrationAnchors(object.list = list(Zheng_seurat, Guo_seurat), 
                                         dims = 1:20,
                                         k.filter = 100) 

Zheng.Guo.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)

# Perform integrated analysis
DefaultAssay(Zheng.Guo.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
Zheng.Guo.combined <- ScaleData(Zheng.Guo.combined, verbose = FALSE)
Zheng.Guo.combined <- RunPCA(Zheng.Guo.combined, npcs = 30, verbose = FALSE)

# t-SNE and Clustering
Zheng.Guo.combined <- RunUMAP(Zheng.Guo.combined, reduction = "pca", dims = 1:20)
Zheng.Guo.combined <- FindNeighbors(Zheng.Guo.combined, reduction = "pca", dims = 1:20)
Zheng.Guo.combined <- FindClusters(Zheng.Guo.combined, resolution = 0.5)

# Visualization
p1 <- DimPlot(Zheng.Guo.combined, reduction = "umap", group.by = "batch")
p2 <- DimPlot(Zheng.Guo.combined, reduction = "umap", label = TRUE, label.size = 5)
plot_grid(p1, p2)

DimPlot(Zheng.Guo.combined, reduction = "umap", split.by = "batch")


###### PART 2: Identify differentially expressed features (cluster biomarkers) ######

library(dplyr)

DefaultAssay(Zheng.Guo.combined) <- "RNA"

# cluster 6: these cells are more abundant in batch 1: HCC
cluster6.markers <- FindMarkers(Zheng.Guo.combined, ident.1 = 6, only.pos = TRUE)
head(cluster6.markers, n = 10)
#cluster6.markers %>% top_n(n=6, wt=avg_logFC)

write.table(cluster6.markers, "5.24.cluster6.markers.txt", sep="\t")

VlnPlot(Zheng.Guo.combined, features = c("SLC4A10", "KLRB1", "KLRG1"), pt.size=0)
VlnPlot(Zheng.Guo.combined, features = c("BCAS1", "PLIN4"), pt.size=0.05)

FeaturePlot(Zheng.Guo.combined, features = c("SLC4A10", "KLRB1"))
FeaturePlot(Zheng.Guo.combined, features = c("KLRG1"))
FeaturePlot(Zheng.Guo.combined, features = c("AGAP3"))

############### Find all markers #################

# find markers for every cluster compared to all remaining cells, report only the positive ones
all.markers <- FindAllMarkers(Zheng.Guo.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

## wt: sort by avg_logFC
all.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC) # select greated avg logFC
#all.markers %>% group_by(cluster) %>% top_n(n = -1, wt = p_val) ## select smallest p value
#max(all.markers$avg_logFC)

write.table(all.markers, "5.28.all.markers.txt" , sep="\t")

## Visualize markers
# c0
VlnPlot(Zheng.Guo.combined, features = c("CRMP1", "IL21"), pt.size=0.5)
FeaturePlot(Zheng.Guo.combined, features = c("CRMP1", "IL21"))

# c1
VlnPlot(Zheng.Guo.combined, features = c("SUSD4", "EDAR"), pt.size=0.5)
FeaturePlot(Zheng.Guo.combined, features = c("SUSD4", "EDAR"))

# c2
VlnPlot(Zheng.Guo.combined, features = c("KLRC1", "VCAM1"), pt.size=0.01)
FeaturePlot(Zheng.Guo.combined, features = c("KLRC1", "VCAM1"))

# c3
VlnPlot(Zheng.Guo.combined, features = c("FGFBP2", "HBB"), pt.size=0)
FeaturePlot(Zheng.Guo.combined, features = c("FGFBP2", "HBB"))

# c4
VlnPlot(Zheng.Guo.combined, features = c("CD177", "CCL22"), pt.size=0)
FeaturePlot(Zheng.Guo.combined, features = c("CD177", "CCL22"))

# c5
VlnPlot(Zheng.Guo.combined, features = c("PI16", "CCR10"), pt.size=0.01)
FeaturePlot(Zheng.Guo.combined, features = c("PI16", "CCR10"))
# c6
VlnPlot(Zheng.Guo.combined, features = c("SLC4A10", "ME1"), pt.size=0)
FeaturePlot(Zheng.Guo.combined, features = c("SLC4A10", "ME1"))
# c7
VlnPlot(Zheng.Guo.combined, features = c("RRM2", "DLGAP5"), pt.size=0)
FeaturePlot(Zheng.Guo.combined, features = c("RRM2", "DLGAP5"))
# c8
VlnPlot(Zheng.Guo.combined, features = c("AQP11", "PAQR6"), pt.size=0.3)
FeaturePlot(Zheng.Guo.combined, features = c("AQP11", "PAQR6"))


# min.cutoff: quantile=9, partition all values to 9 clusters as found
FeaturePlot(Zheng.Guo.combined, features = c("CRMP1", "EDAR", "KLRC1", "HBB", "CD177", "CCR10", "SLC4A10", "RRM2", "PAQR6"), min.cutoff = "q9")

# TODO: ANNOTATE clusters as specific cell types



