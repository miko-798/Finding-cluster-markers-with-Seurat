library(Seurat)
library(cowplot)

pruned_Zheng_seurat@meta.data$cancer_type <- "HCC"
pruned_Guo_seurat@meta.data$cancer_type <- "NSCLC"

dim(pruned_Zheng_seurat)
dim(pruned_Guo_seurat)

# To avoid NA in meta data
# append a unique identifier to the cells in each object before merging
pruned_Zheng_seurat <- RenameCells(object = pruned_Zheng_seurat, add.cell.id = "o1")
pruned_Guo_seurat <- RenameCells(object = pruned_Guo_seurat, add.cell.id = "o2")

pruned_Zheng_seurat@meta.data
pruned_Guo_seurat@meta.data

### Visualize the clusters before batch correction ##

# merge two seurat objects
merged_uncorrected <- merge(x=pruned_Zheng_seurat, y=pruned_Guo_seurat)
which(is.na(merged_uncorrected@meta.data$cancer_type))

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
p1 <- DimPlot(merged_uncorrected, reduction = "umap", group.by = "cancer_type")
p2 <- DimPlot(merged_uncorrected, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

# split.by: name of a metadata column to split plot by (side by side view)
DimPlot(merged_uncorrected, reduction = "umap", split.by = "cancer_type", label = TRUE)


######## Batch Correction: Integrate two datasets ########## 

#anyNA(pruned_Zheng_seurat@assays$RNA@counts)

immune.anchors <- FindIntegrationAnchors(object.list = list(pruned_Zheng_seurat, pruned_Guo_seurat), 
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
p1 <- DimPlot(Zheng.Guo.combined, reduction = "umap", group.by = "cancer_type")
p2 <- DimPlot(Zheng.Guo.combined, reduction = "umap", label = TRUE, label.size = 5)
plot_grid(p1, p2)

DimPlot(Zheng.Guo.combined, reduction = "umap", label = TRUE, label.size = 5, split.by = "cancer_type")


###### PART 2: Identify differentially expressed features (cluster biomarkers) ######

library(dplyr)

DefaultAssay(Zheng.Guo.combined) <- "RNA"

# cluster 6: these cells are more abundant in batch 1: HCC
cluster9.markers <- FindMarkers(Zheng.Guo.combined, ident.1 = 9, only.pos = TRUE)
head(cluster9.markers, n = 10)
#cluster6.markers %>% top_n(n=6, wt=avg_logFC)
#write.table(cluster.markers, "5.24.cluster6.markers.txt", sep="\t")

VlnPlot(Zheng.Guo.combined, features = c("SLC4A10", "ZBTB16"), pt.size=0)
FeaturePlot(Zheng.Guo.combined, features = c("SLC4A10", "ZBTB16"))


############### Find all markers #################

# find markers for every cluster compared to all remaining cells, report only the positive ones
all.markers <- FindAllMarkers(Zheng.Guo.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

## wt: sort by avg_logFC
all_markers_top2 <- all.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC) # select greated avg logFC

all_markers_top5 <- all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) # select greated avg logFC

all_markers_top20 <- all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC) # select greated avg logFC

all_markers_top50 <- all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC) # select greated avg logFC

#all.markers %>% group_by(cluster) %>% top_n(n = -1, wt = p_val) ## select smallest p value
#max(all.markers$avg_logFC)

write.table(all.markers, "6.10.all.markers.txt" , sep="\t")
write.table(all_markers_top2, "6.10.top2.all.markers.txt" , sep="\t")
write.table(all_markers_top5, "6.10.top5.all.markers.txt" , sep="\t")
write.table(all_markers_top20, "6.10.top20.all.markers.txt" , sep="\t")
write.table(all_markers_top50, "6.10.top50.all.markers.txt" , sep="\t")



## Visualize markers: indexed by Seurat's indices
# c0
VlnPlot(Zheng.Guo.combined, features = c("CAPG"), y.max=2.0 , pt.size=0)
FeaturePlot(Zheng.Guo.combined, features = c("CAPG"))

# c1
VlnPlot(Zheng.Guo.combined, features = c("CCR7", "ACTN1"), pt.size=0)
FeaturePlot(Zheng.Guo.combined, features = c("CCR7", "ACTN1"))

# c2
VlnPlot(Zheng.Guo.combined, features = c("PRSS23", "CX3CR1"), pt.size=0)
FeaturePlot(Zheng.Guo.combined, features = c("PRSS23", "CX3CR1"))

# c3
#VlnPlot(Zheng.Guo.combined, features = c("ATP5I", "COX17"), pt.size=0)
#FeaturePlot(Zheng.Guo.combined, features = c("ATP5I", "COX17"))

# c4
VlnPlot(Zheng.Guo.combined, features = c("CCR8", "IL1R2"), pt.size=0)
FeaturePlot(Zheng.Guo.combined, features = c("CCR8", "IL1R2"))

# c5: not ideal
VlnPlot(Zheng.Guo.combined, features = c("EEF1A1", "TPT1"), pt.size=0)
FeaturePlot(Zheng.Guo.combined, features = c("EEF1A1", "TPT1"))

# c6
VlnPlot(Zheng.Guo.combined, features = c("CCR4", "RTKN2"), pt.size=0)
FeaturePlot(Zheng.Guo.combined, features = c("CCR4", "RTKN2"))

# c7
VlnPlot(Zheng.Guo.combined, features = c("CCL3L1", "CCL3L3", "CCL3"), pt.size=0)
FeaturePlot(Zheng.Guo.combined, features = c("CCL3L1", "CCL3L3", "CCL3"))

# c8: results from seurat, not overlapping
VlnPlot(Zheng.Guo.combined, features = c("TVP23A", "NUBP1"), pt.size=0)
FeaturePlot(Zheng.Guo.combined, features = c("TVP23A", "NUBP1"))

# c9: HCC more expression - from seurat, not overlap
VlnPlot(Zheng.Guo.combined, features = c("COLQ", "SLC4A10"), pt.size=0)
FeaturePlot(Zheng.Guo.combined, features = c("COLQ", "SLC4A10"))

# c10 - from seurat, not overlap
#VlnPlot(Zheng.Guo.combined, features = c("CTLA4", "CXCL13"), pt.size=0)
VlnPlot(Zheng.Guo.combined, features = c("CPM", "CD200"), pt.size=0)
FeaturePlot(Zheng.Guo.combined, features = c("CPM", "CD200"))

# c11 - from seurat, not overlap
VlnPlot(Zheng.Guo.combined, features = c("RRM2", "SPC25"), pt.size=0)
FeaturePlot(Zheng.Guo.combined, features = c("RRM2", "SPC25"))

# c12
VlnPlot(Zheng.Guo.combined, features = c("SLC47A1", "LENG9", "MIR142"), pt.size=0)
FeaturePlot(Zheng.Guo.combined, features = c("SLC47A1", "LENG9", "MIR142"))

# c13 - from seurat, not overlap
VlnPlot(Zheng.Guo.combined, features = c("RNF219-AS1", "RNF219"), pt.size=0)
FeaturePlot(Zheng.Guo.combined, features = c("RNF219-AS1", "RNF219"))


# ANNOTATE clusters as specific cell types



