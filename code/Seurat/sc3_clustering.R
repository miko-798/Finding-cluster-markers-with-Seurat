# This file tries to do SC3 clustering on the Zheng dataset

############ SC3 Clustering ###############
# create a SingleCellExperiment object
Zheng_sc <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(Zheng_2),
    logcounts = log2(as.matrix(Zheng_2) + 1)
  )
)

# reassign genesymbol to the matrix
# define feature names in feature_symbol column
Zheng_data[,2]
rowData(Zheng_sc)$feature_symbol <- Zheng_data[,2]
# remove features with duplicated names
Zheng_sc <- Zheng_sc[!duplicated(rowData(Zheng_sc)$feature_symbol), ]

# define spike-ins
# isSpike(sce, "ERCC") <- grepl("ERCC", rowData(sce)$feature_symbol)

plotPCA(Zheng_sc)

Zheng_sc <- sc3(Zheng_sc, ks = 10:12, biology = TRUE, gene_filter = FALSE)

# SC3 results in colData
col_data <- colData(Zheng_sc)
head(col_data[ , grep("sc3_", colnames(col_data))])

plotPCA(Zheng_sc, colour_by = "sc3_12_clusters")

# SC3 results in rowData
row_data <- rowData(Zheng_sc)
head(row_data[ , grep("sc3_", colnames(row_data))])
