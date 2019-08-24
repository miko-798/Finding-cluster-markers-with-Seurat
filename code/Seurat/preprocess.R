# Preprocess the data before using Seurat 
# add gene symbols for rows

library(dplyr)
library(Seurat)
library(tidyverse)

Zheng_data <- read.table("GSE98638_HCC.TCell.S4070.norm.centered.txt", 
                         header = TRUE, sep = "", dec = ".")

Guo_data <- read.table("GSE99254_NSCLC.TCell.S11769.norm.centered.txt", 
                       header = TRUE, sep = "", dec = ".")

# remove duplicate gene symbols 
Zheng_1 <- Zheng_data
Zheng_1 <- Zheng_1[!duplicated(Zheng_1$geneSymbol), ]
dim(Zheng_1) # 14093 x  4072

# remove gene and cell names, add back in later
Zheng_2 <- Zheng_1
Zheng_2[1:2] <- list(NULL) # remove first two cols (indicating genes)
names(Zheng_2) <- NULL # remove first row(header: indicating cells)

# filter data to get rid of negative values, by adding 11
Zheng_2 <- data.matrix(Zheng_2)
Zheng_2 <- Zheng_2 + 11

# same idea for the other dataset
Guo_1 <- Guo_data
Guo_1 <- Guo_1[!duplicated(Guo_1$geneSymbol), ]
dim(Guo_1) # 12395 x 11771

Guo_2 <- Guo_1
Guo_2[1:2] <- list(NULL)
names(Guo_2) <- NULL 

# filter data to get rid of negative values, by adding 13
Guo_2 <- data.matrix(Guo_2)
Guo_2 <- Guo_2 + 13

# set the row names as the gene symbols (for the purpose of finding overlap genes)
rownames(Zheng_2) <- Zheng_1[[2]]
rownames(Guo_2) <- Guo_1[[2]]

# set col names as the patient sample names 
colnames(Zheng_2) <- colnames(Zheng_1)[3:(length(colnames(Zheng_1)))]
colnames(Guo_2) <- colnames(Guo_1)[3:(length(colnames(Guo_1)))]

# GET THE OVERLAPPING GENES 
# prune matrices: only keep rows with the shared genes
overlap_gene <- intersect(rownames(Zheng_2), rownames(Guo_2))
Zheng_3 <- Zheng_2[rownames(Zheng_2) %in% overlap_gene, ]
Guo_3 <- Guo_2[rownames(Guo_2) %in% overlap_gene, ]

dim(Zheng_3)
dim(Guo_3)

# set rownames as gene symbols (for finding markers)
names(Zheng_3) <- rownames(Zheng_3)
names(Guo_3) <- rownames(Guo_3)


#anyNA(rownames(Zheng_3))
#which(is.na(rownames(Zheng_3))) # 30

#anyNA(rownames(Guo_3))
#which(is.na(rownames(Guo_3))) # 63

Zheng_4 <- Zheng_3[-c(30),]
Guo_4 <- Guo_3[-c(63),]
dim(Zheng_4)
dim(Guo_4)

