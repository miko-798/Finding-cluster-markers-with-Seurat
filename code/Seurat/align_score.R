
# This file calculates the alignment score (definition from the Butler et al
# paper) of a particular dataset containing different batches. 


# ctrl + l to clear
# esc to exit, stop the run
# remember to load libraries

library(SingleCellExperiment)
library(SC3)
library(scater)
library(FNN)
library(RDocumentation)
library(rlist)
library("ggplot2")
library(ggfortify)

head(ann)

############ Load data ##################
# Read tabular data into R
Zheng_data <- read.table("GSE98638_HCC.TCell.S4070.norm.centered.txt", 
                         header = TRUE, sep = "", dec = ".")

# Create Zheng_1: remove col and row names
Zheng_1 <- Zheng_data
Zheng_1[1:2] <- list(NULL) # remove first two cols
names(Zheng_1) <- NULL # remove first row(header)

# Create Zheng_2: filter data to get rid of negative values, by adding 11
Zheng_2 <- data.matrix(Zheng_1)
Zheng_2 <- Zheng_2 + 11

# Paper 2
Guo_data <- read.table("GSE99254_NSCLC.TCell.S11769.norm.centered.txt", 
                       header = TRUE, sep = "", dec = ".")
# Create Guo_1: remove col and row names
Guo_1 <- Guo_data
Guo_1[1:2] <- list(NULL) # remove first two cols
names(Guo_1) <- NULL # remove first row(header)

# Create Guo_2: filter data to get rid of negative values, by adding 13
Guo_2 <- data.matrix(Guo_1)
Guo_2 <- Guo_2 + 13

# First test:
# make a copy of Zheng_2: Zheng_2_copy
#Zheng_2_copy <- Zheng_2
#list_dataset <- list(Zheng_2, Zheng_2_copy)

# Second test: Combining Zheng and Guo data

# GET THE OVERLAPPING GENES 

# set the row names as the gene symbols
rownames(Zheng_2) <- Zheng_data[[2]]
rownames(Guo_2) <- Guo_data[[2]]

# set col names as the patient sample names 
colnames(Zheng_2) <- colnames(Zheng_data)[3:(length(colnames(Zheng_data)))]
colnames(Guo_2) <- colnames(Guo_data)[3:(length(colnames(Guo_data)))]

# prune matrices: only keep rows with the shared genes
overlap_gene <- intersect(rownames(Zheng_2), rownames(Guo_2))
Zheng_3 <- Zheng_2[rownames(Zheng_2) %in% overlap_gene, ]
Guo_3 <- Guo_2[rownames(Guo_2) %in% overlap_gene, ]

dim(Zheng_3)
dim(Guo_3)
 

# COMBINE TWO MATRICES

Zheng_Guo_list <- list(Zheng_3, Guo_3)
n_batch <- length(Zheng_Guo_list)
n_batch
cell_per_dataset <- get_min_cells(Zheng_Guo_list)
cell_per_dataset
down_sam_Zheng_Guo <- down_sample(Zheng_Guo_list, cell_per_dataset)

dim(down_sam_Zheng_Guo[[1]])
dim(down_sam_Zheng_Guo[[2]])

#rownames(down_sam_Zheng_Guo[[1]])
#length(rownames(down_sam_Zheng_Guo[[1]]))

# Merge datasets: combine two or more dataset into one (add more cell cols)
# IMPORTANT NOTE: assume same genes in the dataset, only cells are different

# testing
#m1 <- down_sam_Zheng_Guo[[1]][1:5, 1:5]
#m2 <- down_sam_Zheng_Guo[[2]][1:5, 1:5]
#merged <- cbind(m1, m2)

Zheng_Guo <- cbind(down_sam_Zheng_Guo[[1]], down_sam_Zheng_Guo[[2]])
dim(Zheng_Guo)

# calculate alignment score: FUNCTIONS BELOW
percent_neighbors = 0.01
k <-round(percent_neighbors*cell_per_dataset, digits=0)
knn <- find_KNN(Zheng_Guo, k)
x_bar <- find_x_bar(knn)
align_score <- alignment_score(x_bar, k, n_batch, dim(Zheng_Guo)[2])
align_score 

# PCA graph
# Reference: https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html

Zheng_Guo.pca <- prcomp(Zheng_Guo, center = TRUE, scale. = TRUE)
summary(Zheng_Guo.pca)
str(Zheng_Guo.pca)

# CORRECT?
autoplot(Zheng_Guo.pca, colour = Zheng_Guo$batch.info)

## This may not work
d.factanal <- factanal(Zheng_Guo, factors = 3, scores = 'regression')
autoplot(d.factanal, data = Zheng_Guo, colour = 'Gene expression')


########## Functions to calculate ALIGNMENT SCORE  #########

# Function
# get min cells
get_min_cells <- function(list_dataset){
  min_cells <- Inf
  for(i in list_dataset){
    if(dim(i)[2] < min_cells){  # dim[i](2) is the number of cells
      min_cells <- dim(i)[2]
    }
  }
  return(min_cells)
}

# Function
# randomly downsampling: keep "min_cells" number of cols 
# returns a list 
down_sample <- function(list_dataset, min_cells){
  down_samp_dataset <- list()
  for(i in list_dataset){
    i <- i[, sample(ncol(i), min_cells)]
    print(type(i))
    print(dim(i))
    
    down_samp_dataset <- unlist(list(down_samp_dataset, list(i)), recursive=FALSE)
  }
  return(down_samp_dataset)
}

# Function
# Find nearest neighbors for all cells: default is 1% closest neighbors
# Reference: https://www.rdocumentation.org/packages/FNN/versions/1.1.3/topics/get.knn

find_KNN <- function(combined_data, k){
  # get knn (get transpose the matrix)
  KNN <- get.knn(t(combined_data), k, algorithm=c("kd_tree", "cover_tree", "CR", "brute"))
  return(KNN)
}

# Function:
# find x_bar
# for each cell, count the number of NN from its own batch, then average
find_x_bar <- function(KNN){
  # this list stores the number of NN in the same batch for each cell
  same_batch_nn_list <- list()
  
  for (row_num in 1:nrow(KNN$nn.index)){
    num_same_batch_nn <- 0
    for (col_num in 1:ncol(KNN$nn.index)){
      # get the batch number
      neighbor_index <- KNN$nn.index[row_num, col_num]
      cell_index <- row_num
      neighbor_batch <- get_batch(cell_per_dataset, neighbor_index)
      cell_batch <- get_batch(cell_per_dataset, cell_index)
      if (neighbor_batch == cell_batch){
        num_same_batch_nn <- num_same_batch_nn + 1
      }
    }
    same_batch_nn_list <- c(same_batch_nn_list, num_same_batch_nn)
  }
  x_bar <- mean(as.numeric(same_batch_nn_list))
  return(x_bar)
}

# Function
# returns the batch number at index of the knn result nxk matrix
get_batch <- function(min_cells, index){
  batch_num <- floor((index - 1) / min_cells + 1)
  return (batch_num)
}

# Function
# calculate alignment score
alignment_score <- function(x_bar_value, num_k, num_batch, N){
  align_score <- 0
  align_score <- num_batch * (1 - ( (x_bar_value - num_k/N) / (num_k - num_k/N) ) )
  return(align_score)
}


