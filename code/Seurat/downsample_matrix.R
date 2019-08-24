# downsample to 100 cells and 50 genes

Zheng_small <- down_sam_Zheng_Guo[[1]][1:50, 1:100]
dim(Zheng_small)

Guo_small <- down_sam_Zheng_Guo[[2]][1:50, 1:100]
dim(Guo_small)


# entire matrix
Zheng_all <- down_sam_Zheng_Guo[[1]]
dim(Zheng_all)

Guo_all <- down_sam_Zheng_Guo[[2]]
dim(Guo_all)
