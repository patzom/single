###data
library(Seurat)
#load data

tabib.data  <- Read10X("C:/Users/labor/Desktop/patrick/single/single/data/set_4/MD001_DK1/filtered_feature_bc_matrix")


#set up tabib data
tabib.keloid <- CreateSeuratObject(counts=tabib.data , assay="RNA")

#subset samples

tabib.keloid[c("CD3D", "TCL1A", "MS4A1"), 1:300]
