###data

#load data

tabib.data  <- Read10X("C:/Users/labor/Desktop/patrick/single/single/data/set_4/MD001_DK1/filtered_feature_bc_matrix")


#set up tabib data
tabib.keloid <- CreateSeuratObject(counts=test , assay="RNA")

#subset samples

