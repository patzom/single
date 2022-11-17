


#
#D002_computation12_skin and keloid__single integrated.sct
#1 skin samples, 4 keloid samples, 6 skin tabib
#Update Seurat 4, R4


setwd(dir ="D:/Direder/Projekte/D002_10xScarWars/D002_Berechnungen/D002_c12_7s4k")


#define color schemes
{
  rainbow.colors<-palette(rainbow(25))
  library(Polychrome)
  color_UMAP <- as.character(c("#990F26", "#B33E52", "#CC7A88", "#E6B8BF", 
                               "#99600F", "#B3823E",  
                               "#54990F", "#78B33E","#03AC13",
                               "#45b6fe", "#3792cb","#296d98", 
                               "#F9A602", "#967ACC", 
                               "#666666", "#333333", "#999999", "#CCCCCC","#3D0F99"))
  
  color_UMAP_short <- as.character(c("#990F26", "#CC7A88", 
                                     "#99600F",   
                                     "#54990F", "#78B33E",
                                     "#45b6fe", "#3792cb","#296d98", 
                                     "#F9A602", "#967ACC", 
                                     "#666666", "#333333", "#999999", "#CCCCCC","#3D0F99"))
  
  color_UMAP_short_double  <- as.character(c("#990F26","#990F26", "#CC7A88","#CC7A88", 
                                             "#99600F","#99600F",   
                                             "#54990F","#54990F", "#78B33E","#78B33E",
                                             "#45b6fe","#45b6fe", "#3792cb","#3792cb","#296d98","#296d98", 
                                             "#F9A602","#F9A602", "#967ACC","#967ACC", 
                                             "#666666","#666666", "#333333","#333333", "#999999","#999999", "#CCCCCC","#CCCCCC","#3D0F99","#3D0F99"))
  color_UMAP2 <- as.character(c("#1F78C8", "#0000FF","#a6cee3",  "#ff0000", "#33a02c", 
                                "#6A33C2", "#ff7f00", "#FFD700", 
                                "#36648B", 
                                "#00E2E5", "#00FF00", "#778B00", "#BEBE00", 
                                "#8B3B00", "#A52A3C", "#FB6496", "#b2df8a", "#CAB2D6", 
                                "#FDBF6F", "#999999", "#EEE685", "#C8308C", 
                                "#FF83FA", "#C814FA"))
  color_SC <- as.character(c("#338333","#67B932",  "#FFBF00",  "#1F78C8",  "#C8308C", 
                             "#FF83FA", "#F8766D",  "#0000FF",  "#36648B", "#ff0000", "#33a02c","#6A33C2", 
                             "#00E2E5", "#00FF00", "#778B00", "#BEBE00", 
                             "#8B3B00", "#A52A3C", "#FB6496",  "#b2df8a", "#CAB2D6", 
                             "#FDBF6F", "#999999","#C814FA","#FFD700","#ff7f00","#EEE685"))
  color_MAC <- as.character(c("#BA68CB","#00EC9D","#00E2FC","#8530d1"))
  color_DEG <- as.character(c("#54990F","#C51D34"))
  single_col <- as.character("#02025f")
  color_tissue <- as.character(c("#338333","#80CAED","#F8766D"))
  library(Seurat)
  library(dplyr)
  library(magrittr)
  library(ggplot2)
  library(xlsx)
  library(patchwork)
  library(sctransform)
  library(monocle3)
  library(SeuratWrappers)
  library(ggrepel)
  library(tidyr)
  library(clustree)
  library(ggpubr)
  library(ggsignif)
  library(limma)
  library(Matrix.utils)
  library(rstatix)
}

test  <- Read10X("C:/Users/labor/Desktop/patrick/single/single/data/set_4/MD001_DK1/filtered_feature_bc_matrix")
testseu<-CreateSeuratObject(counts=test , assay="RNA")


#############################Skin-Data wrangling############################
###-###-###-###-###-###-###-###-###-###-###-###-###

#load tabib_skin.data
tabib.data<-read.csv("D:/Direder/Projekte/D002_10xScarWars/D002_Data/10x_Tabib_data/Raw Data Tabib/Skin_6Control_rawUMI.csv")
rownames(tabib.data)<-tabib.data$X
tabib.data<-tabib.data[,-1]

tabib.metadata<-read.csv("D:/Direder/Projekte/D002_10xScarWars/D002_Data/10x_Tabib_data/Raw Data Tabib/Skin_6Control_Metadata.csv")  
rownames(tabib.metadata)<-tabib.metadata$X
tabib.metadata<-tabib.metadata[,-1]


#setup Tabib data
tabib.skin<-CreateSeuratObject(counts=tabib.data, assay="RNA", meta.data = tabib.metadata)

#subset samples
skin1.D002.1.sct <- subset(tabib.skin, ident=c("SC18control"))
skin2.D002.1.sct <- subset(tabib.skin, ident=c("SC1control"))
skin3.D002.1.sct <- subset(tabib.skin, ident=c("SC32control"))
skin4.D002.1.sct <- subset(tabib.skin, ident=c("SC33control"))
skin5.D002.1.sct <- subset(tabib.skin, ident=c("SC34control"))
skin6.D002.1.sct <- subset(tabib.skin, ident=c("SC4control"))

#Get original 10x Data, convert to NCBI Symbol, Create Seurat Object
#skin1
skin1.D002.1.sct.data <- GetAssay(skin1.D002.1.sct, assay="RNA")
skin1.D002.1.sct.data<-skin1.D002.1.sct.data@counts
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(skin1.D002.1.sct.data)),paste0(c("D:/Direder/SC-EnsembID/Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(skin1.D002.1.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(skin1.D002.1.sct.data)= as.character(Signlist2$new)


#skin2
skin2.D002.1.sct.data <- GetAssay(skin2.D002.1.sct, assay="RNA")
skin2.D002.1.sct.data<-skin2.D002.1.sct.data@counts
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(skin2.D002.1.sct.data)),paste0(c("D:/Direder/SC-EnsembID/Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(skin2.D002.1.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(skin2.D002.1.sct.data)= as.character(Signlist2$new)


#skin3
skin3.D002.1.sct.data <- GetAssay(skin3.D002.1.sct, assay="RNA")
skin3.D002.1.sct.data<-skin3.D002.1.sct.data@counts
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(skin3.D002.1.sct.data)),paste0(c("D:/Direder/SC-EnsembID/Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(skin3.D002.1.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(skin3.D002.1.sct.data)= as.character(Signlist2$new)



#skin4
skin4.D002.1.sct.data <- GetAssay(skin4.D002.1.sct, assay="RNA")
skin4.D002.1.sct.data<-skin4.D002.1.sct.data@counts
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(skin4.D002.1.sct.data)),paste0(c("D:/Direder/SC-EnsembID/Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(skin4.D002.1.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(skin4.D002.1.sct.data)= as.character(Signlist2$new)


#skin5
skin5.D002.1.sct.data <- GetAssay(skin5.D002.1.sct, assay="RNA")
skin5.D002.1.sct.data<-skin5.D002.1.sct.data@counts
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(skin5.D002.1.sct.data)),paste0(c("D:/Direder/SC-EnsembID/Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(skin5.D002.1.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(skin5.D002.1.sct.data)= as.character(Signlist2$new)


#skin6
skin6.D002.1.sct.data <- GetAssay(skin6.D002.1.sct, assay="RNA")
skin6.D002.1.sct.data<-skin6.D002.1.sct.data@counts
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(skin6.D002.1.sct.data)),paste0(c("D:/Direder/SC-EnsembID/Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(skin6.D002.1.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(skin6.D002.1.sct.data)= as.character(Signlist2$new)



#load skin 7 (MD)
skin7.D002.1.sct.data<-Read10X("D:/Direder/Projekte/D002_10xScarWars/D002_Data/D002_healthy skin_Direder/MD007 DN2 NH human")

#convert to NCBI Symbol
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(skin7.D002.1.sct.data)),paste0(c("D:/Direder/SC-EnsembID/Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(skin7.D002.1.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(skin7.D002.1.sct.data)= as.character(Signlist2$new)

# load keloid
keloid1.D002.1.sct.data<-Read10X("D:/Direder/Projekte/D002_10xScarWars/D002_Data/D002_DK1_Keloid1/filtered feature bc matrix")
keloid2.D002.1.sct.data<-Read10X("D:/Direder/Projekte/D002_10xScarWars/D002_Data/D002_DK3_Keloid3/filtered feature bc matrix")
keloid3L.D002.1.sct.data<-Read10X("D:/Direder/Projekte/D002_10xScarWars/D002_Data/D002_DK4L_Keloid4L/filtered feature bc matrix")
keloid3R.D002.1.sct.data<-Read10X("D:/Direder/Projekte/D002_10xScarWars/D002_Data/D002_DK4R_Keloid4R/filtered feature bc matrix")

#convert to NCBI Symbol
#keloid1
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(keloid1.D002.1.sct.data)),paste0(c("D:/Direder/SC-EnsembID/Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(keloid1.D002.1.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(keloid1.D002.1.sct.data)= as.character(Signlist2$new)

#keloid2
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(keloid2.D002.1.sct.data)),paste0(c("D:/Direder/SC-EnsembID/Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(keloid2.D002.1.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(keloid2.D002.1.sct.data)= as.character(Signlist2$new)

#keloid3L
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(keloid3L.D002.1.sct.data)),paste0(c("D:/Direder/SC-EnsembID/Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(keloid3L.D002.1.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(keloid3L.D002.1.sct.data)= as.character(Signlist2$new)

#keloid3R
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(keloid3R.D002.1.sct.data)),paste0(c("D:/Direder/SC-EnsembID/Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(keloid3R.D002.1.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(keloid3R.D002.1.sct.data)= as.character(Signlist2$new)

#check gene names
#s1<-rownames(skin1.D002.1.sct.data)
#s2<-rownames(skin2.D002.1.sct.data)
#s3<-rownames(skin3.D002.1.sct.data)
#s4<-rownames(skin4.D002.1.sct.data)
#s5<-rownames(skin5.D002.1.sct.data)
#s6<-rownames(skin6.D002.1.sct.data)
#s7<-rownames(skin7.D002.1.sct.data)
#k1<-rownames(keloid1.D002.1.sct.data)
#k2<-rownames(keloid2.D002.1.sct.data)
#k3l<-rownames(keloid3L.D002.1.sct.data)
#k3r<-rownames(keloid3R.D002.1.sct.data)

#write.xlsx(s1,"Lists/features detected/rownamess1.xlsx")
#write.xlsx(s2,"Lists/features detected/rownamess2.xlsx")
#write.xlsx(s3,"Lists/features detected/rownamess3.xlsx")
#write.xlsx(s4,"Lists/features detected/rownamess4.xlsx")
#write.xlsx(s5,"Lists/features detected/rownamess5.xlsx")
#write.xlsx(s6,"Lists/features detected/rownamess6.xlsx")
#write.xlsx(s7,"Lists/features detected/rownamess7.xlsx")
#write.xlsx(k1,"Lists/features detected/rownamesk1.xlsx")
#write.xlsx(k2,"Lists/features detected/rownamesk2.xlsx")
#write.xlsx(k3l,"Lists/features detected/rownamesk3l.xlsx")
#write.xlsx(k3r,"Lists/features detected/rownamesk3r.xlsx")

# keep only Genes detected in all dataset Template
dim(skin7.D002.1.sct.data)
dim(skin1.D002.1.sct.data)
dim(keloid1.D002.1.sct.data)

ROWLIST1 <- skin7.D002.1.sct.data[!duplicated(rownames(skin7.D002.1.sct.data)),]
ROWLIST2 <- skin1.D002.1.sct.data[!duplicated(rownames(skin1.D002.1.sct.data)),]
ROWLIST3 <- keloid1.D002.1.sct.data[!duplicated(rownames(keloid1.D002.1.sct.data)),]

Temp1<-merge(ROWLIST1,ROWLIST2, by.x = rownames(ROWLIST1), by.y=rownames(ROWLIST2),all.x=F, all.y=F)
Temp1 <- Temp1[,1:2]
dim(ROWLIST1)
dim(ROWLIST2)
dim(Temp1)

Temp2<-merge(Temp1,ROWLIST3, by.x = rownames(Temp1), by.y=rownames(ROWLIST3),all.x=F, all.y=F)
Temp2 <- Temp2[,1:2]
dim(Temp1)
dim(ROWLIST3)
dim(Temp2)


which(duplicated(rownames(skin7.D002.1.sct.data)))
dim(skin7.D002.1.sct.data)
dim(ROWLIST1)

#remove duplicates
skin1.D002.1.sct.data <- skin1.D002.1.sct.data[!duplicated(rownames(skin1.D002.1.sct.data)),]
skin2.D002.1.sct.data <- skin2.D002.1.sct.data[!duplicated(rownames(skin2.D002.1.sct.data)),]
skin3.D002.1.sct.data <- skin3.D002.1.sct.data[!duplicated(rownames(skin3.D002.1.sct.data)),]
skin4.D002.1.sct.data <- skin4.D002.1.sct.data[!duplicated(rownames(skin4.D002.1.sct.data)),]
skin5.D002.1.sct.data <- skin5.D002.1.sct.data[!duplicated(rownames(skin5.D002.1.sct.data)),]
skin6.D002.1.sct.data <- skin6.D002.1.sct.data[!duplicated(rownames(skin6.D002.1.sct.data)),]
skin7.D002.1.sct.data <- skin7.D002.1.sct.data[!duplicated(rownames(skin7.D002.1.sct.data)),]
keloid1.D002.1.sct.data <- keloid1.D002.1.sct.data[!duplicated(rownames(keloid1.D002.1.sct.data)),]
keloid2.D002.1.sct.data <- keloid2.D002.1.sct.data[!duplicated(rownames(keloid2.D002.1.sct.data)),]
keloid3L.D002.1.sct.data <- keloid3L.D002.1.sct.data[!duplicated(rownames(keloid3L.D002.1.sct.data)),]
keloid3R.D002.1.sct.data <- keloid3R.D002.1.sct.data[!duplicated(rownames(keloid3R.D002.1.sct.data)),]

# adapt all
Temp3<-join.Matrix(Temp2,skin1.D002.1.sct.data, by.x = rownames(Temp2), by.y=rownames(skin1.D002.1.sct.data), all.x = F, all.y = F)
skin1.D002.1.sct.data <- Temp3[,3:ncol(Temp3)]

Temp3<-join.Matrix(Temp2,skin2.D002.1.sct.data, by.x = rownames(Temp2), by.y=rownames(skin2.D002.1.sct.data), all.x = F, all.y = F)
skin2.D002.1.sct.data <- Temp3[,3:ncol(Temp3)]

Temp3<-join.Matrix(Temp2,skin3.D002.1.sct.data, by.x = rownames(Temp2), by.y=rownames(skin3.D002.1.sct.data), all.x = F, all.y = F)
skin3.D002.1.sct.data <- Temp3[,3:ncol(Temp3)]

Temp3<-join.Matrix(Temp2,skin4.D002.1.sct.data, by.x = rownames(Temp2), by.y=rownames(skin4.D002.1.sct.data), all.x = F, all.y = F)
skin4.D002.1.sct.data <- Temp3[,3:ncol(Temp3)]

Temp3<-join.Matrix(Temp2,skin5.D002.1.sct.data, by.x = rownames(Temp2), by.y=rownames(skin5.D002.1.sct.data), all.x = F, all.y = F)
skin5.D002.1.sct.data <- Temp3[,3:ncol(Temp3)]

Temp3<-join.Matrix(Temp2,skin6.D002.1.sct.data, by.x = rownames(Temp2), by.y=rownames(skin6.D002.1.sct.data), all.x = F, all.y = F)
skin6.D002.1.sct.data <- Temp3[,3:ncol(Temp3)]

Temp3<-join.Matrix(Temp2,skin7.D002.1.sct.data, by.x = rownames(Temp2), by.y=rownames(skin7.D002.1.sct.data), all.x = F, all.y = F)
skin7.D002.1.sct.data <- Temp3[,3:ncol(Temp3)]

Temp3<-join.Matrix(Temp2,keloid1.D002.1.sct.data, by.x = rownames(Temp2), by.y=rownames(keloid1.D002.1.sct.data), all.x = F, all.y = F)
keloid1.D002.1.sct.data <- Temp3[,3:ncol(Temp3)]

Temp3<-join.Matrix(Temp2,keloid2.D002.1.sct.data, by.x = rownames(Temp2), by.y=rownames(keloid2.D002.1.sct.data), all.x = F, all.y = F)
keloid2.D002.1.sct.data <- Temp3[,3:ncol(Temp3)]

Temp3<-join.Matrix(Temp2,keloid3L.D002.1.sct.data, by.x = rownames(Temp2), by.y=rownames(keloid3L.D002.1.sct.data), all.x = F, all.y = F)
keloid3L.D002.1.sct.data <- Temp3[,3:ncol(Temp3)]

Temp3<-join.Matrix(Temp2,keloid3R.D002.1.sct.data, by.x = rownames(Temp2), by.y=rownames(keloid3R.D002.1.sct.data), all.x = F, all.y = F)
keloid3R.D002.1.sct.data <- Temp3[,3:ncol(Temp3)]


#Create Sparcematrix
skin1.D002.1.sct<-CreateSeuratObject(skin1.D002.1.sct.data)
skin2.D002.1.sct<-CreateSeuratObject(skin2.D002.1.sct.data)
skin3.D002.1.sct<-CreateSeuratObject(skin3.D002.1.sct.data)
skin4.D002.1.sct<-CreateSeuratObject(skin4.D002.1.sct.data)
skin5.D002.1.sct<-CreateSeuratObject(skin5.D002.1.sct.data)
skin6.D002.1.sct<-CreateSeuratObject(skin6.D002.1.sct.data)
skin7.D002.1.sct<-CreateSeuratObject(skin7.D002.1.sct.data)
keloid1.D002.1.sct<-CreateSeuratObject(keloid1.D002.1.sct.data)
keloid2.D002.1.sct<-CreateSeuratObject(keloid2.D002.1.sct.data)
keloid3L.D002.1.sct<-CreateSeuratObject(keloid3L.D002.1.sct.data)
keloid3R.D002.1.sct<-CreateSeuratObject(keloid3R.D002.1.sct.data)

#flag-1_sample
skin1.D002.1.sct<-AddMetaData(skin1.D002.1.sct, "skin_1", col.name = "sample")
skin2.D002.1.sct<-AddMetaData(skin2.D002.1.sct, "skin_2", col.name = "sample")
skin3.D002.1.sct<-AddMetaData(skin3.D002.1.sct, "skin_3", col.name = "sample")
skin4.D002.1.sct<-AddMetaData(skin4.D002.1.sct, "skin_4", col.name = "sample")
skin5.D002.1.sct<-AddMetaData(skin5.D002.1.sct, "skin_5", col.name = "sample")
skin6.D002.1.sct<-AddMetaData(skin6.D002.1.sct, "skin_6", col.name = "sample")
skin7.D002.1.sct<-AddMetaData(skin7.D002.1.sct, "skin_7", col.name = "sample")
keloid1.D002.1.sct<-AddMetaData(keloid1.D002.1.sct, "keloid_1", col.name = "sample")
keloid2.D002.1.sct<-AddMetaData(keloid2.D002.1.sct, "keloid_2", col.name = "sample")
keloid3L.D002.1.sct<-AddMetaData(keloid3L.D002.1.sct, "keloid_3L", col.name = "sample")
keloid3R.D002.1.sct<-AddMetaData(keloid3R.D002.1.sct, "keloid_3R", col.name = "sample")

#flag-2_tissue
skin1.D002.1.sct<-AddMetaData(skin1.D002.1.sct, "skin", col.name = "tissue")
skin2.D002.1.sct<-AddMetaData(skin2.D002.1.sct, "skin", col.name = "tissue")
skin3.D002.1.sct<-AddMetaData(skin3.D002.1.sct, "skin", col.name = "tissue")
skin4.D002.1.sct<-AddMetaData(skin4.D002.1.sct, "skin", col.name = "tissue")
skin5.D002.1.sct<-AddMetaData(skin5.D002.1.sct, "skin", col.name = "tissue")
skin6.D002.1.sct<-AddMetaData(skin6.D002.1.sct, "skin", col.name = "tissue")
skin7.D002.1.sct<-AddMetaData(skin7.D002.1.sct, "skin", col.name = "tissue")
keloid1.D002.1.sct<-AddMetaData(keloid1.D002.1.sct, "keloid", col.name = "tissue")
keloid2.D002.1.sct<-AddMetaData(keloid2.D002.1.sct, "keloid", col.name = "tissue")
keloid3L.D002.1.sct<-AddMetaData(keloid3L.D002.1.sct, "keloid", col.name = "tissue")
keloid3R.D002.1.sct<-AddMetaData(keloid3R.D002.1.sct, "keloid", col.name = "tissue")

#flag-3_origin
skin1.D002.1.sct<-AddMetaData(skin1.D002.1.sct, "TT", col.name = "origin")
skin2.D002.1.sct<-AddMetaData(skin2.D002.1.sct, "TT", col.name = "origin")
skin3.D002.1.sct<-AddMetaData(skin3.D002.1.sct, "TT", col.name = "origin")
skin4.D002.1.sct<-AddMetaData(skin4.D002.1.sct, "TT", col.name = "origin")
skin5.D002.1.sct<-AddMetaData(skin5.D002.1.sct, "TT", col.name = "origin")
skin6.D002.1.sct<-AddMetaData(skin6.D002.1.sct, "TT", col.name = "origin")
skin7.D002.1.sct<-AddMetaData(skin7.D002.1.sct, "MD", col.name = "origin")
keloid1.D002.1.sct<-AddMetaData(keloid1.D002.1.sct, "MD", col.name = "origin")
keloid2.D002.1.sct<-AddMetaData(keloid2.D002.1.sct, "MD", col.name = "origin")
keloid3L.D002.1.sct<-AddMetaData(keloid3L.D002.1.sct, "MD", col.name = "origin")
keloid3R.D002.1.sct<-AddMetaData(keloid3R.D002.1.sct, "MD", col.name = "origin")


###skin_tabib###Quality ctrl, scale, normalize
#Quality ctrl, scale, normalize, skin1.D002.1.sct
skin1.D002.1.sct[["percent.mt"]] <- PercentageFeatureSet(skin1.D002.1.sct, pattern = "^MT-")
#skin1.D002.1.sct[["percent.ERY"]] <- PercentageFeatureSet(skin1.D002.1.sct, pattern = "HBB")
#skin1.D002.1.sct <- subset(skin1.D002.1.sct, subset = percent.ERY < 5 )

VlnPlot(skin1.D002.1.sct, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)

pdf("Graphs/Quality control/VlnPlots_nfeature_ncout_percentmt_skin1.D002.1.sct.pdf")
VlnPlot(skin1.D002.1.sct, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
dev.off()

pdf("Graphs/Quality control/FeatureScatter_nfeature_ncount_skin1.D002.1.sct.pdf")
FeatureScatter(skin1.D002.1.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

pdf("Graphs/Quality control/FeatureScatter_percent.mt_ncount_skin1.D002.1.sct.pdf")
FeatureScatter(skin1.D002.1.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")
dev.off()

skin1.D002.1.sct <- SCTransform(skin1.D002.1.sct, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)

#Quality ctrl, scale, normalize, skin2.D002.1.sct
skin2.D002.1.sct[["percent.mt"]] <- PercentageFeatureSet(skin2.D002.1.sct, pattern = "^MT-")
#skin2.D002.1.sct[["percent.ERY"]] <- PercentageFeatureSet(skin2.D002.1.sct, pattern = "HBB")
#skin2.D002.1.sct <- subset(skin2.D002.1.sct, subset = percent.ERY < 5 )

pdf("Graphs/Quality control/VlnPlots_nfeature_ncout_percentmt_skin2.D002.1.sct.pdf")
VlnPlot(skin2.D002.1.sct, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
dev.off()

pdf("Graphs/Quality control/FeatureScatter_nfeature_ncount_skin2.D002.1.sct.pdf")
FeatureScatter(skin2.D002.1.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

pdf("Graphs/Quality control/FeatureScatter_percent.mt_ncount_skin2.D002.1.sct.pdf")
FeatureScatter(skin2.D002.1.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")
dev.off()

skin2.D002.1.sct <- SCTransform(skin2.D002.1.sct, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)

#Quality ctrl, scale, normalize, skin3.D002.1.sct
skin3.D002.1.sct[["percent.mt"]] <- PercentageFeatureSet(skin3.D002.1.sct, pattern = "^MT-")
#skin3.D002.1.sct[["percent.ERY"]] <- PercentageFeatureSet(skin3.D002.1.sct, pattern = "HBB")
#skin3.D002.1.sct <- subset(skin3.D002.1.sct, subset = percent.ERY < 5 )

pdf("Graphs/Quality control/VlnPlots_nfeature_ncout_percentmt_skin3.D002.1.sct.pdf")
VlnPlot(skin3.D002.1.sct, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
dev.off()

pdf("Graphs/Quality control/FeatureScatter_nfeature_ncount_skin3.D002.1.sct.pdf")
FeatureScatter(skin3.D002.1.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

pdf("Graphs/Quality control/FeatureScatter_percent.mt_ncount_skin3.D002.1.sct.pdf")
FeatureScatter(skin3.D002.1.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")
dev.off()

skin3.D002.1.sct <- SCTransform(skin3.D002.1.sct, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)

#Quality ctrl, scale, normalize, skin4.D002.1.sct
skin4.D002.1.sct[["percent.mt"]] <- PercentageFeatureSet(skin4.D002.1.sct, pattern = "^MT-")
#skin4.D002.1.sct[["percent.ERY"]] <- PercentageFeatureSet(skin4.D002.1.sct, pattern = "HBB")
#skin4.D002.1.sct <- subset(skin4.D002.1.sct, subset = percent.ERY < 5 )

pdf("Graphs/Quality control/VlnPlots_nfeature_ncout_percentmt_skin4.D002.1.sct.pdf")
VlnPlot(skin4.D002.1.sct, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
dev.off()

pdf("Graphs/Quality control/FeatureScatter_nfeature_ncount_skin4.D002.1.sct.pdf")
FeatureScatter(skin4.D002.1.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

pdf("Graphs/Quality control/FeatureScatter_percent.mt_ncount_skin4.D002.1.sct.pdf")
FeatureScatter(skin4.D002.1.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")
dev.off()

skin4.D002.1.sct <- SCTransform(skin4.D002.1.sct, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)

#Quality ctrl, scale, normalize, skin5.D002.1.sct
skin5.D002.1.sct[["percent.mt"]] <- PercentageFeatureSet(skin5.D002.1.sct, pattern = "^MT-")
#skin5.D002.1.sct[["percent.ERY"]] <- PercentageFeatureSet(skin5.D002.1.sct, pattern = "HBB")
#skin5.D002.1.sct <- subset(skin5.D002.1.sct, subset = percent.ERY < 5 )

pdf("Graphs/Quality control/VlnPlots_nfeature_ncout_percentmt_skin5.D002.1.sct.pdf")
VlnPlot(skin5.D002.1.sct, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
dev.off()

pdf("Graphs/Quality control/FeatureScatter_nfeature_ncount_skin5.D002.1.sct.pdf")
FeatureScatter(skin5.D002.1.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

pdf("Graphs/Quality control/FeatureScatter_percent.mt_ncount_skin5.D002.1.sct.pdf")
FeatureScatter(skin5.D002.1.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")
dev.off()

skin5.D002.1.sct <- SCTransform(skin5.D002.1.sct, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)

#Quality ctrl, scale, normalize, skin6.D002.1.sct
skin6.D002.1.sct[["percent.mt"]] <- PercentageFeatureSet(skin6.D002.1.sct, pattern = "^MT-")
#skin6.D002.1.sct[["percent.ERY"]] <- PercentageFeatureSet(skin6.D002.1.sct, pattern = "HBB")
#skin6.D002.1.sct <- subset(skin6.D002.1.sct, subset = percent.ERY < 5)

pdf("Graphs/Quality control/VlnPlots_nfeature_ncout_percentmt_skin6.D002.1.sct.pdf")
VlnPlot(skin6.D002.1.sct, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
dev.off()

pdf("Graphs/Quality control/FeatureScatter_nfeature_ncount_skin6.D002.1.sct.pdf")
FeatureScatter(skin6.D002.1.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

pdf("Graphs/Quality control/FeatureScatter_percent.mt_ncount_skin6.D002.1.sct.pdf")
FeatureScatter(skin6.D002.1.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")
dev.off()

skin6.D002.1.sct <- SCTransform(skin6.D002.1.sct, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)

#Quality ctrl, scale, normalize, skin7.D002.1.sct
skin7.D002.1.sct[["percent.mt"]] <- PercentageFeatureSet(skin7.D002.1.sct, pattern = "^MT-")
#skin7.D002.1.sct[["percent.ERY"]] <- PercentageFeatureSet(skin7.D002.1.sct, pattern = "HBB")
#skin7.D002.1.sct <- subset(skin7.D002.1.sct, subset = percent.ERY < 5 )

pdf("Graphs/Quality control/VlnPlots_nfeature_ncout_skin7.D002.1.sct.pdf")
VlnPlot(skin7.D002.1.sct, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

pdf("Graphs/Quality control/FeatureScatter_nfeature_ncount_skin7.D002.1.sct.pdf")
FeatureScatter(skin7.D002.1.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

pdf("Graphs/Quality control/FeatureScatter_percent.mt_ncount_skin7.D002.1.sct.pdf")
FeatureScatter(skin7.D002.1.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")
dev.off()

skin7.D002.1.sct <- SCTransform(skin7.D002.1.sct, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)

###keloid###Quality ctrl, scale, normalize
#Quality ctrl, scale, normalize, keloid1.D002.1.sct
keloid1.D002.1.sct[["percent.mt"]] <- PercentageFeatureSet(keloid1.D002.1.sct, pattern = "^MT-")
#keloid1.D002.1.sct[["percent.ERY"]] <- PercentageFeatureSet(keloid1.D002.1.sct, pattern = "HBB")
#keloid1.D002.1.sct <- subset(keloid1.D002.1.sct, subset = percent.ERY < 5 )

pdf("Graphs/Quality control/VlnPlots_nfeature_ncout_keloid1.D002.1.sct.pdf")
VlnPlot(keloid1.D002.1.sct, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

pdf("Graphs/Quality control/FeatureScatter_nfeature_ncount_keloid1.di.D002.pdf")
FeatureScatter(keloid1.D002.1.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

pdf("Graphs/Quality control/FeatureScatter_percent.mt_ncount_keloid1.D002.1.sct.pdf")
FeatureScatter(keloid1.D002.1.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")
dev.off()

keloid1.D002.1.sct <- SCTransform(keloid1.D002.1.sct,method = "glmGamPoi", vars.to.regress = "percent.mt" ,verbose = F)


#Quality ctrl, scale, normalize, keloid2.D002.1.sct
keloid2.D002.1.sct[["percent.mt"]] <- PercentageFeatureSet(keloid2.D002.1.sct, pattern = "^MT-")
#keloid2.D002.1.sct[["percent.ERY"]] <- PercentageFeatureSet(keloid2.D002.1.sct, pattern = "HBB")
#keloid2.D002.1.sct <- subset(keloid2.D002.1.sct, subset = percent.ERY < 5 )

pdf("Graphs/Quality control/VlnPlots_nfeature_ncout_keloid2.D002.1.sct.pdf")
VlnPlot(keloid2.D002.1.sct, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

pdf("Graphs/Quality control/FeatureScatter_nfeature_ncount_keloid2.D002.1.sct.pdf")
FeatureScatter(keloid2.D002.1.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

pdf("Graphs/Quality control/FeatureScatter_percent.mt_ncount_keloid2.D002.1.sct.pdf")
FeatureScatter(keloid2.D002.1.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")
dev.off()

keloid2.D002.1.sct <- SCTransform(keloid2.D002.1.sct,method = "glmGamPoi", vars.to.regress = "percent.mt",verbose = F)


#Quality ctrl, scale, normalize, keloid3L.D002.1.sct
keloid3L.D002.1.sct[["percent.mt"]] <- PercentageFeatureSet(keloid3L.D002.1.sct, pattern = "^MT-")
#keloid3L.D002.1.sct[["percent.ERY"]] <- PercentageFeatureSet(keloid3L.D002.1.sct, pattern = "HBB")
#keloid3L.D002.1.sct <- subset(keloid3L.D002.1.sct, subset = percent.ERY < 5 )

pdf("Graphs/Quality control/VlnPlots_nfeature_ncout_keloid3L.D002.1.sct.pdf")
VlnPlot(keloid3L.D002.1.sct, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

pdf("Graphs/Quality control/FeatureScatter_nfeature_ncount_keloid3L.D002.1.sct.pdf")
FeatureScatter(keloid3L.D002.1.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

pdf("Graphs/Quality control/FeatureScatter_percent.mt_ncount_keloid3L.D002.1.sct.pdf")
FeatureScatter(keloid3L.D002.1.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")
dev.off()

keloid3L.D002.1.sct <- SCTransform(keloid3L.D002.1.sct, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)


#Quality ctrl, scale, normalize, keloid3R.D002.1.sct
keloid3R.D002.1.sct[["percent.mt"]] <- PercentageFeatureSet(keloid3R.D002.1.sct, pattern = "^MT-")
#keloid3R.D002.1.sct[["percent.ERY"]] <- PercentageFeatureSet(keloid3R.D002.1.sct, pattern = "HBB")
#keloid3R.D002.1.sct <- subset(keloid3R.D002.1.sct, subset = percent.ERY < 5 )

pdf("Graphs/Quality control/VlnPlots_nfeature_ncout_keloid3R.D002.1.sct.pdf")
VlnPlot(keloid3R.D002.1.sct, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

pdf("Graphs/Quality control/FeatureScatter_nfeature_ncount_keloid3R.D002.1.sct.pdf")
FeatureScatter(keloid3R.D002.1.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

pdf("Graphs/Quality control/FeatureScatter_percent.mt_ncount_keloid3R.D002.1.sct.pdf")
FeatureScatter(keloid3R.D002.1.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")
dev.off()

keloid3R.D002.1.sct <- SCTransform(keloid3R.D002.1.sct, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)


##### integrate skin and keloid ####
s.k.total.list<-list(keloid1.D002.1.sct,keloid2.D002.1.sct, keloid3L.D002.1.sct,keloid3R.D002.1.sct,skin1.D002.1.sct, skin2.D002.1.sct, skin3.D002.1.sct, skin4.D002.1.sct, skin5.D002.1.sct, skin6.D002.1.sct, skin7.D002.1.sct)#, skin8.D002.1.sct, skin9.D002.1.sct, skin10.D002.1.sct, skin11.D002.1.sct
s.k.total.features <- SelectIntegrationFeatures(object.list = s.k.total.list, nfeatures = 3000)
s.k.total.list <- PrepSCTIntegration(s.k.total.list, anchor.features = s.k.total.features)
s.k.total.list <- lapply(s.k.total.list, RunPCA, verbose = F, features= s.k.total.features)
s.k.total.anchors<-FindIntegrationAnchors(s.k.total.list,normalization.method = "SCT", anchor.features = s.k.total.features, reduction = "rpca")
s.k.total.x <- IntegrateData(anchorset=s.k.total.anchors, normalization.method = "SCT")
s.k.total <- s.k.total.x
s.k.total <- RunPCA(s.k.total, npcs = 60)
pdf("Graphs/Side/Elbowplot_skin_keloid_7s4k.pdf")
ElbowPlot(s.k.total, ndims = 60)
dev.off()
s.k.total <- RunUMAP(s.k.total, dims = 1:40)
s.k.total <- FindNeighbors(s.k.total, dims = 1:40)
#DefaultAssay(s.k.total)<- "integrated"
s.k.total <- FindClusters(s.k.total, resolution = 0.5)
UMAPPlot(s.k.total, label=T)
DefaultAssay(s.k.total)<- "RNA"
s.k.total <- NormalizeData(s.k.total)
FeaturePlot(s.k.total,c("PMEL","S100B","HBB","TOP2A"), split.by = "tissue")

FeaturePlot(s.k.total,c("FOS","FOSL1","FOSL2","JUN","JUNB","JUND","IL8"),split.by="origin")

#order conditions
s.k.total$tissue <- factor(x = s.k.total$tissue, levels = c("skin", "keloid"))
s.k.total$sample <- factor(x = s.k.total$sample, levels = c("skin_1","skin_2","skin_3", "skin_4", "skin_5", "skin_6", "skin_7", "keloid_1","keloid_2","keloid_3L","keloid_3R"))
s.k.total$origin <- factor(x = s.k.total$origin, levels = c("MD", "VV","TT"))

#UMAPPlots, unnamed clusters
#DefaultAssay(s.k.total)<-"integrated"
UMAPPlot(s.k.total, label=T)+ scale_color_manual(values=color_UMAP_short)+ ggtitle("Basic-UMAP")
UMAPPlot(s.k.total, label=F, split.by="tissue")#+ scale_color_manual(values=color_UMAP_short)+ ggtitle("Basic-UMAP")
UMAPPlot(s.k.total, label=F, split.by="sample", ncol = 3)+ scale_color_manual(values=color_UMAP_short)+ ggtitle("Basic-UMAP")
UMAPPlot(s.k.total, label=F, split.by="origin", ncol = 3)+ scale_color_manual(values=color_UMAP_short)+ ggtitle("Basic-UMAP")


#Characterization of Cluster

###StackedVlnPlot-Code###########################################################
## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}
######################################################################

#Clusteridentification
DotPlot(s.k.total, features=c("PDGFRA","LUM","COL1A1", "DCN","FBLN1","ACTA2","RGS5","KRT10","KRT1","KRT14","KRT5","SELE","PECAM1","VWF","LYVE1","CD3D","CD2","CXCR4","CD79A","MS4A1","CD68","AIF1","FCER1A","S100B","NGFR","PMEL","MLANA","HBB","HBA1","TOP2A","MKI67"), group.by = "seurat_clusters", assay = "RNA") + coord_flip()+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("Basic Cluster identification")

DefaultAssay(s.k.total)<- "RNA"
#Name Cluster
Idents(s.k.total)<-s.k.total$seurat_clusters
s.k.total<-RenameIdents(s.k.total,
                        `0`="FB", `1`="EC", `2`="FB", `3`="FB", `4`="SMC/PC", 
                        `5`="KC", `6`="KC", `7`="EC", `8`="KC", `9`="DC", `10`="TC",
                        `11`="LEC", `12`="SMC/PC", `13`="EC", `14`="MAC", `15`="ERY",
                        `16`="SC", `17`="FB", `18`="FB", `19`="MEL", `20`="FB")
s.k.total$celltype<-Idents(s.k.total)



#Define cluster levels
clusters_ordered<-c(
  "FB","SMC/PC","KC",
  "EC","LEC","TC", 
  "MAC","DC","SC","MEL","ERY")
s.k.total$celltype<- factor(s.k.total$celltype, levels = clusters_ordered)
Idents(s.k.total)<-s.k.total$celltype

#Marker Celltype specific
pdf("Graphs/Main/DotPlot_s.k.total._Clusteridentification.pdf")
DotPlot(s.k.total, features=c("PDGFRA","LUM","COL1A1", "DCN","FBLN1","ACTA2","RGS5","KRT10","KRT1","KRT14","KRT5","SELE","PECAM1","VWF","LYVE1","CD3D","CD2","CXCR4","CD68","AIF1","FCER1A","S100B","NGFR","PMEL","MLANA","HBB","HBA1"), group.by = "celltype", assay = "RNA") + coord_flip()+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("Basic Cluster identification")
dev.off()


pdf("Graphs/Main/DotPlot_s.k.total._SC_Clusteridentification.pdf")
DotPlot(s.k.total, features=c("MPZ","PLP1","PMP22","NCAM1","L1CAM","SCN7A","SOX10","PRX","MAG","DHH","FABP7","GAP43","EGR2"), assay = "RNA") + coord_flip()+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("Schwann cell Cluster identification")
dev.off()

# 27831564
pdf("Graphs/Main/DotPlot_s.k.total._Clusteridentification_extended_ 27831564.pdf")
DotPlot(s.k.total, features=c("SOX2","NCAM1","PROM1","VCAN","NGFR","FN1","PDGFRA","CDH7","B3GAT1","ITGB1"), assay = "RNA", group.by = "celltype.tissue") + coord_flip()+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("Progenitor identification")
dev.off()

DefaultAssay(s.k.total)<-"RNA"

F1<-FeaturePlot(s.k.total, "LUM", min.cutoff = 0)+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("LUM_Fibroblasts")
F2<-FeaturePlot(s.k.total, "DCN", min.cutoff = 0)+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("DCN_Fibroblasts")
F3<-FeaturePlot(s.k.total, "FBLN1", min.cutoff = 0)+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("FBLN1_Fibroblasts")
F4<-FeaturePlot(s.k.total, "ACTA2")+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("ACTA2_Smooth Muscle Cells")
F5<-FeaturePlot(s.k.total, "RGS5")+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("RGS5_Pericytes")
F6<-FeaturePlot(s.k.total, "KRT1")+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("KRT1_Keratinocytes")
F7<-FeaturePlot(s.k.total, "KRT14")+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("KRT14_Keratinocytes")
F8<-FeaturePlot(s.k.total, "ICAM1")+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("ICAM1_Endothelial Cells")
F9<-FeaturePlot(s.k.total, "LYVE1")+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("LYVE1_Lymphatic Endothelial Cells")
F10<-FeaturePlot(s.k.total, "CD3D")+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("CD3D_T Cells")
F11<-FeaturePlot(s.k.total, "CD68")+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("CD68_Macrophages")
F12<-FeaturePlot(s.k.total, "FCER1A")+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("FCER1A_Dendritic Cells")
F13<-FeaturePlot(s.k.total, "S100B")+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("S100B_Schwann Cells")
F14<-FeaturePlot(s.k.total, "MLANA")+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("MLANA_Melanocytes")
F15<-FeaturePlot(s.k.total, "HBB")+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("HBB_Erythrocytes")
pdf("Graphs/Main/FeaturePlots_s.k.total.Clusteridentification.pdf",width=15, height=9)
wrap_plots(F1,F2,F3,F4,F5,F6,F7,F8,F9,F10,F11,F12,F13,F14,F15, nrow = 3)
dev.off()

#UMAPPlots, named clusters
DefaultAssay(s.k.total)<-"integrated"
Idents(s.k.total)<-"celltype"

pdf("Graphs/Main/UMAPPlots_s.k.total.identified.pdf")
UMAPPlot(s.k.total, label=T, label.color = "black") + NoLegend()+ scale_color_manual(values=color_UMAP_short)+ ggtitle("Basic-UMAP")
dev.off()

pdf("Graphs/Main/UMAPPlots_s.k.total.identified.split.by.tissue.pdf" , height=5, width=10)
UMAPPlot(s.k.total, label=F, split.by="tissue")+ scale_color_manual(values=color_UMAP_short)+ ggtitle("Basic-UMAP")
dev.off()


pdf("Graphs/Side/UMAPPlots_s.k.total.identified.split.by.sample.pdf", height=20, width = 20)
UMAPPlot(s.k.total, label=F, split.by="sample", ncol = 3)+ scale_color_manual(values=color_UMAP_short)+ ggtitle("Basic-UMAP")
dev.off()

pdf("Graphs/Side/UMAPPlots_s.k.total.identified.split.by.origin.pdf", height=10, width = 20)
UMAPPlot(s.k.total, label=F, split.by="origin", ncol = 3)+ scale_color_manual(values=color_UMAP_short)+ ggtitle("Basic-UMAP")
dev.off()

pdf("Graphs/Side/UMAPPlots_s.k.total.identified.split.by.celltype.pdf", height=10, width = 20)
UMAPPlot(s.k.total, label=F, split.by="celltype", ncol = 3)+ scale_color_manual(values=color_UMAP_short)+ ggtitle("Basic-UMAP")
dev.off()

#Clustertree
s.k.total<-BuildClusterTree(s.k.total, dims=1:40, assay = "RNA")
DefaultAssay(s.k.total)<-"integrated"
?ape::plot.phylo
pdf("Graphs/Main/Clustertree_s.k.total.pdf", width = 10, height=4)
PlotClusterTree(s.k.total, type="phylo", tip.color=color_UMAP_short ,node.pos=1 ,show.node.label=FALSE, cex=1.5 ,font=4,no.margin=T, rotate.tree=180, cex.main=2)
dev.off()



#Clustermarker + Heatmap total celltype
Idents(s.k.total)<-"celltype"
s.k.total.clustermarker<-FindAllMarkers(s.k.total, assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
s.k.total.clustermarker$Foldchange_UP <- 2^(s.k.total.clustermarker$avg_log2FC)
s.k.total.clustermarker$Foldchange_DOWN <- 2^(-s.k.total.clustermarker$avg_log2FC)
s.k.total.clustermarker$Ratio_pct1_pct2 <- (s.k.total.clustermarker$pct.1)/(s.k.total.clustermarker$pct.2)
write.xlsx(s.k.total.clustermarker, "Lists/Clustermarker_7s.4k_total.xlsx")


#Heatmap Top up and down regulated genes 
top10log.sc <- s.k.total.clustermarker %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
pdf("Graphs/Main/Heatmap_clustermarker_skin_scar_total_avg_logfc.pdf", width= 20, height = 8)
DoHeatmap(s.k.total, assay= "SCT", cells = as.character(WhichCells(s.k.total))[c(T,F)], features = top10log.sc$gene, group.colors = color_UMAP_short)+ scale_color_manual(values=color_UMAP_short)+NoLegend() + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y = element_blank())
dev.off()

#DotPlot Top up and down regulated genes 
s.k.total.clustermarker_subset_UP <- subset(s.k.total.clustermarker, cluster == "SC" & Foldchange_UP >=2 )
top50_UP_scsubs_clustermarker <- s.k.total.clustermarker_subset_UP %>% group_by(cluster) %>% top_n(n=50,wt=Foldchange_UP)
top50_UP_scsubs_clustermarker<-arrange(top50_UP_scsubs_clustermarker,Foldchange_UP)
top50_UP_scsubs_clustermarker <- as.character(top50_UP_scsubs_clustermarker$gene)
top50_UP_scsubs_clustermarker<-top50_UP_scsubs_clustermarker[!duplicated(top50_UP_scsubs_clustermarker)]


pdf("Graphs/Main/DotPlot_Top_50UP_Clustermarker_total_FCmin2.pdf", width = 7, height = 12)
DotPlot(s.k.total, features= top50_UP_scsubs_clustermarker, group.by = "celltype", assay = "RNA") + coord_flip()+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 45)) + ggtitle("Basic Cluster_Top50-SC Clustermarker_FCmin2") + scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()


s.k.total.clustermarker_subset_DOWN <- subset(s.k.total.clustermarker, cluster == "SC" & Foldchange_DOWN >=2 )
top50_DOWN_scsubs_clustermarker <- s.k.total.clustermarker_subset_DOWN %>% group_by(cluster) %>% top_n(n=50,wt=Foldchange_DOWN)
top50_DOWN_scsubs_clustermarker<-arrange(top50_DOWN_scsubs_clustermarker,Foldchange_DOWN)
top50_DOWN_scsubs_clustermarker <- as.character(top50_DOWN_scsubs_clustermarker$gene)
top50_DOWN_scsubs_clustermarker<-top50_DOWN_scsubs_clustermarker[!duplicated(top50_DOWN_scsubs_clustermarker)]

pdf("Graphs/Main/DotPlot_Top_50DOWN_Clustermarker_total_FCmin2_SC.pdf", width = 7, height = 12)
DotPlot(s.k.total, features=top50_DOWN_scsubs_clustermarker, group.by = "celltype", assay = "RNA") + coord_flip()+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 45)) + ggtitle("Basic Cluster_Top50-SC Clustermarker_FCmin2") + scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()




###Pieplot/Barplot cellnumber /tissue/cluster
#subset conditions
DefaultAssay(s.k.total)<-"RNA"
skin.SC_7s4k.subset<-subset(s.k.total, celltype=="SC")

#Data frame frequencies
skin.SC_7s4k.cluster.frequencies<-as.data.frame(table(skin.SC_7s4k.subset$tissue))


pdf("Graphs/Side/PiePlot_SC_7s4k_clusterfrequencies.pdf", width = 10, height = 10)
pie(skin.SC_7s4k.cluster.frequencies$Freq, labels = skin.SC_7s4k.cluster.frequencies$Var1, main="Cell source_schwann cell cluster",  clockwise = T, col= color_tissue)
dev.off()

#Barplot Cluster percentage
skin.SC_7s4k.cluster.frequencies$celltype = "SC"
skin.SC_7s4k.cluster.frequencies<-mutate(skin.SC_7s4k.cluster.frequencies,Percentage= skin.SC_7s4k.cluster.frequencies$Freq/sum(skin.SC_7s4k.cluster.frequencies$Freq)*100)

SC_7s4k.cluster.df <- skin.SC_7s4k.cluster.frequencies
SC_7s4k.cluster.df$tissue <- factor(x = SC_7s4k.cluster.df$celltype, levels = c("SC"))

pdf("Graphs/Main/Barplot_skin_scars_SC_cellfrequencies.pdf", width = 4, height = 6)
ggplot(data=SC_7s4k.cluster.df, aes(x=Var1, y=Freq, fill=Var1)) + ylim(0,1000)+
  geom_bar(stat="identity", position="dodge") +  theme_classic() +
  geom_text_repel(aes(label = Freq, y=Freq+0.5),position=position_dodge(0.9), vjust=0.6, color = 'black', size = 3.5)+ggtitle("Cell source_schwann cell cluster") + xlab("condition")+ylab("total cell count") +scale_fill_manual(values = color_tissue)
dev.off()



#Percentage total
###Pieplot/Barplot cellnumber /tissue/cluster
#subset conditions
skin.subset<-subset(s.k.total, tissue=="skin")
keloid.subset<-subset(s.k.total, tissue=="keloid")

#Data frame frequencies
skin.cluster.frequencies<-as.data.frame(table(skin.subset$celltype))
keloid.cluster.frequencies<-as.data.frame(table(keloid.subset$celltype))

#Barplot Cluster percentage
skin.cluster.frequencies$tissue="skin"
keloid.cluster.frequencies$tissue="keloid"

skin.cluster.frequencies<-mutate(skin.cluster.frequencies,Percentage= skin.cluster.frequencies$Freq/sum(skin.cluster.frequencies$Freq)*100)
keloid.cluster.frequencies<-mutate(keloid.cluster.frequencies,Percentage= keloid.cluster.frequencies$Freq/sum(keloid.cluster.frequencies$Freq)*100)



skin_kel.cluster.df <- rbind(skin.cluster.frequencies,keloid.cluster.frequencies)
skin_kel.cluster.df$tissue <- factor(x = skin_kel.cluster.df$tissue, levels = c("skin", "keloid"))

pdf("Graphs/Main/Barplot_skin_kel_cellfrequencie.pdf", width = 3, height = 5)
ggplot(data=skin_kel.cluster.df, aes(x=tissue, y=Percentage, fill=Var1)) +
  geom_bar(stat="identity", position="stack")+scale_fill_manual(values =color_UMAP_short)+ theme_classic()
dev.off()

pdf("Graphs/Side/PiePlot__skin_kel_cellfrequencie.pdf", width = 10, height = 10)
pie(skin.cluster.frequencies$Freq, labels = skin.cluster.frequencies$Var1, main="skin",  clockwise = T, col = c(color_UMAP_short))
pie(keloid.cluster.frequencies$Freq, labels = keloid.cluster.frequencies$Var1, main="Keloid",  clockwise = T, col = c(color_UMAP_short))
dev.off()


#Percentage total
###Pieplot/Barplot cellnumber /sample/cluster
#subset conditions
skin1.subset<-subset(s.k.total, sample=="skin_1")
skin2.subset<-subset(s.k.total, sample=="skin_2")
skin3.subset<-subset(s.k.total, sample=="skin_3")
skin4.subset<-subset(s.k.total, sample=="skin_4")
skin5.subset<-subset(s.k.total, sample=="skin_5")
skin6.subset<-subset(s.k.total, sample=="skin_6")
skin7.subset<-subset(s.k.total, sample=="skin_7")
keloid1.subset<-subset(s.k.total, sample=="keloid_1")
keloid2.subset<-subset(s.k.total, sample=="keloid_2")
keloid3L.subset<-subset(s.k.total, sample=="keloid_3L")
keloid3R.subset<-subset(s.k.total, sample=="keloid_3R")

#Data frame frequencies
skin1.cluster.frequencies<-as.data.frame(table(skin1.subset$celltype))
skin2.cluster.frequencies<-as.data.frame(table(skin2.subset$celltype))
skin3.cluster.frequencies<-as.data.frame(table(skin3.subset$celltype))
skin4.cluster.frequencies<-as.data.frame(table(skin4.subset$celltype))
skin5.cluster.frequencies<-as.data.frame(table(skin5.subset$celltype))
skin6.cluster.frequencies<-as.data.frame(table(skin6.subset$celltype))
skin7.cluster.frequencies<-as.data.frame(table(skin7.subset$celltype))
keloid1.cluster.frequencies<-as.data.frame(table(keloid1.subset$celltype))
keloid2.cluster.frequencies<-as.data.frame(table(keloid2.subset$celltype))
keloid3L.cluster.frequencies<-as.data.frame(table(keloid3L.subset$celltype))
keloid3R.cluster.frequencies<-as.data.frame(table(keloid3R.subset$celltype))

#Barplot Cluster percentage
skin1.cluster.frequencies$sample="skin_1"
skin2.cluster.frequencies$sample="skin_2"
skin3.cluster.frequencies$sample="skin_3"
skin4.cluster.frequencies$sample="skin_4"
skin5.cluster.frequencies$sample="skin_5"
skin6.cluster.frequencies$sample="skin_6"
skin7.cluster.frequencies$sample="skin_7"
keloid1.cluster.frequencies$sample="keloid_1"
keloid2.cluster.frequencies$sample="keloid_2"
keloid3L.cluster.frequencies$sample="keloid_3L"
keloid3R.cluster.frequencies$sample="keloid_3R"

#Barplot Cluster percentage
skin1.cluster.frequencies$tissue="skin"
skin2.cluster.frequencies$tissue="skin"
skin3.cluster.frequencies$tissue="skin"
skin4.cluster.frequencies$tissue="skin"
skin5.cluster.frequencies$tissue="skin"
skin6.cluster.frequencies$tissue="skin"
skin7.cluster.frequencies$tissue="skin"
keloid1.cluster.frequencies$tissue="keloid"
keloid2.cluster.frequencies$tissue="keloid"
keloid3L.cluster.frequencies$tissue="keloid"
keloid3R.cluster.frequencies$tissue="keloid"

skin1.cluster.frequencies<-mutate(skin1.cluster.frequencies,Percentage= skin1.cluster.frequencies$Freq/sum(skin1.cluster.frequencies$Freq)*100)
skin2.cluster.frequencies<-mutate(skin2.cluster.frequencies,Percentage= skin2.cluster.frequencies$Freq/sum(skin2.cluster.frequencies$Freq)*100)
skin3.cluster.frequencies<-mutate(skin3.cluster.frequencies,Percentage= skin3.cluster.frequencies$Freq/sum(skin3.cluster.frequencies$Freq)*100)
skin4.cluster.frequencies<-mutate(skin4.cluster.frequencies,Percentage= skin4.cluster.frequencies$Freq/sum(skin4.cluster.frequencies$Freq)*100)
skin5.cluster.frequencies<-mutate(skin5.cluster.frequencies,Percentage= skin5.cluster.frequencies$Freq/sum(skin5.cluster.frequencies$Freq)*100)
skin6.cluster.frequencies<-mutate(skin6.cluster.frequencies,Percentage= skin6.cluster.frequencies$Freq/sum(skin6.cluster.frequencies$Freq)*100)
skin7.cluster.frequencies<-mutate(skin7.cluster.frequencies,Percentage= skin7.cluster.frequencies$Freq/sum(skin7.cluster.frequencies$Freq)*100)
keloid1.cluster.frequencies<-mutate(keloid1.cluster.frequencies,Percentage= keloid1.cluster.frequencies$Freq/sum(keloid1.cluster.frequencies$Freq)*100)
keloid2.cluster.frequencies<-mutate(keloid2.cluster.frequencies,Percentage= keloid2.cluster.frequencies$Freq/sum(keloid2.cluster.frequencies$Freq)*100)
keloid3L.cluster.frequencies<-mutate(keloid3L.cluster.frequencies,Percentage= keloid3L.cluster.frequencies$Freq/sum(keloid3L.cluster.frequencies$Freq)*100)
keloid3R.cluster.frequencies<-mutate(keloid3R.cluster.frequencies,Percentage= keloid3R.cluster.frequencies$Freq/sum(keloid3R.cluster.frequencies$Freq)*100)


skin_kel_samp.cluster.df <- rbind(skin1.cluster.frequencies,skin2.cluster.frequencies,skin3.cluster.frequencies,skin4.cluster.frequencies,skin5.cluster.frequencies,skin6.cluster.frequencies,skin7.cluster.frequencies,keloid1.cluster.frequencies,keloid2.cluster.frequencies,keloid3L.cluster.frequencies,keloid3R.cluster.frequencies)
skin_kel_samp.cluster.df$sample <- factor(x = skin_kel_samp.cluster.df$sample, levels = c("skin_1","skin_2","skin_3","skin_4","skin_5","skin_6","skin_7", "keloid_1", "keloid_2", "keloid_3L", "keloid_3R"))
skin_kel_samp.cluster.df$tissue <- factor(x = skin_kel_samp.cluster.df$tissue, levels = c("skin", "keloid"))

pdf("Graphs/Side/Barplot_skin_kel_cell_clusterfrequencie.pdf", width = 10, height = 5)
ggplot(data=skin_kel_samp.cluster.df, aes(x=sample, y=Percentage, fill=Var1)) +
  geom_bar(stat="identity", position="stack")+scale_fill_manual(values =color_UMAP_short)+ theme_classic()
dev.off()

pdf("Graphs/Side/Barplot_skin_scars_cellfrequencies_celltype.pdf", width = 6, height = 6)
ggplot(data=skin_kel_samp.cluster.df, aes(x=sample, y=Freq, fill=Var1)) +
  geom_bar(stat="identity", position="stack")+scale_fill_manual(values =color_UMAP_short) +  theme_classic() +ggtitle("Total Cell distribution") + xlab("condition")+ylab("total cell count")
dev.off()

pdf("Graphs/Side/PiePlot__skin_kel_cell_clusterfrequencie.pdf", width = 10, height = 10)
pie(skin1.cluster.frequencies$Freq, labels = skin1.cluster.frequencies$Var1, main="skin_1",  clockwise = T, col = c(color_UMAP_short))
pie(skin2.cluster.frequencies$Freq, labels = skin2.cluster.frequencies$Var1, main="skin_2",  clockwise = T, col = c(color_UMAP_short))
pie(skin3.cluster.frequencies$Freq, labels = skin3.cluster.frequencies$Var1, main="skin_3",  clockwise = T, col = c(color_UMAP_short))
pie(skin4.cluster.frequencies$Freq, labels = skin4.cluster.frequencies$Var1, main="skin_4",  clockwise = T, col = c(color_UMAP_short))
pie(skin5.cluster.frequencies$Freq, labels = skin5.cluster.frequencies$Var1, main="skin_5",  clockwise = T, col = c(color_UMAP_short))
pie(skin6.cluster.frequencies$Freq, labels = skin6.cluster.frequencies$Var1, main="skin_6",  clockwise = T, col = c(color_UMAP_short))
pie(skin7.cluster.frequencies$Freq, labels = skin7.cluster.frequencies$Var1, main="skin_7",  clockwise = T, col = c(color_UMAP_short))
pie(keloid1.cluster.frequencies$Freq, labels = keloid1.cluster.frequencies$Var1, main="Keloid_1",  clockwise = T, col = c(color_UMAP_short))
pie(keloid2.cluster.frequencies$Freq, labels = keloid2.cluster.frequencies$Var1, main="Keloid_2",  clockwise = T, col = c(color_UMAP_short))
pie(keloid3L.cluster.frequencies$Freq, labels = keloid3L.cluster.frequencies$Var1, main="Keloid_3L",  clockwise = T, col = c(color_UMAP_short))
pie(keloid3R.cluster.frequencies$Freq, labels = keloid3R.cluster.frequencies$Var1, main="Keloid_3R",  clockwise = T, col = c(color_UMAP_short))
dev.off()

#Barplot percentage mean all donors_ self calculate using skin_kel_samp.cluster.df
skin_kel.cluster.percentage.mean <- read.xlsx("D:/Direder/Projekte/D002_10xScarWars/D002_Berechnungen/D002_c12_7s4k/Lists/Percentage_mean.xlsx",1,header = T)
skin_kel.cluster.percentage.mean$tissue <- factor(x = skin_kel.cluster.percentage.mean$tissue, levels = c("skin", "keloid"))
skin_kel.cluster.percentage.mean$Var1 <- factor(x = skin_kel.cluster.percentage.mean$Var1, levels = clusters_ordered)

pdf("Graphs/Main/Barplot_skin_kel_cellfrequencie_mean_percentage.pdf", width = 3, height = 5)
ggplot(data=skin_kel.cluster.percentage.mean, aes(x=tissue, y=Percentage, fill=Var1)) +
  geom_bar(stat="identity", position="stack")+scale_fill_manual(values =color_UMAP_short)+ theme_classic()
dev.off()


#Boxplots Percentage - t-test
pdf("Graphs/Side/BoxPlot_Percentage_comparison.pdf", width = 10, height= 10)
ggplot(skin_kel_samp.cluster.df, aes(x=tissue, y=Percentage, fill=tissue)) +geom_boxplot()+ facet_wrap(~Var1) + theme_classic() +
  labs(title="Percentage comparison",x="", y = "Percentage")+ scale_fill_manual(values=color_tissue)+ stat_compare_means(method = "t.test",method.args = list(var.equal = TRUE),label =  "p.signif", label.x = 1.5)
dev.off()
write.xlsx(skin_kel_samp.cluster.df, "Lists/D002_c12_cellular_Percentage.xlsx")


SC_percentage<-subset(skin_kel_samp.cluster.df,skin_kel_samp.cluster.df$Var1 == "SC")
pdf("Graphs/Side/BoxPlot_SC-Percentage_comparison.pdf", width = 2, height= 5)
ggplot(SC_percentage, aes(x=tissue, y=Percentage, fill=tissue)) +geom_boxplot() + theme_classic() +
  labs(title="SC-Percentage comparison",x="", y = "Percentage")+ ylim(0,3) +scale_fill_manual(values=color_tissue)# + geom_signif(comparisons = list(c("skin", "keloid")))#, test="t.test", test.args = list(var.equal = TRUE),map_signif_level=TRUE
dev.off()


#cellular source general
#Barplot Cluster percentage
skin.keloid.cluster.frequencies<-as.data.frame(table(s.k.total$sample))
skin.keloid.cluster.frequencies$celltype = "SC"
skin.keloid.cluster.frequencies<-mutate(skin.keloid.cluster.frequencies,Percentage= skin.keloid.cluster.frequencies$Freq/sum(skin.keloid.cluster.frequencies$Freq)*100)

pdf("Graphs/Side/Barplot_skin_scars_cellfrequencies.pdf", width = 8, height = 8)
ggplot(data=skin.keloid.cluster.frequencies, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", position="dodge") +  theme_classic() +
  geom_text_repel(aes(label = Freq, y=Freq+0.5),position=position_dodge(0.9), vjust=0.6, color = 'black', size = 3.5)+ggtitle("Cell source_all cells") + xlab("condition")+ylab("total cell count")
dev.off()

#cellular source SC
SC_subset<-subset(s.k.total, celltype=="SC")
skin.keloid.SC.frequencies<-as.data.frame(table(SC_subset$sample))
skin.keloid.SC.frequencies$celltype = "SC"
skin.keloid.SC.frequencies<-mutate(skin.keloid.SC.frequencies,Percentage= skin.keloid.SC.frequencies$Freq/sum(skin.keloid.SC.frequencies$Freq)*100)

pdf("Graphs/Side/Barplot_skin_scars_SC_frequencies.pdf", width = 8, height = 8)
ggplot(data=skin.keloid.SC.frequencies, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", position="dodge") +  theme_classic() +
  geom_text_repel(aes(label = Freq, y=Freq+0.5),position=position_dodge(0.9), vjust=0.6, color = 'black', size = 3.5)+ggtitle("Schwann Cell source_all cells") + xlab("condition")+ylab("total cell count")
dev.off()


#add identity to analyze tissue/health diffgenes
#Set tissue in identy 
Idents(s.k.total) <- s.k.total$celltype
s.k.total$celltype.tissue <- paste(Idents(s.k.total), s.k.total$tissue, sep = "_")
Idents(s.k.total)<-s.k.total$celltype.tissue
clusters_tissue_ordered<-c("FB_skin" ,      "FB_keloid",
                           "SMC/PC_skin",   "SMC/PC_keloid",
                           "KC_skin",       "KC_keloid",
                           "EC_skin",       "EC_keloid",
                           "LEC_skin",      "LEC_keloid",
                           "TC_skin",       "TC_keloid",
                           "MAC_skin",      "MAC_keloid",
                           "DC_skin",       "DC_keloid",
                           "SC_skin",       "SC_keloid", 
                           "MEL_skin",      "MEL_keloid",
                           "ERY_skin",      "ERY_keloid")
s.k.total$celltype.tissue <- factor(s.k.total$celltype.tissue, levels = clusters_tissue_ordered)
Idents(s.k.total)<-s.k.total$celltype.tissue

#Clustermarker + Heatmap total celltype_tissue
Idents(s.k.total)<-"celltype.tissue"
s.k.total.clustermarker<-FindAllMarkers(s.k.total, assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
s.k.total.clustermarker$Foldchange_UP <- 2^(s.k.total.clustermarker$avg_log2FC)
s.k.total.clustermarker$Foldchange_DOWN <- 2^(-s.k.total.clustermarker$avg_log2FC)
s.k.total.clustermarker$Ratio_pct1_pct2 <- (s.k.total.clustermarker$pct.1)/(s.k.total.clustermarker$pct.2)
write.xlsx(s.k.total.clustermarker, "Lists/Clustermarker_7s.4k_c_t.xlsx")
#Idents(s.k.total)<-"celltype"

DEG_SC_skin_vs_AllnonSC<-FindMarkers(s.k.total, ident.1="SC_skin", ident.2=c("EC_keloid","EC_skin","ERY_keloid","ERY_skin","FB_keloid","FB_skin","KC_keloid","KC_skin","LEC_keloid","LEC_skin","MAC_keloid","MAC_skin","DC_keloid","DC_skin","MEL_keloid","MEL_skin","SMC/PC_keloid","SMC/PC_skin","TC_keloid","TC_skin"), assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
DEG_SC_skin_vs_AllnonSC$Foldchange_UP <- 2^(DEG_SC_skin_vs_AllnonSC$avg_log2FC)
DEG_SC_skin_vs_AllnonSC$Foldchange_DOWN <- 2^(-DEG_SC_skin_vs_AllnonSC$avg_log2FC)
DEG_SC_skin_vs_AllnonSC$Ratio_pct1_pct2 <- (DEG_SC_skin_vs_AllnonSC$pct.1)/(DEG_SC_skin_vs_AllnonSC$pct.2)
DEG_SC_skin_vs_AllnonSC$comparison <- "SC-skin_vs_nonSC"
write.xlsx(DEG_SC_skin_vs_AllnonSC, "Lists/DEG_SC_skin_vs_AllnonSC.xlsx")

DEG_SC_keloid_vs_AllnonSC<-FindMarkers(s.k.total, ident.1="SC_keloid", ident.2=c("EC_keloid","EC_skin","ERY_keloid","ERY_skin","FB_keloid","FB_skin","KC_keloid","KC_skin","LEC_keloid","LEC_skin","MAC_keloid","MAC_skin","DC_keloid","DC_skin","MEL_keloid","MEL_skin","SMC/PC_keloid","SMC/PC_skin","TC_keloid","TC_skin"), assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
DEG_SC_keloid_vs_AllnonSC$Foldchange_UP <- 2^(DEG_SC_keloid_vs_AllnonSC$avg_log2FC)
DEG_SC_keloid_vs_AllnonSC$Foldchange_DOWN <- 2^(-DEG_SC_keloid_vs_AllnonSC$avg_log2FC)
DEG_SC_keloid_vs_AllnonSC$Ratio_pct1_pct2 <- (DEG_SC_keloid_vs_AllnonSC$pct.1)/(DEG_SC_keloid_vs_AllnonSC$pct.2)
DEG_SC_keloid_vs_AllnonSC$comparison <- "SC-keloid_vs_nonSC"
write.xlsx(DEG_SC_keloid_vs_AllnonSC, "Lists/DEG_SC_keloid_vs_AllnonSC.xlsx")

#SC_tissue_comparison
SC_tissue_comp_df <- rbind(DEG_SC_skin_vs_AllnonSC,DEG_SC_keloid_vs_AllnonSC)
SC_tissue_comp_df_subset <- subset(SC_tissue_comp_df, Foldchange_UP >=2 | Foldchange_DOWN >=2)
SC_tissue_comp_df_subset$Direction <- ifelse (SC_tissue_comp_df_subset$Foldchange_UP>=2,"UP","DOWN")
SC_tissue_comp_df_subset<- select(SC_tissue_comp_df_subset, Direction, comparison)
row.names(SC_tissue_comp_df_subset)<-NULL
SC_tissue_comp_df_subset$cluster<- as.factor(SC_tissue_comp_df_subset$comparison)
SC_tissue_comp_df_subset$cluster<- as.character(SC_tissue_comp_df_subset$Direction)
SC_tissue_comp_df_subset<- dplyr::count(SC_tissue_comp_df_subset, comparison,Direction) %>% ungroup()
SC_tissue_comp_df_subset$Direction <- factor(SC_tissue_comp_df_subset$Direction, levels= c("UP","DOWN"))
SC_tissue_comp_df_subset$comparison <- factor(SC_tissue_comp_df_subset$comparison, levels= c("SC-skin_vs_nonSC","SC-keloid_vs_nonSC"))
pdf("Graphs/Main/Barplot_nDEG_SC_tissue_comparison.pdf")
ggplot(data=SC_tissue_comp_df_subset, aes(x=comparison, y=n, fill=Direction)) +
  geom_bar(stat="identity", position=position_dodge())+ theme_classic() +
  geom_text_repel(aes(label = n, y=n+0.5),position=position_dodge(0.9), vjust=0.6, color = 'black', size = 4) + ggtitle("DEG_SC-tissue comparison")  +scale_fill_manual(values =color_DEG) + xlab("SC-type")+ylab("total number of DEG")
dev.off()


############################
####Subset Schwann Cells####
############################
Idents(s.k.total)<- s.k.total$celltype
SC_7s4k <- subset(s.k.total,idents = "SC")
DefaultAssay(SC_7s4k)<-"RNA"
SC_7s4k[["percent.mt"]] <- PercentageFeatureSet(SC_7s4k, pattern = "^MT-")
SC_7s4k <- SCTransform(SC_7s4k, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
SC_7s4k <- RunPCA(SC_7s4k, npcs = 50)
ElbowPlot(SC_7s4k, ndims = 25)
SC_7s4k <- RunUMAP(SC_7s4k, dims = 1:14)
SC_7s4k <- FindNeighbors(SC_7s4k, dims = 1:14)
#DefaultAssay(SC_7s4k)<- "SCT"
SC_7s4k <- FindClusters(SC_7s4k, resolution = 0.6)
UMAPPlot(SC_7s4k, label=T, split.by="tissue")
DefaultAssay(SC_7s4k)<- "RNA"
SC_7s4k <- NormalizeData(SC_7s4k)
FeaturePlot(SC_7s4k,c("DCN","TOP2A","SELE","MBP","PLP1","S100B"))

#SC identification
pdf("Graphs/Main/Featureplot_Clusteridentification_SC7s4k.pdf", width = 10, height = 10)
FeaturePlot(SC_7s4k,c("S100B","PLP1","MPZ","NES","IGFBP5","TNFAIP6","MKI67","TOP2A","DCN","LUM","SELE","ICAM1"),ncol = 3, cols = c("lightgrey",single_col), pt.size = 0.5)
dev.off()

pdf("Graphs/Main/Featureplot_Clusteridentification_SC7s4k_R1.pdf", width = 10, height = 10)
FeaturePlot(SC_7s4k,c("S100B","MBP","SCN7A","NES","IGFBP5","TNFAIP6","MKI67","TOP2A","DCN","LUM","SELE","ICAM1"),ncol = 3, cols = c("lightgrey",single_col), pt.size = 0.5, order=T)
dev.off()

pdf("Graphs/Main/Dotplot_Clusteridentification_SC7s4k_R1.pdf")
DotPlot(SC_7s4k,features=c("MBP","PLP1","PMP22","MPZ","NCAM1","L1CAM","SCN7A","POU3F1","GFAP","FOXO4","SOX10","ETV5","ERBB3","ITGA4","TFAP2A","CDH2","NT5E","THY1","ENG"     ,"NES","PAX3","SNAI1","SNAI2","MSX1","TWIST1","TWIST2","MSI1","ATXN1","SOX9"),assay="RNA") + coord_flip()+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 45)) + ggtitle("Mature SC-characterisation") + scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()

pdf("Graphs/Main/Dotplot_SC_Clusteridentification_s.k.total_R2.pdf")
DotPlot(s.k.total,features=c("S100B","NGFR","PLP1","PMP22","MPZ","NCAM1","L1CAM","SCN7A","POU3F1","SOX10","PRX","MAG","EGR2","MLANA","PMEL"     ,"NT5E","THY1","ENG"     ,"NES","PAX3","SNAI1","SNAI2","MSX1","TWIST1","TWIST2","MSI1","ATXN1","SOX9"),assay="RNA") + coord_flip()+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 45)) + ggtitle("SC-Cluster Identification") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(SC_7s4k,features=c("S100B","NGFR","MBP","PLP1","PMP22","MPZ","NCAM1","L1CAM","SCN7A","POU3F1","GFAP","FOXO4","SOX10","ETV5","ERBB3","ITGA4","TFAP2A","CDH2", "PLP1","PMP22","MPZ","NCAM1","L1CAM","SCN7A","POU3F1","SOX10","PRX","MAG","EGR2","MLANA","PMEL"     ,"NT5E","THY1","ENG"     ,"NES","PAX3","SNAI1","SNAI2","MSX1","TWIST1","TWIST2","MSI1","ATXN1","SOX9"),assay="RNA") + coord_flip()+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 45)) + ggtitle("SC-Cluster Identification") + scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()

pdf("Graphs/Main/Featureplot_blend_NGFR_100B_SC7s4k.pdf", width= 10, height = 4)
FeaturePlot(SC_7s4k,c("NGFR","S100B"),blend = T, blend.threshold = 0.5, cols = c("lightgrey","#990F26","#54990F"), pt.size = 1.5, order=T, min.cutoff = 0, max.cutoff = 2)
dev.off()


B1<-FeaturePlot(SC_7s4k,c("CDH19","S100B"),blend = T, blend.threshold = 0.5, cols = c("lightgrey","#990F26","#54990F"), pt.size = 1.5, order=T, min.cutoff = 0, max.cutoff = 3)
B2<-FeaturePlot(SC_7s4k,c("SELE","S100B"),blend = T, blend.threshold = 0.5, cols = c("lightgrey","#990F26","#54990F"), pt.size = 1.5, order=T, min.cutoff = 0, max.cutoff = 3)
B3<-FeaturePlot(SC_7s4k,c("THY1","S100B"),blend = T, blend.threshold = 0.5, cols = c("lightgrey","#990F26","#54990F"), pt.size = 1.5, order=T, min.cutoff = 0, max.cutoff = 3)
B4<-FeaturePlot(SC_7s4k,c("MKI67","S100B"),blend = T, blend.threshold = 0.5, cols = c("lightgrey","#990F26","#54990F"), pt.size = 1.5, order=T, min.cutoff = 0, max.cutoff = 3)

pdf("Graphs/Main/Blendmix_SC7s4k.pdf", width= 10, height = 16)
wrap_plots(B2,B3,B4,B1,nrow = 4)
dev.off()


pdf("Graphs/Main/Featureplot_CCL2_AP1_SC7s4k.pdf", width = 12, height = 10)
FeaturePlot(SC_7s4k,c("CCL2","MMP9","JUN","JUNB","JUND","FOS","FOSL1","FOSL2"), cols = c("lightgrey",single_col), pt.size = 0.5)
dev.off()


#UMAPPlots, unnamed clusters
DefaultAssay(SC_7s4k)<-"RNA"
pdf("Graphs/Side/UMAPPlots_SC_7s4k.labeled.pdf")
UMAPPlot(SC_7s4k, label=T, label.color="black")+ scale_color_manual(values=color_SC)
UMAPPlot(SC_7s4k, label=T, split.by="tissue")+ scale_color_manual(values=color_SC) + NoLegend()   
UMAPPlot(SC_7s4k, label=F, split.by="sample")+ scale_color_manual(values=color_SC)
dev.off()


pdf("Graphs/Side/UMAP_Feature_Revision1_SC_7s4k.identified.split.by.origin.pdf")
FeaturePlot(SC_7s4k,c("MBP","SCN7A"), split.by = "tissue", order=T, cols = c("lightgrey",single_col), pt.size = 0.5)
dev.off()

DotPlot(SC_7s4k,features=c("DCN","TOP2A","SELE","MPZ","PLP1","S100B"))

#Name Cluster
Idents(SC_7s4k)<-SC_7s4k$seurat_clusters
SC_7s4k<-RenameIdents(SC_7s4k,
                      `0`="SC-Repair1", `1`="SC-Repair2", `2`="SC-Repair3", `3`="SC-FB",`4`="SC-Skin",`5`="SC-Promyel",`6`="SC-Prolif",`7`="SC-EC")
SC_7s4k$celltype<-Idents(SC_7s4k)


#Define cluster levels
clusters_ordered<-c("SC-Skin","SC-Promyel","SC-Repair1","SC-Repair2","SC-Repair3","SC-Prolif",
                    "SC-EC","SC-FB")
SC_7s4k$celltype<- factor(SC_7s4k$celltype, levels = clusters_ordered)
Idents(SC_7s4k)<-SC_7s4k$celltype


#Clustermarker total
SC_7s4k.clustermarker<-FindAllMarkers(SC_7s4k,  assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
SC_7s4k.clustermarker$Foldchange_UP <- 2^(SC_7s4k.clustermarker$avg_log2FC)
SC_7s4k.clustermarker$Foldchange_DOWN <- 2^(-SC_7s4k.clustermarker$avg_log2FC)
SC_7s4k.clustermarker$Ratio_pct1_pct2 <- (SC_7s4k.clustermarker$pct.1)/(SC_7s4k.clustermarker$pct.2)

#Barplot Freq up and downregulated Genes (FC>2)
SC_7s4k.clustermarker_subset <- subset(SC_7s4k.clustermarker, Foldchange_UP >=2 | Foldchange_DOWN >=2)
SC_7s4k.clustermarker_subset$Direction <- ifelse (SC_7s4k.clustermarker_subset$Foldchange_UP>=2,"UP","DOWN")
SC_7s4k.clustermarker_subset<- select(SC_7s4k.clustermarker_subset, Direction, cluster)
row.names(SC_7s4k.clustermarker_subset)<-NULL
SC_7s4k.clustermarker_subset$cluster<- as.character(SC_7s4k.clustermarker_subset$cluster)
SC_7s4k.clustermarker_subset<- dplyr::count(SC_7s4k.clustermarker_subset, cluster,Direction) %>% ungroup()
SC_7s4k.clustermarker_subset$Direction <- factor(x = SC_7s4k.clustermarker_subset$Direction, levels = c("UP","DOWN"))
SC_7s4k.clustermarker_subset$cluster <- factor(x = SC_7s4k.clustermarker_subset$cluster, levels = c("SC-Skin","SC-Promyel","SC-Repair","SC-Prolif","SC-Myel","SC-EC","SC-FB"))

pdf("Graphs/Main/Barplot_nDEG_SC_all repair.pdf", width = 5, height = 5)
ggplot(data=SC_7s4k.clustermarker_subset, aes(x=cluster, y=n, fill=Direction)) + theme_classic() +
  geom_bar(stat="identity", position=position_dodge())+ ggtitle("DEG-SC") +scale_fill_manual(values =color_DEG) + xlab("SC-type")+ylab("total number of DEG")+theme(axis.text.x = element_text(angle = 90))
dev.off()

#DEG SC tissue total
Idents(SC_7s4k)<-"tissue"
SC_7s4k.tissue.DEG<-FindAllMarkers(SC_7s4k,  assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
SC_7s4k.tissue.DEG$Foldchange_UP <- 2^(SC_7s4k.tissue.DEG$avg_log2FC)
SC_7s4k.tissue.DEG$Foldchange_DOWN <- 2^(-SC_7s4k.tissue.DEG$avg_log2FC)
SC_7s4k.tissue.DEG$Ratio_pct1_pct2 <- (SC_7s4k.tissue.DEG$pct.1)/(SC_7s4k.tissue.DEG$pct.2)

#Barplot Freq up and downregulated Genes (FC>2)
SC_7s4k.tissue.DEG_subset <- subset(SC_7s4k.tissue.DEG, Foldchange_UP >=2 | Foldchange_DOWN >=2)
SC_7s4k.tissue.DEG_subset$Direction <- ifelse (SC_7s4k.tissue.DEG_subset$Foldchange_UP>=2,"UP","DOWN")
SC_7s4k.tissue.DEG_subset<- select(SC_7s4k.tissue.DEG_subset, Direction, cluster)
row.names(SC_7s4k.tissue.DEG_subset)<-NULL
SC_7s4k.tissue.DEG_subset$cluster<- as.character(SC_7s4k.tissue.DEG_subset$cluster)
SC_7s4k.tissue.DEG_subset<- dplyr::count(SC_7s4k.tissue.DEG_subset, cluster,Direction) %>% ungroup()
SC_7s4k.tissue.DEG_subset$Direction <- factor(x = SC_7s4k.tissue.DEG_subset$Direction, levels = c("UP","DOWN"))
SC_7s4k.tissue.DEG_subset$cluster <- factor(x = SC_7s4k.tissue.DEG_subset$cluster, levels = c("keloid","skin"))

pdf("Graphs/Main/Barplot_nDEG_SC_tissue.pdf")
ggplot(data=SC_7s4k.tissue.DEG_subset, aes(x=cluster, y=n, fill=Direction)) + theme_classic() +
  geom_bar(stat="identity", position=position_dodge())+ ggtitle("DEG_SC-tissue") +scale_fill_manual(values =color_DEG) + xlab("SC-type")+ylab("total number of DEG")
dev.off()
Idents(SC_7s4k)<-"celltype"

#Name Cluster
Idents(SC_7s4k)<-SC_7s4k$seurat_clusters
SC_7s4k<-RenameIdents(SC_7s4k,
                      `0`="SC-Keloid", `1`="SC-Keloid", `2`="SC-Keloid", `3`="SC-FB",`4`="SC-Skin",`5`="SC-Promyel",`6`="SC-Prolif",`7`="SC-EC")
SC_7s4k$celltype<-Idents(SC_7s4k)


#Define cluster levels
clusters_ordered<-c("SC-Skin","SC-Promyel","SC-Keloid","SC-Prolif",
                    "SC-EC","SC-FB")
SC_7s4k$celltype<- factor(SC_7s4k$celltype, levels = clusters_ordered)
Idents(SC_7s4k)<-SC_7s4k$celltype
save(s.k.total, SC_7s4k,file="Object1.RData")

DotPlot(SC_7s4k,features=c("MYF5","FBLN1","RGS5","ITGB1","THY1","CD164","B4GALNT1","CD151","ACTA2","L1CAM","ITGA6","NCAM1","MCAM","SOX2","SOX9"))

#UMAPPlots, named clusters
DefaultAssay(SC_7s4k)<-"integrated"
pdf("Graphs/Main/UMAPPlots_SC_7s4k.identified.pdf")
UMAPPlot(SC_7s4k, label=T)+ scale_color_manual(values=color_SC)+ ggtitle("SC-UMAP")
dev.off()

pdf("Graphs/Side/UMAPPlots_SC_7s4k.identified.split.by.tissue.pdf", width= 20, height = 10)
UMAPPlot(SC_7s4k, label=F, split.by="tissue")+ scale_color_manual(values=color_SC)+ ggtitle("SC-UMAP")
dev.off()

pdf("Graphs/Side/UMAPPlots_SC_7s4k.identified.split.by.celltype.pdf", width= 20, height = 4)
UMAPPlot(SC_7s4k, label=F, split.by="celltype")+ scale_color_manual(values=color_SC)+ ggtitle("SC-UMAP")
dev.off()

pdf("Graphs/Side/UMAPPlots_SC_7s4k.identified.group.by.tissue.pdf", width= 20, height = 10)
UMAPPlot(SC_7s4k, label=F, group.by="tissue")+ scale_color_manual(values=color_tissue)+ ggtitle("SC-UMAP")
dev.off()

pdf("Graphs/Side/UMAPPlots_SC_7s4k.identified.group.by.sample.pdf")
UMAPPlot(SC_7s4k, label=F, group.by="sample")+scale_fill_brewer(palette="Dark2")+ ggtitle("SC-UMAP")
dev.off()

pdf("Graphs/Side/UMAPPlots_SC_7s4k.identified.split.by.origin.pdf")
UMAPPlot(SC_7s4k, label=F, split.by="origin")+scale_fill_brewer(palette="Dark2")+ ggtitle("SC-UMAP")
dev.off()

#Qualitycheck subset
pdf("Graphs/Main/QC_pctMT_SC_7s4k.pdf")
VlnPlot(SC_7s4k, features = "percent.mt", cols = c(color_SC), pt.size = 0)+ geom_hline(yintercept=5, color="red")
dev.off()
pdf("Graphs/Main/QC_nFeature_SC_7s4k.pdf")
VlnPlot(SC_7s4k, features = "nFeature_RNA", cols = c(color_SC), pt.size = 0)+ geom_hline(yintercept=300, color="red")+ geom_hline(yintercept=2500, color="red")
dev.off()
pdf("Graphs/Main/QC_nCount_SC_7s4k.pdf")
VlnPlot(SC_7s4k, features = "nCount_RNA", cols = c(color_SC), pt.size = 0)
dev.off()

pdf("Graphs/Main/QC_Featurescatter_SC7s4k.pdf", width=5, height=5)
FeatureScatter(SC_7s4k, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = c(color_SC))
FeatureScatter(SC_7s4k, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",cols = c(color_SC))
dev.off()


#Clustertree
SC_7s4k<-BuildClusterTree(SC_7s4k, dims=1:14)
DefaultAssay(SC_7s4k)<-"integrated"
pdf("Graphs/Main/Clustertree_SC_7s4k.pdf", width = 10, height=3)
PlotClusterTree(SC_7s4k, tip.color=color_SC, type="phylo", node.pos=1 ,show.node.label=FALSE, cex=1.5 ,font=4,no.margin=T, rotate.tree=180, cex.main=2)
dev.off()

# confirmation top50 SC-marker check
SC_only_clustermarker <- subset(s.k.total.clustermarker,s.k.total.clustermarker$cluster== c("SC_skin","SC_keloid"))
SC_only_clustermarker <- SC_only_clustermarker%>% top_n(n=50,wt=avg_log2FC)
SC_only_clustermarker<-arrange(SC_only_clustermarker,cluster,Foldchange_UP)
SC_only_clustermarker <- as.character(SC_only_clustermarker$gene)
SC_only_clustermarker<-SC_only_clustermarker[!duplicated(SC_only_clustermarker)]
pdf("Graphs/Main/DotPlot_Top_50_Clustermarker_SC_FCmin2_SC.pdf", width = 7, height = 12)
DotPlot(SC_7s4k, features=SC_only_clustermarker, group.by = "celltype",assay="RNA") + coord_flip()+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 45)) + ggtitle("SC_Top50 Clustermarker") + scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()

#Clustermarker total
SC_7s4k.clustermarker<-FindAllMarkers(SC_7s4k,  assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
SC_7s4k.clustermarker$Foldchange_UP <- 2^(SC_7s4k.clustermarker$avg_log2FC)
SC_7s4k.clustermarker$Foldchange_DOWN <- 2^(-SC_7s4k.clustermarker$avg_log2FC)
SC_7s4k.clustermarker$Ratio_pct1_pct2 <- (SC_7s4k.clustermarker$pct.1)/(SC_7s4k.clustermarker$pct.2)
write.xlsx(SC_7s4k.clustermarker, "Lists/Clustermarker_SC_7s4k.xlsx")

#Barplot Freq up and downregulated Genes (FC>2)
SC_7s4k.clustermarker_subset <- subset(SC_7s4k.clustermarker, Foldchange_UP >=2 | Foldchange_DOWN >=2)
SC_7s4k.clustermarker_subset$Direction <- ifelse (SC_7s4k.clustermarker_subset$Foldchange_UP>=2,"UP","DOWN")
SC_7s4k.clustermarker_subset<- select(SC_7s4k.clustermarker_subset, Direction, cluster)
row.names(SC_7s4k.clustermarker_subset)<-NULL
SC_7s4k.clustermarker_subset$cluster<- as.character(SC_7s4k.clustermarker_subset$cluster)
SC_7s4k.clustermarker_subset<- dplyr::count(SC_7s4k.clustermarker_subset, cluster,Direction) %>% ungroup()
SC_7s4k.clustermarker_subset$Direction <- factor(x = SC_7s4k.clustermarker_subset$Direction, levels = c("UP","DOWN"))
SC_7s4k.clustermarker_subset$cluster <- factor(x = SC_7s4k.clustermarker_subset$cluster, levels = c("SC-Skin","SC-Promyel","SC-Repair","SC-Prolif","SC-Myel","SC-EC","SC-FB"))

pdf("Graphs/Main/Barplot_nDEG_SC_only.pdf")
ggplot(data=SC_7s4k.clustermarker_subset, aes(x=cluster, y=n, fill=Direction)) + theme_classic() +
  geom_bar(stat="identity", position=position_dodge())+ ggtitle("DEG-SC") +scale_fill_manual(values =color_DEG) + xlab("SC-type")+ylab("total number of DEG")
dev.off()


#Heatmap
top10log.sub <- SC_7s4k.clustermarker %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
pdf("Graphs/Main/Heatmap_clustermarker_SC_7s4k_avg_logfc.pdf", 
    width=13, height=8)
DoHeatmap(SC_7s4k,assay = "SCT", cells = as.character(WhichCells(SC_7s4k))[c(T,F)], features = top10log.sub$gene, group.colors = color_SC, size = 3.5)+NoLegend()+ ggtitle("Heatmap_SC-Keloid")# + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y = element_blank())
dev.off()


#DotPlot Top up and down regulated genes 
#up
SC_7s4k.clustermarker_subset_UP <- subset(SC_7s4k.clustermarker, Foldchange_UP >=2 )
top10_UP_scsubs_clustermarker <- SC_7s4k.clustermarker_subset_UP %>% group_by(cluster) %>% top_n(n=10,wt=Foldchange_UP)
top10_UP_scsubs_clustermarker<-arrange(top10_UP_scsubs_clustermarker,cluster,Foldchange_UP)
top10_UP_scsubs_clustermarker <- as.character(top10_UP_scsubs_clustermarker$gene)
top10_UP_scsubs_clustermarker<-top10_UP_scsubs_clustermarker[!duplicated(top10_UP_scsubs_clustermarker)]

top10_UP_scsubs_clustermarker[34]="MKI67"
pdf("Graphs/Main/DotPlot_Top_10UP_Clustermarker_SC_FCmin2_SC.pdf", width = 7, height = 12)
DotPlot(SC_7s4k, features=top10_UP_scsubs_clustermarker, group.by = "celltype",assay="RNA") + coord_flip()+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 45)) + ggtitle("SC Cluster_Top-SC DEG") + scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()


#down
SC_7s4k.clustermarker_subset_DOWN <- subset(SC_7s4k.clustermarker, Foldchange_DOWN >=2 )
top10_DOWN_scsubs_clustermarker <- SC_7s4k.clustermarker_subset_DOWN %>% group_by(cluster) %>% top_n(n=10,wt=Foldchange_DOWN)
top10_DOWN_scsubs_clustermarker<-arrange(top10_DOWN_scsubs_clustermarker,cluster,Foldchange_DOWN)
top10_DOWN_scsubs_clustermarker <- as.character(top10_DOWN_scsubs_clustermarker$gene)
top10_DOWN_scsubs_clustermarker<-top10_DOWN_scsubs_clustermarker[!duplicated(top10_DOWN_scsubs_clustermarker)]

pdf("Graphs/Main/DotPlot_Top_10DOWN_Clustermarker_SC_FCmin2_SC.pdf", width = 7, height = 12)
DotPlot(SC_7s4k, features=top10_DOWN_scsubs_clustermarker, group.by = "celltype",assay="RNA") + coord_flip()+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 45)) + ggtitle("SC Cluster_Top10-SC DEG") + scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()



#Percentage total
###Pieplot/Barplot cellnumber /tissue/cluster
#subset conditions
skin_SC.subset<-subset(SC_7s4k, tissue=="skin")
keloid_SC.subset<-subset(SC_7s4k, tissue=="keloid")

#Data frame frequencies
skin_SC.cluster.frequencies<-as.data.frame(table(skin_SC.subset$celltype))
keloid_SC.cluster.frequencies<-as.data.frame(table(keloid_SC.subset$celltype))

#Barplot Cluster percentage
skin_SC.cluster.frequencies$tissue="skin"
keloid_SC.cluster.frequencies$tissue="keloid"

skin_SC.cluster.frequencies<-mutate(skin_SC.cluster.frequencies,Percentage= skin_SC.cluster.frequencies$Freq/sum(skin_SC.cluster.frequencies$Freq)*100)
keloid_SC.cluster.frequencies<-mutate(keloid_SC.cluster.frequencies,Percentage= keloid_SC.cluster.frequencies$Freq/sum(keloid_SC.cluster.frequencies$Freq)*100)



skin_kel_SC.cluster.df <- rbind(skin_SC.cluster.frequencies,keloid_SC.cluster.frequencies)
skin_kel_SC.cluster.df$tissue <- factor(x = skin_kel_SC.cluster.df$tissue, levels = c("skin", "keloid"))

pdf("Graphs/Main/Barplot_SC_subset_skin_kel_cellfrequencie.pdf", width = 4, height = 6)
ggplot(data=skin_kel_SC.cluster.df, aes(x=tissue, y=Freq, fill=Var1)) +
  geom_bar(stat="identity", position="stack")+scale_fill_manual(values =color_SC)+
  ylim(0,1000) +  theme_classic() +ggtitle("Cell source_schwann cell cluster") + xlab("condition")+ylab("total cell count") 
dev.off()

pdf("Graphs/Main/Barplot_SC_subset_skin_kel_cellpercentage.pdf", width = 4, height = 6)
ggplot(data=skin_kel_SC.cluster.df, aes(x=tissue, y=Percentage, fill=Var1)) +
  geom_bar(stat="identity", position="stack")+scale_fill_manual(values =color_SC) +  theme_classic() +ggtitle("Cell source_schwann cell cluster") + xlab("condition")+ylab("total cell count") 
dev.off()

pdf("Graphs/Side/PiePlot_SC_subset_skin_kel_cellfrequencie.pdf", width = 10, height = 10)
pie(skin_SC.cluster.frequencies$Freq, labels = skin_SC.cluster.frequencies$Var1, main="skin",  clockwise = T, col = c(color_SC))
pie(keloid_SC.cluster.frequencies$Freq, labels = keloid_SC.cluster.frequencies$Var1, main="Keloid",  clockwise = T, col = c(color_SC))
dev.off()


#Percentage total
###Pieplot/Barplot cellnumber /sample/cluster
#subset conditions
SC_skin1.subset<-subset(SC_7s4k, sample=="skin_1")
SC_skin2.subset<-subset(SC_7s4k, sample=="skin_2")
SC_skin4.subset<-subset(SC_7s4k, sample=="skin_4")
SC_skin3.subset<-subset(SC_7s4k, sample=="skin_3")
SC_skin5.subset<-subset(SC_7s4k, sample=="skin_5")
SC_skin6.subset<-subset(SC_7s4k, sample=="skin_6")
SC_skin7.subset<-subset(SC_7s4k, sample=="skin_7")
SC_keloid1.subset<-subset(SC_7s4k, sample=="keloid_1")
SC_keloid2.subset<-subset(SC_7s4k, sample=="keloid_2")
SC_keloid3L.subset<-subset(SC_7s4k, sample=="keloid_3L")
SC_keloid3R.subset<-subset(SC_7s4k, sample=="keloid_3R")

#Data frame frequencies
SC_skin1.cluster.frequencies<-as.data.frame(table(SC_skin1.subset$celltype))
SC_skin2.cluster.frequencies<-as.data.frame(table(SC_skin2.subset$celltype))
SC_skin3.cluster.frequencies<-as.data.frame(table(SC_skin3.subset$celltype))
SC_skin4.cluster.frequencies<-as.data.frame(table(SC_skin4.subset$celltype))
SC_skin5.cluster.frequencies<-as.data.frame(table(SC_skin5.subset$celltype))
SC_skin6.cluster.frequencies<-as.data.frame(table(SC_skin6.subset$celltype))
SC_skin7.cluster.frequencies<-as.data.frame(table(SC_skin7.subset$celltype))
SC_keloid1.cluster.frequencies<-as.data.frame(table(SC_keloid1.subset$celltype))
SC_keloid2.cluster.frequencies<-as.data.frame(table(SC_keloid2.subset$celltype))
SC_keloid3L.cluster.frequencies<-as.data.frame(table(SC_keloid3L.subset$celltype))
SC_keloid3R.cluster.frequencies<-as.data.frame(table(SC_keloid3R.subset$celltype))

#Barplot Cluster percentage
SC_skin1.cluster.frequencies$sample="skin_1"
SC_skin2.cluster.frequencies$sample="skin_2"
SC_skin3.cluster.frequencies$sample="skin_3"
SC_skin4.cluster.frequencies$sample="skin_4"
SC_skin5.cluster.frequencies$sample="skin_5"
SC_skin6.cluster.frequencies$sample="skin_6"
SC_skin7.cluster.frequencies$sample="skin_7"
SC_keloid1.cluster.frequencies$sample="keloid_1"
SC_keloid2.cluster.frequencies$sample="keloid_2"
SC_keloid3L.cluster.frequencies$sample="keloid_3L"
SC_keloid3R.cluster.frequencies$sample="keloid_3R"

SC_skin1.cluster.frequencies<-mutate(SC_skin1.cluster.frequencies,Percentage= SC_skin1.cluster.frequencies$Freq/sum(SC_skin1.cluster.frequencies$Freq)*100)
SC_skin2.cluster.frequencies<-mutate(SC_skin2.cluster.frequencies,Percentage= SC_skin2.cluster.frequencies$Freq/sum(SC_skin2.cluster.frequencies$Freq)*100)
SC_skin3.cluster.frequencies<-mutate(SC_skin3.cluster.frequencies,Percentage= SC_skin3.cluster.frequencies$Freq/sum(SC_skin3.cluster.frequencies$Freq)*100)
SC_skin4.cluster.frequencies<-mutate(SC_skin4.cluster.frequencies,Percentage= SC_skin4.cluster.frequencies$Freq/sum(SC_skin4.cluster.frequencies$Freq)*100)
SC_skin5.cluster.frequencies<-mutate(SC_skin5.cluster.frequencies,Percentage= SC_skin5.cluster.frequencies$Freq/sum(SC_skin5.cluster.frequencies$Freq)*100)
SC_skin6.cluster.frequencies<-mutate(SC_skin6.cluster.frequencies,Percentage= SC_skin6.cluster.frequencies$Freq/sum(SC_skin6.cluster.frequencies$Freq)*100)
SC_skin7.cluster.frequencies<-mutate(SC_skin7.cluster.frequencies,Percentage= SC_skin7.cluster.frequencies$Freq/sum(SC_skin7.cluster.frequencies$Freq)*100)
SC_keloid1.cluster.frequencies<-mutate(SC_keloid1.cluster.frequencies,Percentage= SC_keloid1.cluster.frequencies$Freq/sum(SC_keloid1.cluster.frequencies$Freq)*100)
SC_keloid2.cluster.frequencies<-mutate(SC_keloid2.cluster.frequencies,Percentage= SC_keloid2.cluster.frequencies$Freq/sum(SC_keloid2.cluster.frequencies$Freq)*100)
SC_keloid3L.cluster.frequencies<-mutate(SC_keloid3L.cluster.frequencies,Percentage= SC_keloid3L.cluster.frequencies$Freq/sum(SC_keloid3L.cluster.frequencies$Freq)*100)
SC_keloid3R.cluster.frequencies<-mutate(SC_keloid3R.cluster.frequencies,Percentage= SC_keloid3R.cluster.frequencies$Freq/sum(SC_keloid3R.cluster.frequencies$Freq)*100)


SC_skin_kel_samp.cluster.df <- rbind(SC_skin1.cluster.frequencies,SC_skin2.cluster.frequencies,SC_skin3.cluster.frequencies,SC_skin4.cluster.frequencies,SC_skin5.cluster.frequencies,SC_skin6.cluster.frequencies,SC_skin7.cluster.frequencies,SC_keloid1.cluster.frequencies,SC_keloid2.cluster.frequencies,SC_keloid3L.cluster.frequencies,SC_keloid3R.cluster.frequencies)
SC_skin_kel_samp.cluster.df$sample <- factor(x = SC_skin_kel_samp.cluster.df$sample, levels = c("skin_1","skin_2","skin_3","skin_4","skin_5","skin_6","skin_7", "keloid_1", "keloid_2", "keloid_3L", "keloid_3R"))

pdf("Graphs/Side/Barplot_SC_skin_kel_cell_clusterpercentage.pdf", width = 10, height = 5)
ggplot(data=SC_skin_kel_samp.cluster.df, aes(x=sample, y=Percentage, fill=Var1)) +
  geom_bar(stat="identity", position="stack")+scale_fill_manual(values =color_SC)+ theme_classic()
dev.off()

pdf("Graphs/Side/Barplot_SC_skin_kel_cell_clusterfrequencie.pdf", width = 7, height = 7)
ggplot(data=SC_skin_kel_samp.cluster.df, aes(x=sample, y=Freq, fill=Var1)) +
  geom_bar(stat="identity", position="stack") +  theme_classic() +
  ggtitle("Cell source_schwann cell cluster") + xlab("condition")+ylab("total cell count") +scale_fill_manual(values =color_SC)
dev.off()

pdf("Graphs/Side/PiePlot_SC_skin_kel_cell_clusterfrequencie.pdf", width = 10, height = 10)
pie(SC_skin1.cluster.frequencies$Freq, labels = SC_skin1.cluster.frequencies$Var1, main="skin_1",  clockwise = T, col = c(color_SC))
pie(SC_skin2.cluster.frequencies$Freq, labels = SC_skin2.cluster.frequencies$Var1, main="skin_2",  clockwise = T, col = c(color_SC))
pie(SC_skin4.cluster.frequencies$Freq, labels = SC_skin4.cluster.frequencies$Var1, main="skin_4",  clockwise = T, col = c(color_SC))
pie(SC_skin5.cluster.frequencies$Freq, labels = SC_skin5.cluster.frequencies$Var1, main="skin_5",  clockwise = T, col = c(color_SC))
pie(SC_skin6.cluster.frequencies$Freq, labels = SC_skin6.cluster.frequencies$Var1, main="skin_6",  clockwise = T, col = c(color_SC))
pie(SC_skin7.cluster.frequencies$Freq, labels = SC_skin7.cluster.frequencies$Var1, main="skin_7",  clockwise = T, col = c(color_SC))
pie(SC_keloid1.cluster.frequencies$Freq, labels = SC_keloid1.cluster.frequencies$Var1, main="Keloid_1",  clockwise = T, col = c(color_SC))
pie(SC_keloid2.cluster.frequencies$Freq, labels = SC_keloid2.cluster.frequencies$Var1, main="Keloid_2",  clockwise = T, col = c(color_SC))
pie(SC_keloid3L.cluster.frequencies$Freq, labels = SC_keloid3L.cluster.frequencies$Var1, main="Keloid_3L",  clockwise = T, col = c(color_SC))
pie(SC_keloid3R.cluster.frequencies$Freq, labels = SC_keloid3R.cluster.frequencies$Var1, main="Keloid_3R",  clockwise = T, col = c(color_SC))
dev.off()


#Paper-Marker and other stuff
DefaultAssay(SC_7s4k)<-"RNA"
SC_marker_general<- c("NGFR","SOX10","CDH19","EGR1","GAP43","NCAM1","S100B","EGR2")
pdf("Graphs/Side/StackedVlnplot_SC_marker_general_SCsubset.pdf", width = 20, height =20)
StackedVlnPlot(SC_7s4k, features=SC_marker_general)
dev.off()

SC_marker_Toma_etal._32349983<- c("NGFR","SOX10","CDH2","L1CAM","EDNRB","EMP1","SEMA3E","POU3F1","EGR2","MAG","MBP","PMP22","MPZ","PLP1","AIF1","PDGFRA","PCOLCE2","DPP4","DPT","LY6C1","COMP","ETVL","WIFL","SOX9","OSR2","MEOX1","SLC2A1","CASQ2","MSLN","MKI67","TOP2A")
pdf("Graphs/Side/DotPlot_SC_marker_Toma_etal._32349983_SCsubset.pdf")
DotPlot(s.k.total, features= SC_marker_Toma_etal._32349983 , group.by = "celltype") + coord_flip()+scale_color_gradient(low="grey", "high"="red")+theme(axis.text.x = element_text(angle = 45))
DotPlot(SC_7s4k, features=SC_marker_Toma_etal._32349983, group.by = "celltype") + coord_flip()+scale_color_gradient(low="grey", "high"="red")+theme(axis.text.x = element_text(angle = 45))
dev.off()

SC_Marker_Tamara <-c("SOX10", "S100B", "NGFR", "GFAP", "ERBB3","TYRP1","MLANA","PMEL","MITF","DCT","L1CAM","MBP","MPZ","NCAM1","PLP1","GLDN","NRXN1","TJP1","SCN7A","HLA-DRA","JUN","SOX2","ZEB2","SHH","KIT","GDNF","LIF","CLCF1","BTC","CCL2","UCN2")
pdf("Graphs/Side/DotPlot_SC_Marker_Tamara_SCsubset.pdf")
DotPlot(s.k.total, features= SC_Marker_Tamara , group.by = "celltype") + coord_flip()+scale_color_gradient(low="grey", "high"="red")+theme(axis.text.x = element_text(angle = 45))
DotPlot(SC_7s4k, features=SC_Marker_Tamara, group.by = "celltype") + coord_flip()+scale_color_gradient(low="grey", "high"="red")+theme(axis.text.x = element_text(angle = 45))
dev.off()

NatureSC_marker <- c("ITGA4","TFAP2A","CDH2","ERBB3","L1CAM","NGFR","SOX10","CDH19","FABP7","DHH","MPZ","GAP43","PMP22","PLP1","EGR1","GFAP","NCAM1","S100B","POU3F1","FOXO4", "EGR2")
pdf("Graphs/Side/DotPlot_SC_NatureSC_marker.pdf")
DotPlot(s.k.total, features= NatureSC_marker , group.by = "celltype") + coord_flip()+scale_color_gradient(low="grey", "high"="red")+theme(axis.text.x = element_text(angle = 45))
DotPlot(SC_7s4k, features=NatureSC_marker, group.by = "celltype") + coord_flip()+scale_color_gradient(low="grey", "high"="red")+theme(axis.text.x = element_text(angle = 45))
dev.off()


Idents(s.k.total)<-"celltype.tissue"
pdf("Graphs/Main/D002.1_CCL2.pdf", height = 10, width = 10)
VlnPlot(s.k.total,c("CCL2","CCR2","CCR4","MMP12"), pt.size=0, cols=color_UMAP_short_double, ncol = 2)
dev.off()
pdf("Graphs/Main/DotPlot_SC_Weiss_Gene_Expression.pdf", height = 10, width = 10)
DotPlot(s.k.total, features = c("TNFSF9","IL6","IL1R1","CXCL1","CXCL2","IL8","ICAM1","IL1B","HLA-C","B2M","HLA-B","HLA-A","CD40","HLA-DRB1","HLA-DRB5","HLA-DRA","IL6R","TNFRSF14","IL1A","CD58","PDGFB","PDGFA","CCL7","CCL2","LILRB2","LILRB4","CD86","HAVCR2","CCL5","CCL4","CCL3","TNF","TNFSF14","CD80","TNFRSF4","IFNG","ICOSLG","TIMD4","CSF2","CCL8","CD70","CD276","VEGFA","LIF","CXCL6","SEMA3A","CD274","PDCD1LG2","CXCL12","VTCN1","TNFSF18"))+coord_flip()+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 45)) + ggtitle("SC_Genes _Weiss_paper") + scale_y_discrete(guide = guide_axis(angle = 90))+ theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short_double))
dev.off()
Idents(s.k.total)<-"celltype"
pdf("Graphs/Main/VlnPlot_SC_HLA_check.pdf", height = 20, width = 17)
VlnPlot(SC_7s4k,c("HLA-A","HLA-B","HLA-C","HLA-DRA","HLA-DRB1","HLA-E","HLA-G","HLA-DMA","HLA-DMB","HLA-DRB5","HLA-DPB1","HLA-DQA1","HLA-DQB1","HLA-F"), ncol=3, cols = c(color_SC))
dev.off()


##############
##Pseudotime##
##############
library(monocle3)

#### Create a Monocle CDS Object
# Project PC dimensions to whole data set
SC_convert <- ProjectDim(SC_7s4k, reduction = "pca")

# Create an expression matrix
expression_matrix <- SC_convert@assays$RNA@counts

# Get cell metadata
cell_metadata <- SC_convert@meta.data
if (all.equal(colnames(expression_matrix), rownames(cell_metadata))) {
  print(sprintf("Cell identifiers match"))
} else {
  print(sprintf("Cell identifier mismatch - %i cells in expression matrix, %i cells in metadata",
                ncol(expression_matrix), nrow(cell_metadata)))
  print("If the counts are equal, sort differences will throw this error")
}

# get gene annotations
gene_annotation <- data.frame(gene_short_name = rownames(SC_convert@assays$RNA), row.names = rownames(SC_convert@assays$RNA))
if (all.equal(rownames(expression_matrix), rownames(gene_annotation))) {
  print(sprintf("Gene identifiers all match"))
} else {
  print(sprintf("Gene identifier mismatch - %i genes in expression matrix, %i gene in gene annotation",
                nrow(expression_matrix), nrow(gene_annotation)))
  print("If the counts are equal, sort differences will throw this error")
}

# Seurat-derived CDS
SC_7s4k.cds <- new_cell_data_set(expression_matrix,
                                 cell_metadata = cell_metadata,
                                 gene_metadata = gene_annotation)

# Transfer Seurat embeddings
# Note that these may be calculated on the Integrated object, not the counts
#   and thus will involve fewer genes
reducedDim(SC_7s4k.cds, type = "PCA") <- SC_7s4k@reductions$pca@cell.embeddings 
SC_7s4k.cds@preprocess_aux$prop_var_expl <- SC_7s4k@reductions$pca@stdev
plot_pc_variance_explained(SC_7s4k.cds)

# Transfer Seurat UMAP embeddings
SC_7s4k.cds@int_colData@listData$reducedDims$UMAP <- SC_7s4k@reductions$umap@cell.embeddings
#    plot_cells(SC_7s4k.cds)

# Copy cluster info from Seurat
SC_7s4k.cds@clusters$UMAP_so$clusters <- SC_7s4k@meta.data$gt_tp_cell_type_integrated_.0.9

SC_7s4k.cds <- cluster_cells(SC_7s4k.cds, reduction_method = "UMAP")

# Fix from https://gitmemory.com/cole-trapnell-lab
rownames(SC_7s4k.cds@principal_graph_aux$UMAP$dp_mst) <- NULL
colnames(SC_7s4k.cds@int_colData@listData$reducedDims$UMAP) <- NULL

DimPlot(SC_7s4k, reduction = "umap")

SC_7s4k.cds = cluster_cells(SC_7s4k.cds)

p1 <- plot_cells(SC_7s4k.cds,group_label_size = 3.5, cell_size = 2) + ggtitle("SC_7s4k_ _ cluster")
p2 <- plot_cells(SC_7s4k.cds, color_cells_by="celltype", group_cells_by="partition",group_label_size = 3.5, cell_size = 2) + ggtitle("SC_7s4k_partition")
wrap_plots(p1, p2)

SC_7s4k.cds <- learn_graph(SC_7s4k.cds,learn_graph_control = list(prune_graph=F, ncenter=20))

pdf("Graphs/Main/Pseudo_Principal_graph_SC_7s4k.pdf")
plot_cells(SC_7s4k.cds, label_groups_by_cluster = T,color_cells_by="celltype", label_leaves = FALSE, label_branch_points = FALSE, cell_size=2, label_cell_groups = F)+ scale_color_manual(values=color_SC) + ggtitle("Principal graph")
dev.off()

SC_7s4k.cds<- order_cells(SC_7s4k.cds)

pdf("Graphs/Main/Pseudo_Myel_start_SC_7s4k.pdf")
plot_cells(SC_7s4k.cds, label_groups_by_cluster = T,color_cells_by="pseudotime", label_leaves = FALSE, label_branch_points = FALSE, cell_size=2, label_cell_groups = F) + ggtitle("Pseudotime_ SC-Skin as root")
dev.off()

plot_cells(SC_7s4k.cds,
           genes=c("S100B"),           
           label_cell_groups=F,
           show_trajectory_graph=T, cell_size=1.5)


Mod_cds <- SC_7s4k.cds[,grepl("SC", colData(SC_7s4k.cds)$celltype, ignore.case=TRUE)]
plot_cells(Mod_cds, color_cells_by="celltype")
pr_graph_test_res<- graph_test (Mod_cds, neighbor_graph = "principal_graph", cores = 1)
pr_deg_ids <- row.names(subset(pr_graph_test_res,q_value < 0.05))
gene_module_df <- find_gene_modules(Mod_cds[pr_deg_ids,], resolution=1e-2)#, resolution=1e-1
plot_cells(SC_7s4k.cds, genes=gene_module_df, 
           show_trajectory_graph=FALSE, 
           label_cell_groups=FALSE)

#Adapt Excel-output
X1 <- subset(pr_graph_test_res, q_value < 0.05)
X1 <- cbind(rownames(X1),X1)
X1 <- subset(X1, select=c("rownames(X1)","gene_short_name","morans_I"))
write.xlsx(X1,"Lists/pr_graph_SC_7s4k.xlsx")


#Analyse Pseudotime
IGFBP_genes <- c("IGFBP3", "IGFBP5", "IGFBP7")
IGFBP_lineage_cds <- SC_7s4k.cds[rowData(SC_7s4k.cds)$gene_short_name %in% IGFBP_genes,]#

pdf("Graphs/Main/Pseudo_IGFBP_SC_7s4k.pdf", width = 8, height=10)
plot_genes_in_pseudotime(IGFBP_lineage_cds, color_cells_by = "celltype", min_expr=0.5, cell_size=1.5)  + scale_color_manual(values=color_SC) + ggtitle("IGFBP expression along Pseudotime")
dev.off()

AP1_genes <- c("JUN","JUNB","JUND","FOS","FOSL1","FOSL2")
AP1_lineage_cds <- SC_7s4k.cds[rowData(SC_7s4k.cds)$gene_short_name %in% AP1_genes,]

pdf("Graphs/Main/Pseudo_AP1_SC_7s4k.pdf", width = 7, height=10)
plot_genes_in_pseudotime(AP1_lineage_cds, color_cells_by = "celltype", min_expr=0.5, cell_size=1.5)+ scale_color_manual(values=color_SC) + ggtitle("AP-1 member expression along Pseudotime")
dev.off()

SC_genes <- c("S100B","PMP22","NGFR","NRXN","NES","NOV","MBP","MPZ","PLP1","NCAM1","EGR2","EGR1")
SC_lineage_cds <- SC_7s4k.cds[rowData(SC_7s4k.cds)$gene_short_name %in% SC_genes,]

pdf("Graphs/Main/Pseudo_SC_SC_7s4k.pdf", width = 8, height=20)
plot_genes_in_pseudotime(SC_lineage_cds, color_cells_by = "celltype", min_expr=0.5, cell_size=1.5)+ scale_color_manual(values=color_SC) + ggtitle("AP-1 member expression along Pseudotime")
dev.off()

SC_genes_2 <- c("MBP","CD9","MPZ","PLP1","AHNAK","CAV1","GPM6B","PMP22","ZEB2","NRXN1","SPARCL1","PTN","PTPRZ1","PDGFA","CCN3","TGFBI","TNC","CALR")
SC2_lineage_cds <- SC_7s4k.cds[rowData(SC_7s4k.cds)$gene_short_name %in% SC_genes_2,colData(SC_7s4k.cds)$celltype %in% c("SC-Skin","SC-Promyel","SC-Repair")]

pdf("Graphs/Main/Pseudo_SC2_SC_7s4k.pdf", width = 10, height=10)
plot_genes_in_pseudotime(SC2_lineage_cds, color_cells_by = "celltype", min_expr=0.5, cell_size=1.5, ncol = 3)+ scale_color_manual(values=color_SC) + ggtitle("SCs dedifferentiate along Pseudotime")+ylim(0.5,6)
dev.off()

MAC_genes <- c("EDNRB","LGALS3","RARRES2","GFRA3","NRXN1","EMP2","FXYD1","CTSC","PLP1","ZEB2")
MAC_lineage_cds <- SC_7s4k.cds[rowData(SC_7s4k.cds)$gene_short_name %in% MAC_genes,]

pdf("Graphs/Main/Pseudo_SC_SC_7s4k.pdf", width = 8, height=20)
plot_genes_in_pseudotime(MAC_lineage_cds, color_cells_by = "celltype", min_expr=0.1, cell_size=1.5)+ scale_color_manual(values=color_SC) + ggtitle("AP-1 member expression along Pseudotime")
dev.off()

# analyse SC-Skin to sC-Repair dedifferentiation
FeaturePlot(s.k.total,"NES")
#Identify pseudotime relevant genes
Idents(SC_7s4k)<-"celltype"
SC_Skin_vs_Promyel<-FindMarkers(SC_7s4k, ident.1="SC-Skin", ident.2=c("SC-Promyel"), assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
SC_Skin_vs_Promyel$Foldchange_UP <- 2^(SC_Skin_vs_Promyel$avg_log2FC)
SC_Skin_vs_Promyel$Foldchange_DOWN <- 2^(-SC_Skin_vs_Promyel$avg_log2FC)
SC_Skin_vs_Promyel$Ratio_pct1_pct2 <- (SC_Skin_vs_Promyel$pct.1)/(SC_Skin_vs_Promyel$pct.2)
SC_Skin_vs_Promyel$comparison <- "SC-skin vs promyel"
write.xlsx(SC_Skin_vs_Promyel, "Lists/Pseudotime_relevant/SC_Skin_vs_Promyel.xlsx")

SC_Promyel_vs_Skin_repair<-FindMarkers(SC_7s4k, ident.1="SC-Promyel", ident.2=c("SC-Repair","SC-Skin"), assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
SC_Promyel_vs_Skin_repair$Foldchange_UP <- 2^(SC_Promyel_vs_Skin_repair$avg_log2FC)
SC_Promyel_vs_Skin_repair$Foldchange_DOWN <- 2^(-SC_Promyel_vs_Skin_repair$avg_log2FC)
SC_Promyel_vs_Skin_repair$Ratio_pct1_pct2 <- (SC_Promyel_vs_Skin_repair$pct.1)/(SC_Promyel_vs_Skin_repair$pct.2)
SC_Promyel_vs_Skin_repair$comparison <- "SC_Promyel_vs_Skin_repair"
write.xlsx(SC_Promyel_vs_Skin_repair, "Lists/Pseudotime_relevant/SC_Promyel_vs_Skin_repair.xlsx")

SC_Repair_vs_SC_END<-FindMarkers(SC_7s4k, ident.1="SC-Repair", ident.2=c("SC-EC","SC-FB","SC-Prolif"), assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
SC_Repair_vs_SC_END$Foldchange_UP <- 2^(SC_Repair_vs_SC_END$avg_log2FC)
SC_Repair_vs_SC_END$Foldchange_DOWN <- 2^(-SC_Repair_vs_SC_END$avg_log2FC)
SC_Repair_vs_SC_END$Ratio_pct1_pct2 <- (SC_Repair_vs_SC_END$pct.1)/(SC_Repair_vs_SC_END$pct.2)
SC_Repair_vs_SC_END$comparison <- "SC_Repair_vs_SC_END"
write.xlsx(SC_Repair_vs_SC_END, "Lists/Pseudotime_relevant/SC_Repair_vs_SC_END.xlsx")

SC_Repair_vs_SC_START<-FindMarkers(SC_7s4k, ident.1="SC-Repair", ident.2=c("SC-Skin","SC-Promyel"), assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
SC_Repair_vs_SC_START$Foldchange_UP <- 2^(SC_Repair_vs_SC_START$avg_log2FC)
SC_Repair_vs_SC_START$Foldchange_DOWN <- 2^(-SC_Repair_vs_SC_START$avg_log2FC)
SC_Repair_vs_SC_START$Ratio_pct1_pct2 <- (SC_Repair_vs_SC_START$pct.1)/(SC_Repair_vs_SC_START$pct.2)
SC_Repair_vs_SC_START$comparison <- "SC_Repair_vs_SC_START"
write.xlsx(SC_Repair_vs_SC_START, "Lists/Pseudotime_relevant/SC_Repair_vs_SC_START.xlsx")

#relevant Genes by Pseudotime
Gene_extract<- read.xlsx("D:/Direder/Projekte/D002_10xScarWars/D002_Berechnungen/D002_c12_7s4k/Lists/Pseudotime_relevant/Gene_extrakt.xlsx",2, header = T)
Gene_extract1 <- Gene_extract$T0
Gene_extract2 <- Gene_extract$T5
Gene_extract3 <- Gene_extract$T10
Gene_extract4 <- Gene_extract$T15
Gene_extract5 <- Gene_extract$Tmix

Test_lineage_cds <- SC_7s4k.cds[rowData(SC_7s4k.cds)$gene_short_name %in% Gene_extract1,colData(SC_7s4k.cds)$celltype %in% c("SC-Skin","SC-Promyel","SC-Repair")]
pdf("Graphs/Main/Pseudo_T0_SC_7s4k.pdf", width = 20, height=15)
plot_genes_in_pseudotime(Test_lineage_cds, color_cells_by = "celltype", min_expr=0.5, cell_size=1,ncol = 10)+ scale_color_manual(values=color_SC) + ggtitle("AP-1 member expression along Pseudotime")
dev.off()

Test_lineage_cds <- SC_7s4k.cds[rowData(SC_7s4k.cds)$gene_short_name %in% Gene_extract2,colData(SC_7s4k.cds)$celltype %in% c("SC-Skin","SC-Promyel","SC-Repair")]
pdf("Graphs/Main/Pseudo_T5_SC_7s4k.pdf", width = 20, height=10)
plot_genes_in_pseudotime(Test_lineage_cds, color_cells_by = "celltype", min_expr=0.5, cell_size=1,ncol = 10)+ scale_color_manual(values=color_SC) + ggtitle("AP-1 member expression along Pseudotime")
dev.off()

Test_lineage_cds <- SC_7s4k.cds[rowData(SC_7s4k.cds)$gene_short_name %in% Gene_extract3,colData(SC_7s4k.cds)$celltype %in% c("SC-Skin","SC-Promyel","SC-Repair")]
pdf("Graphs/Main/Pseudo_T10_SC_7s4k.pdf", width = 20, height=10)
plot_genes_in_pseudotime(Test_lineage_cds, color_cells_by = "celltype", min_expr=0.5, cell_size=1,ncol = 10)+ scale_color_manual(values=color_SC) + ggtitle("AP-1 member expression along Pseudotime")
dev.off()

Test_lineage_cds <- SC_7s4k.cds[rowData(SC_7s4k.cds)$gene_short_name %in% Gene_extract4,colData(SC_7s4k.cds)$celltype %in% c("SC-Skin","SC-Promyel","SC-Repair")]
pdf("Graphs/Main/Pseudo_T15_SC_7s4k.pdf", width = 20, height=10)
plot_genes_in_pseudotime(Test_lineage_cds, color_cells_by = "celltype", min_expr=0.5, cell_size=1,ncol = 10)+ scale_color_manual(values=color_SC) + ggtitle("AP-1 member expression along Pseudotime")
dev.off()

Test_lineage_cds <- SC_7s4k.cds[rowData(SC_7s4k.cds)$gene_short_name %in% Gene_extract5,colData(SC_7s4k.cds)$celltype %in% c("SC-Skin","SC-Promyel","SC-Repair")]
pdf("Graphs/Main/Pseudo_Tmix_SC_7s4k.pdf", width = 20, height=10)
plot_genes_in_pseudotime(Test_lineage_cds, color_cells_by = "celltype", min_expr=0.5, cell_size=1,ncol = 10)+ scale_color_manual(values=color_SC) + ggtitle("AP-1 member expression along Pseudotime")
dev.off()

#relevant Genes by Alphabet
Gene_extract<- read.xlsx("D:/Direder/Projekte/D002_10xScarWars/D002_Berechnungen/D002_c12_7s4k/Lists/Pseudotime_relevant/Gene_extrakt.xlsx",3, header = F)
Gene_extract1 <- Gene_extract$X1
Gene_extract2 <- Gene_extract$X2


Test_lineage_cds <- SC_7s4k.cds[rowData(SC_7s4k.cds)$gene_short_name %in% Gene_extract1,colData(SC_7s4k.cds)$celltype %in% c("SC-Skin","SC-Promyel","SC-Repair")]
pdf("Graphs/Main/Alph1_SC_7s4k.pdf", width = 20, height=15)
plot_genes_in_pseudotime(Test_lineage_cds, color_cells_by = "celltype", min_expr=0.5, cell_size=1,ncol = 10)+ scale_color_manual(values=color_SC) + ggtitle("AP-1 member expression along Pseudotime")
dev.off()

Test_lineage_cds <- SC_7s4k.cds[rowData(SC_7s4k.cds)$gene_short_name %in% Gene_extract2,colData(SC_7s4k.cds)$celltype %in% c("SC-Skin","SC-Promyel","SC-Repair")]
pdf("Graphs/Main/Alph2_SC_7s4k.pdf", width = 20, height=15)
plot_genes_in_pseudotime(Test_lineage_cds, color_cells_by = "celltype", min_expr=0.5, cell_size=1,ncol = 10)+ scale_color_manual(values=color_SC) + ggtitle("AP-1 member expression along Pseudotime")
dev.off()



# other Pseudotimes
SC_7s4k.cds<- order_cells(SC_7s4k.cds)

pdf("Graphs/Main/Pseudo_Prolif_start_SC_7s4k.pdf")
plot_cells(SC_7s4k.cds, label_groups_by_cluster = T,color_cells_by="pseudotime", label_leaves = FALSE, label_branch_points = FALSE, cell_size=2, label_cell_groups = F) + ggtitle("Pseudotime_Proliferating SCs as root")
dev.off()

SC_7s4k.cds<- order_cells(SC_7s4k.cds)

pdf("Graphs/Main/Pseudo_EC_start_SC_7s4k.pdf")
plot_cells(SC_7s4k.cds, label_groups_by_cluster = T,color_cells_by="pseudotime", label_leaves = FALSE, label_branch_points = FALSE, cell_size=2, label_cell_groups = F) + ggtitle("Pseudotime_ SC-EC as root")
dev.off()

SC_7s4k.cds<- order_cells(SC_7s4k.cds)

pdf("Graphs/Main/Pseudo_FB_start_SC_7s4k.pdf")
plot_cells(SC_7s4k.cds, label_groups_by_cluster = T,color_cells_by="pseudotime", label_leaves = FALSE, label_branch_points = FALSE, cell_size=2, label_cell_groups = F) + ggtitle("Pseudotime_ SC-FB as root")
dev.off()

plot_cells(SC_7s4k.cds, genes=c("IBSP"), label_cell_groups=FALSE,cell_size = 3, label_leaves=FALSE)



####validation SC-EC / SC-FB
Idents(s.k.total)<- s.k.total$celltype
EC_FB_KC_total<- subset(s.k.total,idents = c("FB","EC","LEC","KC"))
DefaultAssay(EC_FB_KC_total)<-"RNA"
EC_FB_KC_total[["percent.mt"]] <- PercentageFeatureSet(EC_FB_KC_total, pattern = "^MT-")
EC_FB_KC_total <- SCTransform(EC_FB_KC_total, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
EC_FB_KC_total[['integrated']] <- NULL

EC_FB_SC <- SC_7s4k
DefaultAssay(s.k.total)<-"RNA"
EC_FB_SC[["percent.mt"]] <- PercentageFeatureSet(EC_FB_SC, pattern = "^MT-")
EC_FB_SC <- SCTransform(EC_FB_SC, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
EC_FB_SC[['integrated']] <- NULL


EC_FB_KC_compare.list<-list(EC_FB_KC_total,EC_FB_SC)
EC_FB_KC_compare.features <- SelectIntegrationFeatures(object.list = EC_FB_KC_compare.list, nfeatures = 3000)
EC_FB_KC_compare.list <- PrepSCTIntegration(EC_FB_KC_compare.list, anchor.features = EC_FB_KC_compare.features)
EC_FB_KC_compare.list <- lapply(EC_FB_KC_compare.list, RunPCA, verbose = F, features= EC_FB_KC_compare.features)
EC_FB_KC_compare.anchors<-FindIntegrationAnchors(EC_FB_KC_compare.list,normalization.method = "SCT", anchor.features = EC_FB_KC_compare.features, reduction = "rpca")
EC_FB_KC_compare<-IntegrateData(anchorset=EC_FB_KC_compare.anchors, normalization.method = "SCT")
EC_FB_KC_compare <- RunPCA(EC_FB_KC_compare)
ElbowPlot(EC_FB_KC_compare, ndims = 20)
EC_FB_KC_compare <- RunUMAP(EC_FB_KC_compare, dims = 1:20)
EC_FB_KC_compare <- FindNeighbors(EC_FB_KC_compare, dims = 1:20)
#DefaultAssay(EC_FB_KC_compare)<- "integrated"
#EC_FB_KC_compare <- FindClusters(EC_FB_KC_compare, resolution = 0.45)
UMAPPlot(EC_FB_KC_compare, label=T, group.by="celltype")
DefaultAssay(EC_FB_KC_compare)<- "RNA"
EC_FB_KC_compare <- NormalizeData(EC_FB_KC_compare)
FeaturePlot(EC_FB_KC_compare,c("DCN","TOP2A","SELE","MPZ","PLP1","S100B"))

#FB/EC
Idents(EC_FB_KC_compare)<- EC_FB_KC_compare$celltype
SC_FB_vs_FB_DEG <- FindMarkers(EC_FB_KC_compare, ident.1 = "SC-FB", ident.2=c("FB"))
SC_FB_vs_FB_DEG$Foldchange_UP <- 2^(SC_FB_vs_FB_DEG$avg_log2FC)
SC_FB_vs_FB_DEG$Foldchange_DOWN <- 2^(-SC_FB_vs_FB_DEG$avg_log2FC)
SC_FB_vs_FB_DEG$Ratio_pct1_pct2 <- (SC_FB_vs_FB_DEG$pct.1)/(SC_FB_vs_FB_DEG$pct.2)
SC_FB_vs_FB_DEG$celltype <- "SC-FB"
SC_FB_vs_FB_DEG$comparison <- "SC-FB_vs_FB"
write.xlsx(SC_FB_vs_FB_DEG, "Lists/DEG_SCFB_vs_FB.xlsx")

Idents(EC_FB_KC_compare)<- EC_FB_KC_compare$celltype
SC_EC_vs_EC_DEG <- FindMarkers(EC_FB_KC_compare, ident.1 = "SC-EC", ident.2=c("EC","LEC"))
SC_EC_vs_EC_DEG$Foldchange_UP <- 2^(SC_EC_vs_EC_DEG$avg_log2FC)
SC_EC_vs_EC_DEG$Foldchange_DOWN <- 2^(-SC_EC_vs_EC_DEG$avg_log2FC)
SC_EC_vs_EC_DEG$Ratio_pct1_pct2 <- (SC_EC_vs_EC_DEG$pct.1)/(SC_EC_vs_EC_DEG$pct.2)
SC_EC_vs_EC_DEG$celltype <- "SC-EC"
SC_EC_vs_EC_DEG$comparison <- "SC-EC_vs_EC"
write.xlsx(SC_EC_vs_EC_DEG, "Lists/DEG_SCEC_vs_EC.xlsx")

#KC
Idents(EC_FB_KC_compare)<- EC_FB_KC_compare$celltype
SC_FB_vs_KC_DEG <- FindMarkers(EC_FB_KC_compare, ident.1 = "SC-FB", ident.2=c("KC"))
SC_FB_vs_KC_DEG$Foldchange_UP <- 2^(SC_FB_vs_KC_DEG$avg_log2FC)
SC_FB_vs_KC_DEG$Foldchange_DOWN <- 2^(-SC_FB_vs_KC_DEG$avg_log2FC)
SC_FB_vs_KC_DEG$Ratio_pct1_pct2 <- (SC_FB_vs_KC_DEG$pct.1)/(SC_FB_vs_KC_DEG$pct.2)
SC_FB_vs_KC_DEG$celltype <- "SC-FB"
SC_FB_vs_KC_DEG$comparison <- "SC-FB_vs_KC"
write.xlsx(SC_FB_vs_KC_DEG, "Lists/DEG_SCFB_vs_KC.xlsx")

Idents(EC_FB_KC_compare)<- EC_FB_KC_compare$celltype
SC_EC_vs_KC_DEG <- FindMarkers(EC_FB_KC_compare, ident.1 = "SC-EC", ident.2=c("KC"))
SC_EC_vs_KC_DEG$Foldchange_UP <- 2^(SC_EC_vs_KC_DEG$avg_log2FC)
SC_EC_vs_KC_DEG$Foldchange_DOWN <- 2^(-SC_EC_vs_KC_DEG$avg_log2FC)
SC_EC_vs_KC_DEG$Ratio_pct1_pct2 <- (SC_EC_vs_KC_DEG$pct.1)/(SC_EC_vs_KC_DEG$pct.2)
SC_EC_vs_KC_DEG$celltype <- "SC-EC"
SC_EC_vs_KC_DEG$comparison <- "SC-EC_vs_KC"
write.xlsx(SC_EC_vs_KC_DEG, "Lists/DEG_SCEC_vs_KC.xlsx")

#SC
Idents(EC_FB_KC_compare)<- EC_FB_KC_compare$celltype
SC_FB_vs_SC_DEG <- FindMarkers(EC_FB_KC_compare, ident.1 = "SC-FB", ident.2=c("SC-Skin","SC-Promyel","SC-Repair","SC-Prolif","SC-EC"))
SC_FB_vs_SC_DEG$Foldchange_UP <- 2^(SC_FB_vs_SC_DEG$avg_log2FC)
SC_FB_vs_SC_DEG$Foldchange_DOWN <- 2^(-SC_FB_vs_SC_DEG$avg_log2FC)
SC_FB_vs_SC_DEG$Ratio_pct1_pct2 <- (SC_FB_vs_SC_DEG$pct.1)/(SC_FB_vs_SC_DEG$pct.2)
SC_FB_vs_SC_DEG$celltype <- "SC-FB"
SC_FB_vs_SC_DEG$comparison <- "SC-FB_vs_SC"
write.xlsx(SC_FB_vs_SC_DEG, "Lists/DEG_SCFB_vs_SC.xlsx")

Idents(EC_FB_KC_compare)<- EC_FB_KC_compare$celltype
SC_EC_vs_SC_DEG <- FindMarkers(EC_FB_KC_compare, ident.1 = "SC-EC", ident.2=c("SC-Skin","SC-Promyel","SC-Repair","SC-Prolif","SC-FB"))
SC_EC_vs_SC_DEG$Foldchange_UP <- 2^(SC_EC_vs_SC_DEG$avg_log2FC)
SC_EC_vs_SC_DEG$Foldchange_DOWN <- 2^(-SC_EC_vs_SC_DEG$avg_log2FC)
SC_EC_vs_SC_DEG$Ratio_pct1_pct2 <- (SC_EC_vs_SC_DEG$pct.1)/(SC_EC_vs_SC_DEG$pct.2)
SC_EC_vs_SC_DEG$celltype <- "SC-EC"
SC_EC_vs_SC_DEG$comparison <- "SC-EC_vs_SC"
write.xlsx(SC_EC_vs_SC_DEG, "Lists/DEG_SCEC_vs_SC.xlsx")

#combi
Idents(EC_FB_KC_compare)<- EC_FB_KC_compare$celltype
SC_FB_vs_SCFB_DEG <- FindMarkers(EC_FB_KC_compare, ident.1 = "SC-FB", ident.2=c("SC-Skin","SC-Promyel","SC-Repair","SC-Prolif","SC-EC","FB"))
SC_FB_vs_SCFB_DEG$Foldchange_UP <- 2^(SC_FB_vs_SCFB_DEG$avg_log2FC)
SC_FB_vs_SCFB_DEG$Foldchange_DOWN <- 2^(-SC_FB_vs_SCFB_DEG$avg_log2FC)
SC_FB_vs_SCFB_DEG$Ratio_pct1_pct2 <- (SC_FB_vs_SCFB_DEG$pct.1)/(SC_FB_vs_SCFB_DEG$pct.2)
SC_FB_vs_SCFB_DEG$celltype <- "SC-FB"
SC_FB_vs_SCFB_DEG$comparison <- "vs_SC+FB"
write.xlsx(SC_FB_vs_SCFB_DEG, "Lists/DEG_SCFB_vs_SC+FB.xlsx")

Idents(EC_FB_KC_compare)<- EC_FB_KC_compare$celltype
SC_EC_vs_SCEC_DEG <- FindMarkers(EC_FB_KC_compare, ident.1 = "SC-EC", ident.2=c("SC-Skin","SC-Promyel","SC-Repair","SC-Prolif","SC-FB","EC","LEC"))
SC_EC_vs_SCEC_DEG$Foldchange_UP <- 2^(SC_EC_vs_SCEC_DEG$avg_log2FC)
SC_EC_vs_SCEC_DEG$Foldchange_DOWN <- 2^(-SC_EC_vs_SCEC_DEG$avg_log2FC)
SC_EC_vs_SCEC_DEG$Ratio_pct1_pct2 <- (SC_EC_vs_SCEC_DEG$pct.1)/(SC_EC_vs_SCEC_DEG$pct.2)
SC_EC_vs_SCEC_DEG$celltype <- "SC-EC"
SC_EC_vs_SCEC_DEG$comparison <- "vs_SC+EC"
write.xlsx(SC_EC_vs_SCEC_DEG, "Lists/DEG_SCEC_vs_SC+EC.xlsx")

#Barplot Freq up and downregulated Genes (FC>2)
#EC
EC_SC_comp_df <- rbind(SC_EC_vs_EC_DEG,SC_EC_vs_SC_DEG,SC_EC_vs_KC_DEG)
EC_SC_comp_df_subset <- subset(EC_SC_comp_df, Foldchange_UP >=2 | Foldchange_DOWN >=2)
EC_SC_comp_df_subset$Direction <- ifelse (EC_SC_comp_df_subset$Foldchange_UP>=2,"UP","DOWN")
EC_SC_comp_df_subset<- select(EC_SC_comp_df_subset, Direction, comparison)
row.names(EC_SC_comp_df_subset)<-NULL
EC_SC_comp_df_subset$cluster<- as.factor(EC_SC_comp_df_subset$comparison)
EC_SC_comp_df_subset$cluster<- as.character(EC_SC_comp_df_subset$Direction)
EC_SC_comp_df_subset<- dplyr::count(EC_SC_comp_df_subset, comparison,Direction) %>% ungroup()
EC_SC_comp_df_subset$Direction <- factor(EC_SC_comp_df_subset$Direction, levels= c("UP","DOWN"))
EC_SC_comp_df_subset$comparison <- factor(EC_SC_comp_df_subset$comparison, levels= c("SC-EC_vs_SC","SC-EC_vs_EC","SC-EC_vs_KC"))
pdf("Graphs/Main/Barplot_nDEG_SC-EC_only.pdf")
ggplot(data=EC_SC_comp_df_subset, aes(x=comparison, y=n, fill=Direction)) +
  geom_bar(stat="identity", position=position_dodge())+ theme_classic() +
  geom_text_repel(aes(label = n, y=n+0.5),position=position_dodge(0.9), vjust=0.6, color = 'black', size = 4) +ggtitle("Comparison DEG_SC-EC") +scale_fill_manual(values =color_DEG) + xlab("SC-type")+ylab("total number of DEG")
dev.off()
#FB
FB_SC_comp_df <- rbind(SC_FB_vs_FB_DEG,SC_FB_vs_SC_DEG,SC_FB_vs_KC_DEG)
FB_SC_comp_df_subset <- subset(FB_SC_comp_df, Foldchange_UP >=2 | Foldchange_DOWN >=2)
FB_SC_comp_df_subset$Direction <- ifelse (FB_SC_comp_df_subset$Foldchange_UP>=2,"UP","DOWN")
FB_SC_comp_df_subset<- select(FB_SC_comp_df_subset, Direction, comparison)
row.names(FB_SC_comp_df_subset)<-NULL
FB_SC_comp_df_subset$cluster<- as.factor(FB_SC_comp_df_subset$comparison)
FB_SC_comp_df_subset$cluster<- as.character(FB_SC_comp_df_subset$Direction)
FB_SC_comp_df_subset<- dplyr::count(FB_SC_comp_df_subset, comparison,Direction) %>% ungroup()
FB_SC_comp_df_subset$Direction <- factor(FB_SC_comp_df_subset$Direction, levels= c("UP","DOWN"))
FB_SC_comp_df_subset$comparison <- factor(FB_SC_comp_df_subset$comparison, levels= c("SC-FB_vs_SC","SC-FB_vs_FB","SC-FB_vs_KC"))
pdf("Graphs/Main/Barplot_nDEG_SC-FB_only.pdf")
ggplot(data=FB_SC_comp_df_subset, aes(x=comparison, y=n, fill=Direction)) +
  geom_bar(stat="identity", position=position_dodge())+ theme_classic() +
  geom_text_repel(aes(label = n, y=n+0.5),position=position_dodge(0.9), vjust=0.6, color = 'black', size = 4) + ggtitle("Comparison DEG_SC-FB")  +scale_fill_manual(values =color_DEG) + xlab("SC-type")+ylab("total number of DEG")
dev.off()


### staining corresponding Markercheck
DefaultAssay(s.k.total)<-"RNA"
DefaultAssay(SC_7s4k)<-"RNA"

p1<- VlnPlot(s.k.total,"S100B", pt.size = 0, cols=color_tissue, split.by="tissue")
p1<- p1 + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p2<- VlnPlot(SC_7s4k,"S100B", pt.size = 0, cols=color_SC)
pdf("Graphs/Main/SC_VlnPlot_S100B.pdf")
wrap_plots(p1,p2, nrow=2)
dev.off()


p1<- VlnPlot(s.k.total,"NGFR", pt.size = 0, cols=color_tissue, split.by="tissue")
p1<- p1 + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p2<- VlnPlot(SC_7s4k,"NGFR", pt.size = 0, cols=color_SC)
pdf("Graphs/Main/SC_VlnPlot_NGFR.pdf")
wrap_plots(p1,p2, nrow=2)
dev.off()


p1<- VlnPlot(s.k.total,"NF1", pt.size = 0, cols=color_tissue, split.by="tissue")
p1<- p1 + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p2<- VlnPlot(SC_7s4k,"NF1", pt.size = 0, cols=color_SC)
pdf("Graphs/Main/SC_VlnPlot_NF1.pdf")
wrap_plots(p1,p2, nrow=2)
dev.off()


p1<- VlnPlot(s.k.total,"VIM", pt.size = 0, cols=color_tissue, split.by="tissue")
p1<- p1 + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p2<- VlnPlot(SC_7s4k,"VIM", pt.size = 0, cols=color_SC)
pdf("Graphs/Main/SC_VlnPlot_VIM.pdf")
wrap_plots(p1,p2, nrow=2)
dev.off()


p1<- VlnPlot(s.k.total,"NES", pt.size = 0, cols=color_tissue, split.by="tissue")
p1<- p1 + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p2<- VlnPlot(SC_7s4k,"NES", pt.size = 0, cols=color_SC)
pdf("Graphs/Main/SC_VlnPlot_NES.pdf")
wrap_plots(p1,p2, nrow=2)
dev.off()

p1<- VlnPlot(s.k.total,"SOX10", pt.size = 0, cols=color_tissue, split.by="tissue")
p1<- p1 + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p2<- VlnPlot(SC_7s4k,"SOX10", pt.size = 0, cols=color_SC)
pdf("Graphs/Main/SC_VlnPlot_SOX10.pdf")
wrap_plots(p1,p2, nrow=2)
dev.off()

p1<- VlnPlot(s.k.total,"ACTA2", pt.size = 0, cols=color_tissue, split.by="tissue")
p1<- p1 + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p2<- VlnPlot(SC_7s4k,"ACTA2", pt.size = 0, cols=color_SC)
pdf("Graphs/Main/SC_VlnPlot_ACTA2.pdf")
wrap_plots(p1,p2, nrow=2)
dev.off()

p1<- VlnPlot(s.k.total,"CDH19", pt.size = 0, cols=color_tissue, split.by="tissue")
p1<- p1 + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p2<- VlnPlot(SC_7s4k,"CDH19", pt.size = 0, cols=color_SC)
pdf("Graphs/Main/SC_VlnPlot_CDH19.pdf")
wrap_plots(p1,p2, nrow=2)
dev.off()

p1<- VlnPlot(s.k.total,"PECAM1", pt.size = 0, cols=color_tissue, split.by="tissue")
p1<- p1 + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p2<- VlnPlot(SC_7s4k,"PECAM1", pt.size = 0, cols=color_SC)
pdf("Graphs/Main/SC_VlnPlot_PECAM1.pdf")
wrap_plots(p1,p2, nrow=2)
dev.off()

p1<- VlnPlot(s.k.total,"MBP", pt.size = 0, cols=color_tissue, split.by="tissue")
p1<- p1 + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p2<- VlnPlot(SC_7s4k,"MBP", pt.size = 0, cols=color_SC)
pdf("Graphs/Main/SC_VlnPlot_MBP.pdf")
wrap_plots(p1,p2, nrow=2)
dev.off()

p1<- VlnPlot(s.k.total,"COL1A1", pt.size = 0, cols=color_tissue, split.by="tissue")
p1<- p1 + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p2<- VlnPlot(SC_7s4k,"COL1A1", pt.size = 0, cols=color_SC)
pdf("Graphs/Main/SC_VlnPlot_COL1A1.pdf")
wrap_plots(p1,p2, nrow=2)
dev.off()

p1<- VlnPlot(s.k.total,"IGFBP5", pt.size = 0, cols=color_tissue, split.by="tissue")
p1<- p1 + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p2<- VlnPlot(SC_7s4k,"IGFBP5", pt.size = 0, cols=color_SC)
pdf("Graphs/Main/SC_VlnPlot_IGFBP5.pdf")
wrap_plots(p1,p2, nrow=2)
dev.off()

p1<- VlnPlot(s.k.total,"NCAM1", pt.size = 0, cols=color_tissue, split.by="tissue")
p1<- p1 + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p2<- VlnPlot(SC_7s4k,"NCAM1", pt.size = 0, cols=color_SC)
pdf("Graphs/Main/SC_VlnPlot_NCAM1.pdf")
wrap_plots(p1,p2, nrow=2)
dev.off()

p1<- VlnPlot(s.k.total,"MKI67", pt.size = 0, cols=color_tissue, split.by="tissue")
p1<- p1 + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p2<- VlnPlot(SC_7s4k,"MKI67", pt.size = 0, cols=color_SC)
pdf("Graphs/Main/SC_VlnPlot_MKI67.pdf")
wrap_plots(p1,p2, nrow=2)
dev.off()

p1<- VlnPlot(s.k.total,"THY1", pt.size = 0, cols=color_tissue, split.by="tissue")
p1<- p1 + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p2<- VlnPlot(SC_7s4k,"THY1", pt.size = 0, cols=color_SC)
pdf("Graphs/Main/SC_VlnPlot_THY1.pdf")
wrap_plots(p1,p2, nrow=2)
dev.off()

p1<- VlnPlot(s.k.total,"JUN", pt.size = 0, cols=color_tissue, split.by="tissue")
p1<- p1 + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p2<- VlnPlot(SC_7s4k,"JUN", pt.size = 0, cols=color_SC)
pdf("Graphs/Main/SC_VlnPlot_JUN.pdf")
wrap_plots(p1,p2, nrow=2)
dev.off()

p1<- VlnPlot(s.k.total,c("IGFBP2","IGFBP3","IGFBP4","IGFBP5","IGFBP6","IGFBP7"), pt.size = 0, cols=color_UMAP_short_double, ncol = 6)
p1<- p1 & theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p2<- VlnPlot(SC_7s4k,c("IGFBP2","IGFBP3","IGFBP4","IGFBP5","IGFBP6","IGFBP7"), pt.size = 0, cols=color_SC, ncol = 6)
pdf("Graphs/Main/SC_VlnPlot_IGFBPs.pdf")
wrap_plots(p1,p2, nrow=2)
dev.off()

p1<- VlnPlot(s.k.total,c("CCL2"), pt.size = 0, cols=color_tissue, split.by="tissue")+stat_summary(fun = "mean", geom = "crossbar")
p1<- p1 & theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p2<- VlnPlot(SC_7s4k,c("CCL2"), pt.size = 0, cols=color_SC)+stat_summary(fun = "mean", geom = "crossbar")
pdf("Graphs/Main/SC_VlnPlot_CCL2+mean.pdf")
wrap_plots(p1,p2, nrow=2)
dev.off()

p1<- VlnPlot(s.k.total,c("JUNB"), pt.size = 0, cols=color_tissue, split.by="tissue")
p1<- p1 & theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p2<- VlnPlot(SC_7s4k,c("JUNB"), pt.size = 0, cols=color_SC)
pdf("Graphs/Main/SC_VlnPlot_JUNB.pdf")
wrap_plots(p1,p2, nrow=2)
dev.off()

p1<- VlnPlot(s.k.total,c("JUND"), pt.size = 0, cols=color_tissue, split.by="tissue")
p1<- p1 & theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p2<- VlnPlot(SC_7s4k,c("JUND"), pt.size = 0, cols=color_SC)
pdf("Graphs/Main/SC_VlnPlot_JUND.pdf")
wrap_plots(p1,p2, nrow=2)
dev.off()

p1<- VlnPlot(s.k.total,c("FOSL1"), pt.size = 0, cols=color_tissue, split.by="tissue")
p1<- p1 & theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p2<- VlnPlot(SC_7s4k,c("FOSL1"), pt.size = 0, cols=color_SC)
pdf("Graphs/Main/SC_VlnPlot_FOSL1.pdf")
wrap_plots(p1,p2, nrow=2)
dev.off()

p1<- VlnPlot(s.k.total,c("FOSL2"), pt.size = 0, cols=color_tissue, split.by="tissue")
p1<- p1 & theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p2<- VlnPlot(SC_7s4k,c("FOSL2"), pt.size = 0, cols=color_SC)
pdf("Graphs/Main/SC_VlnPlot_FOSL2.pdf")
wrap_plots(p1,p2, nrow=2)
dev.off()

p1<- VlnPlot(s.k.total,c("FOS"), pt.size = 0, cols=color_tissue, split.by="tissue")
p1<- p1 & theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p2<- VlnPlot(SC_7s4k,c("FOS"), pt.size = 0, cols=color_SC)
pdf("Graphs/Main/SC_VlnPlot_FOS.pdf")
wrap_plots(p1,p2, nrow=2)
dev.off()

p1<- VlnPlot(s.k.total,c("TNFAIP6"), pt.size = 0, cols=color_tissue, split.by="tissue")+stat_summary(fun = "mean", geom = "crossbar")
p1<- p1 & theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p2<- VlnPlot(SC_7s4k,c("TNFAIP6"), pt.size = 0, cols=color_SC)+stat_summary(fun = "mean", geom = "crossbar")
pdf("Graphs/Main/SC_VlnPlot_TNFAIP6+mean.pdf")
wrap_plots(p1,p2, nrow=2)
dev.off()

p1<- VlnPlot(s.k.total,c("CCN3"), pt.size = 0, cols=color_tissue, split.by="tissue")+stat_summary(fun = "mean", geom = "crossbar")
p1<- p1 & theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p2<- VlnPlot(SC_7s4k,c("CCN3"), pt.size = 0, cols=color_SC)+stat_summary(fun = "mean", geom = "crossbar")
pdf("Graphs/Main/SC_VlnPlot_CCN3+mean.pdf")
wrap_plots(p1,p2, nrow=2)
dev.off()

p1<- VlnPlot(s.k.total,c("IGFBP5"), pt.size = 0, cols=color_tissue, split.by="tissue")+stat_summary(fun = "mean", geom = "crossbar")
p1<- p1 & theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p2<- VlnPlot(SC_7s4k,c("IGFBP5"), pt.size = 0, cols=color_SC)+stat_summary(fun = "mean", geom = "crossbar")
pdf("Graphs/Main/SC_VlnPlot_IGFBP5+mean.pdf")
wrap_plots(p1,p2, nrow=2)
dev.off()

p1<- VlnPlot(s.k.total,c("CRYAB"), pt.size = 0, cols=color_tissue, split.by="tissue")
p1<- p1 & theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p2<- VlnPlot(SC_7s4k,c("CRYAB"), pt.size = 0, cols=color_SC)
pdf("Graphs/Main/SC_VlnPlot_CRYAB.pdf")
wrap_plots(p1,p2, nrow=2)
dev.off()

p1<- VlnPlot(s.k.total,c("IGFBP3"), pt.size = 0, cols=color_tissue, split.by="tissue")+stat_summary(fun = "mean", geom = "crossbar")
p1<- p1 & theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p2<- VlnPlot(SC_7s4k,c("IGFBP3"), pt.size = 0, cols=color_SC)+stat_summary(fun = "mean", geom = "crossbar")
pdf("Graphs/Main/SC_VlnPlot_IGFBP3+mean.pdf")
wrap_plots(p1,p2, nrow=2)
dev.off()

p1<- VlnPlot(s.k.total,c("IGFBP7"), pt.size = 0, cols=color_tissue, split.by="tissue")+stat_summary(fun = "mean", geom = "crossbar")
p1<- p1 & theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p2<- VlnPlot(SC_7s4k,c("IGFBP7"), pt.size = 0, cols=color_SC)+stat_summary(fun = "mean", geom = "crossbar")
pdf("Graphs/Main/SC_VlnPlot_IGFBP7+mean.pdf")
wrap_plots(p1,p2, nrow=2)
dev.off()

# Interaction Plot SCs
pdf("Graphs/Main/SC_VlnPlot_Interactionplot.pdf", width = 15, height = 5)
VlnPlot(SC_7s4k,features=c("TNFAIP6","CCN3","IGFBP5"), group.by="tissue",pt.size = 0, cols=color_tissue)& ylim(0,6)
dev.off()

### Gene check
## Matrix gene based on PMID:26163349
# Website:http://matrisomeproject.mit.edu/other-resources/human-matrisome/

Idents(s.k.total)<-s.k.total$celltype.tissue

#Core matrisome
ECM_Glycoproteins <- c("ABI3BP","ADIPOQ","AEBP1","AGRN","AMBN","AMELX","AMELY","BGLAP","BMPER","BSPH1","CDCP2","CILP","CILP2","COCH","COLQ","COMP","CRELD1","CRELD2","CRIM1","CRISPLD1","CRISPLD2","CTGF","CTHRC1","CYR61","DDX26B","DMBT1","DMP1","DPT","DSPP","ECM1","ECM2","EDIL3","EFEMP1","EFEMP2","EGFLAM","ELN","ELSPBP1","EMID1","EMILIN1","EMILIN2","EMILIN3","EYS","FBLN1","FBLN2","FBLN5","FBLN7","FBN1","FBN2","FBN3","FGA","FGB","FGG","FGL1","FGL2","FN1","FNDC1","FNDC7","FNDC8","FRAS1","GAS6","GLDN","HMCN1","HMCN2","IBSP","IGFALS","IGFBP1","IGFBP2","IGFBP3","IGFBP4","IGFBP5","IGFBP6","IGFBP7","IGFBPL1","IGSF10","ANOS1","KCP","LAMA1","LAMA2","LAMA3","LAMA4","LAMA5","LAMB1","LAMB2","LAMB3","LAMB4","LAMC1","LAMC2","LAMC3","LGI1","LGI2","LGI3","LGI4","LRG1","LTBP1","LTBP2","LTBP3","LTBP4","MATN1","MATN2","MATN3","MATN4","MEPE","MFAP1","MFAP2","MFAP3","MFAP4","MFAP5","MFGE8","MGP","MMRN1","MMRN2","MXRA5","NDNF","NELL1","NELL2","NID1","NID2","NOV","NPNT","NTN1","NTN3","NTN4","NTN5","NTNG1","NTNG2","OIT3","OTOG","OTOL1","PAPLN","PCOLCE","PCOLCE2","POMZP3","POSTN","PXDN","PXDNL","RELN","RSPO1","RSPO2","RSPO3","RSPO4","SBSPON","SLIT1","SLIT2","SLIT3","SMOC1","SMOC2","SNED1","SPARC","SPARCL1","SPON1","SPON2","SPP1","SRPX","SRPX2","SSPO","SVEP1","TECTA","TECTB","TGFBI","THBS1","THBS2","THBS3","THBS4","THSD4","TINAG","TINAGL1","TNC","TNFAIP6","TNN","TNR","TNXB","TSKU","TSPEAR","VIT","VTN","VWA1","VWA2","VWA3A","VWA3B","VWA5A","VWA5B1","VWA5B2","VWA7","VWA9","VWCE","VWDE","VWF","WISP1","WISP2","WISP3","ZP1","ZP2","ZP3","ZP4","ZPLD1")
Collagens <- c("COL1A1","COL1A2","COL2A1","COL3A1","COL4A1","COL4A2","COL4A3","COL4A4","COL4A5","COL4A6","COL5A1","COL5A2","COL5A3","COL6A1","COL6A2","COL6A3","COL6A5","COL6A6","COL7A1","COL8A1","COL8A2","COL9A1","COL9A2","COL9A3","COL10A1","COL11A1","COL11A2","COL12A1","COL13A1","COL14A1","COL15A1","COL16A1","COL17A1","COL18A1","COL19A1","COL20A1","COL21A1","COL22A1","COL23A1","COL24A1","COL25A1","COL26A1","COL27A1","COL28A1")
Proteoglycans <- c("ACAN","ASPN","BCAN","BGN","CHAD","CHADL","DCN","EPYC","ESM1","FMOD","HAPLN1","HAPLN2","HAPLN3","HAPLN4","HSPG2","IMPG1","IMPG2","KERA","LUM","NCAN","NYX","OGN","OMD","OPTC","PODN","PODNL1","PRELP","PRG2","PRG3","PRG4","SPOCK1","SPOCK2","SPOCK3","SRGN","VCAN")

#Matrisome Associated
ECM_affiliated_Proteins<- c("ANXA1","ANXA10","ANXA11","ANXA13","ANXA2","ANXA3","ANXA4","ANXA5","ANXA6","ANXA7","ANXA8","ANXA8L1","ANXA9","C1QA","C1QB","C1QC","C1QL1","C1QL2","C1QL3","C1QL4","C1QTNF1","C1QTNF2","C1QTNF3","C1QTNF4","C1QTNF5","C1QTNF6","C1QTNF7","C1QTNF8","C1QTNF9","CD209","CLC","CLEC10A","CLEC11A","CLEC12A","CLEC12B","CLEC14A","CLEC17A","CLEC18A","CLEC18B","CLEC18C","CLEC19A","CLEC1A","CLEC1B","CLEC2A","CLEC2B","CLEC2D","CLEC2L","CLEC3A","CLEC3B","CLEC4A","CLEC4C","CLEC4D","CLEC4E","CLEC4F","CLEC4G","CLEC4M","CLEC5A","CLEC6A","CLEC7A","CLEC9A","COLEC10","COLEC11","COLEC12","CSPG4","CSPG5","ELFN1","ELFN2","EMCN","FCN1","FCN2","FCN3","FREM1","FREM2","FREM3","GPC1","GPC2","GPC3","GPC4","GPC5","GPC6","GREM1","GRIFIN","HPX","LGALSL","ITLN1","ITLN2","LGALS1","LGALS12","LGALS13","LGALS14","LGALS16","LGALS2","LGALS3","LGALS4","LGALS7","LGALS8","LGALS9","LGALS9B","LGALS9C","LMAN1","LMAN1L","MBL2","MUC1","MUC12","MUC13","MUC15","MUC16","MUC17","MUC19","MUC2","MUC20","MUC21","MUC22","MUC3A","MUC4","MUC5AC","MUC5B","MUC6","MUC7","MUC8","MUCL1","OVGP1","PARM1","PLXDC1","PLXDC2","PLXNA1","PLXNA2","PLXNA3","PLXNA4","PLXNB1","PLXNB2","PLXNB3","PLXNC1","PLXND1","PROL1","REG1A","REG1B","REG3A","REG3G","REG4","SDC1","SDC2","SDC3","SDC4","SEMA3A","SEMA3B","SEMA3C","SEMA3D","SEMA3E","SEMA3F","SEMA3G","SEMA4A","SEMA4B","SEMA4C","SEMA4D","SEMA4F","SEMA4G","SEMA5A","SEMA5B","SEMA6A","SEMA6B","SEMA6C","SEMA6D","SEMA7A","SFTA2","SFTA3","SFTPA1","SFTPA2","SFTPB","SFTPC","SFTPD")
ECM_Regulators <- c("A2M","A2ML1","ADAM10","ADAM11","ADAM12","ADAM15","ADAM17","ADAM18","ADAM19","ADAM2","ADAM20","ADAM21","ADAM22","ADAM23","ADAM28","ADAM29","ADAM30","ADAM32","ADAM33","ADAM7","ADAM8","ADAM9","ADAMDEC1","ADAMTS1","ADAMTS10","ADAMTS12","ADAMTS13","ADAMTS14","ADAMTS15","ADAMTS16","ADAMTS17","ADAMTS18","ADAMTS19","ADAMTS2","ADAMTS20","ADAMTS3","ADAMTS4","ADAMTS5","ADAMTS6","ADAMTS7","ADAMTS8","ADAMTS9","ADAMTSL1","ADAMTSL2","ADAMTSL3","ADAMTSL4","ADAMTSL5","AGT","AMBP","ASTL","BMP1","C17orf58","CD109","CELA1","CELA2A","CELA2B","CELA3A","CELA3B","CPAMD8","CPN2","CST1","CST11","CST2","CST3","CST4","CST5","CST6","CST7","CST8","CST9","CST9L","CSTA","CSTB","CSTL1","CTSA","CTSB","CTSC","CTSD","CTSE","CTSF","CTSG","CTSH","CTSK","CTSL","CTSO","CTSS","CTSV","CTSW","CTSZ","EGLN1","EGLN2","EGLN3","ELANE","F10","F12","F13A1","F13B","F2","F7","F9","FAM20A","FAM20B","FAM20C","HABP2","HMSD","HPSE","HPSE2","HRG","HTRA1","HTRA3","HTRA4","HYAL1","HYAL2","HYAL3","HYAL4","ITIH1","ITIH2","ITIH3","ITIH4","ITIH5","ITIH6","KAZALD1","KNG1","KY","P3H1","P3H2","P3H3","LOX","LOXL1","LOXL2","LOXL3","LOXL4","LPA","MASP1","MASP2","MEP1A","MEP1B","MMP1","MMP10","MMP11","MMP12","MMP13","MMP14","MMP15","MMP16","MMP17","MMP19","MMP2","MMP20","MMP21","MMP23B","MMP24","MMP25","MMP26","MMP27","MMP28","MMP3","MMP7","MMP8","MMP9","NGLY1","OGFOD1","OGFOD2","P4HA1","P4HA2","P4HA3","P4HTM","PAMR1","PAPPA","PAPPA2","PCSK5","PCSK6","PI3","PLAT","PLAU","PLG","PLOD1","PLOD2","PLOD3","PRSS1","PRSS12","PRSS2","PRSS3","PZP","SERPINA1","SERPINA10","SERPINA11","SERPINA12","SERPINA2","SERPINA3","SERPINA4","SERPINA5","SERPINA6","SERPINA7","SERPINA9","SERPINB1","SERPINB10","SERPINB11","SERPINB12","SERPINB13","SERPINB2","SERPINB3","SERPINB4","SERPINB5","SERPINB6","SERPINB7","SERPINB8","SERPINB9","SERPINC1","SERPIND1","SERPINE1","SERPINE2","SERPINE3","SERPINF1","SERPINF2","SERPING1","SERPINH1","SERPINI1","SERPINI2","SLPI","SPAM1","ST14","SULF1","SULF2","TGM1","TGM2","TGM3","TGM4","TGM5","TGM6","TGM7","TIMP1","TIMP2","TIMP3","TIMP4","TLL1","TLL2","TMPRSS15")
Secreted_Factors <- c("AMH","ANGPT1","ANGPT2","ANGPT4","ANGPTL1","ANGPTL2","ANGPTL3","ANGPTL4","ANGPTL5","ANGPTL6","ANGPTL7","AREG","ARTN","BDNF","BMP10","BMP15","BMP2","BMP3","BMP4","BMP5","BMP6","BMP7","BMP8A","BMP8B","BRINP2","BRINP3","BTC","C1QTNF9B","CBLN1","CBLN2","CBLN3","CBLN4","CCBE1","CCL1","CCL11","CCL13","CCL14","CCL15","CCL16","CCL17","CCL18","CCL19","CCL2","CCL20","CCL21","CCL22","CCL23","CCL24","CCL25","CCL26","CCL27","CCL28","CCL3","CCL3L3","CCL4","CCL4L1","CCL4L2","CCL5","CCL7","CCL8","CFC1","CFC1B","CHRD","CHRDL1","CHRDL2","CLCF1","CNTF","CRHBP","CRLF1","CRLF3","CRNN","CSF1","CSF2","CSF3","CSH1","CSH2","CSHL1","CTF1","CX3CL1","CXCL1","CXCL10","CXCL11","CXCL12","CXCL13","CXCL14","CXCL2","CXCL3","CXCL5","CXCL6","CXCL8","CXCL9","DHH","EBI3","EDA","EGF","EGFL6","EGFL7","EGFL8","EPGN","EPO","EREG","FASLG","FGF1","FGF10","FGF11","FGF12","FGF13","FGF14","FGF16","FGF17","FGF18","FGF19","FGF2","FGF20","FGF21","FGF22","FGF23","FGF3","FGF4","FGF5","FGF6","FGF7","FGF8","FGF9","FGFBP1","FGFBP2","FGFBP3","FIGF","FLG","FLG2","FLT3LG","FRZB","FST","FSTL1","FSTL3","GDF1","GDF10","GDF11","GDF15","GDF2","GDF3","GDF5","GDF6","GDF7","GDF9","GDNF","GH1","GH2","HBEGF","HCFC1","HCFC2","HGF","HGFAC","HHIP","HRNR","IFNA1","IFNA10","IFNA13","IFNA14","IFNA16","IFNA17","IFNA2","IFNA21","IFNA4","IFNA5","IFNA6","IFNA7","IFNA8","IFNB1","IFNE","IFNG","IFNK","IFNW1","IGF1","IGF2","IHH","IL10","IL11","IL12A","IL12B","IL13","IL15","IL16","IL17A","IL17B","IL17C","IL17D","IL17F","IL18","IL19","IL1A","IL1B","IL1F10","IL36RN","IL36A","IL37","IL36B","IL1F9","IL1RN","IL2","IL20","IL22","IL23A","IL24","IL25","IL26","IL3","IL34","IL4","IL5","IL6","IL7","IL9","INHA","INHBA","INHBB","INHBC","INHBE","INS","INS-IGF2","INSL3","INSL5","INSL6","ISM1","ISM2","KITLG","LEFTY1","LEFTY2","LEP","LIF","LTA","LTB","MDK","MEGF10","MEGF11","MEGF6","MEGF8","MEGF9","MST1","MSTN","NGF","NODAL","NRG1","NRG2","NRG3","NRG4","NRTN","NTF3","NTF4","OSM","PDGFA","PDGFB","PDGFC","PDGFD","PF4","PF4V1","PGF","PIK3IP1","PPBP","PRL","PSPN","PTN","RPTN","S100A1","S100A10","S100A11","S100A12","S100A13","S100A14","S100A16","S100A2","S100A3","S100A4","S100A5","S100A6","S100A7","S100A7A","S100A7L2","S100A8","S100A9","S100B","S100G","S100P","S100Z","SCUBE1","SCUBE2","SCUBE3","SFRP1","SFRP2","SFRP4","SFRP5","SHH","TCHH","TCHHL1","TDGF1","TGFA","TGFB1","TGFB2","TGFB3","THPO","TNF","TNFSF10","TNFSF11","TNFSF12","TNFSF13","TNFSF13B","TNFSF14","TNFSF15","TNFSF18","TNFSF4","TNFSF8","TNFSF9","TPO","VEGFA","VEGFB","VEGFC","VWC2","VWC2L","WFIKKN1","WFIKKN2","WIF1","WNT1","WNT10A","WNT10B","WNT11","WNT16","WNT2","WNT2B","WNT3","WNT3A","WNT4","WNT5A","WNT5B","WNT6","WNT7A","WNT7B","WNT8A","WNT8B","WNT9A","WNT9B","XCL1","XCL2","ZFP91")



pdf("Graphs/Side/Dotplot_ECM_Glycoproteins.pdf", width = 10, height = 40)
DotPlot(s.k.total,features=ECM_Glycoproteins)+scale_color_gradient(low="grey", "high"=single_col) + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short_double))+coord_flip() +theme(axis.text.x = element_text(angle = 90)) + ggtitle("ECM_Glycoproteins")
dev.off()

pdf("Graphs/Side/Dotplot_Collagens.pdf", width = 10, height = 10)
DotPlot(s.k.total,features=Collagens)+scale_color_gradient(low="grey", "high"=single_col) + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short_double))+coord_flip() +theme(axis.text.x = element_text(angle = 90)) + ggtitle("Collagens")
dev.off()

pdf("Graphs/Side/Dotplot_Proteoglycans.pdf", width = 10, height = 8)
DotPlot(s.k.total,features=Proteoglycans)+scale_color_gradient(low="grey", "high"=single_col) + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short_double))+coord_flip() +theme(axis.text.x = element_text(angle = 90)) + ggtitle("Proteoglycans")
dev.off()

pdf("Graphs/Side/Dotplot_ECM_affiliated_Proteins.pdf", width = 10, height = 40)
DotPlot(s.k.total,features=ECM_affiliated_Proteins)+scale_color_gradient(low="grey", "high"=single_col) + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short_double))+coord_flip() +theme(axis.text.x = element_text(angle = 90)) + ggtitle("ECM_affiliated_Proteins")
dev.off()

pdf("Graphs/Side/Dotplot_ECM_Regulators.pdf", width = 10, height = 40)
DotPlot(s.k.total,features=ECM_Regulators)+scale_color_gradient(low="grey", "high"=single_col) + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short_double))+coord_flip() +theme(axis.text.x = element_text(angle = 90)) + ggtitle("ECM_Regulators")
dev.off()

pdf("Graphs/Side/Dotplot_Secreted_Factors.pdf", width = 10, height = 60)
DotPlot(s.k.total,features=Secreted_Factors)+scale_color_gradient(low="grey", "high"=single_col) + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short_double))+coord_flip() +theme(axis.text.x = element_text(angle = 90)) + ggtitle("Secreted_Factors")
dev.off()

#Top Matrisome long list
Idents(s.k.total)<-s.k.total$celltype.tissue
SC_matrisome<-c("AGRN","CTHRC1","ELN","FBLN2","FGL2","FN1","GAS6","GLDN","HMCN1","HMCN2","IBSP","IGFBP2","IGFBP3","IGFBP5","IGFBP6","IGFBP7","LAMA2","LAMA4","LAMC1","LGI4","LTBP3","LTBP4","MATN2","MFAP1","MFGE8","MXRA5","NID1","NID2","NPNT","POSTN","PXDN","SPARC","SPARCL1","SRPX2","TGFBI","THBS1","THBS2","TNC","TNFAIP6","VWA1","ANXA2","ANXA4","ANXA3","ANXA5","ANXA6","ANXA7","LGALS3","LMAN1","SDC3","SEMA3C","CCL2","CRLF1","FGFBP2","FSTL3","HBEGF","HCFC2","PDGFA","PDGFB","PDGFC","PIK3IP1","PTN","S100A10","S100A11","S100A13","S100A16","S100A4","S100A6","S100B","CTSF","CTSA","CTSC","CTSD","CTSL","EGLN2","ITIH5","P3H1","P3H2","LOXL1","LOXL2","LOXL3","MMP1","MMP10","MMP11","MMP12","MMP13","MMP14","MMP15","MMP16","MMP17","MMP19","MMP2","MMP20","MMP21","MMP23B","MMP24","MMP25","MMP26","MMP27","MMP28","MMP3","MMP7","MMP8","MMP9","OGFOD1","P4HA2","PLAT","PLOD1","SERPINB1","SERPINB6","SERPINE1","SERPINE2","SERPING1","SERPINH1","TIMP1","TIMP2","TIMP3","TIMP4","COL1A1","COL1A2","COL2A1","COL3A1","COL4A1","COL4A2","COL4A3","COL4A4","COL4A5","COL4A6","COL5A1","COL5A2","COL5A3","COL6A1","COL6A2","COL6A3","COL6A5","COL6A6","COL7A1","COL8A1","COL8A2","COL9A1","COL9A2","COL9A3","COL10A1","COL11A1","COL11A2","COL12A1","COL13A1","COL14A1","COL15A1","COL16A1","COL17A1","COL18A1","COL19A1","COL20A1","COL21A1","COL22A1","COL23A1","COL24A1","COL25A1","COL26A1","COL27A1","COL28A1","CHADL","PRELP")
SC_matrisome <- sort(SC_matrisome)
#save to order alphabetical 
#write.xlsx(SC_matrisome,"Lists/Matrisome_list.xlsx")
#load ordered list
SC_matrisome<- read.xlsx("Lists/Matrisome_list.xlsx",1,header = T)
SC_matrisome<-SC_matrisome[,-1]

pdf("Graphs/Main/Dotplot_SC_matrisome_genes.pdf", width = 10, height = 35)
DotPlot(s.k.total,features=SC_matrisome)+scale_color_gradient(low="grey", "high"=single_col) + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short_double))+coord_flip() +theme(axis.text.x = element_text(angle = 90)) + ggtitle("SC_Matrisome_genes")
dev.off()

#Short list
SC_matrisome_shortlist <- c("COL1A1","COL3A1","COL4A1","COL4A2","COL5A1","COL7A1","COL8A1","COL11A1","COL12A1","COL18A1","ELN","IBSP","IGFBP2","IGFBP3","IGFBP5","ITIH5","MMP15","CCN3","TGFBI","TNC")
pdf("Graphs/Main/Dotplot_SC_short_listmatrisome.pdf", width = 8, height = 8)
DotPlot(s.k.total,features=SC_matrisome_shortlist)+scale_color_gradient(low="grey", "high"=single_col) + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short_double))+coord_flip() +theme(axis.text.x = element_text(angle = 90)) + ggtitle("SC_Matrisome_genes")
dev.off()

#Vln_Matrix_tissue
Idents(SC_7s4k)<-"tissue"
SC_matrisome_Vlnlist <- c("CCN3","COL1A1","COL3A1","COL4A1","COL4A2","COL5A1","COL5A2","COL7A1","COL8A1","COL12A1","COL18A1","ELN","IGFBP3","IGFBP5","TNC","TGFBI")
pdf("Graphs/Main/Vlnplot_SC_matrisome.pdf", width=15, height=5)
VlnPlot(SC_7s4k,y=10,features=SC_matrisome_Vlnlist, ncol = 8, pt.size = 0,cols=color_tissue)& stat_summary(fun = "mean", geom = "crossbar")& stat_compare_means(comparisons = list(c("skin","keloid")) ,method = "t.test", label = "p.signif")
dev.off()

Idents(s.k.total)<-"tissue"
pdf("Graphs/Main/Vlnplot_total_matrisome.pdf", width=15, height=5)
VlnPlot(s.k.total,features=SC_matrisome_Vlnlist, ncol = 8, pt.size = 0,cols=color_tissue)& stat_summary(fun = "mean", geom = "crossbar")& stat_compare_means(comparisons = list(c("skin","keloid")) ,method = "t.test", label = "p.signif")
dev.off()

### Gene check
## Inflammation genes based on PMID:17940599
# Website:https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0001035#s5
Inflamm_List <- read.xlsx("D:/Direder/Projekte/D002_10xScarWars/D002_Berechnungen/D002_c6_skin&normscar&keloid_3s3n4k_si_SCT/3s4k/UpdateSeurat4R4/Inflammation Lists 17940599/Inflamm_genelist_adapted.xlsx",1, header = T)

Adhesion_Extravasation_Migration <- Inflamm_List[which(Inflamm_List$Primary.Pathway =="Adhesion-Extravasation-Migration"),]
Apoptosis_Signaling<- Inflamm_List[which(Inflamm_List$Primary.Pathway =="Apoptosis Signaling"),]
Calcium_Signaling<- Inflamm_List[which(Inflamm_List$Primary.Pathway =="Calcium Signaling"),]
Complement_Cascade<- Inflamm_List[which(Inflamm_List$Primary.Pathway =="Complement Cascase"),]
Cytokine_Signaling<- Inflamm_List[which(Inflamm_List$Primary.Pathway =="Cytokine signaling"),]
Eicosanoid_Signaling<- Inflamm_List[which(Inflamm_List$Primary.Pathway =="Eicosanoid Signaling"),]
Glucocorticoid_PPAR_Signaling<- Inflamm_List[which(Inflamm_List$Primary.Pathway =="Glucocorticoid/PPAR signaling"),]
G_Protein_Coupled_Receptor_Signaling<- Inflamm_List[which(Inflamm_List$Primary.Pathway =="G-Protein Coupled Receptor Signaling"),]
Innate_pathogen_detection<- Inflamm_List[which(Inflamm_List$Primary.Pathway =="Innate pathogen detection"),]
Leukocyte_signaling<- Inflamm_List[which(Inflamm_List$Primary.Pathway =="Leukocyte signaling"),]
MAPK_signaling<- Inflamm_List[which(Inflamm_List$Primary.Pathway =="MAPK signaling"),]
Natural_Killer_Cell_Signaling<- Inflamm_List[which(Inflamm_List$Primary.Pathway =="Natural Killer Cell Signaling"),]
NF_kB_signaling<- Inflamm_List[which(Inflamm_List$Primary.Pathway =="NF-kB signaling"),]
Phagocytosis_Ag_presentation<- Inflamm_List[which(Inflamm_List$Primary.Pathway =="Phagocytosis-Ag presentation"),]
PI3K_AKT_Signaling<- Inflamm_List[which(Inflamm_List$Primary.Pathway =="PI3K/AKT Signaling"),]
ROS_Glutathione_Cytotoxic_granules<- Inflamm_List[which(Inflamm_List$Primary.Pathway =="ROS/Glutathione/Cytotoxic granules"),]
TNF_Superfamily_Signaling<- Inflamm_List[which(Inflamm_List$Primary.Pathway =="TNF Superfamily Signaling"),]

pdf("Graphs/Side/Dotplot_Adhesion_Extravasation_Migration.pdf", width = 10, height = 30)
DotPlot(s.k.total,features=Adhesion_Extravasation_Migration$geneID)+scale_color_gradient(low="grey", "high"=single_col) + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short_double))+coord_flip() +theme(axis.text.x = element_text(angle = 90)) + ggtitle("Adhesion_Extravasation_Migration")
dev.off()
pdf("Graphs/Side/Dotplot_Apoptosis_Signaling.pdf", width = 10, height = 20)
DotPlot(s.k.total,features=Apoptosis_Signaling$geneID)+scale_color_gradient(low="grey", "high"=single_col) + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short_double))+coord_flip() +theme(axis.text.x = element_text(angle = 90)) + ggtitle("Apoptosis_Signaling")
dev.off()
pdf("Graphs/Side/Dotplot_Calcium_Signaling.pdf", width = 10, height = 5)
DotPlot(s.k.total,features=Calcium_Signaling$geneID)+scale_color_gradient(low="grey", "high"=single_col) + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short_double))+coord_flip() +theme(axis.text.x = element_text(angle = 90)) + ggtitle("Calcium_Signaling")
dev.off()
pdf("Graphs/Side/Dotplot_Complement_Cascade.pdf", width = 10, height = 10)
DotPlot(s.k.total,features=Complement_Cascade$geneID)+scale_color_gradient(low="grey", "high"=single_col) + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short_double))+coord_flip() +theme(axis.text.x = element_text(angle = 90)) + ggtitle("Complement_Cascade")
dev.off()
pdf("Graphs/Side/Dotplot_Cytokine_Signaling.pdf", width = 10, height = 30)
DotPlot(s.k.total,features=Cytokine_Signaling$geneID)+scale_color_gradient(low="grey", "high"=single_col) + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short_double))+coord_flip() +theme(axis.text.x = element_text(angle = 90)) + ggtitle("Cytokine_Signaling")
dev.off()
pdf("Graphs/Side/Dotplot_Eicosanoid_Signaling.pdf", width = 10, height = 10)
DotPlot(s.k.total,features=Eicosanoid_Signaling$geneID)+scale_color_gradient(low="grey", "high"=single_col) + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short_double))+coord_flip() +theme(axis.text.x = element_text(angle = 90)) + ggtitle("Eicosanoid_Signaling")
dev.off()
pdf("Graphs/Side/Dotplot_Glucocorticoid_PPAR_Signaling.pdf", width = 10, height = 5)
DotPlot(s.k.total,features=Glucocorticoid_PPAR_Signaling$geneID)+scale_color_gradient(low="grey", "high"=single_col) + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short_double))+coord_flip() +theme(axis.text.x = element_text(angle = 90)) + ggtitle("Glucocorticoid_PPAR_Signaling")
dev.off()
pdf("Graphs/Side/Dotplot_G_Protein_Coupled_Receptor_Signaling.pdf", width = 10, height = 12)
DotPlot(s.k.total,features=G_Protein_Coupled_Receptor_Signaling$geneID)+scale_color_gradient(low="grey", "high"=single_col) + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short_double))+coord_flip() +theme(axis.text.x = element_text(angle = 90)) + ggtitle("G_Protein_Coupled_Receptor_Signaling")
dev.off()
pdf("Graphs/Side/Dotplot_Innate_pathogen_detection.pdf", width = 10, height = 10)
DotPlot(s.k.total,features=Innate_pathogen_detection$geneID)+scale_color_gradient(low="grey", "high"=single_col) + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short_double))+coord_flip() +theme(axis.text.x = element_text(angle = 90)) + ggtitle("Innate_pathogen_detection")
dev.off()
pdf("Graphs/Side/Dotplot_Leukocyte_signaling.pdf", width = 10, height = 25)
DotPlot(s.k.total,features=Leukocyte_signaling$geneID)+scale_color_gradient(low="grey", "high"=single_col) + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short_double))+coord_flip() +theme(axis.text.x = element_text(angle = 90)) + ggtitle("Leukocyte_signaling")
dev.off()
pdf("Graphs/Side/Dotplot_MAPK_signaling.pdf", width = 10, height = 30)
DotPlot(s.k.total,features=MAPK_signaling$geneID)+scale_color_gradient(low="grey", "high"=single_col) + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short_double))+coord_flip() +theme(axis.text.x = element_text(angle = 90)) + ggtitle("MAPK_signaling")
dev.off()
pdf("Graphs/Side/Dotplot_Natural_Killer_Cell_Signaling.pdf", width = 10, height = 8)
DotPlot(s.k.total,features=Natural_Killer_Cell_Signaling$geneID)+scale_color_gradient(low="grey", "high"=single_col) + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short_double))+coord_flip() +theme(axis.text.x = element_text(angle = 90)) + ggtitle("Natural_Killer_Cell_Signaling")
dev.off()
pdf("Graphs/Side/Dotplot_NF_kB_signaling.pdf", width = 10, height = 8)
DotPlot(s.k.total,features=NF_kB_signaling$geneID)+scale_color_gradient(low="grey", "high"=single_col) + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short_double))+coord_flip() +theme(axis.text.x = element_text(angle = 90)) + ggtitle("NF_kB_signaling")
dev.off()
pdf("Graphs/Side/Dotplot_Phagocytosis_Ag_presentation.pdf", width = 10, height = 10)
DotPlot(s.k.total,features=Phagocytosis_Ag_presentation$geneID)+scale_color_gradient(low="grey", "high"=single_col) + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short_double))+coord_flip() +theme(axis.text.x = element_text(angle = 90)) + ggtitle("Phagocytosis_Ag_presentation")
dev.off()
pdf("Graphs/Side/Dotplot_PI3K_AKT_Signaling.pdf", width = 10, height = 10)
DotPlot(s.k.total,features=PI3K_AKT_Signaling$geneID)+scale_color_gradient(low="grey", "high"=single_col) + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short_double))+coord_flip() +theme(axis.text.x = element_text(angle = 90)) + ggtitle("PI3K_AKT_Signaling")
dev.off()
pdf("Graphs/Side/Dotplot_ROS_Glutathione_Cytotoxic_granules.pdf", width = 10, height = 7)
DotPlot(s.k.total,features=ROS_Glutathione_Cytotoxic_granules$geneID)+scale_color_gradient(low="grey", "high"=single_col) + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short_double))+coord_flip() +theme(axis.text.x = element_text(angle = 90)) + ggtitle("ROS_Glutathione_Cytotoxic_granules")
dev.off()
pdf("Graphs/Side/Dotplot_TNF_Superfamily_Signaling.pdf", width = 10, height = 9)
DotPlot(s.k.total,features=TNF_Superfamily_Signaling$geneID)+scale_color_gradient(low="grey", "high"=single_col) + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short_double))+coord_flip() +theme(axis.text.x = element_text(angle = 90)) + ggtitle("TNF_Superfamily_Signaling")
dev.off()

pdf("Graphs/Side/Dotplot_Inflammation_long.pdf", width = 10, height = 12)
DotPlot(s.k.total,features=c("CD59","ITGB1","VCL","CD9","MYL6","MMP14","CCL2","MYH9","ALCAM","CD99","ILF2","SOCS2","HMGB1","STAT3","SOCS3","MYDGF","IFNGR2","IL1RAP","IRF1","IFNGR1","PRKACB","PRKAR1A","PSMA1","PSMB5","PSME1","PSME2","XBP1","HLA-A2","HLA-B","HLA-C","MGST3","PTGES3","PTGDS","NFKBIA","ATF3","DDIT3","ATF4","HINT1","MEFC2","RAC1","HSPB1","YWHAZ","TNFRSF1A","TNFSF13B","TNFSF12","TNFAIP6","TNFAIP3","TRAF1","TRAF2","HSP90B1","PRDX","CAT","PRDX2","SOD1","SOD2","PRDX4","LMNA","BAD","BIRC3","BAX","DAP","CYCS","SLC3A2","EDNRB","MAL","NR4A1","NR4A2","NR3C1","CITED2","NCAM1","HLA-E"))+scale_color_gradient(low="grey", "high"=single_col) + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short_double))+coord_flip() +theme(axis.text.x = element_text(angle = 90)) + ggtitle("TNF_Superfamily_Signaling")
dev.off()

#Vln_inflamm_tissue
Idents(SC_7s4k)<-"tissue"
SC_inflamm_Vlnlist <- c("TPSG1","TPSD1","TPSB2","TPSAB1","IL1B","IL1A","TNF","PTGDS","PTGES","CYSLTR1","LTC4S","LTB4R","CXCL8","IL6","IFNG","ALDH7A1","MAOB","HNMT","AOC1","HRH4","HRH3","HDC","HRH1","CSF2","CD59","SOD1","C1R","SOD2")
SC_inflamm_Vlnlist <- c("CXCL8","IL6","IL1B","IL1A","IL6","TNF","IFNG")
SC_inflamm_Vlnlist <- c("PTGDS","PTGES","LTB4R","MAOB","HNMT","HLA-E","CD59","SOD1","C1R","SOD2","SERPING1","NFKBIA","HLA-C","HLA-A","HSPB1","HSP90AA1")
SC_inflamm_Vlnlist <- c("CRP","IL1B","IL6","CXCL8","TNF","")
pdf("Graphs/Main/Vlnplot_SC_inflamm.pdf", width=15, height=5)
VlnPlot(SC_7s4k,y=5,features=SC_inflamm_Vlnlist, ncol = 8, pt.size = 0,cols=color_tissue)& stat_summary(fun = "mean", geom = "crossbar")& stat_compare_means(comparisons = list(c("skin","keloid")) ,method = "t.test", label = "p.signif")
dev.off()

Idents(SC_7s4k)<-"tissue"
SC_inflamm_Vlnlist <- c("CXCL8","IFNG","IL1B","IL1A","IL6","TNF","LTC4S","CSF2")
pdf("Graphs/Main/Vlnplot_SC_inflamm.pdf", width=15, height=2.5)
VlnPlot(SC_7s4k,y=8,features=SC_inflamm_Vlnlist, ncol = 8, pt.size = 0,cols=color_tissue)& stat_summary(fun = "mean", geom = "crossbar")& stat_compare_means(comparisons = list(c("skin","keloid")) ,method = "t.test", label = "p.signif")
dev.off()

Idents(s.k.total)<-"tissue"
SC_inflamm_Vlnlist <- c("CXCL8","IFNG","IL1B","IL1A","IL6","TNF","LTC4S","CSF2")
pdf("Graphs/Main/Vlnplot_total_inflamm.pdf", width=15, height=2.5)
VlnPlot(s.k.total,features=SC_inflamm_Vlnlist, ncol = 8, pt.size = 0,cols=color_tissue)& stat_summary(fun = "mean", geom = "crossbar")& stat_compare_means(comparisons = list(c("skin","keloid")) ,method = "t.test", label = "p.signif")
dev.off()

####ECM/Inflamm Figure Main
Inflamm_tiss<- read.xlsx("D:/Direder/Projekte/D002_10xScarWars/D002_Berechnungen/D002_c12_7s4k/Lists/Inflamm_list.xlsx",1, header = T)


pdf("Graphs/Main/Dotplot_Inflamm_short.pdf", width = 10, height = 14)
DotPlot(s.k.total,features=Inflamm_tiss)+scale_color_gradient(low="grey", "high"=single_col) + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short_double))+coord_flip() +theme(axis.text.x = element_text(angle = 90)) + ggtitle("TNF_Superfamily_Signaling")
dev.off()

#Genes included in GO-Calculation
SC_7s4k.clustermarker_subset_UP <- subset(SC_7s4k.clustermarker, Foldchange_UP >=2 )

#SC-EC
GO_SC_EC <- subset(SC_7s4k.clustermarker_subset_UP, SC_7s4k.clustermarker_subset_UP$cluster == "SC-EC")
GO_SC_EC<-arrange(GO_SC_EC,Foldchange_UP)
GO_SC_EC <- as.character(GO_SC_EC$gene)
pdf("Graphs/Main/DotPlot_GO_SC_EC_FCmin2_SC.pdf", width = 25, height = 7)
DotPlot(SC_7s4k, features=GO_SC_EC, group.by = "celltype",assay="RNA") +scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("SC-EC_DEG_FC>=2") 
dev.off()

#SC-FB
GO_SC_FB <- subset(SC_7s4k.clustermarker_subset_UP, SC_7s4k.clustermarker_subset_UP$cluster == "SC-FB")
GO_SC_FB<-arrange(GO_SC_FB,Foldchange_UP)
GO_SC_FB <- as.character(GO_SC_FB$gene)
pdf("Graphs/Main/DotPlot_GO_SC_FB_FCmin2_SC.pdf", width = 25, height = 7)
DotPlot(SC_7s4k, features=GO_SC_FB, group.by = "celltype",assay="RNA") +scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("SC-FB_DEG_FC>=2") 
dev.off()

#SC-Skin
SC_7s4k.clustermarker_subset_UP_FCmin3 <- subset(SC_7s4k.clustermarker, Foldchange_UP >=3 )
GO_SC_Skin <- subset(SC_7s4k.clustermarker_subset_UP_FCmin3, SC_7s4k.clustermarker_subset_UP_FCmin3$cluster == "SC-Skin")
GO_SC_Skin<-arrange(GO_SC_Skin,Foldchange_UP)
GO_SC_Skin <- as.character(GO_SC_Skin$gene)
#GO_SC_Skin<-GO_SC_Skin[!duplicated(GO_SC_Skin)]
pdf("Graphs/Main/DotPlot_GO_SC_Skin_FCmin3_SC.pdf", width = 30, height = 7)
DotPlot(SC_7s4k, features=GO_SC_Skin, group.by = "celltype",assay="RNA") +scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("SC-Skin_DEG_FC>=3") 
dev.off()

#SC-Prolif
GO_SC_Prolif <- subset(SC_7s4k.clustermarker_subset_UP, SC_7s4k.clustermarker_subset_UP$cluster == "SC-Prolif")
GO_SC_Prolif<-arrange(GO_SC_Prolif,Foldchange_UP)
GO_SC_Prolif <- as.character(GO_SC_Prolif$gene)
pdf("Graphs/Main/DotPlot_GO_SC_Prolif_FCmin2_SC.pdf", width = 25, height = 7)
DotPlot(SC_7s4k, features=GO_SC_Prolif, group.by = "celltype",assay="RNA") +scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("SC-Prolif_DEG_FC>=2") 
dev.off()

#SC-Promyel
GO_SC_Promyel <- subset(SC_7s4k.clustermarker_subset_UP, SC_7s4k.clustermarker_subset_UP$cluster == "SC-Promyel")
GO_SC_Promyel<-arrange(GO_SC_Promyel,Foldchange_UP)
GO_SC_Promyel <- as.character(GO_SC_Promyel$gene)
pdf("Graphs/Main/DotPlot_GO_SC_Promyel_FCmin2_SC.pdf", width = 25, height = 7)
DotPlot(SC_7s4k, features=GO_SC_Promyel, group.by = "celltype",assay="RNA") +scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("SC-Promyel_DEG_FC>=2") 
dev.off()

#SC-Repair
SC_7s4k.clustermarker_subset_UP_FC1.5 <- subset(SC_7s4k.clustermarker, Foldchange_UP >=1.5 )
GO_SC_Repair <- subset(SC_7s4k.clustermarker_subset_UP_FC1.5, SC_7s4k.clustermarker_subset_UP_FC1.5$cluster == "SC-Repair")
GO_SC_Repair<-arrange(GO_SC_Repair,Foldchange_UP)
GO_SC_Repair <- as.character(GO_SC_Repair$gene)
pdf("Graphs/Main/DotPlot_GO_SC_Repair_FCmin1.5_SC.pdf", width = 25, height = 7)
DotPlot(SC_7s4k, features=GO_SC_Repair, group.by = "celltype",assay="RNA") +scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("SC-Repair_DEG_FC>=1.5") 
dev.off()

###GO-included Genes
####Template_adapt celltype
library(reshape2)
GO_df<- read.xlsx("D:/Direder/Projekte/D002_10xScarWars/D002_Berechnungen/D002_c12_7s4k/Metascape/SC-EC-UP/D002_c12_SC-EC_FCmin2_UP.xlsx",2, header = F)

head(GO_df)
rownames(GO_df)= as.character(GO_df$X1)
GO_df <- GO_df[,-1 ]
GO_names<- c("GO:1","GO:2","GO:3","GO:4","GO:5","GO:6","GO:7","GO:8","GO:9","GO:10","GO:11","GO:12","GO:13","GO:14","GO:15","GO:16","GO:17","GO:18","GO:19","GO:20")
colnames(GO_df) <- GO_names
GO_df<- as.matrix(GO_df)
GO_df<- melt(GO_df)
head(GO_df)

pdf("Metascape/D002_c12_SC_GO-template.pdf", width=12, height=12)
ggplot(GO_df, aes(x = Var2, y = Var1, fill=value)) +
  geom_tile() + scale_fill_manual(breaks=c("0","1"),values=c("#fdf6c0","#900603") , guide= "colorbar")   + 
  theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),axis.text.y=element_text(size=9),plot.title=element_text(size=11))+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.ontop = TRUE,
        panel.background = element_rect(fill = "transparent"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())+coord_flip()+
  scale_y_discrete(guide = guide_axis(angle = -90))+
  geom_vline(xintercept=seq(1.5, length(GO_df$Var2)),colour = "#636363", size=0.02) + 
  geom_hline(yintercept=seq(1.5, length(GO_df$Var1)),colour = "#636363", size=0.02)+ theme(legend.position = "none")
dev.off()

#Genes included in GO-Calculation
SC_7s4k.clustermarker_subset_UP <- subset(SC_7s4k.clustermarker, Foldchange_UP >=2 )

#SC-celltype
GO_SC_celltype <- subset(SC_7s4k.clustermarker_subset_UP, SC_7s4k.clustermarker_subset_UP$cluster == "SC-EC")
GO_SC_celltype<-arrange(GO_SC_celltype,gene)
GO_SC_celltype <- as.character(GO_SC_celltype$gene)
pdf("Graphs/Main/DotPlot_SC_Template_GO_included.pdf", width = 18, height = 3)
DotPlot(SC_7s4k, features=GO_SC_celltype, group.by = "celltype",assay="RNA") +scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("SC-repair Top-Clustermarker avglogFC=>1.5") 
dev.off()

############################
####Subset Macrophages####
############################
Idents(s.k.total)<- s.k.total$celltype
MAC_7s4k <- subset(s.k.total,idents = "MAC")
Idents(MAC_7s4k)<- MAC_7s4k$sample
MAC_7s4k <- subset(MAC_7s4k,idents = c("skin_1","skin_2","skin_3","skin_4","skin_5","skin_6","keloid_1","keloid_2","keloid_3L","keloid_3R"))
DefaultAssay(MAC_7s4k)<-"RNA"
MAC_7s4k[["percent.mt"]] <- PercentageFeatureSet(MAC_7s4k, pattern = "^MT-")
MAC_7s4k <- SCTransform(MAC_7s4k, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
MAC_7s4k <- RunPCA(MAC_7s4k, npcs = 50)
ElbowPlot(MAC_7s4k, ndims = 50)
MAC_7s4k <- RunUMAP(MAC_7s4k, dims = 1:30)
MAC_7s4k <- FindNeighbors(MAC_7s4k, dims = 1:30)
#DefaultAssay(MAC_7s4k)<- "integrated"
MAC_7s4k <- FindClusters(MAC_7s4k, resolution = 0.2)
UMAPPlot(MAC_7s4k, label=T, split.by="tissue")
VlnPlot(MAC_7s4k,c("CD68","CD40","MRC1","SELE","DCN"))
DefaultAssay(MAC_7s4k)<- "RNA"
MAC_7s4k <- NormalizeData(MAC_7s4k)
StackedVlnPlot(MAC_7s4k,features=c("CD68","CD40","MRC1","SELE","DCN","LUM","COL1A1","AIFM2","ACTA2"))
MAC_7s4k <- subset(MAC_7s4k,idents = c("0","1","2","3"))
DefaultAssay(MAC_7s4k)<-"RNA"
MAC_7s4k[["percent.mt"]] <- PercentageFeatureSet(MAC_7s4k, pattern = "^MT-")
MAC_7s4k <- SCTransform(MAC_7s4k, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
MAC_7s4k <- RunPCA(MAC_7s4k, npcs = 50)
ElbowPlot(MAC_7s4k, ndims = 50)
MAC_7s4k <- RunUMAP(MAC_7s4k, dims = 1:30)
MAC_7s4k <- FindNeighbors(MAC_7s4k, dims = 1:30)
#DefaultAssay(MAC_7s4k)<- "integrated"
MAC_7s4k <- FindClusters(MAC_7s4k, resolution = 0.2)
UMAPPlot(MAC_7s4k, label=T)
DefaultAssay(MAC_7s4k)<- "RNA"
MAC_7s4k <- NormalizeData(MAC_7s4k)

FeaturePlot(MAC_7s4k,features= c("CD68", "CD80","CD86","CD40","TNF","IL1B","CXCL2","HLA-DMA","HLA-DOA","HLA-DPA1","HLA-DQA1","HLA-DRA","MRC1","CD163","CSF1R","IL10","LUM","DCN","COL1A1"))+ coord_flip()
FeaturePlot(MAC_7s4k,c("CD68","DCN"),blend = T, blend.threshold = 0.5, cols = c("lightgrey","#990F26","#54990F"), pt.size = 1.5, order=T, min.cutoff = 0, max.cutoff = 2)
FeaturePlot(MAC_7s4k,c("CD68","MRC1"),blend = T, blend.threshold = 0.5, cols = c("lightgrey","#990F26","#54990F"), pt.size = 1.5, order=T, min.cutoff = 0, max.cutoff = 2)
FeaturePlot(MAC_7s4k,c("MRC1","IL1B"),blend = T, blend.threshold = 0.5, cols = c("lightgrey","#990F26","#54990F"), pt.size = 1.5, order=T, min.cutoff = 0, max.cutoff = 2)
FeaturePlot(MAC_7s4k,c("MRC1","CD163"),blend = T, blend.threshold = 0.5, cols = c("lightgrey","#990F26","#54990F"), pt.size = 1.5, order=T, min.cutoff = 0, max.cutoff = 2)

DefaultAssay(MAC_7s4k)<-"RNA"
M1<-FeaturePlot(MAC_7s4k, "CD68", min.cutoff = 0, order=T)+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("CD68_allMAC")
M2<-FeaturePlot(MAC_7s4k, "IL1B", min.cutoff = 0, order=T)+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("IL1B_M1-MAC")
M3<-FeaturePlot(MAC_7s4k, "CXCL2", min.cutoff = 0, order=T)+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("CXCL2_M1-MAC")
M4<-FeaturePlot(MAC_7s4k, "MRC1", order=T)+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("MRC1_M2-MAC")
M5<-FeaturePlot(MAC_7s4k, "CD163", order=T)+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("CD163_M2-MAC")
M6<-FeaturePlot(MAC_7s4k, "LUM", order=T)+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("LUM_MAC_to_FB")
pdf("Graphs/MAC/FeaturePlots_MAC.Clusteridentification.total.pdf",width=12, height=15)
wrap_plots(M1,M2,M3,M4,M5,M6, ncol = 2)
dev.off()



FeaturePlot(MAC_7s4k,c("IL1B","IL6","TNF"), split.by = "tissue")

#Name Cluster
Idents(MAC_7s4k)<-MAC_7s4k$seurat_clusters
MAC_7s4k<-RenameIdents(MAC_7s4k,
                       `0`="MAC-M2", `1`="MAC-M1/M2", `2`="MAC_to_FB", `3`="MAC-M1")
MAC_7s4k$celltype<-Idents(MAC_7s4k)


#VlnPlot(s.k.total,c("CCL2","CCL3","CCL5","TNF","MMP9","TNFAIP6","TIMP1"), split.by="tissue", pt.size=0)
#VlnPlot(SC_7s4k,c("CCL2","CCL3","CCL5","TNF","MMP9","TNFAIP6"))
#VlnPlot(SC_7s4k,c("IGF1","NRG1","NRG2","NRG3","MAPK1","MAPK3","TIMP1"))
#DotPlot(SC_7s4k,features=c("S100A9","ADAMTS4","S100A8","TAGLN","KLF5","MMP9","NEFM","MMP14","TIMP1","CD44","CDH1","IL1RN","LCN2","IL15","IFNB1","ATF5"))

#Define cluster levels
clusters_ordered<-c("MAC-M1", "MAC-M2","MAC-M1/M2","MAC_to_FB")
MAC_7s4k$celltype<- factor(MAC_7s4k$celltype, levels = clusters_ordered)
Idents(MAC_7s4k)<-MAC_7s4k$celltype
pdf("Graphs/MAC/UMAP_MAC.pdf")
UMAPPlot(MAC_7s4k, label=T)+ scale_color_manual(values=color_MAC)
dev.off()
pdf("Graphs/MAC/UMAP_MAC_split.by.tissue.pdf", width = 7, height=5)
UMAPPlot(MAC_7s4k, label=T, split.by="tissue")+ scale_color_manual(values=color_MAC)
dev.off()
pdf("Graphs/MAC/VlnPlot_MAC.pdf")
VlnPlot(MAC_7s4k, c("CD68","AIF1","PLAUR","MMP9","MRC1"), cols = color_MAC)
dev.off()

##############
##Pseudotime##
##############
library(monocle3)

#### Create a Monocle CDS Object
# Project PC dimensions to whole data set
MAC_convert <- ProjectDim(MAC_7s4k, reduction = "pca")

# Create an expression matrix
expression_matrix <- MAC_convert@assays$RNA@counts

# Get cell metadata
cell_metadata <- MAC_convert@meta.data
if (all.equal(colnames(expression_matrix), rownames(cell_metadata))) {
  print(sprintf("Cell identifiers match"))
} else {
  print(sprintf("Cell identifier mismatch - %i cells in expression matrix, %i cells in metadata",
                ncol(expression_matrix), nrow(cell_metadata)))
  print("If the counts are equal, sort differences will throw this error")
}

# get gene annotations
gene_annotation <- data.frame(gene_short_name = rownames(MAC_convert@assays$RNA), row.names = rownames(MAC_convert@assays$RNA))
if (all.equal(rownames(expression_matrix), rownames(gene_annotation))) {
  print(sprintf("Gene identifiers all match"))
} else {
  print(sprintf("Gene identifier mismatch - %i genes in expression matrix, %i gene in gene annotation",
                nrow(expression_matrix), nrow(gene_annotation)))
  print("If the counts are equal, sort differences will throw this error")
}

# Seurat-derived CDS
MAC_7s4k.cds <- new_cell_data_set(expression_matrix,
                                 cell_metadata = cell_metadata,
                                 gene_metadata = gene_annotation)

# Transfer Seurat embeddings
# Note that these may be calculated on the Integrated object, not the counts
#   and thus will involve fewer genes
reducedDim(MAC_7s4k.cds, type = "PCA") <- MAC_7s4k@reductions$pca@cell.embeddings 
MAC_7s4k.cds@preprocess_aux$prop_var_expl <- MAC_7s4k@reductions$pca@stdev
plot_pc_variance_explained(MAC_7s4k.cds)

# Transfer Seurat UMAP embeddings
MAC_7s4k.cds@int_colData@listData$reducedDims$UMAP <- MAC_7s4k@reductions$umap@cell.embeddings
#    plot_cells(MAC_7s4k.cds)

# Copy cluster info from Seurat
MAC_7s4k.cds@clusters$UMAP_so$clusters <- MAC_7s4k@meta.data$gt_tp_cell_type_integrated_.0.9

MAC_7s4k.cds <- cluster_cells(MAC_7s4k.cds, reduction_method = "UMAP")

# Fix from https://gitmemory.com/cole-trapnell-lab
rownames(MAC_7s4k.cds@principal_graph_aux$UMAP$dp_mst) <- NULL
colnames(MAC_7s4k.cds@int_colData@listData$reducedDims$UMAP) <- NULL

DimPlot(MAC_7s4k, reduction = "umap")

MAC_7s4k.cds = cluster_cells(MAC_7s4k.cds, resolution = 0.002, partition_qval= 1)

p1 <- plot_cells(MAC_7s4k.cds,group_label_size = 3.5, cell_size = 2) + ggtitle("MAC_7s4k_ _ cluster")
p2 <- plot_cells(MAC_7s4k.cds, color_cells_by="celltype", group_cells_by="partition",group_label_size = 3.5, cell_size = 2) + ggtitle("MAC_7s4k_partition")
wrap_plots(p1, p2)

MAC_7s4k.cds <- learn_graph(MAC_7s4k.cds,learn_graph_control = list(prune_graph=F, ncenter=9))

pdf("Graphs/MAC/Pseudo_Principal_graph_MAC_7s4k.pdf")
plot_cells(MAC_7s4k.cds, label_groups_by_cluster = T,color_cells_by="celltype", label_leaves = FALSE, label_branch_points = FALSE, cell_size=2, label_cell_groups = F)+ scale_color_manual(values=color_MAC) + ggtitle("Principal graph")
dev.off()

MAC_7s4k.cds<- order_cells(MAC_7s4k.cds)

pdf("Graphs/MAC/Pseudo_M1_start_MAC_7s4k.pdf")
plot_cells(MAC_7s4k.cds, label_groups_by_cluster = T,color_cells_by="pseudotime", label_leaves = FALSE, label_branch_points = FALSE, cell_size=2, label_cell_groups = F) + ggtitle("Pseudotime_ MAC-Skin as root")
dev.off()

plot_cells(MAC_7s4k.cds,
           genes=c("CD68"),           
           label_cell_groups=F,
           show_trajectory_graph=T, cell_size=1.5)


Mod_cds <- MAC_7s4k.cds[,grepl("MAC", colData(MAC_7s4k.cds)$celltype, ignore.case=TRUE)]
plot_cells(Mod_cds, color_cells_by="celltype")
pr_graph_test_res<- graph_test (Mod_cds, neighbor_graph = "principal_graph", cores = 1)
pr_deg_ids <- row.names(subset(pr_graph_test_res,q_value < 0.05))
gene_module_df <- find_gene_modules(Mod_cds[pr_deg_ids,], resolution=1e-2)#, resolution=1e-1
plot_cells(MAC_7s4k.cds, genes=gene_module_df, 
           show_trajectory_graph=FALSE, 
           label_cell_groups=FALSE)

#Adapt Excel-output
X1 <- subset(pr_graph_test_res, q_value < 0.05)
X1 <- cbind(rownames(X1),X1)
X1 <- subset(X1, select=c("rownames(X1)","gene_short_name","morans_I"))
write.xlsx(X1,"Lists/pr_graph_MAC_7s4k.xlsx")


#Analyse Pseudotime
Test_genes <- c("MRC1", "CD68", "CXCL2","TNF","IL1B","KLF4","THY1")
Test_lineage_cds <- MAC_7s4k.cds[rowData(MAC_7s4k.cds)$gene_short_name %in% Test_genes,colData(MAC_7s4k.cds)$celltype %in% c("MAC-M1","MAC-M1/M2","MAC_to_FB")]

pdf("Graphs/MAC/Pseudo_Test_MAC_7s4k.pdf", width = 8, height=10)
plot_genes_in_pseudotime(Test_lineage_cds, color_cells_by = "celltype", min_expr=0.1, cell_size=1.5)  + scale_color_manual(values=color_MAC) + ggtitle("Test expression along Pseudotime")
dev.off()

Clustermarker_FC2<- read.xlsx("D:/Direder/Projekte/D002_10xScarWars/D002_Berechnungen/D002_c12_7s4k/Lists/Clustermarker_MAC_7s4k.xlsx",2, header = T)
Clustermarker_M1_FC2 <- Clustermarker_FC2$M1
Clustermarker_M1.M2_FC2 <- Clustermarker_FC2$M1.M2
Clustermarker_M2_FC2 <- Clustermarker_FC2$M2.FC2
Clustermarker_MAC_to_FB_FC2 <- Clustermarker_FC2$MAC_to_FB.FC2

#M1
M1_lineage_cds <- MAC_7s4k.cds[rowData(MAC_7s4k.cds)$gene_short_name %in% Clustermarker_M1_FC2,colData(MAC_7s4k.cds)$celltype %in% c("MAC-M1","MAC-M1/M2","MAC-M2")]

pdf("Graphs/MAC/Pseudo_M1_MAC_7s4k.pdf", width = 20, height=15)
plot_genes_in_pseudotime(M1_lineage_cds, color_cells_by = "celltype", min_expr=0.5, cell_size=1,ncol = 10)  + scale_color_manual(values=color_MAC) + ggtitle("M1-genes expression along Pseudotime")
dev.off()

#M1.M2
M1.M2_lineage_cds <- MAC_7s4k.cds[rowData(MAC_7s4k.cds)$gene_short_name %in% Clustermarker_M1.M2_FC2]

pdf("Graphs/MAC/Pseudo_M1.M2_MAC_7s4k.pdf", width = 20, height=15)
plot_genes_in_pseudotime(M1.M2_lineage_cds, color_cells_by = "celltype", min_expr=0.5, cell_size=1,ncol = 10)  + scale_color_manual(values=color_MAC) + ggtitle("M1-genes expression along Pseudotime")
dev.off()

#M2
M2_lineage_cds <- MAC_7s4k.cds[rowData(MAC_7s4k.cds)$gene_short_name %in% Clustermarker_M2_FC2_1.5,colData(MAC_7s4k.cds)$celltype %in% c("MAC-M1","MAC-M1/M2","MAC-M2")]

pdf("Graphs/MAC/Pseudo_M2_MAC_7s4k.pdf", width = 20, height=15)
plot_genes_in_pseudotime(M2_lineage_cds, color_cells_by = "celltype", min_expr=0.5, cell_size=1,ncol = 10)  + scale_color_manual(values=color_MAC) + ggtitle("M1-genes expression along Pseudotime")
dev.off()

#MAC_to_FB
MAC_to_FB_lineage_cds <- MAC_7s4k.cds[rowData(MAC_7s4k.cds)$gene_short_name %in% Clustermarker_MAC_to_FB_FC2,colData(MAC_7s4k.cds)$celltype %in% c("MAC-M1","MAC-M1/M2","MAC_to_FB")]

pdf("Graphs/MAC/Pseudo_MAC_to_FB_MAC_7s4k.pdf", width = 20, height=15)
plot_genes_in_pseudotime(MAC_to_FB_lineage_cds, color_cells_by = "celltype", min_expr=0.5, cell_size=1,ncol = 10)  + scale_color_manual(values=color_MAC) + ggtitle("M1-genes expression along Pseudotime")
dev.off()


#Clustermarker + Heatmap total
MAC_7s4k.clustermarker<-FindAllMarkers(MAC_7s4k,  assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
MAC_7s4k.clustermarker$Foldchange_UP <- 2^(MAC_7s4k.clustermarker$avg_log2FC)
MAC_7s4k.clustermarker$Foldchange_DOWN <- 2^(-MAC_7s4k.clustermarker$avg_log2FC)
MAC_7s4k.clustermarker$Ratio_pct1_pct2 <- (MAC_7s4k.clustermarker$pct.1)/(MAC_7s4k.clustermarker$pct.2)
write.xlsx(MAC_7s4k.clustermarker, "Lists/Clustermarker_MAC_7s4k.xlsx")


#DotPlot Top up and down regulated genes 
#up
MAC_7s4k.clustermarker_subset_UP <- subset(MAC_7s4k.clustermarker, Foldchange_UP >=2 )
top10_UP_MACsubs_clustermarker <- MAC_7s4k.clustermarker_subset_UP %>% group_by(cluster) %>% top_n(n=10,wt=Foldchange_UP)
top10_UP_MACsubs_clustermarker<-arrange(top10_UP_MACsubs_clustermarker,cluster,Foldchange_UP)
top10_UP_MACsubs_clustermarker <- as.character(top10_UP_MACsubs_clustermarker$gene)
top10_UP_MACsubs_clustermarker<-top10_UP_MACsubs_clustermarker[!duplicated(top10_UP_MACsubs_clustermarker)]

VlnPlot(MAC_7s4k,c("CCL2","CCR2","CCR4","MMP12"))

pdf("Graphs/MAC/DotPlot_Top_10UP_Clustermarker_MAC_FCmin2_MAC.pdf", width = 7, height = 12)
DotPlot(MAC_7s4k, features=top10_UP_MACsubs_clustermarker, group.by = "celltype",assay="RNA") + coord_flip()+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 45)) + ggtitle("MAC Cluster_Top-MAC DEG") + scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()

#down
MAC_7s4k.clustermarker_subset_DOWN <- subset(MAC_7s4k.clustermarker, Foldchange_DOWN >=2 )
top10_DOWN_MACsubs_clustermarker <- MAC_7s4k.clustermarker_subset_DOWN %>% group_by(cluster) %>% top_n(n=10,wt=Foldchange_DOWN)
top10_DOWN_MACsubs_clustermarker<-arrange(top10_DOWN_MACsubs_clustermarker,cluster,Foldchange_DOWN)
top10_DOWN_MACsubs_clustermarker <- as.character(top10_DOWN_MACsubs_clustermarker$gene)
top10_DOWN_MACsubs_clustermarker<-top10_DOWN_MACsubs_clustermarker[!duplicated(top10_DOWN_MACsubs_clustermarker)]

pdf("Graphs/MAC/DotPlot_Top_10DOWN_Clustermarker_MAC_FCmin2_MAC.pdf", width = 7, height = 12)
DotPlot(MAC_7s4k, features=top10_DOWN_MACsubs_clustermarker, group.by = "celltype",assay="RNA") + coord_flip()+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 45)) + ggtitle("MAC Cluster_Top10-MAC DEG") + scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()

#Percentage total
###Pieplot/Barplot cellnumber /tissue/cluster
#subset conditions
skin_MAC.subset<-subset(MAC_7s4k, tissue=="skin")
keloid_MAC.subset<-subset(MAC_7s4k, tissue=="keloid")

#Data frame frequencies
skin_MAC.cluster.frequencies<-as.data.frame(table(skin_MAC.subset$celltype))
keloid_MAC.cluster.frequencies<-as.data.frame(table(keloid_MAC.subset$celltype))

#Barplot Cluster percentage
skin_MAC.cluster.frequencies$tissue="skin"
keloid_MAC.cluster.frequencies$tissue="keloid"

skin_MAC.cluster.frequencies<-mutate(skin_MAC.cluster.frequencies,Percentage= skin_MAC.cluster.frequencies$Freq/sum(skin_MAC.cluster.frequencies$Freq)*100)
keloid_MAC.cluster.frequencies<-mutate(keloid_MAC.cluster.frequencies,Percentage= keloid_MAC.cluster.frequencies$Freq/sum(keloid_MAC.cluster.frequencies$Freq)*100)

skin_kel_MAC.cluster.df <- rbind(skin_MAC.cluster.frequencies,keloid_MAC.cluster.frequencies)
skin_kel_MAC.cluster.df$tissue <- factor(x = skin_kel_MAC.cluster.df$tissue, levels = c("skin", "keloid"))

pdf("Graphs/MAC/Barplot_MAC_subset_skin_kel_cellfrequencie.pdf", width = 4, height = 6)
ggplot(data=skin_kel_MAC.cluster.df, aes(x=tissue, y=Freq, fill=Var1)) +
  geom_bar(stat="identity", position="stack")+scale_fill_manual(values =color_MAC)+
  ylim(0,1000) +  theme_classic() +ggtitle("Cell source_macrophage cluster") + xlab("condition")+ylab("total cell count") 
dev.off()

pdf("Graphs/MAC/Barplot_MAC_subset_skin_kel_cellpercentage.pdf", width = 4, height = 6)
ggplot(data=skin_kel_MAC.cluster.df, aes(x=tissue, y=Percentage, fill=Var1)) +
  geom_bar(stat="identity", position="stack")+scale_fill_manual(values =color_MAC) +  theme_classic() +ggtitle("Cell source_macrophage cluster") + xlab("condition")+ylab("total cell count") 
dev.off()

pdf("Graphs/MAC/PiePlot_MAC_subset_skin_kel_cellfrequencie.pdf", width = 10, height = 10)
pie(skin_MAC.cluster.frequencies$Freq, labels = skin_MAC.cluster.frequencies$Var1, main="skin",  clockwise = T, col = c(color_MAC))
pie(keloid_MAC.cluster.frequencies$Freq, labels = keloid_MAC.cluster.frequencies$Var1, main="Keloid",  clockwise = T, col = c(color_MAC))
dev.off()

#Percentage total
###Pieplot/Barplot cellnumber /sample/cluster
#subset conditions
MAC_skin1.subset<-subset(MAC_7s4k, sample=="skin_1")
MAC_skin2.subset<-subset(MAC_7s4k, sample=="skin_2")
MAC_skin4.subset<-subset(MAC_7s4k, sample=="skin_4")
MAC_skin3.subset<-subset(MAC_7s4k, sample=="skin_3")
MAC_skin5.subset<-subset(MAC_7s4k, sample=="skin_5")
MAC_skin6.subset<-subset(MAC_7s4k, sample=="skin_6")
#MAC_skin7.subset<-subset(MAC_7s4k, sample=="skin_7")
MAC_keloid1.subset<-subset(MAC_7s4k, sample=="keloid_1")
MAC_keloid2.subset<-subset(MAC_7s4k, sample=="keloid_2")
MAC_keloid3L.subset<-subset(MAC_7s4k, sample=="keloid_3L")
MAC_keloid3R.subset<-subset(MAC_7s4k, sample=="keloid_3R")

#Data frame frequencies
MAC_skin1.cluster.frequencies<-as.data.frame(table(MAC_skin1.subset$celltype))
MAC_skin2.cluster.frequencies<-as.data.frame(table(MAC_skin2.subset$celltype))
MAC_skin3.cluster.frequencies<-as.data.frame(table(MAC_skin3.subset$celltype))
MAC_skin4.cluster.frequencies<-as.data.frame(table(MAC_skin4.subset$celltype))
MAC_skin5.cluster.frequencies<-as.data.frame(table(MAC_skin5.subset$celltype))
MAC_skin6.cluster.frequencies<-as.data.frame(table(MAC_skin6.subset$celltype))
#MAC_skin7.cluster.frequencies<-as.data.frame(table(MAC_skin7.subset$celltype))
MAC_keloid1.cluster.frequencies<-as.data.frame(table(MAC_keloid1.subset$celltype))
MAC_keloid2.cluster.frequencies<-as.data.frame(table(MAC_keloid2.subset$celltype))
MAC_keloid3L.cluster.frequencies<-as.data.frame(table(MAC_keloid3L.subset$celltype))
MAC_keloid3R.cluster.frequencies<-as.data.frame(table(MAC_keloid3R.subset$celltype))

#Barplot Cluster percentage
MAC_skin1.cluster.frequencies$sample="skin_1"
MAC_skin2.cluster.frequencies$sample="skin_2"
MAC_skin3.cluster.frequencies$sample="skin_3"
MAC_skin4.cluster.frequencies$sample="skin_4"
MAC_skin5.cluster.frequencies$sample="skin_5"
MAC_skin6.cluster.frequencies$sample="skin_6"
#MAC_skin7.cluster.frequencies$sample="skin_7"
MAC_keloid1.cluster.frequencies$sample="keloid_1"
MAC_keloid2.cluster.frequencies$sample="keloid_2"
MAC_keloid3L.cluster.frequencies$sample="keloid_3L"
MAC_keloid3R.cluster.frequencies$sample="keloid_3R"

MAC_skin1.cluster.frequencies<-mutate(MAC_skin1.cluster.frequencies,Percentage= MAC_skin1.cluster.frequencies$Freq/sum(MAC_skin1.cluster.frequencies$Freq)*100)
MAC_skin2.cluster.frequencies<-mutate(MAC_skin2.cluster.frequencies,Percentage= MAC_skin2.cluster.frequencies$Freq/sum(MAC_skin2.cluster.frequencies$Freq)*100)
MAC_skin3.cluster.frequencies<-mutate(MAC_skin3.cluster.frequencies,Percentage= MAC_skin3.cluster.frequencies$Freq/sum(MAC_skin3.cluster.frequencies$Freq)*100)
MAC_skin4.cluster.frequencies<-mutate(MAC_skin4.cluster.frequencies,Percentage= MAC_skin4.cluster.frequencies$Freq/sum(MAC_skin4.cluster.frequencies$Freq)*100)
MAC_skin5.cluster.frequencies<-mutate(MAC_skin5.cluster.frequencies,Percentage= MAC_skin5.cluster.frequencies$Freq/sum(MAC_skin5.cluster.frequencies$Freq)*100)
MAC_skin6.cluster.frequencies<-mutate(MAC_skin6.cluster.frequencies,Percentage= MAC_skin6.cluster.frequencies$Freq/sum(MAC_skin6.cluster.frequencies$Freq)*100)
#MAC_skin7.cluster.frequencies<-mutate(MAC_skin7.cluster.frequencies,Percentage= MAC_skin7.cluster.frequencies$Freq/sum(MAC_skin7.cluster.frequencies$Freq)*100)
MAC_keloid1.cluster.frequencies<-mutate(MAC_keloid1.cluster.frequencies,Percentage= MAC_keloid1.cluster.frequencies$Freq/sum(MAC_keloid1.cluster.frequencies$Freq)*100)
MAC_keloid2.cluster.frequencies<-mutate(MAC_keloid2.cluster.frequencies,Percentage= MAC_keloid2.cluster.frequencies$Freq/sum(MAC_keloid2.cluster.frequencies$Freq)*100)
MAC_keloid3L.cluster.frequencies<-mutate(MAC_keloid3L.cluster.frequencies,Percentage= MAC_keloid3L.cluster.frequencies$Freq/sum(MAC_keloid3L.cluster.frequencies$Freq)*100)
MAC_keloid3R.cluster.frequencies<-mutate(MAC_keloid3R.cluster.frequencies,Percentage= MAC_keloid3R.cluster.frequencies$Freq/sum(MAC_keloid3R.cluster.frequencies$Freq)*100)


MAC_skin_kel_samp.cluster.df <- rbind(MAC_skin1.cluster.frequencies,MAC_skin2.cluster.frequencies,MAC_skin3.cluster.frequencies,MAC_skin4.cluster.frequencies,MAC_skin5.cluster.frequencies,MAC_skin6.cluster.frequencies,MAC_keloid1.cluster.frequencies,MAC_keloid2.cluster.frequencies,MAC_keloid3L.cluster.frequencies,MAC_keloid3R.cluster.frequencies)#,MAC_skin7.cluster.frequencies
MAC_skin_kel_samp.cluster.df$sample <- factor(x = MAC_skin_kel_samp.cluster.df$sample, levels = c("skin_1","skin_2","skin_3","skin_4","skin_5","skin_6", "keloid_1", "keloid_2", "keloid_3L", "keloid_3R"))#"skin_7",

pdf("Graphs/MAC/Barplot_MAC_skin_kel_cell_clusterpercentage.pdf", width = 10, height = 5)
ggplot(data=MAC_skin_kel_samp.cluster.df, aes(x=sample, y=Percentage, fill=Var1)) +
  geom_bar(stat="identity", position="stack")+scale_fill_manual(values =color_MAC)+ theme_classic()
dev.off()

pdf("Graphs/MAC/Barplot_MAC_skin_kel_cell_clusterfrequencie.pdf", width = 7, height = 7)
ggplot(data=MAC_skin_kel_samp.cluster.df, aes(x=sample, y=Freq, fill=Var1)) +
  geom_bar(stat="identity", position="stack") +  theme_classic() +
  ggtitle("Cell source_macrophage cluster") + xlab("condition")+ylab("total cell count") +scale_fill_manual(values =color_MAC)
dev.off()

pdf("Graphs/MAC/PiePlot_MAC_skin_kel_cell_clusterfrequencie.pdf", width = 10, height = 10)
pie(MAC_skin1.cluster.frequencies$Freq, labels = MAC_skin1.cluster.frequencies$Var1, main="skin_1",  clockwise = T, col = c(color_MAC))
pie(MAC_skin2.cluster.frequencies$Freq, labels = MAC_skin2.cluster.frequencies$Var1, main="skin_2",  clockwise = T, col = c(color_MAC))
pie(MAC_skin4.cluster.frequencies$Freq, labels = MAC_skin4.cluster.frequencies$Var1, main="skin_4",  clockwise = T, col = c(color_MAC))
pie(MAC_skin5.cluster.frequencies$Freq, labels = MAC_skin5.cluster.frequencies$Var1, main="skin_5",  clockwise = T, col = c(color_MAC))
pie(MAC_skin6.cluster.frequencies$Freq, labels = MAC_skin6.cluster.frequencies$Var1, main="skin_6",  clockwise = T, col = c(color_MAC))
#pie(MAC_skin7.cluster.frequencies$Freq, labels = MAC_skin7.cluster.frequencies$Var1, main="skin_7",  clockwise = T, col = c(color_MAC))
pie(MAC_keloid1.cluster.frequencies$Freq, labels = MAC_keloid1.cluster.frequencies$Var1, main="Keloid_1",  clockwise = T, col = c(color_MAC))
pie(MAC_keloid2.cluster.frequencies$Freq, labels = MAC_keloid2.cluster.frequencies$Var1, main="Keloid_2",  clockwise = T, col = c(color_MAC))
pie(MAC_keloid3L.cluster.frequencies$Freq, labels = MAC_keloid3L.cluster.frequencies$Var1, main="Keloid_3L",  clockwise = T, col = c(color_MAC))
pie(MAC_keloid3R.cluster.frequencies$Freq, labels = MAC_keloid3R.cluster.frequencies$Var1, main="Keloid_3R",  clockwise = T, col = c(color_MAC))
dev.off()

#Identify pseudotime relevant genes
Idents(MAC_7s4k)<-"celltype"
MAC_M1_vs_M1.M2<-FindMarkers(MAC_7s4k, ident.1="MAC-M1", ident.2=c("MAC-M1/M2"), assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
MAC_M1_vs_M1.M2$Foldchange_UP <- 2^(MAC_M1_vs_M1.M2$avg_log2FC)
MAC_M1_vs_M1.M2$Foldchange_DOWN <- 2^(-MAC_M1_vs_M1.M2$avg_log2FC)
MAC_M1_vs_M1.M2$Ratio_pct1_pct2 <- (MAC_M1_vs_M1.M2$pct.1)/(MAC_M1_vs_M1.M2$pct.2)
MAC_M1_vs_M1.M2$comparison <- "MAC-M1 vs M1.M2"
write.xlsx(MAC_M1_vs_M1.M2, "Lists/Pseudotime_relevant/MAC_M1_vs_M1.M2.xlsx")

SC_M1.M2_vs_Skin_repair<-FindMarkers(MAC_7s4k, ident.1="MAC_M1.M2", ident.2=c("SC-Repair","MAC-M1"), assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
SC_M1.M2_vs_Skin_repair$Foldchange_UP <- 2^(SC_M1.M2_vs_Skin_repair$avg_log2FC)
SC_M1.M2_vs_Skin_repair$Foldchange_DOWN <- 2^(-SC_M1.M2_vs_Skin_repair$avg_log2FC)
SC_M1.M2_vs_Skin_repair$Ratio_pct1_pct2 <- (SC_M1.M2_vs_Skin_repair$pct.1)/(SC_M1.M2_vs_Skin_repair$pct.2)
SC_M1.M2_vs_Skin_repair$comparison <- "SC_M1.M2_vs_Skin_repair"
write.xlsx(SC_M1.M2_vs_Skin_repair, "Lists/Pseudotime_relevant/SC_M1.M2_vs_Skin_repair.xlsx")


####validation MAC_to_FB
Idents(s.k.total)<- s.k.total$celltype
FB_KC_total<- subset(s.k.total,idents = c("FB","KC"))
DefaultAssay(FB_KC_total)<-"RNA"
FB_KC_total[["percent.mt"]] <- PercentageFeatureSet(FB_KC_total, pattern = "^MT-")
FB_KC_total <- SCTransform(FB_KC_total, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
FB_KC_total[['integrated']] <- NULL

FB_MAC <- MAC_7s4k
DefaultAssay(FB_MAC)<-"RNA"
FB_MAC[["percent.mt"]] <- PercentageFeatureSet(FB_MAC, pattern = "^MT-")
FB_MAC <- SCTransform(FB_MAC, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
FB_MAC[['integrated']] <- NULL


FB_KC_compare.list<-list(FB_KC_total,FB_MAC)
FB_KC_compare.features <- SelectIntegrationFeatures(object.list = FB_KC_compare.list, nfeatures = 3000)
FB_KC_compare.list <- PrepSCTIntegration(FB_KC_compare.list, anchor.features = FB_KC_compare.features)
FB_KC_compare.list <- lapply(FB_KC_compare.list, RunPCA, verbose = F, features= FB_KC_compare.features)
FB_KC_compare.anchors<-FindIntegrationAnchors(FB_KC_compare.list,normalization.method = "SCT", anchor.features = FB_KC_compare.features, reduction = "rpca")
FB_KC_compare<-IntegrateData(anchorset=FB_KC_compare.anchors, normalization.method = "SCT")
FB_KC_compare <- RunPCA(FB_KC_compare)
ElbowPlot(FB_KC_compare, ndims = 50)
FB_KC_compare <- RunUMAP(FB_KC_compare, dims = 1:50)
FB_KC_compare <- FindNeighbors(FB_KC_compare, dims = 1:20)
#DefaultAssay(FB_KC_compare)<- "integrated"
#FB_KC_compare <- FindClusters(FB_KC_compare, resolution = 0.45)
UMAPPlot(FB_KC_compare, label=T, group.by="celltype")
DefaultAssay(FB_KC_compare)<- "RNA"
FB_KC_compare <- NormalizeData(FB_KC_compare)
FeaturePlot(FB_KC_compare,c("DCN","MRC1","MPZ","PLP1","CD68"))

#FB
Idents(FB_KC_compare)<- FB_KC_compare$celltype
MAC_to_FB_vs_FB_DEG <- FindMarkers(FB_KC_compare, ident.1 = "MAC_to_FB", ident.2=c("FB"))
MAC_to_FB_vs_FB_DEG$Foldchange_UP <- 2^(MAC_to_FB_vs_FB_DEG$avg_log2FC)
MAC_to_FB_vs_FB_DEG$Foldchange_DOWN <- 2^(-MAC_to_FB_vs_FB_DEG$avg_log2FC)
MAC_to_FB_vs_FB_DEG$Ratio_pct1_pct2 <- (MAC_to_FB_vs_FB_DEG$pct.1)/(MAC_to_FB_vs_FB_DEG$pct.2)
MAC_to_FB_vs_FB_DEG$celltype <- "MAC_to_FB"
MAC_to_FB_vs_FB_DEG$comparison <- "MAC_to_FB_vs_FB"
write.xlsx(MAC_to_FB_vs_FB_DEG, "Lists/DEG_MAC_to_FB_vs_FB.xlsx")


#KC
Idents(FB_KC_compare)<- FB_KC_compare$celltype
MAC_to_FB_vs_KC_DEG <- FindMarkers(FB_KC_compare, ident.1 = "MAC_to_FB", ident.2=c("KC"))
MAC_to_FB_vs_KC_DEG$Foldchange_UP <- 2^(MAC_to_FB_vs_KC_DEG$avg_log2FC)
MAC_to_FB_vs_KC_DEG$Foldchange_DOWN <- 2^(-MAC_to_FB_vs_KC_DEG$avg_log2FC)
MAC_to_FB_vs_KC_DEG$Ratio_pct1_pct2 <- (MAC_to_FB_vs_KC_DEG$pct.1)/(MAC_to_FB_vs_KC_DEG$pct.2)
MAC_to_FB_vs_KC_DEG$celltype <- "MAC_to_FB"
MAC_to_FB_vs_KC_DEG$comparison <- "MAC_to_FB_vs_KC"
write.xlsx(MAC_to_FB_vs_KC_DEG, "Lists/DEG_MAC_to_FB_vs_KC.xlsx")


#MAC
Idents(FB_KC_compare)<- FB_KC_compare$celltype
MAC_to_FB_vs_MAC_DEG <- FindMarkers(FB_KC_compare, ident.1 = "MAC_to_FB", ident.2=c("MAC-M1","MAC-M1/M2","MAC-M2"))
MAC_to_FB_vs_MAC_DEG$Foldchange_UP <- 2^(MAC_to_FB_vs_MAC_DEG$avg_log2FC)
MAC_to_FB_vs_MAC_DEG$Foldchange_DOWN <- 2^(-MAC_to_FB_vs_MAC_DEG$avg_log2FC)
MAC_to_FB_vs_MAC_DEG$Ratio_pct1_pct2 <- (MAC_to_FB_vs_MAC_DEG$pct.1)/(MAC_to_FB_vs_MAC_DEG$pct.2)
MAC_to_FB_vs_MAC_DEG$celltype <- "MAC_to_FB"
MAC_to_FB_vs_MAC_DEG$comparison <- "MAC_to_FB_vs_MAC"
write.xlsx(MAC_to_FB_vs_MAC_DEG, "Lists/DEG_MAC_to_FB_vs_MAC.xlsx")


#combi
Idents(FB_KC_compare)<- FB_KC_compare$celltype
MAC_to_FB_vs_MACFB_DEG <- FindMarkers(FB_KC_compare, ident.1 = "MAC_to_FB", ident.2=c("MAC-M1","MAC-M1/M2","MAC-M2","FB"))
MAC_to_FB_vs_MACFB_DEG$Foldchange_UP <- 2^(MAC_to_FB_vs_MACFB_DEG$avg_log2FC)
MAC_to_FB_vs_MACFB_DEG$Foldchange_DOWN <- 2^(-MAC_to_FB_vs_MACFB_DEG$avg_log2FC)
MAC_to_FB_vs_MACFB_DEG$Ratio_pct1_pct2 <- (MAC_to_FB_vs_MACFB_DEG$pct.1)/(MAC_to_FB_vs_MACFB_DEG$pct.2)
MAC_to_FB_vs_MACFB_DEG$celltype <- "MAC_to_FB"
MAC_to_FB_vs_MACFB_DEG$comparison <- "vs_MAC+FB"
write.xlsx(MAC_to_FB_vs_MACFB_DEG, "Lists/DEG_MAC_to_FB_vs_MAC+FB.xlsx")


#Barplot Freq up and downregulated Genes (FC>2)
#FB
FB_MAC_comp_df <- rbind(MAC_to_FB_vs_FB_DEG,MAC_to_FB_vs_MAC_DEG,MAC_to_FB_vs_KC_DEG)
FB_MAC_comp_df_subset <- subset(FB_MAC_comp_df, Foldchange_UP >=2 | Foldchange_DOWN >=2)
FB_MAC_comp_df_subset$Direction <- ifelse (FB_MAC_comp_df_subset$Foldchange_UP>=2,"UP","DOWN")
FB_MAC_comp_df_subset<- select(FB_MAC_comp_df_subset, Direction, comparison)
row.names(FB_MAC_comp_df_subset)<-NULL
FB_MAC_comp_df_subset$cluster<- as.factor(FB_MAC_comp_df_subset$comparison)
FB_MAC_comp_df_subset$cluster<- as.character(FB_MAC_comp_df_subset$Direction)
FB_MAC_comp_df_subset<- dplyr::count(FB_MAC_comp_df_subset, comparison,Direction) %>% ungroup()
FB_MAC_comp_df_subset$Direction <- factor(FB_MAC_comp_df_subset$Direction, levels= c("UP","DOWN"))
FB_MAC_comp_df_subset$comparison <- factor(FB_MAC_comp_df_subset$comparison, levels= c("MAC_to_FB_vs_MAC","MAC_to_FB_vs_FB","MAC_to_FB_vs_KC"))
pdf("Graphs/MAC/Barplot_nDEG_MAC-FB_only.pdf")
ggplot(data=FB_MAC_comp_df_subset, aes(x=comparison, y=n, fill=Direction)) +
  geom_bar(stat="identity", position=position_dodge())+ theme_classic() +
  geom_text_repel(aes(label = n, y=n+0.5),position=position_dodge(0.9), vjust=0.6, color = 'black', size = 4) +ggtitle("Comparison DEG_MAC-FB") +scale_fill_manual(values =color_DEG) +ylab("total number of DEG")
dev.off()

#Qualitycheck subset
pdf("Graphs/MAC/QC_pctMT_MAC_7s4k.pdf")
VlnPlot(MAC_7s4k, features = "percent.mt", cols = c(color_MAC))
dev.off()
pdf("Graphs/MAC/QC_nFeature_SC_7s4k.pdf")
VlnPlot(MAC_7s4k, features = "nFeature_RNA", cols = c(color_MAC))
dev.off()
pdf("Graphs/MAC/QC_nCount_SC_7s4k.pdf")
VlnPlot(MAC_7s4k, features = "nCount_RNA", cols = c(color_MAC))
dev.off()

p1<- VlnPlot(s.k.total,c("MMP9"), pt.size = 0, cols=color_tissue, split.by="tissue")+stat_summary(fun = "mean", geom = "crossbar")
p1<- p1 & theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p2<- VlnPlot(MAC_7s4k,c("MMP9"), pt.size = 0, cols=color_MAC)+stat_summary(fun = "mean", geom = "crossbar")
pdf("Graphs/Main/MAC_VlnPlot_MMP9+mean.pdf")
wrap_plots(p1,p2, nrow=2)
dev.off()

p1<- VlnPlot(s.k.total,c("CCL3"), pt.size = 0, cols=color_tissue, split.by="tissue")+stat_summary(fun = "mean", geom = "crossbar")
p1<- p1 & theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p2<- VlnPlot(MAC_7s4k,c("CCL3"), pt.size = 0, cols=color_MAC)+stat_summary(fun = "mean", geom = "crossbar")
pdf("Graphs/Main/MAC_VlnPlot_CCL3+mean.pdf")
wrap_plots(p1,p2, nrow=2)
dev.off()

Idents(s.k.total)="celltype"
p1<- VlnPlot(s.k.total,c("GAS6"), pt.size = 0, cols=color_tissue, split.by="tissue")+stat_summary(fun = "mean", geom = "crossbar")
p1<- p1 & theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p2<- VlnPlot(MAC_7s4k,c("GAS6"),y=5, pt.size = 0, cols=color_MAC)+ ggtitle("MAC-subset")+stat_summary(fun = "mean", geom = "crossbar")& stat_compare_means(comparisons = list(c("skin","keloid")) ,method = "t.test", label = "p.signif")
pdf("Graphs/Main/MAC_VlnPlot_GAS6+mean.pdf")
wrap_plots(p1,p2, nrow=2)
dev.off()
Idents(MAC_7s4k)="celltype"

p1<- VlnPlot(s.k.total,c("TNF"), pt.size = 0, cols=color_tissue, split.by="tissue")+stat_summary(fun = "mean", geom = "crossbar")
p1<- p1 & theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p2<- VlnPlot(MAC_7s4k,c("TNF"), pt.size = 0, cols=color_MAC)+stat_summary(fun = "mean", geom = "crossbar")
pdf("Graphs/Main/MAC_VlnPlot_TNF+mean.pdf")
wrap_plots(p1,p2, nrow=2)
dev.off()

p1<- VlnPlot(s.k.total,c("TIMP1"), pt.size = 0, cols=color_tissue, split.by="tissue")
p1<- p1 & theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p2<- VlnPlot(MAC_7s4k,c("TIMP1"), pt.size = 0, cols=color_MAC)
pdf("Graphs/Main/MAC_VlnPlot_TIMP1.pdf")
wrap_plots(p1,p2, nrow=2)
dev.off()

# Interaction Plot MACs
pdf("Graphs/MAC/MAC_VlnPlot_Interactionplot.pdf")
VlnPlot(MAC_7s4k,features=c("CCL3","TNF","MMP9"), group.by="tissue",pt.size = 0, cols=color_tissue)& ylim(0,6)
dev.off()

pdf("Graphs/Main/Total_VlnPlot_Interactionplot.pdf")
VlnPlot(s.k.total,features=c("CCL2","GAS6"), group.by="tissue",pt.size = 0, cols=color_tissue)& ylim(0,6)
dev.off()



#ECM and MAC
pdf("Graphs/MAC/Dotplot_MAC_ECM_Glycoproteins.pdf", width = 10, height = 30)
DotPlot(MAC_7s4k,features=ECM_Glycoproteins)+scale_color_gradient(low="grey", "high"=single_col) + theme(axis.text.x =  element_text(angle = 45, hjust = 1))+coord_flip() +theme(axis.text.x = element_text(angle = 90)) + ggtitle("ECM_Glycoproteins")
dev.off()

pdf("Graphs/MAC/Dotplot_MAC_Collagens.pdf", width = 10, height = 30)
DotPlot(MAC_7s4k,features=Collagens)+scale_color_gradient(low="grey", "high"=single_col) + theme(axis.text.x =  element_text(angle = 45, hjust = 1))+coord_flip() +theme(axis.text.x = element_text(angle = 90)) + ggtitle("Collagens")
dev.off()

pdf("Graphs/MAC/Dotplot_MAC_Proteoglycans.pdf", width = 10, height = 30)
DotPlot(MAC_7s4k,features=Proteoglycans)+scale_color_gradient(low="grey", "high"=single_col) + theme(axis.text.x =  element_text(angle = 45, hjust = 1))+coord_flip() +theme(axis.text.x = element_text(angle = 90)) + ggtitle("Proteoglycans")
dev.off()

pdf("Graphs/MAC/Dotplot_MAC_ECM_affiliated_Proteins.pdf", width = 10, height = 30)
DotPlot(MAC_7s4k,features=ECM_affiliated_Proteins)+scale_color_gradient(low="grey", "high"=single_col) + theme(axis.text.x =  element_text(angle = 45, hjust = 1))+coord_flip() +theme(axis.text.x = element_text(angle = 90)) + ggtitle("ECM_affiliated_Proteins")
dev.off()

pdf("Graphs/MAC/Dotplot_MAC_ECM_Regulators.pdf", width = 10, height = 35)
DotPlot(MAC_7s4k,features=ECM_Regulators)+scale_color_gradient(low="grey", "high"=single_col) + theme(axis.text.x =  element_text(angle = 45, hjust = 1))+coord_flip() +theme(axis.text.x = element_text(angle = 90)) + ggtitle("ECM_Regulators")
dev.off()

pdf("Graphs/MAC/Dotplot_MAC_Secreted_Factors.pdf", width = 10, height = 50)
DotPlot(MAC_7s4k,features=Secreted_Factors)+scale_color_gradient(low="grey", "high"=single_col) + theme(axis.text.x =  element_text(angle = 45, hjust = 1))+coord_flip() +theme(axis.text.x = element_text(angle = 90)) + ggtitle("Secreted_Factors")
dev.off()

pdf("Graphs/MAC/Dotplot_MAC_ECMassozGenes.pdf", width = 6, height = 20)
DotPlot(MAC_7s4k,features=c("CCL2","CCL3","CCL3L3","CCL4","CCL4L1","CCL4L2","CCL8","CCL11","CCL13","CCL17","CCL18","CCL19","CCL20","CCL22","CXCL1","CXCL12","CXCL14","CXCL2","CXCL3","CXCL8","EGFL7","IGF1","IL10","IL1B","IL1RN","S100A10","S100A11","S100A4","S100A6","S100A9","S100A8","TGFB1","TNF","VEGFA","SRGN","DCN","A2M","CST3","CSTB","CTSA","CTSB","CTSC","CTSD","CTSH","CTSL","CTSS","CTSZ","F13A1","MMP19","MMP9","PLAU","SERPINB1","SERPINB6","SERPINB9","TIMP1","TIMP2","EMILIN2","FGL2","GAS6","IGFBP7","SPARC","TGFBI","COL1A1","COL1A2","COL3A1","ANXA1","ANXA11","ANXA2","ANXA4","ANXA5","ANXA7","C1QA","C1QB","C1QC","CLEC10A","COLEC12","LGALS1","LGALS2","LGALS3"))+scale_color_gradient(low="grey", "high"=single_col) + theme(axis.text.x =  element_text(angle = 45, hjust = 1))+coord_flip() +theme(axis.text.x = element_text(angle = 90)) + ggtitle("Secreted_Factors")
dev.off()


###Macrophages classification
#M1 MHC-II high
StackedVlnPlot(MAC_7s4k,features=c("HLA-DMA","HLA-DOA","HLA-DPA1","HLA-DQA1","HLA-DRA","CD80","CD86","IFNAR1","TLR4"),group.by="celltype")
#M2a MHC-II low
StackedVlnPlot(MAC_7s4k,features=c("MRC1","CD86","IL1R1","CD163"),group.by="celltype")
#M2b MHC-II low
StackedVlnPlot(MAC_7s4k,features=c("IL10RB","IL6R","IL12RB1"),group.by="celltype")
#M2c 
StackedVlnPlot(MAC_7s4k,features=c("MRC1","TLR2","TLR8","CD163"),group.by="celltype")
#M2d 
StackedVlnPlot(MAC_7s4k,features=c("IL10RA","VEGFA","IL12RB2"),group.by="celltype")
#signal 
StackedVlnPlot(MAC_7s4k,features=c("CD68","CD163","CD14","FCGR3A","PLAC8","TIMD4","LYVE1","FOLR2","IGF1","IL1B","IFIT3","IRF7","IFIT1"),group.by="celltype")






###ELISA Auswertung
#CCL2
CCL2_data<- read.xlsx("D:/Direder/Projekte/D002_10xScarWars/D002_ELISA/050521/D002_Auswertung_CCL2.xlsx",3, header = F)
CCL2_data$X1 <- factor(x = CCL2_data$X1, levels = c("Skin", "Keloid"))

#Boxplots Percentage - t-test
pdf("Graphs/Main/BoxPlot_ELISA_CCL2.pdf", width = 5, height= 5)
ggplot(CCL2_data, aes(x= X1, y=X3, fill=X1)) +geom_boxplot() + theme_classic() +
  labs(title="CCL2-ELISA",x="", y = "pg CCL2/mg Protein")+ scale_fill_manual(values=color_tissue)+
  geom_signif(comparisons = list(c("Skin","Keloid")), test="t.test", test.args = list(var.equal = TRUE),
              map_signif_level=TRUE)
dev.off()

#Barplot Mean
pdf("Graphs/Main/BarPlot_ELISA_CCL2_Mean.pdf", width = 10, height= 10)
ggplot(CCL2_data, aes(x= X1, y=X4, fill=X1)) +geom_bar(stat = "identity") + theme_classic() +
  labs(title="CCL2-ELISA",x="", y = "pg CCL2/mg Protein")+ scale_fill_manual(values=color_tissue)
dev.off()

#MMP9
MMP9_data<- read.xlsx("D:/Direder/Projekte/D002_10xScarWars/D002_ELISA/050521/D002_Auswertung_MMP9.xlsx",3, header = F)
MMP9_data$X1 <- factor(x = MMP9_data$X1, levels = c("Skin", "Keloid"))

#Boxplots Percentage - t-test
pdf("Graphs/Main/BoxPlot_ELISA_MMP9.pdf", width = 5, height= 5)
ggplot(MMP9_data, aes(x= X1, y=X3, fill=X1)) +geom_boxplot() + theme_classic() + geom_jitter() +
  labs(title="MMP9-ELISA",x="", y = "ng MMP9/mg Protein")+ scale_fill_manual(values=color_tissue)+
  geom_signif(comparisons = list(c("Skin","Keloid")), test="t.test", test.args = list(var.equal = TRUE),
              map_signif_level=TRUE)
dev.off()

#Barplot Mean
pdf("Graphs/Main/BarPlot_ELISA_MMP9_Mean.pdf", width = 10, height= 10)
ggplot(MMP9_data, aes(x= X1, y=X4, fill=X1)) +geom_bar(stat = "identity") + theme_classic() +
  labs(title="MMP9-ELISA",x="", y = "ng MMP9/mg Protein")+ scale_fill_manual(values=color_tissue)
dev.off()



#_____________________________________________________________________
#####      Workspace

###########

############

VlnPlot(s.k.total,c("PLP1","MPZ","MATN2","PCSK2","CD9","NRXN1","PRKCDBP"), pt.size=0, cols=color_UMAP_short_double)
VlnPlot(SC_7s4k,c("PLP1","MPZ","MATN2","PCSK2","CD9","NRXN1","PRKCDBP"), pt.size=0, cols=color_SC)
VlnPlot(s.k.total,c("GFRA3","PPAP2A","DKK3","NDRG2","SCN7A"), pt.size=0, cols=color_UMAP_short_double)
VlnPlot(SC_7s4k,c("GFRA3","PPAP2A","DKK3","NDRG2","SCN7A"), pt.size=0, cols=color_SC)
VlnPlot(s.k.total,c("SELM","SCCPDH","PMP2","SMIM5","SGCE","APOD","NRN1"), pt.size=0, cols=color_UMAP_short_double)
VlnPlot(SC_7s4k,c("SELM","SCCPDH","PMP2","SMIM5","SGCE","APOD","NRN1"), pt.size=0, cols=color_SC)
VlnPlot(s.k.total,c("ART3","CCL2","TUBA1A","TUBA1B","COMT","CLU","TSPAN15","TMEM66","SOX2","PLLP"), pt.size=0, cols=color_UMAP_short_double)
VlnPlot(SC_7s4k,c("ART3","CCL2","TUBA1A","TUBA1B","COMT","CLU","TSPAN15","TMEM66","SOX2","PLLP"), pt.size=0, cols=color_SC)
VlnPlot(s.k.total,c("HSPA12A","TAX1BP3","MYOT","SNCA","HES1","SORBS2","XKR4","TMEM176B","MAL","CNP"), pt.size=0, cols=color_UMAP_short_double)
VlnPlot(SC_7s4k,c("HSPA12A","TAX1BP3","MYOT","SNCA","HES1","SORBS2","XKR4","TMEM176B","MAL","CNP"), pt.size=0, cols=color_SC)
VlnPlot(s.k.total,c("CBR1","GATM","CKB","ITM2C","CD59","CADM2","GPR155","RGCC","C8orf4","ARPC1B","COL28A1","FGL2"), pt.size=0, cols=color_UMAP_short_double)
VlnPlot(SC_7s4k,c("CBR1","GATM","CKB","ITM2C","CD59","CADM2","GPR155","RGCC","C8orf4","ARPC1B","COL28A1","FGL2"), pt.size=0, cols=color_SC)
VlnPlot(s.k.total,c("MAGED2","PLEKHB1","APLP2","RAPGEF5","C1orf63","OLFML2A","VWA1","FXYD1","IGFBP7","LTBP4","CAPS","U2AF1"), pt.size=0, cols=color_UMAP_short_double)
VlnPlot(SC_7s4k,c("MAGED2","PLEKHB1","APLP2","RAPGEF5","C1orf63","OLFML2A","VWA1","FXYD1","IGFBP7","LTBP4","CAPS","U2AF1"), pt.size=0, cols=color_SC)
VlnPlot(s.k.total,c("DDIT3","DYNLRB1","LGI4","SEPP1","PCDH20","DNAJB1","CNPY2","CYR61","IQGAP2","FEZ1","JUNB","JUN","SEP15"), pt.size=0, cols=color_UMAP_short_double)
VlnPlot(SC_7s4k,c("DDIT3","DYNLRB1","LGI4","SEPP1","PCDH20","DNAJB1","CNPY2","CYR61","IQGAP2","FEZ1","JUNB","JUN","SEP15"), pt.size=0, cols=color_SC)
DotPlot(SC_7s4k,features = c("IGFBP5","NOV","PMEPA1","PTPRZ1","TGFBI","NES","ITGB1","PTN","COL7A1","S100A6","TNFRSF12A","CDH6","COL18A1","TNC","IGFBP3","COL8A1","DAG1","TNFAIP6","CDH2","COL4A2","ENC1","SEMA3C","PDGFA","COL4A1","FRMD4A","S100A10","CCND1","CAVIN3","S100A16","EMP3","HSBP1","SH3PXD2A","COL5A3","LSM7","ZEB2","PLPP1","RHOB","PPP1R14B","PPFIBP1","MCAM","QKI","ARPC5","SEM1","HBEGF","PLS3","CBX5","TMSB10","SCARB2","HMGB1","S100A11","BEX3","VCL","SPPL2A","C4orf48","CALU","SPARCL1","CALM1","DYNLT1","FAM96B","ANKRD10","ANXA5","NDUFA4","MAP4","ATOX1","SERPINE2","MYL6","TXNDC17","RCN2","GNG5"))
VlnPlot(s.k.total,c("IGFBP3","IGFBP5","IGFBP7"), pt.size=0, cols=color_UMAP_short_double)

pdf("Graphs/MAC/M2_Mac.pdf", height=10, width = 10)
FeaturePlot(s.k.total, c())&scale_color_gradient(low="grey", "high"=single_col)
FeaturePlot(s.k.total, c())&scale_color_gradient(low="grey", "high"=single_col) 
dev.off()
pdf("Graphs/MAC/M1_Mac.pdf", height =5, width = 5)
FeaturePlot(s.k.total, c())&scale_color_gradient(low="grey", "high"=single_col) 
FeaturePlot(s.k.total, c())&scale_color_gradient(low="grey", "high"=single_col)
FeaturePlot(s.k.total, c())&scale_color_gradient(low="grey", "high"=single_col)
FeaturePlot(s.k.total, c())&scale_color_gradient(low="grey", "high"=single_col)
FeaturePlot(s.k.total, c())&scale_color_gradient(low="grey", "high"=single_col)
dev.off()
######################




Test1 <- SC_7s4k
Test1 <- RunVelocity(object = Test1, deltaT = 1, kCells = 25, fit.quantile = 0.02)
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = Test1)))
names(x = ident.colors) <- levels(x = Test1)
cell.colors <- ident.colors[Idents(object = Test1)]
names(x = cell.colors) <- colnames(x = Test1)
show.velocity.on.embedding.cor(emb = Embeddings(object = Test1, reduction = "umap"), vel = Tool(object = Test1, 
                                                                                                slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)





DotPlot(s.k.total, features = c("RPL39"
                                ,"ATP5F1E"
                                ,"RPL41"
                                ,"COL1A1"
                                ,"PMEPA1"
                                ,"IGFBP5"
                                ,"ITGB1"
                                ,"IGFBP3"
                                ,"H3F3A"
                                ,"SEM1"
                                ,"ELOB"
                                ,"TNC"
), split.by="tissue")

DotPlot(s.k.total, features = c("COL5A3"
                                ,"DAG1"
                                ,"CDH2"
                                ,"ENC1"
                                ,"NOV"
                                ,"COL7A1"
                                ,"PTPRZ1"
                                ,"NES"
                                ,"TGFBI","IGFBP3"
), split.by="tissue")                           
DotPlot(s.k.total, features = c("S100A16"
                                ,"COL18A1"
                                ,"CCND1"
                                ,"CBX5"
                                ,"PHLDA2"
                                ,"TNFAIP6"
                                ,"SERPINE2"
                                ,"PLPP1"
                                ,"CDH6"
                                ,"TNFRSF12A","ELN","COL8A1","PDGFA","PTN"
                                ,"IGFBP5","PMEPA1"), split.by="tissue")






VlnPlot(s.k.total,c("IGFBP5"
                    ,"NOV"
                    ,"PMEPA1"
                    ,"PTPRZ1"
                    ,"TGFBI"
                    ,"NES"
                    ,"ITGB1"
                    ,"PTN"
                    ,"COL7A1"
                    ,"S100A6"
                    ,"TNFRSF12A"
                    ,"GAS7"
                    ,"CDH6"
                    ,"COL18A1"), pt.size = 0)


p1<- VlnPlot(s.k.total,"JUN", pt.size = 0, cols=color_tissue, split.by="tissue")
p1<- p1 + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p2<- VlnPlot(SC_7s4k,"JUN", pt.size = 0, cols=color_SC)
p3<- VlnPlot(s.k.total,"JUNB", pt.size = 0, cols=color_tissue, split.by="tissue")
p3<- p3 + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p4<- VlnPlot(SC_7s4k,"JUNB", pt.size = 0, cols=color_SC)
p5<- VlnPlot(s.k.total,"JUND", pt.size = 0, cols=color_tissue, split.by="tissue")
p5<- p5 + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p6<- VlnPlot(SC_7s4k,"JUND", pt.size = 0, cols=color_SC)
p7<- VlnPlot(s.k.total,"FOS", pt.size = 0, cols=color_tissue, split.by="tissue")
p7<- p7 + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p8<- VlnPlot(SC_7s4k,"FOS", pt.size = 0, cols=color_SC)
p9<- VlnPlot(s.k.total,"FOSB", pt.size = 0, cols=color_tissue, split.by="tissue")
p9<- p9 + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p10<- VlnPlot(SC_7s4k,"FOSB", pt.size = 0, cols=color_SC)
p11<- VlnPlot(s.k.total,"FOSL1", pt.size = 0, cols=color_tissue, split.by="tissue")
p11<- p11 + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p12<- VlnPlot(SC_7s4k,"FOSL1", pt.size = 0, cols=color_SC)
p13<- VlnPlot(s.k.total,"FOSL2", pt.size = 0, cols=color_tissue, split.by="tissue")
p13<- p13 + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short))
p14<- VlnPlot(SC_7s4k,"FOSL2", pt.size = 0, cols=color_SC)
pdf("Graphs/Side/JUN.pdf", width = 10, height = 30)
wrap_plots(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14, nrow=7)
dev.off()


############################
####Subset Dendritic cells####
############################
Idents(s.k.total)<- s.k.total$celltype
DC_7s4k <- subset(s.k.total,idents = "DC")
DefaultAssay(DC_7s4k)<-"RNA"
DC_7s4k[["percent.mt"]] <- PercentageFeatureSet(DC_7s4k, pattern = "^MT-")
DC_7s4k <- SCTransform(DC_7s4k, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
DC_7s4k <- RunPCA(DC_7s4k, npcs = 50)
ElbowPlot(DC_7s4k, ndims = 50)
DC_7s4k <- RunUMAP(DC_7s4k, dims = 1:30)
DC_7s4k <- FindNeighbors(DC_7s4k, dims = 1:30)
#DefaultAssay(DC_7s4k)<- "integrated"
DC_7s4k <- FindClusters(DC_7s4k, resolution = 0.2)
UMAPPlot(DC_7s4k, label=T)
DefaultAssay(DC_7s4k)<- "RNA"
DC_7s4k <- NormalizeData(DC_7s4k)
VlnPlot(DC_7s4k,c("MAP1A","POSTN","CD163","MRC1","MSR1","CD80","CD86","FCGR1A","FCGR3A","FCGR2A"), split.by = "tissue")
VlnPlot(DC_7s4k,c("FCN1","SPP1","FABP4","CD68"), split.by = "tissue")
FeaturePlot(DC_7s4k,"HPGDS")
UMAPPlot(s.k.total)
FeaturePlot(s.k.total,"TPSB2", split.by = "tissue")
VlnPlot(SC_7s4k,"NINJ1", pt.size = 0, split.by = "tissue")


FeaturePlot(SC_7s4k,"CDH19")
#____________________________________________________________________
#######testarea start
Idents(s.k.total) <- s.k.total$celltype
DefaultAssay(SC_7s4k)<- "RNA"
UMAPPlot(SC_7s4k, split.by="sample")

FeaturePlot(SC_7s4k,c("THY1","S100B"),blend = T, blend.threshold = 0.0, cols = c("lightgrey","#990F26","#54990F"), pt.size = 1.5, order=T)


FeaturePlot(s.k.total,c("WNT3A","SNAIL2"))
FeaturePlot(SC_7s4k,c("WNT3A"),order = T, pt.size = 1, split.by="sample")

VlnPlot(SC_7s4k, "ZEB1", split.by = "tissue", pt.size = 0)
DotPlot(SC_7s4k, features=c("BEX3"), group.by = "seurat_clusters") + coord_flip()+scale_color_gradient(low="grey", "high"="red")+theme(axis.text.x = element_text(angle = 45))
DotPlot(SC_7s4k, features=NatureSC_marker, group.by = "seurat_clusters") + coord_flip()+scale_color_gradient(low="grey", "high"="red")+theme(axis.text.x = element_text(angle = 45))
StackedVlnPlot(SC_7s4k, c("S100B","IGFBP5","NOV","PTPRZ1","PMEPA1","NES","TGFBI","CDH19","S100A4","CRYAB","NTM","ITGB1","COL7A1"), pt.size = 0)
FeaturePlot(s.k.total, "TRIM28", split.by= "tissue", order=T, min.cutoff = 2)

Idents(SC_7s4k)<- SC_7s4k$celltype
Test1 <- subset(SC_7s4k,idents = "SC-intermed")
DefaultAssay(Test1)<-"RNA"
Test1[["percent.mt"]] <- PercentageFeatureSet(Test1, pattern = "^MT-")
Test1 <- SCTransform(Test1, vars.to.regress = "percent.mt",verbose = F)
Test1 <- RunPCA(Test1, npcs = 50)
ElbowPlot(Test1, ndims = 20)
Test1 <- RunUMAP(Test1, dims = 1:10)
Test1 <- FindNeighbors(Test1, dims = 1:10)
#DefaultAssay(Test1)<- "integrated"
Test1 <- FindClusters(Test1, resolution = 0.9)
UMAPPlot(Test1, label=T, split.by="sample")
DefaultAssay(Test1)<- "RNA"
Test1 <- NormalizeData(Test1)
FeaturePlot(Test1,c("XIST","GAP43","NCAM1","L1CAM","GFAP","RAN","A5E3","EGR1"))

FeaturePlot(Test1,c("GAP43","EGR1","SPARC","OMD","IBSP", "MAP5"))
FeaturePlot(Test1,c("NGFR","GAP43","NCAM1","L1CAM","GFAP","RAN","A5E3","EGR1"))
FeaturePlot(Test1,c("AAK1","ITGA4","CDH2","GFAP","S100B","CDH19","NOTCH1","SOX10","POU3F1"))
FeaturePlot(Test1,c("EGR1","IBSP","GAP43","L1CAM","CALB2","TPM2","CSRP2","GAPDH"))
FeaturePlot(Test1,c("IBSP","CALB2","IGFBP2","PLCG2","XIST","ITIH5","NRXN1","SEC61G"))
FeaturePlot(Test1,c("XIST","POU3F1","SOX2","IGFBP2","POSTN"))

table(Test1$celltype)


#Clustermarker + Heatmap total
Test1.clustermarker<-FindAllMarkers(Test1,  assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
write.xlsx(Test1.clustermarker, "Lists/Clustermarker_Test1.xlsx")


DotPlot(Test1, features= SC_marker_Toma_etal._32349983 , group.by = "seurat_clusters") + coord_flip()+scale_color_gradient(low="grey", "high"="red")+theme(axis.text.x = element_text(angle = 45))


#### Create a Monocle CDS Object
# Project PC dimensions to whole data set
SC_convert <- ProjectDim(s.k.total, reduction = "pca")

# Create an expression matrix
expression_matrix <- SC_convert@assays$RNA@counts

# Get cell metadata
cell_metadata <- SC_convert@meta.data
if (all.equal(colnames(expression_matrix), rownames(cell_metadata))) {
  print(sprintf("Cell identifiers match"))
} else {
  print(sprintf("Cell identifier mismatch - %i cells in expression matrix, %i cells in metadata",
                ncol(expression_matrix), nrow(cell_metadata)))
  print("If the counts are equal, sort differences will throw this error")
}

# get gene annotations
gene_annotation <- data.frame(gene_short_name = rownames(SC_convert@assays$RNA), row.names = rownames(SC_convert@assays$RNA))
if (all.equal(rownames(expression_matrix), rownames(gene_annotation))) {
  print(sprintf("Gene identifiers all match"))
} else {
  print(sprintf("Gene identifier mismatch - %i genes in expression matrix, %i gene in gene annotation",
                nrow(expression_matrix), nrow(gene_annotation)))
  print("If the counts are equal, sort differences will throw this error")
}

# Seurat-derived CDS
Test.cds <- new_cell_data_set(expression_matrix,
                              cell_metadata = cell_metadata,
                              gene_metadata = gene_annotation)

# Transfer Seurat embeddings
# Note that these may be calculated on the Integrated object, not the counts
#   and thus will involve fewer genes
reducedDim(Test.cds, type = "PCA") <- s.k.total@reductions$pca@cell.embeddings 
Test.cds@preprocess_aux$prop_var_expl <- s.k.total@reductions$pca@stdev
plot_pc_variance_explained(Test.cds)

# Transfer Seurat UMAP embeddings
Test.cds@int_colData@listData$reducedDims$UMAP <- s.k.total@reductions$umap@cell.embeddings
#    plot_cells()

# Copy cluster info from Seurat
Test.cds@clusters$UMAP_so$clusters <- s.k.total@meta.data$gt_tp_cell_type_integrated_.0.9

Test.cds <- cluster_cells(Test.cds, reduction_method = "UMAP")

# Fix from https://gitmemory.com/cole-trapnell-lab
rownames(Test.cds@principal_graph_aux$UMAP$dp_mst) <- NULL
colnames(Test.cds@int_colData@listData$reducedDims$UMAP) <- NULL


cds_3d <- reduce_dimension(Test.cds, max_components = 3)
cds_3d <- cluster_cells(cds_3d, resolution = 0.00002)
#cds_3d@clusters@listData$UMAP$clusters<-s.k.total@active.ident
cds_3d <- learn_graph(cds_3d, close_loop = F, use_partition = T, learn_graph_control = list(prune_graph=F, ncenter=14))
t1<-plot_cells(cds_3d, color_cells_by = "cluster")
t2<-plot_cells(cds_3d, color_cells_by = "celltype")
wrap_plots(t1,t2)
plot_cells(cds_3d, genes=c("S100B"), show_trajectory_graph = F)
cds_3d <- order_cells(cds_3d)
plot_cells_3d(cds_3d, color_cells_by = "cluster", cell_size=25)
plot_cells_3d(cds_3d, genes=c("S100B"))
cds_3d_subset <- choose_cells(cds_3d)

plot_cells_3d(cds_3d_subset, color_cells_by = "cluster", cell_size=25, show_trajectory_graph = F)
#cds_3d_subset <- reduce_dimension(cds_3d_subset, max_components = 3)
cds_3d_subset <- cluster_cells(cds_3d_subset, resolution = 0.009)
cds_3d_subset <- learn_graph(cds_3d_subset, close_loop = F, use_partition = T, learn_graph_control = list(prune_graph=F, ncenter=15))
plot_cells_3d(cds_3d_subset, color_cells_by = "cluster", cell_size=30, show_trajectory_graph = F)
plot_cells_3d(cds_3d_subset, genes=c("MKI67"), show_trajectory_graph = F)
t1<-plot_cells(cds_3d_subset, color_cells_by = "cluster")
t2<-plot_cells(cds_3d_subset, color_cells_by = "celltype")
wrap_plots(t1,t2)

cds_3d_1 <- reduce_dimension(SC_7s4k.cds, max_components = 3)
cds_3d_1 <- cluster_cells(cds_3d_1, resolution = 0.00002, partition_qval = 1)
#cds_3d_1@clusters@listData$UMAP$clusters<-s.k.total@active.ident
cds_3d_1 <- learn_graph(cds_3d_1, close_loop = F, use_partition = T, learn_graph_control = list(prune_graph=F, ncenter=30))
t1<-plot_cells(cds_3d_1, color_cells_by = "cluster")
t2<-plot_cells(cds_3d_1, color_cells_by = "celltype")
wrap_plots(t1,t2)
plot_cells(cds_3d_1, genes=c("S100B"), show_trajectory_graph = F)
cds_3d_1 <- order_cells(cds_3d_1)
plot_cells_3d(cds_3d_1, color_cells_by = "cluster", cell_size=25)
plot_cells_3d(cds_3d_1, genes=c("CCL2"))
cds_3d_1_subset <- choose_cells(cds_3d_1)

plot_cells_3d(cds_3d_1, color_cells_by = "tissue", cell_size=25, show_trajectory_graph = F)
plot_cells_3d(cds_3d_1_subset, genes=c("S100B"), show_trajectory_graph = F)
cds_3d_1<- order_cells(cds_3d_1)

plot_cells_3d(cds_3d_subset, genes=c("MKI67"), show_trajectory_graph = F, cell_size=100)
plot_cells_3d(cds_3d_1, show_trajectory_graph = T,color_cells_by = "pseudotime")
plot_cells_3d(cds_3d_1, genes=c("MPZ"), show_trajectory_graph = T, cell_size=100)


cds_subs <- choose_graph_segments(SC_7s4k.cds)
cds_subs <- order_cells(cds_subs)
gene_fits <- fit_models(cds_3d_1, model_formula_str = "~celltype")
genes_of_interest <- c("PLP1","MKI67","TOP2A","IBSP")
SC_lineate_cds <- cds_3d_1[rowData(cds_3d_1)$gene_short_name%in%genes_of_interest,
                           colData(cds_3d_1)$cell.type %in% c("SC-Prolif","SC-X1","SC-X2","SC-Myel")]
gene_fits <- fit_models(cds_3d_1, model_formula_str = "~celltype")
plot_genes_in_pseudotime(SC_lineate_cds, color_cells_by = "celltype", min_expr = 0.5)


#test-reintegration
UMAPPlot(s.k.total)
Idents(s.k.total)<- s.k.total$celltype
Test1 <- subset(s.k.total,idents = c("FB1","FB2","FB3","SMC/PC","KC1","KC2","MAC","DC","MEL","EC1","EC2","LEC","TC","ERY"))
DefaultAssay(Test1)<-"SCT"
Test1[["percent.mt"]] <- PercentageFeatureSet(Test1, pattern = "^MT-")
Test1 <- SCTransform(Test1, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
Test1[['integrated']] <- NULL

Test2 <- SC_7s4k
DefaultAssay(Test2)<-"RNA"
Test2[["percent.mt"]] <- PercentageFeatureSet(Test2, pattern = "^MT-")
Test2 <- SCTransform(Test2, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
Test2[['integrated']] <- NULL

Test3.list<-list(Test1,Test2)
Test3.features <- SelectIntegrationFeatures(object.list = Test3.list, nfeatures = 3000)
Test3.list <- PrepSCTIntegration(Test3.list, anchor.features = Test3.features)
Test3.list <- lapply(Test3.list, RunPCA, verbose = F, features= Test3.features)
Test3.anchors<-FindIntegrationAnchors(Test3.list,normalization.method = "SCT", anchor.features = Test3.features, reduction = "rpca")
Test3<-IntegrateData(anchorset=Test3.anchors, normalization.method = "SCT")
Test3 <- RunPCA(Test3, npcs = 60)
pdf("Graphs/Elbowplot_skin_keloid_7s4k.pdf")
ElbowPlot(Test3, ndims = 50)
dev.off()
Test3 <- RunUMAP(Test3, dims = 1:35)
Test3 <- FindNeighbors(Test3, dims = 1:35)
#DefaultAssay(Test3)<- "integrated"
Test3 <- FindClusters(Test3, resolution = 0.3)
UMAPPlot(Test3, label=T, group.by="celltype")
DefaultAssay(Test3)<- "RNA"
Test3 <- NormalizeData(Test3)
FeaturePlot(Test3,"PMEL", split.by = "tissue")


Test1<-FindMarkers(s.k.total, ident.1="10", ident.2=c("0","1","2","3","4","5","6","7","8","9","11","12"), assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
head(Test1)
write.xlsx(Test1,"Test1.xlsx")

Idents(s.k.total)<-"celltype.tissue"
VlnPlot(s.k.total,features = "GAS6", pt.size=0)



VlnPlot(s.k.total,features=c("GAS6"),pt.size=0, split.by = "celltype")
VlnPlot(MAC_7s4k,features=c("GAS6"), group.by = "celltype")

vln_df = data.frame(GAS6 = s.k.total[["RNA"]]@data["GAS6",], cluster = s.k.total$celltype.tissue)
ggplot(vln_df, aes(x = cluster, y = GAS6)) + geom_violin(aes(fill = cluster), trim=TRUE, scale = "width")
vln_df.2 <- vln_df
noise <- rnorm(n = length(x = vln_df.2[, "GAS6"])) / 10
vln_df.2$GAS6 <- vln_df.2$GAS6  + noise

# violin plot with noise
ggplot(vln_df.2, aes(x = cluster, y = GAS6, fill = cluster)) + geom_violin( adjust =1,trim=TRUE, scale = "width")

DotPlot(SC_7s4k,features=c("SHH","BTC","BDNF","GDNF"))


DotPlot(s.k.total,features=c("ANGT","APOA","ASPN","COL2A1","COL5A3","COL11A1","CRTAP","ECM2","GT251","HABP2","LAMA3","MXRA5","AMD","PGCA","PGS1","PLOD2","PODN","S10A7","TIMP1"), split.by="tissue")


######testarea end
#___________________________________________________________________

FeaturePlot(SC_KEL_VS_NF,features="ST3GAL5")
Idents(s.k.total)="celltype"
FeaturePlot(s.k.total,features=c("TGFB1","TGFB2","TGFB3","TGFBR1"), split.by="tissue")
FeaturePlot(SC_7s4k,c("TGFB1","TGFB2","TGFB3","TGFBR1"), order=T)

DotPlot(s.k.total, features=c("SOX2","NCAM1","PROM1","VCAN","NGFR","FN1","PDGFRA","CDH7","B3GAT1"), group.by="celltype.tissue")+coord_flip()
DotPlot(SC_7s4k, features=c("SOX2","NCAM1","PROM1","VCAN","NGFR","FN1","PDGFRA","CDH7","B3GAT1","ITGB1"))+coord_flip()

        
############################
####Subset Schwann Cells####
############################
Idents(test1)<- test1$celltype
SC_test1 <- subset(test1,idents = "SC")
DefaultAssay(SC_test1)<-"RNA"
SC_test1[["percent.mt"]] <- PercentageFeatureSet(SC_test1, pattern = "^MT-")
SC_test1 <- SCTransform(SC_test1, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
SC_test1 <- RunPCA(SC_test1, npcs = 50)
ElbowPlot(SC_test1, ndims = 25)
SC_test1 <- RunUMAP(SC_test1, dims = 1:14)
SC_test1 <- FindNeighbors(SC_test1, dims = 1:14)
#DefaultAssay(SC_test1)<- "SCT"
SC_test1 <- FindClusters(SC_test1, resolution = 0.2)
UMAPPlot(SC_test1, label=T)
DefaultAssay(SC_test1)<- "RNA"
SC_test1 <- NormalizeData(SC_test1)
FeaturePlot(SC_test1,c("DCN","TOP2A","SELE","MBP","PLP1","S100B"))


f1<-FeaturePlot(s.k.total,"PPIC", order=T)
f2<-FeaturePlot(SC_7s4k,"PPIC", order=T)
wrap_plots(f1,f2)

DotPlot(SC_7s4k,features = c("ETV5","ERBB3","L1CAM","ITGA4","TFAP2A","CDH2","CDH19"))
############################
####Subset Melanocytes####
############################
Idents(s.k.total)<- s.k.total$celltype
MEL_7s4k <- subset(s.k.total,idents = "MEL")
DefaultAssay(MEL_7s4k)<-"RNA"
MEL_7s4k[["percent.mt"]] <- PercentageFeatureSet(MEL_7s4k, pattern = "^MT-")
MEL_7s4k <- SCTransform(MEL_7s4k, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
MEL_7s4k <- RunPCA(MEL_7s4k, npcs = 50)
ElbowPlot(MEL_7s4k, ndims = 25)
MEL_7s4k <- RunUMAP(MEL_7s4k, dims = 1:10)
MEL_7s4k <- FindNeighbors(MEL_7s4k, dims = 1:10)
#DefaultAssay(MEL_7s4k)<- "MELT"
MEL_7s4k <- FindClusters(MEL_7s4k, resolution = 0.6)
UMAPPlot(MEL_7s4k, label=T, split.by="tissue")
DefaultAssay(MEL_7s4k)<- "RNA"
MEL_7s4k <- NormalizeData(MEL_7s4k)
FeaturePlot(MEL_7s4k,c("DCN","TOP2A","SELE","MBP","PLP1","S100B","NGFR"))

pdf("Graphs/Main/Featureblend_DCN_S100B_SC_7s4k.pdf", width=15, height=7)
FeaturePlot(SC_7s4k,c("NGFR","S100B"),blend = T, blend.threshold = 0.5, cols = c("lightgrey","#990F26","#54990F"), pt.size = 1.5, order=T, min.cutoff = 0, max.cutoff = 1)
FeaturePlot(SC_7s4k,c("DCN","S100B"),blend = T, blend.threshold = 0.5, cols = c("lightgrey","#990F26","#54990F"), pt.size = 1.5, order=T, min.cutoff = 0, max.cutoff = 3)
FeaturePlot(SC_7s4k,c("DCN","NGFR"),blend = T, blend.threshold = 0.5, cols = c("lightgrey","#990F26","#54990F"), pt.size = 1.5, order=T, min.cutoff = 0, max.cutoff = 0.5)
dev.off()

pdf("Graphs/Main/Featureblend_THY1_S100B_SC_7s4k.pdf", width=15, height=7)
FeaturePlot(SC_7s4k,c("NGFR","S100B"),blend = T, blend.threshold = 0.5, cols = c("lightgrey","#990F26","#54990F"), pt.size = 1.5, order=T, min.cutoff = 0, max.cutoff = 1)
FeaturePlot(SC_7s4k,c("THY1","S100B"),blend = T, blend.threshold = 0.5, cols = c("lightgrey","#990F26","#54990F"), pt.size = 1.5, order=T, min.cutoff = 0, max.cutoff = 1)
FeaturePlot(SC_7s4k,c("THY1","NGFR"),blend = T, blend.threshold = 0.5, cols = c("lightgrey","#990F26","#54990F"), pt.size = 1.5, order=T, min.cutoff = 0, max.cutoff = 0.5)
dev.off()

pdf("Graphs/Main/Featureblend_SELE_S100B_SC_7s4k.pdf", width=15, height=7)
FeaturePlot(SC_7s4k,c("NGFR","S100B"),blend = T, blend.threshold = 0.5, cols = c("lightgrey","#990F26","#54990F"), pt.size = 1.5, order=T, min.cutoff = 0, max.cutoff = 1)
FeaturePlot(SC_7s4k,c("SELE","S100B"),blend = T, blend.threshold = 0.5, cols = c("lightgrey","#990F26","#54990F"), pt.size = 1.5, order=T, min.cutoff = 0, max.cutoff = 1)
FeaturePlot(SC_7s4k,c("SELE","NGFR"),blend = T, blend.threshold = 0.5, cols = c("lightgrey","#990F26","#54990F"), pt.size = 1.5, order=T, min.cutoff = 0, max.cutoff = 0.5)
dev.off()

pdf("Graphs/Main/Featureblend_ICAM1_S100B_SC_7s4k.pdf", width=15, height=7)
FeaturePlot(SC_7s4k,c("NGFR","S100B"),blend = T, blend.threshold = 0.5, cols = c("lightgrey","#990F26","#54990F"), pt.size = 1.5, order=T, min.cutoff = 0, max.cutoff = 1)
FeaturePlot(SC_7s4k,c("ICAM1","S100B"),blend = T, blend.threshold = 0.5, cols = c("lightgrey","#990F26","#54990F"), pt.size = 1.5, order=T, min.cutoff = 0, max.cutoff = 1)
FeaturePlot(SC_7s4k,c("ICAM1","NGFR"),blend = T, blend.threshold = 0.5, cols = c("lightgrey","#990F26","#54990F"), pt.size = 1.5, order=T, min.cutoff = 0, max.cutoff = 0.5)
dev.off()


FeaturePlot(SC_7s4k,"SOX10",split.by = "tissue", order=T)

#genes: Response to oxidative stress
ros_gene<- c("ABCC2","ABL1","ADA","ADAM9","ADIPOQ","ADNP2","ADPRHL2","AGAP3","AIF1","AKR1C3","AKT1","ALDH3B1","ALS2","ANGPTL7","ANXA1","APEX1","APOD","APOE","APP",
             "APTX","AQP1","AQP8","AREG","ARG1","ATOX1","ATP7A","ATRN","AXL","BAD","BAK1","BCL2","BNIP3","BTK","CA3","CAT","CCL19","CCR7","CCS","CD36","CD38","CDK1",
             "CDK2","CHRNA4","CHUK","CLN8","COL1A1","COQ7","CPEB2","CRYAB","CRYGD","CYCS","CYGB","CYP11A1","CYP1B1","CYP2E1","DDIAS","DGKK","DHCR24","DHRS2","DPEP1",
             "DUOX1","DUOX2","DUSP1","ECT2","EDN1","EGFR","EGLN1","EPX","ERCC1","ERCC2","ERCC3","ERCC6","ERCC8","ERO1L","ETFDH","ETS1","ETV5","EZH2","FABP1","FER","FGF8",
             "FKBP1B","FOS","FOSL1","FOXO1","FOXO3","FXN","G6PD","GAB1","GATM","GCLC","GCLM","GJB2","GLRX2","GNAO1","GPX1","GPX2","GPX3","GPX4","GPX5","GPX6","GPX7","GPX8",
             "GSR","GSS","GSTP1","HAO1","HBB","HDAC6","HMOX1","HP","HTRA2","HYAL1","HYAL2","IDH1","IL18BP","IL6","IPCEF1","JAK2","JUN","KCNC2","KDM6B","KLF2","KLF4",
             "KPNA4","KRT1","LCK","LCN2","LIAS","LIG1","LONP1","LPO","LRRK2","MAP2K1","MAP3K5","MAPK7","MAPK8","MB","MBL2","MDM2","MGMT","MGST1","MICB","MMP14","MMP3",
             "MPO","MPV17","MSRA","MSRB1","MSRB2","MSRB3","MT3","MTF1","NAPRT","NDUFA12","NDUFA6","NDUFB4","NDUFS2","NDUFS8","NEIL1","NET1","NFE2L1","NFE2L2","NFKB1",
             "NGF","NR4A2","NR4A3","NUDT1","OGG1","OLR1","OSER1","OXR1","OXSR1","P4HB","PARK2","PARK7","PARP1","PAX2","PCGF2","PDCD10","PDGFRB","PDK2","PDLIM1","PINK1",
             "PKD2","PLEKHA1","PLK3","PNKP","PNPT1","PON2","PPARGC1A","PPARGC1B","PPIF","PPP2CB","PRDX1","PRDX2","PRDX3","PRDX5","PRDX6","PRKAA1","PRKD1","PRKRA","PRNP",
             "PSEN1","PSIP1","PSMB5","PTGS1","PTGS2","PTK2B","PTPRK","PTPRN","PXDN","PXDNL","PXN","PYCR1","RAD52","RBM11","RCAN2","RELA","RGS14","RHOB","ROMO1","RRM2B",
             "S100A7","SCARA3","SCGB1A1","SDC1","SELK","SEPP1","SETX","SFTPC","SGK2","SIRT1","SIRT2","SLC23A2","SLC25A24","SLC7A11","SNCA","SOD1","SOD2","SOD3","SPHK1",
             "SRC","SRXN1","STAR","STAT1","STAT6","STC2","STK24","STK25","STK26","STX2","STX4","TACR1","TAT","THBS1",
             "TMEM161A","TNFAIP3","TOR1A","TP53INP1","TPM1","TPO","TRAF2","TRPA1","TRPC6","TRPM2","TXN","TXN2","TXNDC2","TXNDC8","TXNIP","TXNL1","TXNRD1","TXNRD2","UCN", 
             "UCP2","UCP3","VIMP","VKORC1L1","VNN1","VRK2","WRN","XPA","ZNF277","ZNF580")

Idents(s.k.total)="celltype.tissue"

pdf("Graphs/Side/Dotplot_ros_gene.pdf", width = 10, height = 60)
DotPlot(s.k.total,features=ros_gene)+scale_color_gradient(low="grey", "high"=single_col) + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short_double))+coord_flip() +theme(axis.text.x = element_text(angle = 90)) + ggtitle("Response to Oxidative Stress")
dev.off()
pdf("Graphs/Side/Dotplot_GLUT1.pdf")
DotPlot(s.k.total,features="SLC2A1")+scale_color_gradient(low="grey", "high"=single_col) + theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_UMAP_short_double))+coord_flip() +theme(axis.text.x = element_text(angle = 90)) + ggtitle("Glut1")
dev.off()
