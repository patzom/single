###data
library(Seurat)
#load data

#Keloid (sternum) patient sample - bulk cells after tissue dissociation and death cell removal
tabib.data1  <- Read10X("C:/Users/labor/Desktop/patrick/single/single/data/set_4/MD001_DK1/filtered_feature_bc_matrix")

#Keloid (earlobe) patient sample - bulk cells after tissue dissociation and death cell removal.
tabib.data2  <- Read10X("C:/Users/labor/Desktop/patrick/single/single/data/set_16/MD004_DK3_Direder_4/filtered_feature_bc_matrix")

#Keloid (earlobe) patient sample - bulk cells after tissue dissociation and death cell removal,
tabib.data3  <- Read10X("C:/Users/labor/Desktop/patrick/single/single/data/set_20/MD005_DK4L_Direder_5/filtered_feature_bc_matrix")

#Keloid (earlobe) patient sample - bulk cells after tissue dissociation and death cell removal,
tabib.data4  <- Read10X("C:/Users/labor/Desktop/patrick/single/single/data/MD005_DK4R_Direder_5/filtered_feature_bc_matrix")



#set up tabib data
tabib.keloid1 <- CreateSeuratObject(counts=tabib.data1 , assay="RNA")
tabib.keloid2 <- CreateSeuratObject(counts=tabib.data2 , assay="RNA")
tabib.keloid3 <- CreateSeuratObject(counts=tabib.data3 , assay="RNA")
tabib.keloid4 <- CreateSeuratObject(counts=tabib.data4 , assay="RNA")



tabib.keloid <- merge(tabib.keloid1, y = c(tabib.keloid2, tabib.keloid3, tabib.keloid4), add.cell.ids = c("k1", "k2", "k3", "k4"), project = "")
tabib.keloid
#subset samples


tabib.keloid[["percent.mt"]] <- PercentageFeatureSet(tabib.keloid, pattern = "^MT-")

head(tabib.keloid@meta.data, 5)

VlnPlot(tabib.keloid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)



plot1 <- FeatureScatter(tabib.keloid, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(tabib.keloid, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


tabib.keloid <- subset(tabib.keloid, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


### normalize data
tabib.keloid <- NormalizeData(tabib.keloid, normalization.method = "LogNormalize", scale.factor = 10000)


tabib.keloid <- FindVariableFeatures(tabib.keloid, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(tabib.keloid), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(tabib.keloid)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#scaling the data
all.genes <- rownames(tabib.keloid)
tabib.keloid <- ScaleData(tabib.keloid, features = all.genes)


tabib.keloid <- RunPCA(tabib.keloid, features = VariableFeatures(object = tabib.keloid))

#PCA
print(tabib.keloid[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(tabib.keloid, dims = 1:2, reduction = "pca")
DimPlot(tabib.keloid, reduction = "pca")

#heatmap
DimHeatmap(tabib.keloid, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(tabib.keloid, dims = 1:15, cells = 500, balanced = TRUE)


###
tabib.keloid <- JackStraw(tabib.keloid, num.replicate = 100)
tabib.keloid <- ScoreJackStraw(tabib.keloid, dims = 1:20)


JackStrawPlot(tabib.keloid, dims = 1:15)
ElbowPlot(tabib.keloid)


########## Cluster the cells based on elbow plot

tabib.keloid <- FindNeighbors(tabib.keloid, dims = 1:10)
tabib.keloid <- FindClusters(tabib.keloid, resolution = 0.5)
head(Idents(tabib.keloid), 5)




tabib.keloid <- RunUMAP(tabib.keloid, dims = 1:10)
DimPlot(tabib.keloid, reduction = "umap")

# find all markers of cluster 2
cluster2.markers <- FindMarkers(tabib.keloid, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)


# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(tabib.keloid, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)


# find markers for every cluster compared to all remaining cells, report only the positive
# ones
keloid.markers <- FindAllMarkers(tabib.keloid, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
keloid.markers%>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)


cluster0.markers <- FindMarkers(tabib.keloid, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)


keloid.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(tabib.keloid, features = top10$gene) + NoLegend()

#####
ew.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                    "NK", "DC", "Platelet")

names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()