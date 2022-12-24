#0==============================================================================

#1==============================================================================
#load the data
library(Seurat) # please update to Seurat V4
library(tidyverse)
#Sharir seuratobject
pool.data<- ReadMtx(
  mtx = "/Users/shibatashouta/Desktop/a large pool data/GSM3767568control/GSM3767568_control_matrix.mtx", features = "/Users/shibatashouta/Desktop/a large pool data/GSM3767568control/GSM3767568_control_genes.tsv",cells = "/Users/shibatashouta/Desktop/a large pool data/GSM3767568control/GSM3767568_control_barcodes.tsv"
)
pool.meta <- read.table("/Users/shibatashouta/Desktop/a large pool data/GSE131204/GSE131204_cell_info_8594x25.tsv", 
                        sep = "\t", header = TRUE)
control.meta.data = filter(pool.meta, condition == "control")
control.meta.data.df = control.meta.data[, c(-1, -2)]
rownames(control.meta.data.df) = control.meta.data[,2] 
#Quality check
pool <- CreateSeuratObject(counts = pool.data, project = "pool_control", min.cells = 3, min.features = 3)
pool <- AddMetaData(pool, control.meta.data.df)
pool[["percent.mt"]] <- PercentageFeatureSet(pool, pattern = "^mt-")
VlnPlot(pool, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#Filter 
pool <- subset(pool, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 5)
pool#3464
pool$population%>%table()

#Krivanek seuratobject
atlas.data = read.table( "/Users/shibatashouta/Desktop/dental cell type atlas/counts_SS2_mouse_incisor.txt",
                         header = TRUE, sep = " " ,skip = 0 )
atlas1 = CreateSeuratObject(atlas.data, project = "celltypeatlas", min.cells = 3, min.features = 3)
atlas1.meta = read.table("/Users/shibatashouta/Desktop/dental cell type atlas/annotation_mouse_incisor.txt", header = FALSE, sep = "\t")
atlas1.meta.ch = atlas1.meta[,-1]
cellname = 1:2889
atlas1 = AddMetaData(atlas1, cellname, "cellname")
cellname = atlas1.meta.ch
atlas1 = AddMetaData(atlas1, cellname, "cellname")
atlas1[["percent.mt"]] <- PercentageFeatureSet(atlas1, pattern = "^mt-")
VlnPlot(atlas1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#Filter
atlas1 = subset(atlas1, subset = nFeature_RNA > 200 & nFeature_RNA < 7500)
atlas1$cellname%>%table()
#simply integrated two data
dat1.combined <- merge(pool, y = atlas1, add.cell.ids = c("pool", "atlas"), project = "pool_atlas_combined")
dat1.combined$orig.ident%>%table()
saveRDS(dat1.combined, "/Users/shibatashouta/Desktop/dental cell type atlas/pool_atlas_combined")
dat1.combined = readRDS("/Users/shibatashouta/Desktop/dental cell type atlas/pool_atlas_combined")
dat1.combined$orig.ident%>%table()

#間違った統合方法


# split the dataset into a list of two seurat objects (stim and CTRL)
dat1.list <- SplitObject(dat1.combined, split.by = "orig.ident")
# normalize and identify variable features for each dataset independently
dat1.list <- lapply(X = dat1.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = dat1.list)
#perform integration
incisor.anchors <- FindIntegrationAnchors(object.list = dat1.list, anchor.features = features)
# this command creates an 'integrated' data assay
incisor.combined <- IntegrateData(anchorset = incisor.anchors)
DefaultAssay(incisor.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
incisor.combined <- ScaleData(incisor.combined, verbose = FALSE)
incisor.combined <- RunPCA(incisor.combined, npcs = 30, verbose = FALSE)
incisor.combined <- RunUMAP(incisor.combined, reduction = "pca", dims = 1:30)
incisor.combined <- FindNeighbors(incisor.combined, reduction = "pca", dims = 1:30)
incisor.combined <- FindClusters(incisor.combined, resolution = 0.5)
saveRDS(incisor.combined, "/Users/shibatashouta/Desktop/dental cell type atlas/pool_atlas_combined.after clustering")
#save point
sessionInfo()
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SeuratData)
library(ggplot2)
library(patchwork)
incisor.combined = readRDS("/Users/shibatashouta/Desktop/dental cell type atlas/pool_atlas_combined.after clustering")
head(dat1,5)
dat1 = incisor.combined
#louvain法
incisor.combined <- FindNeighbors(incisor.combined, dims = 1:30)
incisor.combined <- FindClusters(incisor.combined, resolution = 0.5) #0.4のほうがええかも
# Visualization
p1 <- DimPlot(incisor.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(incisor.combined, reduction = "umap", label = TRUE, repel = TRUE)
p2
p1 + p2
DimPlot(incisor.combined, reduction = "umap", split.by = "orig.ident",label = TRUE, repel = TRUE)
p1 <- DimPlot(incisor.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(incisor.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
#2==============================================================================
#cell labeling
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SeuratData)
library(ggplot2)
library(patchwork)
incisor.combined = readRDS("/Users/shibatashouta/Desktop/dental cell type atlas/pool_atlas_combined.after clustering")
incisor.combined #6489
#Sharir氏らのメタデータを用いて細胞周期にある細胞を取り除く regressing out cell cycle effect
incisor.combined$population%>%table()
incisor.combined$cellname%>%table()
dat1 = subset(incisor.combined, population == "ctr_M_G1", invert = TRUE)
dat1 = subset(dat1, population == "ctr_G2_M", invert = TRUE)
dat1 = subset(dat1, population == "ctr_S", invert = TRUE)

dat1 <- ScaleData(dat1, verbose = FALSE)
dat1 <- RunPCA(dat1, npcs = 30, verbose = FALSE)
dat1 <- RunUMAP(dat1, reduction = "pca", dims = 1:30)
dat1 <- FindNeighbors(dat1, reduction = "pca", dims = 1:30)
dat1 <- FindClusters(dat1, resolution = 0.5)
# Visualization
p1 <- DimPlot(dat1, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(dat1, reduction = "umap", label = TRUE, repel = TRUE)
p1
p2
p1 + p2
DimPlot(dat1, reduction = "umap",label = TRUE, repel = TRUE)
dat1$orig.ident%>%table()
saveRDS(dat1,"/Users/shibatashouta/Desktop/atlas and pool integration/no_cycling_cell_pool_atlas_integration" )
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SeuratData)
library(ggplot2)
library(patchwork)
#savepoint
dat1 = readRDS("/Users/shibatashouta/Desktop/atlas and pool integration/no_cycling_cell_pool_atlas_integration")
dat1$orig.ident%>%table()
DimPlot(dat1, reduction = "umap",label = TRUE, repel = TRUE)
dat1$seurat_clusters%>%table()
dat1$cellname%>%table()
dat1$population%>%table()
DimPlot(dat1, reduction = "umap", label = TRUE)
Idents(dat1)%>%table()
FeaturePlot(dat1, features = c("Vegfa","Ngfr"))



cluster0 = subset(dat1, idents = 0)
cluster0$population%>%table()
cluster0$cellname%>%table()
cluster0


cluster1$population%>%table()
cluster1$cellname%>%table()

cluster2 = subset(dat1,idents =2)
cluster2$population%>%table()
cluster2$cellname%>%table()

cluster3 = subset(dat1, idents = 3)
cluster3$population%>%table()
cluster3$cellname%>%table()

cluster4 = subset(dat1,idents = 4)
cluster4$population%>%table()
cluster4$cellname%>%table()

cluster5 = subset(dat1, idents = 5)
cluster5$population%>%table()
cluster5$cellname%>%table()

cluster6 = subset(dat1, idents = 6)
cluster6$population%>%table()
cluster6$cellname%>%table()

cluster7 = subset(dat1, idents = 7)
cluster7$population%>%table()
cluster7$cellname%>%table()
#dental follicle

cluster8 = subset(dat1, idents = 8)
cluster8$population%>%table()
cluster8$cellname%>%table()

cluster9 = subset(dat1, idents = 9)
cluster9$population%>%table()
cluster9$cellname%>%table()

cluster10 = subset(dat1, idents = 10)
cluster10$population%>%table()
cluster10$cellname%>%table()

cluster11 = subset(dat1, idents = 11)
cluster11$population%>%table()
cluster11$cellname%>%table()

cluster12 = subset(dat1, idents = 12)
cluster12$population%>%table()
cluster12$cellname%>%table()

cluster13 = subset(dat1, idents = 13)
cluster13$population%>%table()
cluster13$cellname%>%table()

cluster14 = subset(dat1, idents = 14)
cluster14$population%>%table()
cluster14$cellname%>%table()

cluster15 = subset(dat1, idents = 15)
cluster15$population%>%table()
cluster15$cellname%>%table()

cluster16 = subset(dat1, idents = 16)
cluster16$population%>%table()
cluster16$cellname%>%table()

cluster17 = subset(dat1, idents = 17)
cluster17$population%>%table()
cluster17$cellname%>%table()

cluster18 = subset(dat1, idents = 18)
cluster18$population%>%table()
cluster18$cellname%>%table()
dat1

====================================================================================-
  #単純な細胞集団を取り除く
  cluster9= subset(dat1, ident = "9")
Idents(cluster9) = "Endothelial"
celltype = Idents(cluster9)
cluster9 = AddMetaData(cluster9, celltype, "celltype")
cluster9$celltype%>%table()
Endothelial = saveRDS(cluster9, file ="/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/Endothelial")
Endothelial = readRDS("/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/Endothelial")

cluster10 = subset(dat1, ident = "10")
Idents(cluster10) = "Lymphocytes"
Lymphocytes = saveRDS(cluster9, file ="/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/Lymphocytes")
Lymphocytes = readRDS("/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/Lymphocytes")

cluster17= subset(dat1, ident = "17")
Idents(cluster17) = "Alveolar_osteo"
celltype = Idents(cluster17)
cluster17 = AddMetaData(cluster17, celltype, "celltype")
cluster17$celltype%>%table()
Alveolar_osteo = saveRDS(cluster17, file ="/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/Alveolar_osteo")
Alveolar_osteo = readRDS("/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/Alveolar_osteo")

cluster16= subset(dat1, ident = "16")
Idents(cluster16) = "Perivascular"
celltype = Idents(cluster16)
cluster16 = AddMetaData(cluster16, celltype, "celltype")
cluster16$celltype%>%table()
Perivascular = saveRDS(cluster16, file ="/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/Perivascular")
Perivascular = readRDS("/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/Perivascular")

cluster15= subset(dat1, ident = "15")
Idents(cluster15) = "Innate_leukocytes"
celltype = Idents(cluster15)
cluster15 = AddMetaData(cluster15, celltype, "celltype")
cluster15$celltype%>%table()
Innate_leukocytes = saveRDS(cluster15, file ="/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/Innate_leukocytes")
Innate_leukocytes = readRDS("/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/Innate_leukocytes")

cluster11= subset(dat1, ident = "11")
Idents(cluster11) = "Glia"
celltype = Idents(cluster11)
cluster11 = AddMetaData(cluster11, celltype, "celltype")
cluster11$celltype%>%table()
Glia = saveRDS(cluster11, file ="/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/Glia")
Glia = readRDS("/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/Glia")

cluster18= subset(dat1, ident = "18")
Idents(cluster18) = "Ameloblasts"
celltype = Idents(cluster18)
cluster18 = AddMetaData(cluster18, celltype, "celltype")
cluster18$celltype%>%table()
Ameloblasts = saveRDS(cluster18, file ="/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/Ameloblasts")
Ameloblasts = readRDS("/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/Ameloblasts")

cluster14は名無しの細胞集団だったので今回は判断が出来ないとし取り除いた。
==========================================================================================================================
  #Immune cell 系
  Immunocell = subset(dat1, ident = c("0","10"))
#louvain法
Immunocell <- FindNeighbors(Immunocell, dims = 1:30)
Immunocell <- FindClusters(Immunocell, resolution = 0.3) 
#Umap法
Immunocell <- RunUMAP(Immunocell, dims = 1:30)
DimPlot(Immunocell, reduction = "umap",label = TRUE, repel = TRUE)
IM_0 = subset(Immunocell, ident = "0")
IM_1 = subset(Immunocell, ident = "1")
IM_2 = subset(Immunocell, ident = "2")
IM_3 = subset(Immunocell, ident = "3")
IM_4 = subset(Immunocell, ident = "4")
IM_0$cellname%>%table()
IM_1$cellname%>%table()
IM_2$cellname%>%table()
new.cluster.ids <- c("Lyve1 Macrophage","Macrophage","Lymphocytes")
names(new.cluster.ids) <- levels(Immunocell)
Immunocell <- RenameIdents(Immunocell, new.cluster.ids)
UMAPPlot(Immunocell, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
Idents(Immunocell)%>%table()
VlnPlot(Immunocell, features = "Cd163", assay = "RNA")#Lyve1 macrophageはCd163陽性 https://www.nature.com/articles/s41598-022-08987-3
celltype = Idents(Immunocell)
Immunocell = AddMetaData(Immunocell, celltype, 'celltype')
Macrophages_Lyve1_Macrophages = saveRDS(Immunocell, file = "/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/Macrophages_Lyve1_Macrophages")
Macrophages_Lyve1_Macrophages = readRDS("/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/Macrophages_Lyve1_Macrophages")

DF = subset(dat1, ident ="6")
#louvain法
DF <- FindNeighbors(DF, dims = 1:30)
DF <- FindClusters(DF, resolution = 0.5) 
#Umap法
DF <- RunUMAP(DF, dims = 1:30)
DimPlot(DF, reduction = "umap",label = TRUE, repel = TRUE)
DF0 = subset(DF, ident = "0")
DF1 = subset(DF, ident = "1")
DF2 = subset(DF, ident = "2")
DF0$cellname%>%table()
DF1$cellname%>%table()
DF2$cellname%>%table()
new.cluster.ids <- c("Dental_follicle_2","Dental_follicle_2","Dental_follicle_1")
names(new.cluster.ids) <- levels(DF)
DF <- RenameIdents(DF, new.cluster.ids)
UMAPPlot(DF, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
Idents(DF)%>%table()
celltype = Idents(DF)
DF = AddMetaData(DF, celltype, 'celltype')
Dental_follicle = saveRDS(DF, file = "/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/Dental_follicle")
Dental_follicle = readRDS("/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/Dental_follicle")

#=================================================================================================================================================
上皮系
epithelium = subset(dat1, idents = c("1","2","7","12"))
#louvain法
epithelium <- FindNeighbors(epithelium, dims = 1:30)
epithelium <- FindClusters(epithelium, resolution = 0.5) 
#Umap法
epithelium <- RunUMAP(epithelium, dims = 1:30)
DimPlot(epithelium, reduction = "umap",label = TRUE, repel = TRUE)

epi0 = subset(epithelium, ident = "0")
epi1 = subset(epithelium, ident = "1")
epi2 = subset(epithelium, ident = "2")
epi3 = subset(epithelium, ident = "3")
epi4 = subset(epithelium, ident = "4")
epi5 = subset(epithelium, ident = "5")
epi6 = subset(epithelium, ident = "6")

epi0$population%>%table()
epi0$cellname%>%table()
epi1$population%>%table()
epi1$cellname%>%table()
epi2$population%>%table()
epi2$cellname%>%table()
epi3$population%>%table()
epi3$cellname%>%table()
epi4$population%>%table()
epi4$cellname%>%table()
epi5$population%>%table()
epi5$cellname%>%table()
epi6$population%>%table()
epi6$cellname%>%table()
epithelium

#epi0はpreAMB主
#louvain法
epi0<- FindNeighbors(epi0, dims = 1:30)
epi0 <- FindClusters(epi0, resolution = 0.3) 
#Umap法
epi0 <- RunUMAP(epi0, dims = 1:30)
DimPlot(epi0, reduction = "umap",label = TRUE, repel = TRUE)
epi0$population%>%table()
new.cluster.ids <- c("ctr_pre_AMB_0","ctr_pre_AMB_1")
names(new.cluster.ids) <- levels(epi0)
epi0<- RenameIdents(epi0, new.cluster.ids)
UMAPPlot(epi0, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
#epi1はOSRとSIが主
#louvain法
epi1 <- FindNeighbors(epi1, dims = 1:30)
epi1 <- FindClusters(epi1, resolution = 0.3) 
#Umap法
epi1 <- RunUMAP(epi1, dims = 1:30)
DimPlot(epi1, reduction = "umap",label = TRUE, repel = TRUE)
epi1_0 = subset(epi1, ident = "0")
epi1_1 = subset(epi1, ident = "1")
epi1_0$population%>%table()
epi1_0$cellname%>%table()
epi1_1$population%>%table()
epi1_1$cellname%>%table()
new.cluster.ids <- c("ctr_OSR","ctr_SI")
names(new.cluster.ids) <- levels(epi1)
epi1<- RenameIdents(epi1, new.cluster.ids)
UMAPPlot(epi1, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
#epi2はOEE_IEEとupperIEEが主
#louvain法
epi2 <- FindNeighbors(epi2, dims = 1:30)
epi2 <- FindClusters(epi2, resolution = 0.3) 
#Umap法
epi2 <- RunUMAP(epi2, dims = 1:30)
DimPlot(epi2, reduction = "umap",label = TRUE, repel = TRUE)
epi2_0 = subset(epi2, ident = "0")
epi2_1 = subset(epi2, ident = "1")
new.cluster.ids <- c("ctr_upper_IEE","ctr_OEE_IEE")
names(new.cluster.ids) <- levels(epi2)
epi2<- RenameIdents(epi2, new.cluster.ids)
UMAPPlot(epi2, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

#ISR/SIとOEE1が主
#louvain法
epi3 <- FindNeighbors(epi3, dims = 1:30)
epi3 <- FindClusters(epi3, resolution = 0.3) 
#Umap法
epi3 <- RunUMAP(epi3, dims = 1:30)
DimPlot(epi3, reduction = "umap",label = TRUE, repel = TRUE)
epi3_0 = subset(epi3, ident = "0")
epi3_1 = subset(epi3, ident = "1")
epi3_0$population%>%table()
epi3_0$cellname%>%table()
epi3_1$population%>%table()
epi3_1$cellname%>%table()
new.cluster.ids <- c("ctr_OEE_1","ctr_ISR_SI")
names(new.cluster.ids) <- levels(epi3)
epi3<- RenameIdents(epi3, new.cluster.ids)
UMAPPlot(epi3, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
#ctr_OEE2/VEE
#louvain法
epi4 <- FindNeighbors(epi4, dims = 1:30)
epi4 <- FindClusters(epi4, resolution = 0.3) 
#Umap法
epi4 <- RunUMAP(epi4, dims = 1:30)
DimPlot(epi4, reduction = "umap",label = TRUE, repel = TRUE)
epi4_0 = subset(epi4, ident = "0")
epi4_1 = subset(epi4, ident = "1")
epi4_2 = subset(epi4, ident = "2")
epi4_0$population%>%table()
epi4_0$cellname%>%table()
epi4_1$population%>%table()
epi4_1$cellname%>%table()
epi4_2$population%>%table()
epi4_2$cellname%>%table()
new.cluster.ids <- c("ctr_OEE_2","ctr_VEE","OEEnew")
names(new.cluster.ids) <- levels(epi4)
epi4<- RenameIdents(epi4, new.cluster.ids)
UMAPPlot(epi4, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
#ctr_DEEx
#louvain法
epi5 <- FindNeighbors(epi5, dims = 1:30)
epi5 <- FindClusters(epi5, resolution = 0.3) 
#Umap法
epi5 <- RunUMAP(epi5, dims = 1:30)
DimPlot(epi5, reduction = "umap",label = TRUE, repel = TRUE)
epi5_0 = subset(epi5, ident = "0")
epi5_1 = subset(epi5, ident = "1")
Idents(epi5) = "ctr_DEEx"
epi5 <- RunUMAP(epi5, dims = 1:30, label = TRUE)
DimPlot(epi5, reduction = "umap", label = TRUE)
#AMB
#louvain法
epi6 <- FindNeighbors(epi6, dims = 1:30)
epi6 <- FindClusters(epi6, resolution = 0.3) 
#Umap法
epi6 <- RunUMAP(epi6, dims = 1:30)
DimPlot(epi6, reduction = "umap",label = TRUE, repel = TRUE)
epi6_0 = subset(epi6, ident = "0")
epi6_1 = subset(epi6, ident = "1")
epi6_0$population%>%table()
epi6_1$population%>%table()
new.cluster.ids <- c("ctr_AMB_prox","ctr_AMB_dist")
names(new.cluster.ids) <- levels(epi6)
epi6<- RenameIdents(epi6, new.cluster.ids)
UMAPPlot(epi6, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
epi.subcluster = merge(epi0, y = c(epi1, epi2, epi3, epi4, epi5,epi6),
                       add.cell.ids = c("epi0","epi1","epi2","epi3","epi4","epi5","epi6"), project = "epi.subcluster")
epi.subcluster
Idents(epi.subcluster)%>%table()
celltype = c(
  Idents(epi0),
  Idents(epi1),        
  Idents(epi2),
  Idents(epi3),
  Idents(epi4),
  Idents(epi5),
  Idents(epi6)
)
epi.subcluster = AddMetaData(epi.subcluster, celltype, 'celltype')
epi.subcluster$celltype%>%table()
Idents(epi.subcluster)%>%table()
saveRDS(epi.subcluster, file = "/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/epi.subcluster")
#=============================================================================================================================================
歯髄系
pulp = subset(dat1, ident = c("3","4","5","8","13"))
#louvain法
pulp <- FindNeighbors(pulp, dims = 1:30)
pulp <- FindClusters(pulp, resolution = 0.3) 
#Umap法
pulp <- RunUMAP(pulp, dims = 1:30)
DimPlot(pulp, reduction = "umap",label = TRUE, repel = TRUE)
pulp0= subset(pulp, idents = "0")
pulp0$population%>%table()
pulp0$cellname%>%table()
pulp1= subset(pulp, idents = "1")
pulp1$population%>%table()
pulp1$cellname%>%table()
pulp2= subset(pulp, idents = "2")
pulp2$population%>%table()
pulp2$cellname%>%table()
pulp3= subset(pulp, idents = "3")
pulp3$population%>%table()
pulp3$cellname%>%table()
pulp4= subset(pulp, idents = "4")
pulp4$population%>%table()
pulp4$cellname%>%table()
pulp5= subset(pulp, idents = "5")
pulp5$population%>%table()
pulp5$cellname%>%table()

#louvain法
pulp0 <- FindNeighbors(pulp0, dims = 1:30)
pulp0 <- FindClusters(pulp0, resolution = 0.3) 
#Umap法
pulp0 <- RunUMAP(pulp0, dims = 1:30)
DimPlot(pulp0, reduction = "umap",label = TRUE, repel = TRUE)
pulp0_0 = subset(pulp0, ident = "0")
pulp0_1 = subset(pulp0, ident = "1")
pulp0_0$cellname%>%table()
pulp0_1$cellname%>%table() 
new.cluster.ids <- c("Apical_pulp_0","Apical_pulp_1")
names(new.cluster.ids) <- levels(pulp0)
pulp0<- RenameIdents(pulp0, new.cluster.ids)
UMAPPlot(pulp0, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

#louvain法
pulp2 <- FindNeighbors(pulp2, dims = 1:30)
pulp2 <- FindClusters(pulp2, resolution = 0.8) 
#Umap法
pulp2 <- RunUMAP(pulp2, dims = 1:30)
DimPlot(pulp2, reduction = "umap",label = TRUE, repel = TRUE)
pulp2_0 = subset(pulp2, ident = "0")
pulp2_1 = subset(pulp2, ident = "1")
pulp2_2 = subset(pulp2, ident = "2")
pulp2_3 = subset(pulp2, ident = "3")
pulp2_0$cellname%>%table()
pulp2_0$population%>%table()
pulp2_1$cellname%>%table() 
pulp2_1$population%>%table() 
pulp2_2$cellname%>%table() 
pulp2_2$population%>%table() 
pulp2_3$cellname%>%table() 
pulp2_3$population%>%table() 
new.cluster.ids <- c("Apical_pulp_2","Pre-odontoblasts_0","Apical_pulp_2","Pre-odontoblasts_1")

names(new.cluster.ids) <- levels(pulp2)
pulp2<- RenameIdents(pulp2, new.cluster.ids)
UMAPPlot(pulp2, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
pulp$cellname%>%table

Idents(pulp4) = "Distal_pulp"
Idents(pulp1) = "Maturing_pulp"

#louvain法
pulp4 <- FindNeighbors(pulp4, dims = 1:30)
pulp4 <- FindClusters(pulp4, resolution = 0.8) 
#Umap法
pulp4 <- RunUMAP(pulp4, dims = 1:30)
DimPlot(pulp4, reduction = "umap",label = TRUE, repel = TRUE)


pulp3と５は名無しが多かったり細胞腫に統一性が見られないためノイズとみなす

pulp.subcluster = merge(pulp0, y = c(pulp1, pulp2, pulp4),
                        add.cell.ids = c("pulp0","pulp1","pulp2","pulp4"), project = "pulp.subcluster")
Idents(pulp.subcluster)%>%table()

celltype = c(
  Idents(pulp0),
  Idents(pulp1),        
  Idents(pulp2),
  Idents(pulp4))

pulp.subcluster = AddMetaData(pulp.subcluster, celltype, 'celltype')
pulp.subcluster$celltype%>%table()
saveRDS(pulp.subcluster,"/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/pulp.subcluster")
#===========================================================================================================================================================

Macrophages_Lyve1_Macrophages = readRDS("/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/Macrophages_Lyve1_Macrophages")
Ameloblasts = readRDS("/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/Ameloblasts")
Glia = readRDS("/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/Glia")
Innate_Leukocytes = readRDS("/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/Innate_leukocytes")
Perivascular = readRDS("/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/Perivascular")
Alveolar_osteo = readRDS("/Users/shibatashouta/Desktop/atlas and pool integration/Alveolar_osteo")
Endothelial = readRDS("/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/Endothelial")
epi.subcluster =  readRDS("/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/epi.subcluster")
pulp.subcluster = readRDS("/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/pulp.subcluster")
Lymphocytes = readRDS("/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/Lymphocytes")



#3==============================================================================
#Nichenet analysis
library(devtools)
require(devtools)
install_github("saeyslab/nichenetr", force = TRUE)
library(igraph)
library(interp)

library(nichenetr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(cowplot)
library(ggpubr)

epi.subcluster =  readRDS("/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/epi.subcluster")
pulp.subcluster = readRDS("/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/pulp.subcluster")
Macrophages_Lyve1_Macrophages = readRDS("/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/Macrophages_Lyve1_Macrophages")
Ameloblasts = readRDS("/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/Ameloblasts")
Glia = readRDS("/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/Glia")
Innate_Leukocytes = readRDS("/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/Innate_leukocytes")
Perivascular = readRDS("/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/Perivascular")
Alveolar_osteo = readRDS("/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/Alveolar_osteo")
Endothelial = readRDS("/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/Endothelial")
Lymphocytes = readRDS("/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/Lymphocytes")
Dental_follicle = readRDS("/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/Dental_follicle")
dat1 <- merge(epi.subcluster, y = c(pulp.subcluster,Macrophages_Lyve1_Macrophages,Lymphocytes,Ameloblasts,Glia,Innate_Leukocytes,Perivascular,Alveolar_osteo,Endothelial,Dental_follicle), 
              add.cell.ids = c("epi.subcluster","pulp.subcluster","Macrophages_Lyve1_Macrophages","Lymphocytes","Ameloblasts","Glia","Innate_Leukocytes","Perivascular","Alveolar_osteo","Endothelial","Dental_follicle")
              , project = "comprehensive")
seuratObj = dat1
DotPlot(dat1, features = c("F2r","Dspp","Amelx","Mmp20","Klk4","Enam","Amtn","Odam"), assay = "RNA")
#remove cells whose idents = Ameloblasts.it is because these cells do not have any gene features, which codes enamel matrix protein and I regard these cells as noise cells.
dat1 = subset(dat1, ident = "Ameloblasts", invert = TRUE)
saveRDS(dat1,"/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/incisor.big")
dat1$celltype%>%table()

library(nichenetr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(cowplot)
library(ggpubr)
#次回はここから
dat1 = readRDS("/Users/shibatashouta/Desktop/comprehensive nichenet/celltyping/incisor.big")
seuratObj = dat1
dat1$celltype%>%table()
DotPlot(dat1, features = "Gas1", assay = "RNA")
dat1
DotPlot(dat1, features = c("Nts","Ptn","Twist2","Peg3","Prrx1","Rps19","Pcp4","Dlx6as1","Adamts20","Tmtc2","Dlx1as","Vit","Frzb","Meg3","Gli1","Cd24a"),cols = "RdYlBu",assay = "RNA")
DotPlot(dat1, features = c("Crabp1","E130309F12Rik","1500012F01Rik","Moxd1","Rps24","Rpl3","Rpl29","Klkb1","B4galt1","Mpped1","Mycn","Smpd3","Trpm1","Extl1","Grik1"),cols = "RdYlBu",assay = "RNA")

#loading geneset concerning dental formation from "http://www.informatics.jax.org/go/term/GO:0042476"
geneset_dental_formation = read.csv("/Users/shibatashouta/Desktop/geneset_dental_formation.csv" )
geneset_GO = geneset_dental_formation$Symbol
x <- c("Cldn10","Pthlh","Igfbpl1","Sparcl1","Krt17","Tacstd2","Enpp2","Hes1","Krt15"
       ,"Dcn","Igfbp5","Sfrp5","Igfbp2","Tgfb2","Sfrp2"
       ,"Mmp2","Alpl","Acta2","Calr","Fgf3","Fgf3r","Fgf9","Fgf20","Jag1","Notch1","Notch2","Dlx5",
       "Gdf5","Sox2","Fn1","Itgb1","Sox21","Hapln1","Barx1","Lhx6","Lhx7",
       "Smad1","Sox5","Sox6","Mef2c","Tcf7","Mef2a","Mef2b","Mef2d",
       "Smad1","Smad2","Smad3","Smad4","Smad5","Smad6","Smad7","Smad8",
       "Fam83h","Gli","Yank1","Fos","Jun","Junb","Smoc2","Sparc","Tubb3","Foxd1",
       "Syt6","Notum","Sall1","Fzd1","Sfrp1","Rspo1","Trabd2b","Wif1","Cebpb","Fos","Nme1",
       "Mgfe8","Sox2","Ccnd1","Ccnd2","S1pr1","Sox9","Satb1","Id1","Twist1","Cdkn1a",
       "Gnao1","Eno1","Efnb1","Calm1","Siah2","Atp6v0a1","Kdelr2","Gtpbp1","Polr2c","Sort1",
       "Calb1","Cst6","Foxo1","Bglap","Ibsp",
       "Anxa2","Hmgb2","Hmgb1","Spry2","Spry4","Trpm7","Twist1","Tac1","Foxi3","Foxo1")
#追加したかったもの Foxi3
orignal_geneset =c(paste0(geneset_GO),x)
orignal_geneset


#Read in NicheNet’s ligand-target prior model, ligand-receptor network and weighted integrated networks:
ligand_target_matrix = readRDS("/Users/shibatashouta/Desktop/ligand_target_matrix.rds")#heavy file
ligand_target_matrix[1:5,1:5] 
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
head(lr_network)
weighted_networks = readRDS("/Users/shibatashouta/Desktop/weighted_networks.rds")
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct( from ,to ), by = c("from","to"))
head(weighted_networks$lr_sig)
head(weighted_networks$gr)
#change gene name into mouse symbol
lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()
ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]
weighted_networks_lr = weighted_networks_lr %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()




#NicheNet analysis on Seurat object

#Define sender and receiver
Idents(dat1)%>%table() #Idents are available
dat1$celltype%>%unique()
#choose receiver 
receiver = "Pre-odontoblasts_1"
#assay_oi
expressed_genes_receiver = get_expressed_genes(receiver , dat1, pct = 0.10, assay_oi = "RNA" ) #assay = "RNA"
expressed_genes_receiver
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

class(background_expressed_genes)
Idents(dat1)%>%unique()
dat1$celltype%>%unique()
#all senders
sender_celltypes = c(
  "ctr_OEE_IEE","ctr_upper_IEE","ctr_pre_AMB_0","ctr_pre_AMB_1", "ctr_AMB_prox","ctr_AMB_dist","ctr_SI", "ctr_ISR_SI","ctr_OSR"
  ,"ctr_DEEx","ctr_OEE_1","ctr_OEE_2","ctr_VEE","OEE",
  "Apical_pulp_0" ,"Apical_pulp_1","Apical_pulp_2","Maturing_pulp","Distal_pulp","Pre-odontoblasts_0","Pre-odontoblasts_1","Dental_follicle1","Dental_follicle2",
  "Lyve1 Macrophage","Macrophage","Lymphocytes" ,"Glia","Innate_Leukocytes","Perivascular","Alveolar_osteo","Endothelial" 
)

sender_celltypes_seuratObj = subset(seuratObj, idents = c(sender_celltypes))
#assay = "RNA"
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, dat1, 0.10, assay = "RNA")
# lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

#2 Define a gene set of interest
Idents(dat1)%>%unique()
#----------------------------------------
#pattern1  geneset_dental_formation
geneset_oi = geneset_dental_formation$Symbol
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

#---------------------------------------
#patern2
DE_table_receiver = FindMarkers(dat1, ident.1 = receiver, ident.2 = c("Pre-odontoblasts_0") 
                                ,min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA")
DE_table_receiver1 = DE_table_receiver%>% rownames_to_column("gene")
geneset_oi = DE_table_receiver1 %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
geneset_oi
length(geneset_oi)
#-----------------------------------------
#pattern3
DE_table_receiver = FindMarkers(dat1, ident.1 = receiver, ident.2 = c("Pre-odontoblasts_0") 
                                ,min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA")
DE_table_receiver1 = DE_table_receiver%>% rownames_to_column("gene")
geneset_oi = DE_table_receiver1 %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
dental_gene = intersect(geneset_oi, orignal_geneset) #Extract only genes associated with tooth development
geneset_oi = dental_gene
length(geneset_oi)
#------------------------------------------


#3Define a set of potential ligands: these are ligands that are expressed by the “sender/niche” 
#cell population and bind a (putative) receptor expressed by the “receiver/target” population
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)  #sender達が発現していたシグナル分子？
write.csv(expressed_ligands) #どんなものが有るか確認 #174個あった
expressed_receptors = intersect(receptors,expressed_genes_receiver)#受け取る細胞発現していた受容体遺伝子？
potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>%unique() 
write.csv(potential_ligands)

ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes,
                                              ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson))) 

best_upstream_ligands = ligand_activities %>% top_n(40,pearson) %>%arrange(-pearson) %>%pull(test_ligand) %>% unique()
best_upstream_ligands
#assay = "RNA"
DotPlot(dat1, features = best_upstream_ligands %>% rev(), cols = "RdYlBu", assay = "RNA") + RotatedAxis()
DotPlot(sender_celltypes_seuratObj, features = best_upstream_ligands %>% rev(), cols = "RdYlBu",assay = "RNA") + RotatedAxis()#senderのみ
p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  # geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(40, pearson) %>% pull(pearson))), color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()


#5) infer receptors and top-predicted target genes of ligands that are top-ranked in the ligand activity analysis
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links, geneset = geneset_oi,
                                                                 ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names()
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names()
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
class(vis_ligand_target)

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
target_gene_list = unique(p_ligand_target_network$data$x) #重複なくtarget geneを抽出する

#receptors of top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")

#Receptors of top-ranked ligands, but after considering only bona fide ligand-receptor interactions documented in literature and publicly available databases
lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()

lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))

lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
lr_network_top_matrix_strict = lr_network_top_df_strict %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)

dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))

vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()
p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)")
ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized CAF-ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)")

#図一覧

p0 = p_ligand_pearson
P1 = p_ligand_target_network
ligand_pearson_matrix
P2 = DotPlot(sender_celltypes_seuratObj, features = best_upstream_ligands %>% rev(), cols = "RdYlBu", assay = "RNA") + RotatedAxis() #senderのみのパターン
P3 = p_ligand_receptor_network
P4 = p_ligand_receptor_network_strict

combined_plot = cowplot::plot_grid( P3, P4, P1,P2,rel_heights = c(4,5), nrow = 2, align = "hv")

#==================================================================================================
#書き出し
combined_plot
length(expressed_genes_sender)#5000-10000の間（500や20000ではない）であることが望ましい。
p_hist_lig_activity
match_gene = intersect(target_gene_list, orignal_geneset)
write.csv(match_gene, file = "/match_gene.csv")
ligand_pearson_matrix
write.csv(vis_ligand_target, file = "/vis_ligand_target.csv")

write.csv(ligand_pearson_matrix, file = "/ligand_pearson_matrix.csv")
#差次的遺伝子の数
length(geneset_oi)