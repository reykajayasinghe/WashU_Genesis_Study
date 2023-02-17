library(Seurat)
library(ggplot2)
library(reshape2)
library(viridis)
library(SingleR)
library(celldex)
library(dplyr)
library(tidyr)
library(ggrepel)

data_dir="/diskmnt/Datasets/mmy_scratch/Dipersio/genesis/data/"

P1<-readRDS(paste0(data_dir,"MGI1150_TWJJ-012-037-P37/MGI1150_TWJJ-012-037-P37_processed.rds"))
P2<-readRDS(paste0(data_dir,"MGI1150_TWJJ-012-049-B49/MGI1150_TWJJ-012-049-B49_processed.rds"))
P3<-readRDS(paste0(data_dir,"MGI1150_TWJJ-1104-A1104/MGI1150_TWJJ-1104-A1104_processed.rds"))
P4<-readRDS(paste0(data_dir,"MGI1150_TWJJ-012-015-P15/MGI1150_TWJJ-012-015-P15_processed.rds"))
P5<-readRDS(paste0(data_dir,"MGI1150_TWJJ-012-045-P45/MGI1150_TWJJ-012-045-P45_processed.rds"))
P6<-readRDS(paste0(data_dir,"MGI1150_TWJJ-012-055-P55/MGI1150_TWJJ-012-055-P55_processed.rds"))
P7<-readRDS(paste0(data_dir,"MGI1150_TWJJ-012-036-B36/MGI1150_TWJJ-012-036-B36_processed.rds"))
P8<-readRDS(paste0(data_dir,"MGI1150_TWJJ-1102-A1102/MGI1150_TWJJ-1102-A1102_processed.rds"))

P1$experiment<-"Placebo"
P2$experiment<-"BL-8040"
P3$experiment<-"Plerixafor"
P4$experiment<-"Placebo"
P5$experiment<-"Placebo"
P6$experiment<-"Placebo"
P7$experiment<-"BL-8040"
P8$experiment<-"Plerixafor"

P1$sample<-"012-037-P37"
P2$sample<-"012-049-B49"
P3$sample<-"1104-A1104"
P4$sample<-"012-015-P15"
P5$sample<-"012-045-P45"
P6$sample<-"012-055-P55"
P7$sample<-"012-036-B36"
P8$sample<-"1102-A1102"

###Check initial stats of data
pdf(paste("after_QC_in_sample_all.pdf", sep=""), width=15, height=9)
print(VlnPlot(object = P1, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), ncol = 3))
print(VlnPlot(object = P2, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), ncol = 3))
print(VlnPlot(object = P3, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), ncol = 3))
print(VlnPlot(object = P4, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), ncol = 3))
print(VlnPlot(object = P5, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), ncol = 3))
print(VlnPlot(object = P6, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), ncol = 3))
print(VlnPlot(object = P7, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), ncol = 3))
print(VlnPlot(object = P8, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), ncol = 3))
dev.off()

################################
###INTEGRATION 2 MERGING DATA###
################################

##### Merge datasets and process them
combined = merge(P1,y=c(P2,P3,P4,P5,P6,P7,P8),add.cell.ids = c("P01","P02","P03","P04","P05","P06","P07","P08"),project = "Merged")
combined <- SCTransform(combined, vars.to.regress = c("nCount_RNA","percent.mito"),return.only.var.genes = F)
combined <- RunPCA(combined, npcs = 30)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30,force.recalc=TRUE)
combined <- FindClusters(combined, resolution = 0.5)
combined <- RunUMAP(combined, dims = 1:30)

###Check DimPlot and PCA to make sure no batch effect
pdf(paste("batch_effect_check_merge.pdf", sep=""), width=12, height=12)
print(DimPlot(combined, reduction = "pca",group.by="experiment")&coord_equal())
print(DimPlot(combined, reduction = "pca",group.by="sample")&coord_equal())
print(DimPlot(combined, reduction = "umap",group.by="experiment")&coord_equal())
print(DimPlot(combined, reduction = "umap",group.by="sample")&coord_equal())
print(DimPlot(combined, reduction = "umap",group.by="cell_type",label=TRUE)&coord_equal())
dev.off()

saveRDS(combined,paste("all_integrated_merge.rds",sep=""))

# ##################################
# ###INTEGRATION 2 ANCHORING DATA### - NOT NEEDED NO BATCH EFFECT NOTED WITH MERGEING
# ##################################

# ifnb.list=c(P2,P3,P4,P5,P6,P7,P8,P9)

# # normalize and identify variable features for each dataset independently
# ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
# 	DefaultAssay(x)<-"RNA"
#     x <- NormalizeData(x)
#     x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
# })

# all_genes <- lapply(ifnb.list, row.names) %>% Reduce(intersect, .) 

# # select features that are repeatedly variable across datasets for integration run PCA on each
# # dataset using these features
# features <- SelectIntegrationFeatures(object.list = ifnb.list)
# ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
#     x <- ScaleData(x, features = features, verbose = FALSE)
#     x <- RunPCA(x, features = features, verbose = FALSE)
# })

# # normalize and identify variable features for each dataset independently
# ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
#     x <- NormalizeData(x)
#     x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
# })

# immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features, reduction = "rpca")
# immune.combined <- IntegrateData(anchorset = immune.anchors,features.to.integrate=all_genes)
# DefaultAssay(immune.combined) <- "integrated"
# # Run the standard workflow for visualization and clustering
# immune.combined <- ScaleData(immune.combined, verbose = FALSE)
# immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
# immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
# immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
# immune.combined <- FindClusters(immune.combined, resolution = 0.5)


# ###Check DimPlot and PCA to make sure no batch effect
# pdf(paste("batch_effect_check_anchor.pdf", sep=""), width=12, height=12)
# print(DimPlot(combined, reduction = "pca",group.by="experiment")&coord_equal())
# print(DimPlot(combined, reduction = "pca",group.by="sample")&coord_equal())
# print(DimPlot(combined, reduction = "umap",group.by="experiment")&coord_equal())
# print(DimPlot(combined, reduction = "umap",group.by="sample")&coord_equal())
# dev.off()

# saveRDS(immune.combined,paste("all_integrated_anchor.rds",sep=""))

#run SingleR

#Convert seurat object to single-cell format for singleR 
#https://satijalab.org/seurat/archive/v3.1/conversion_vignette.html
Seurat_Object_Diet <- DietSeurat(combined, graphs = "umap") #https://github.com/satijalab/seurat/issues/4633
obj.sce <- as.SingleCellExperiment(Seurat_Object_Diet)

#Reference SingleR
#https://bioconductor.org/packages/3.13/data/experiment/vignettes/celldex/inst/doc/userguide.html

ref <- MonacoImmuneData()
#singleR usage
#http://bioconductor.org/books/devel/SingleRBook/using-the-classic-mode.html#annotating-the-test-dataset
pred <- SingleR(test = obj.sce, ref = ref, 
    labels = ref$label.fine, assay.type.test=1)
#colnames(pred)
ref_pred<-as.data.frame(pred)

#Extract columns of interest to add to seurat object
ref_predictions<- ref_pred %>%
  select(pruned.labels)
colnames(ref_predictions) <- "MonacoImmuneData"

##Add predictions to seurat object
#Add metadata to seurat object
combined <- AddMetaData(
  object = combined,
  metadata = ref_predictions)


ref <- NovershternHematopoieticData()
pred <- SingleR(test = obj.sce, ref = ref, 
    labels = ref$label.fine, assay.type.test=1)
ref_pred<-as.data.frame(pred)
ref_predictions<- ref_pred %>%
  select(pruned.labels)
colnames(ref_predictions) <- "NovershternHematopoieticData"
combined <- AddMetaData(
  object = combined,
  metadata = ref_predictions)

ref <- DatabaseImmuneCellExpressionData()
pred <- SingleR(test = obj.sce, ref = ref, 
    labels = ref$label.fine, assay.type.test=1)
ref_pred<-as.data.frame(pred)
ref_predictions<- ref_pred %>%
  select(pruned.labels)
colnames(ref_predictions) <- "DatabaseImmuneCellExpressionData"
combined <- AddMetaData(
  object = combined,
  metadata = ref_predictions)

ref <- HumanPrimaryCellAtlasData()
pred <- SingleR(test = obj.sce, ref = ref, 
    labels = ref$label.fine, assay.type.test=1)
ref_pred<-as.data.frame(pred)
ref_predictions<- ref_pred %>%
  select(pruned.labels)
colnames(ref_predictions) <- "HumanPrimaryCellAtlasData"
combined <- AddMetaData(
  object = combined,
  metadata = ref_predictions)

saveRDS(combined,paste("all_integrated_merge_celltypes.rds",sep=""))

##DEGs derived from publication
#https://doi.org/10.1038/s41467-019-10291-0
#A comprehensive single cell transcriptional landscape of human hematopoietic progenitors
GOI<-c("CD34","CD164","BST2","CD37","KIT","ICAM3","CD44","EVI2B",
	"CD63","ITGA2B",
	"CSF2RB",
	"PROM1","CD34","ATP1B3","GYPC","CD47","IL2RG",
	"CD79A","MME","IL7R",
	"CSF3R","CD99","ITGB2","CD53","CD48","LAIR1","FLT3","CD74","SELL")

Meg: "ITGA2B", "PF4", "VWF"
E: "CA1", "HBB", "KLF1", "TFR2"
DC: "CCR2", "IRF8", "MPEG1"
G: "ELANE", "MPO", "LYZ", "CSF1R", "CTSG", "PRTN3", "AZU1" 
Ly1: "RGS1", "NPTX2", "DDIT4", "ID2"
Ly2: "DNTT", "RAG1", "RAG2" 
HSC: "CRHBP", "HLF", "DUSP1", "PCDH9"
And for the Lin−CD34/CD164 data set: 
E: "KLF1", "CA1"
Meg: "ITGA2B", "PLEK"
BEM: "CLC", "CPA3", "HDC"
Ly: "DNTT", "CD79A", "VPREB1" DC: "IRF8", "SPIB", "IGKC"
M: "LYZ", "MS4A6A", "ANXA2" 
N: "ELANE"
HSC: "HLF", "ADGRG6", "CRHBP", "PCDH9"

GOI_supp<-c("ITGA2B", "PF4", "VWF","PLEK",
"CA1", "HBB", "KLF1", "TFR2",
"CCR2", "IRF8", "MPEG1", "SPIB", "IGKC",
"ELANE", "MPO", "LYZ", "CSF1R", "CTSG", "PRTN3", "AZU1" ,"MS4A6A", "ANXA2",
"RGS1", "NPTX2", "DDIT4", "ID2", "DNTT", "RAG1", "RAG2" ,"CD79A", "VPREB1",
"CRHBP", "HLF", "DUSP1", "PCDH9","ADGRG6","PCDH9",
"CLC", "CPA3", "HDC")
#And for the Lin−CD34/CD164 data set: 
#E: "KLF1", "CA1"
#Meg: "ITGA2B", 

#Plerixafor mobilizies unique population of CD34+ HSC Subsets
#E2-2 = TCF4
#CD303=CLEC4C
#ILT7=LILRA4
GOI_pDC_genes<-c("TCF4","SPIB","IRF7","CLEC4C","ILT7")

pdf(paste("batch_effect_check_merge.pdf", sep=""), width=24, height=12)
print(DimPlot(combined, reduction = "pca",group.by="experiment")&coord_equal())
print(DimPlot(combined, reduction = "pca",group.by="sample")&coord_equal())
print(DimPlot(combined, reduction = "umap",group.by="experiment")&coord_equal())
print(DimPlot(combined, reduction = "umap",group.by="HumanPrimaryCellAtlasData")&coord_equal())
print(DimPlot(combined, reduction = "umap",group.by="MonacoImmuneData")&coord_equal())
print(DimPlot(combined, reduction = "umap",group.by="seurat_clusters",label=TRUE)&coord_equal())
print(DotPlot(combined, features=unique(GOI_pDC_genes),group.by="seurat_clusters")&coord_flip())
print(DotPlot(combined, features=unique(GOI),group.by="seurat_clusters")&coord_flip())
print(DotPlot(combined, features=unique(GOI_supp),group.by="seurat_clusters")&coord_flip())
dev.off()




library(Seurat)
library(dplyr)
library(tibble)
library(stringr)
library(reshape2)
library(ggplot2)

sobj<-readRDS("~rjayasin/mmy_scratch/Dipersio/genesis/analysis/all_integrated_merge_celltypes.rds")

Idents(sobj)<-sobj@meta.data$celltype_experiment
all.markers <- FindAllMarkers(object = sobj,logfc.threshold=0.4,only.pos=TRUE)
write.csv(all.markers,file="genesis_cluster_DEGs_experiment_celltype.csv",quote=FALSE)

####################################
#Analysis of HSC-I, HSC-II, HSC-III#
####################################

Idents(sobj)<-sobj@meta.data$experiment
all.markers <- FindAllMarkers(object = sobj,logfc.threshold=0.4,only.pos=TRUE)
write.csv(all.markers,file="genesis_DEGs_experimentonly.csv",quote=FALSE)


table(sobj@meta.data$sample,sobj@meta.data$cell_type)

Idents(sobj)<-sobj@meta.data$cell_type
HSC_only<-subset(sobj,idents=c("HSC-I","HSC-II","HSC-III"))

HSC_only<- SCTransform(HSC_only, vars.to.regress = c("percent.mito","nCount_RNA"),return.only.var.genes = F)
  #These are now standard steps in the 
  #Seurat workflow for visualization and clustering
HSC_only <- RunPCA(HSC_only,verbose = FALSE)
HSC_only <- RunUMAP(HSC_only, dims = 1:30, reduction="pca")
HSC_only <- FindNeighbors(HSC_only, dims = 1:30, reduction="pca",force.recalc = TRUE)
HSC_only <- FindClusters(HSC_only, resolution = 0.5)
HSC_only <- RunUMAP(HSC_only, dims = 1:30)

saveRDS(HSC_only,"HSC.rds")

# Downsample the number of cells per identity class
HSC_only@meta.data$sample_celltype<-paste0(HSC_only@meta.data$sample,"_",HSC_only@meta.data$cell_type)
Idents(HSC_only)<-HSC_only@meta.data$sample_celltype
HSC_subset<-subset(x = HSC_only, downsample = 100)
HSC_subset<- SCTransform(HSC_subset, vars.to.regress = c("percent.mito","nCount_RNA"),return.only.var.genes = F)
  #These are now standard steps in the 
  #Seurat workflow for visualization and clustering
HSC_subset <- RunPCA(HSC_subset,verbose = FALSE)
HSC_subset <- RunUMAP(HSC_subset, dims = 1:30, reduction="pca")
HSC_subset <- FindNeighbors(HSC_subset, dims = 1:30, reduction="pca",force.recalc = TRUE)
HSC_subset <- FindClusters(HSC_subset, resolution = 0.5)
HSC_subset <- RunUMAP(HSC_subset, dims = 1:30)

pdf(paste0("hsc.pdf"),height=8,width=10)
DimPlot(HSC_only,group.by=c("cell_type","sample","seurat_clusters"),label=TRUE)+coord_equal()
FeaturePlot(HSC_only,features=c("CD74"),split.by="experiment",order=TRUE)
VlnPlot(HSC_only,features=c("CD74"),group.by="cell_type",split.by="experiment")
VlnPlot(HSC_only,features=c("CD74"),group.by="cell_type",split.by="seurat_clusters")
DimPlot(HSC_subset,group.by=c("cell_type","sample","seurat_clusters"),label=TRUE)+coord_equal()
FeaturePlot(HSC_subset,features=c("CD74"),split.by="experiment",order=TRUE)
VlnPlot(HSC_subset,features=c("CD74"),group.by="cell_type",split.by="seurat_clusters")
VlnPlot(HSC_subset,features=c("CD74"),group.by="cell_type",split.by="experiment")
DotPlot()
dev.off()

Idents(HSC_only)<-HSC_only@meta.data$cell_type
all.markers <- FindAllMarkers(object = HSC_only,logfc.threshold=0.4,only.pos=TRUE)
write.csv(all.markers,file="genesis_DEGs_HSConly.csv",quote=FALSE


#downsample all samples to smallest sample to look at distribution
#012-015-P15 012-036-B36 012-037-P37 012-045-P45 012-049-B49 012-055-P55
#       4578        4687        7727        6200       20062       15494
# 1102-A1102  1104-A1104
#       8834        2086
Idents(sobj)<-sobj@meta.data$sample
sample_subset<-subset(x = sobj, downsample = 2086)


#######################################
##CALCULATE PROPORTIONS OF ALL SAMPLES#
#USING DOWNSAMPLED OBJECT##############
#######################################

#Make proportion table - extract proprotions and counts for plotting
data<-table(sample_subset@meta.data$sample,sample_subset@meta.data$cell_type)
#The value of each cell divided by the sum of the row cells - row is experiment in this case
data_prop<-prop.table(data,1)
data_adj<-data_prop*100
plot_data<-melt(data_adj)
plot_data$key<-paste0(plot_data$Var1,"_",plot_data$Var2)
colnames(plot_data)<-c("Sample","Seurat_Cluster","Proportion","key")
og_data<-melt(data)
og_data$key<-paste0(og_data$Var1,"_",og_data$Var2)
og_data$Var1<-NULL
og_data$Var2<-NULL
colnames(og_data)<-c("count","key")
merged_data<-merge(plot_data, og_data, by.x="key", by.y="key")

prop_plot<-ggplot(merged_data,aes(x=factor(Sample), y=factor(Seurat_Cluster), size=Proportion, fill=Sample,label=count)) +
    geom_point(alpha=0.8, shape=21, color="black") + geom_text(size=4,color="black")+
    scale_size(range = c(.1, 20), name="Proportion of Cells") +
    scale_fill_viridis(discrete=TRUE, guide=FALSE, option="D") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    theme_classic()+
    ylab("Predicted CellType") +
    xlab("Sample")


sobj@meta.data$sample<-factor(sobj@meta.data$sample,levels=c("012−015−P15","012−037−P37","012−045−P45","012−055−P55","012−036−B36","012−049−B49","1102−A1102","1104−A1104"))

a<-DimPlot(sample_subset,group.by=c("cell_type"),split.by=c("sample"),ncol=4)+coord_equal()

pdf(paste0("subset.pdf"),height=8,width=8)
print(a)
print(prop_plot)
dev.off()


############################
###BATCH 2 DATA INTEGRATION#
############################

library(Seurat)
library(ggplot2)
library(reshape2)
library(viridis)
library(SingleR)
library(celldex)
library(dplyr)
library(tidyr)
library(ggrepel)


sobj<-readRDS("all_integrated_merge_celltypes.rds")

P1<-readRDS(paste0(data_dir,"MGI2482_TWJG-1100-A1100/MGI2482_TWJG-1100-A1100_processed.rds"))
P2<-readRDS(paste0(data_dir,"MGI2482_TWJG-DIV29-A029/MGI2482_TWJG-DIV29-A029_processed.rds"))
P3<-readRDS(paste0(data_dir,"MGI2482_TWJG-012-052-B012-052/MGI2482_TWJG-012-052-B012-052_processed.rds"))
P4<-readRDS(paste0(data_dir,"MGI2482_TWJG-BL010-BL010/MGI2482_TWJG-BL010-BL010_processed.rds"))
P5<-readRDS(paste0(data_dir,"MGI2482_TWJG-S-005-G005/MGI2482_TWJG-S-005-G005_processed.rds"))
P6<-readRDS(paste0(data_dir,"MGI2482_TWJG-012-059-B012-059/MGI2482_TWJG-012-059-B012-059_processed.rds"))
P7<-readRDS(paste0(data_dir,"MGI2482_TWJG-BL011-BL011/MGI2482_TWJG-BL011-BL011_processed.rds"))
P8<-readRDS(paste0(data_dir,"MGI2482_TWJG-S-009-G009/MGI2482_TWJG-S-009-G009_processed.rds"))
P9<-readRDS(paste0(data_dir,"MGI2482_TWJG-D026-A026/MGI2482_TWJG-D026-A026_processed.rds"))
P10<-readRDS(paste0(data_dir,"MGI2482_TWJG-1096-A1096/MGI2482_TWJG-1096-A1096_processed.rds"))

P1$experiment<-"Plerixafor"
P2$experiment<-"Healthy_Plerixafor"
P3$experiment<-"BL-8040"
P4$experiment<-"Healthy_BL-8040"
P5$experiment<-"Healthy_Placebo"
P6$experiment<-"BL-8040"
P7$experiment<-"Healthy_BL-8040"
P8$experiment<-"Healthy_Placebo"
P9$experiment<-"Healthy_Plerixafor"
P10$experiment<-"Plerixafor"

P1$sample<-"1100-A1100"
P2$sample<-"DIV29-A029"
P3$sample<-"012-052-B012-052"
P4$sample<-"BL010-BL010"
P5$sample<-"S-005-G005"
P6$sample<-"012-059-B012-059"
P7$sample<-"BL011-BL011"
P8$sample<-"S-009-G009"
P9$sample<-"D026-A026"
P10$sample<-"1096-A1096"

###Check initial stats of data
pdf(paste("after_QC_in_sample_all.pdf", sep=""), width=15, height=9)
print(VlnPlot(object = P1, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), ncol = 3))
print(VlnPlot(object = P2, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), ncol = 3))
print(VlnPlot(object = P3, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), ncol = 3))
print(VlnPlot(object = P4, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), ncol = 3))
print(VlnPlot(object = P5, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), ncol = 3))
print(VlnPlot(object = P6, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), ncol = 3))
print(VlnPlot(object = P7, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), ncol = 3))
print(VlnPlot(object = P8, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), ncol = 3))
print(VlnPlot(object = P9, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), ncol = 3))
print(VlnPlot(object = P10, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), ncol = 3))
dev.off()

################################
###INTEGRATION 2 MERGING DATA###
###FROM BOTH BATCHES############
################################

##### Merge datasets and process them
combined = merge(sobj,y=c(P1,P2,P3,P4,P5,P6,P7,P8,P9,P10),cell_ids=c("B1","B2_P09","B2_P10","B2_P11","B2_P12","B2_P13","B2_P14","B2_P15","B2_P16","B2_P17","B2_P18"),project = "Merged")
combined <- SCTransform(combined, vars.to.regress = c("nCount_RNA","percent.mito"),return.only.var.genes = F)
combined <- RunPCA(combined, npcs = 30)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30,force.recalc=TRUE)
combined <- FindClusters(combined, resolution = 0.5)
combined <- RunUMAP(combined, dims = 1:30)

###Check DimPlot and PCA to make sure no batch effect
pdf(paste("batch_effect_check_merge.pdf", sep=""), width=12, height=12)
print(DimPlot(combined, reduction = "pca",group.by="experiment")&coord_equal())
print(DimPlot(combined, reduction = "pca",group.by="sample")&coord_equal())
print(DimPlot(combined, reduction = "umap",group.by="experiment")&coord_equal())
print(DimPlot(combined, reduction = "umap",group.by="sample")&coord_equal())
print(DimPlot(combined, reduction = "umap",group.by="cell_type",raster=FALSE,label=TRUE)&coord_equal())
print(DimPlot(combined, reduction = "umap",group.by="seurat_clusters",label=TRUE)&coord_equal())
print(prop_plot)
dev.off()

saveRDS(combined,paste("all_integrated_merge_012022.rds",sep=""))


write.table(a,file=paste0("cluster_celltype.txt"), sep="\t", quote = FALSE, row.names = TRUE)

library(Seurat)
library(ggplot2)
library(reshape2)
library(viridis)
library(SingleR)
library(celldex)
library(dplyr)
library(tidyr)
library(ggrepel)
combined<-readRDS("all_integrated_merge_012022.rds")

###############################################
#Goal - linking facs populations to scrna data#
###############################################

#CD45 = PTPRC
#CD123 = IL3RA
#CD38 = CD38
#CD90 = THY1
#CD49f = ITGA6
#CD10 =  MME
#CD303 =  CLEC4C

goi<-c("IL3RA","CD38","THY1","ITGA6","MME","CLEC4C")
genes<-c("TYMS","TUBA1B","HIST1H4C","IRF8","IRF7","LGALS1","LYZ","MPO","CRIP1","HIST1H1D","HIST1H1C","HIST1H1E","AVP","FTH1","HLA-E","KLF2","LMNA","CLEC3B","EGR1","JUN","ZFP36","HSPA1A","HSPA6","HSPA1B","ACY3","LSP1","SAMHD1","DNTT","LTB","JCHAIN","TPSAB1","PRG2","HDC","HBD","APOC1","MYC","HBB","BLVRB","CA1","C1QTNF4","SMIM24","CD52","FCER1A","PLEK","AL157895.1","MT-ATP8","HIST1H2AI")

pdf(paste0("hsc1.pdf"),height=12,width=12)
a<-DotPlot(combined,assay="SCT",features=goi,group.by="seurat_clusters",cols=c("paleturquoise1", "violetred4"))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
b<-DotPlot(combined,assay="SCT",features=genes,group.by="seurat_clusters",cols=c("paleturquoise1", "violetred4"))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
#b<-DotPlot(combined,assay="SCT", features=goi,group.by="cell_type",cols=c("paleturquoise1", "violetred4"))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
c<-FeaturePlot(combined, features=goi,order=TRUE,cols=c("paleturquoise1", "violetred4"))&coord_equal()
print(a)
#print(b)
print(c)
dev.off()

##########################
#SINGLE R ANNOTATION######
##########################

Seurat_Object_Diet <- DietSeurat(combined, graphs = "umap") #https://github.com/satijalab/seurat/issues/4633
obj.sce <- as.SingleCellExperiment(Seurat_Object_Diet)

#Reference SingleR
#https://bioconductor.org/packages/3.13/data/experiment/vignettes/celldex/inst/doc/userguide.html

ref <- MonacoImmuneData()
#singleR usage
#http://bioconductor.org/books/devel/SingleRBook/using-the-classic-mode.html#annotating-the-test-dataset
pred <- SingleR(test = obj.sce, ref = ref, 
    labels = ref$label.fine, assay.type.test=1)
#colnames(pred)
ref_pred<-as.data.frame(pred)

#Extract columns of interest to add to seurat object
ref_predictions<- ref_pred %>%
  select(pruned.labels)
colnames(ref_predictions) <- "MonacoImmuneData"

##Add predictions to seurat object
#Add metadata to seurat object
combined <- AddMetaData(
  object = combined,
  metadata = ref_predictions)

ref <- NovershternHematopoieticData()
pred <- SingleR(test = obj.sce, ref = ref, 
    labels = ref$label.fine, assay.type.test=1)
ref_pred<-as.data.frame(pred)
ref_predictions<- ref_pred %>%
  select(pruned.labels)
colnames(ref_predictions) <- "NovershternHematopoieticData"
combined <- AddMetaData(
  object = combined,
  metadata = ref_predictions)

ref <- DatabaseImmuneCellExpressionData()
pred <- SingleR(test = obj.sce, ref = ref, 
    labels = ref$label.fine, assay.type.test=1)
ref_pred<-as.data.frame(pred)
ref_predictions<- ref_pred %>%
  select(pruned.labels)
colnames(ref_predictions) <- "DatabaseImmuneCellExpressionData"
combined <- AddMetaData(
  object = combined,
  metadata = ref_predictions)

ref <- HumanPrimaryCellAtlasData()
pred <- SingleR(test = obj.sce, ref = ref, 
    labels = ref$label.fine, assay.type.test=1)
ref_pred<-as.data.frame(pred)
ref_predictions<- ref_pred %>%
  select(pruned.labels)
colnames(ref_predictions) <- "HumanPrimaryCellAtlasData"
combined <- AddMetaData(
  object = combined,
  metadata = ref_predictions)

saveRDS(combined,paste("all_integrated_merge_012022_v1.rds",sep=""))

all.markers <- FindAllMarkers(object = combined,only.pos=TRUE,min.pct=0.3)
write.table(all.markers,file="markers_012022.txt",quote=FALSE,sep="\t")

pdf(paste("batch_effect_check_merge.pdf", sep=""), width=24, height=12)
print(DimPlot(combined, reduction = "pca",group.by="experiment")&coord_equal())
print(DimPlot(combined, reduction = "pca",group.by="sample")&coord_equal())
print(DimPlot(combined, reduction = "umap",group.by="experiment")&coord_equal())
print(DimPlot(combined, reduction = "umap",group.by="HumanPrimaryCellAtlasData")&coord_equal())
print(DimPlot(combined, reduction = "umap",group.by="MonacoImmuneData")&coord_equal())
print(DimPlot(combined, reduction = "umap",group.by="NovershternHematopoieticData")&coord_equal())
print(DimPlot(combined, reduction = "umap",group.by="DatabaseImmuneCellExpressionData")&coord_equal())
print(DimPlot(combined, reduction = "umap",group.by="seurat_clusters",label=TRUE)&coord_equal())
b<-DotPlot(combined,assay="SCT",features=hsc_genes,group.by="seurat_clusters",cols=c("paleturquoise1", "violetred4"))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
print(b)
c<-DotPlot(combined,assay="SCT",features=genes,group.by="seurat_clusters",cols=c("paleturquoise1", "violetred4"))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
print(c)
d<-DotPlot(combined,assay="SCT",features=GOI,group.by="seurat_clusters",cols=c("paleturquoise1", "violetred4"))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
print(d)
dev.off()


########################
#ADD HSC SIGNATURE######
########################

HSC_score<-list(c("CD164","BST2","CD37","KIT","ICAM3","CD44","EVI2B"))
combined <- AddModuleScore(object = combined,assay="SCT",features = HSC_score,name = 'HSC_score')


##DEGs derived from publication
#https://doi.org/10.1038/s41467-019-10291-0
#A comprehensive single cell transcriptional landscape of human hematopoietic progenitors
GOI<-c("CD34","CD164","BST2","CD37","KIT","ICAM3","CD44","EVI2B",
  "CD63","ITGA2B",
  "CSF2RB",
  "PROM1","CD34","ATP1B3","GYPC","CD47","IL2RG",
  "CD79A","MME","IL7R",
  "CSF3R","CD99","ITGB2","CD53","CD48","LAIR1","FLT3","CD74","SELL")

Meg: "ITGA2B", "PF4", "VWF"
E: "CA1", "HBB", "KLF1", "TFR2"
DC: "CCR2", "IRF8", "MPEG1"
G: "ELANE", "MPO", "LYZ", "CSF1R", "CTSG", "PRTN3", "AZU1" 
Ly1: "RGS1", "NPTX2", "DDIT4", "ID2"
Ly2: "DNTT", "RAG1", "RAG2" 
HSC: "CRHBP", "HLF", "DUSP1", "PCDH9"
And for the Lin−CD34/CD164 data set: 
E: "KLF1", "CA1"
Meg: "ITGA2B", "PLEK"
BEM: "CLC", "CPA3", "HDC"
Ly: "DNTT", "CD79A", "VPREB1" DC: "IRF8", "SPIB", "IGKC"
M: "LYZ", "MS4A6A", "ANXA2" 
N: "ELANE"
HSC: "HLF", "ADGRG6", "CRHBP", "PCDH9"

GOI_supp<-c("ITGA2B", "PF4", "VWF","PLEK",
"CA1", "HBB", "KLF1", "TFR2",
"CCR2", "IRF8", "MPEG1", "SPIB", "IGKC",
"ELANE", "MPO", "LYZ", "CSF1R", "CTSG", "PRTN3", "AZU1" ,"MS4A6A", "ANXA2",
"RGS1", "NPTX2", "DDIT4", "ID2", "DNTT", "RAG1", "RAG2" ,"CD79A", "VPREB1",
"CRHBP", "HLF", "DUSP1","ADGRG6","PCDH9",
"CLC", "CPA3", "HDC")

GOI_hayetal_others<-c("TCL1A","IRF4","DNTT","VPREB3","CYGB","ADA","LTB","HMBG1","MS4A1","LEF1","GZMA","SH2D1A","SCT","TGFB1","SAMHD1","LGALS1",
  "ACY3","NEGR1","AVP","CRHBP","HLF", "DUSP1","ADGRG6","PCDH9","TMEM154","MYOZ3","C1QTNF4","SPINK2","AZU1","ELANE","CMTM5","PF4","FCER1A","TIMP3","SLC40A1","PKIG","APOC1","BLVRB","CLC","HDC",
  "CSF3R","SMIM24","TOP2A","MKI67")

pdf("score.pdf",height=8,width=8,useDingbats=FALSE)
FeaturePlot(object = combined, features = c("HSC_score1"), cols=c("paleturquoise1","magenta4"),order=TRUE, ncol=1,reduction = "umap",pt.size = 2)&coord_equal()
VlnPlot(object = combined, features = c("HSC_score1"),group.by="seurat_clusters", pt.size = 0)
#FeaturePlot(object = combined, features = c("cstem_score1"), cols=c("paleturquoise1","magenta4"),order=TRUE, ncol=1,reduction = "umap",pt.size = 2)&coord_equal()
#VlnPlot(object = combined, features = c("cstem_score1"),group.by="seurat_clusters",pt.size = 0)
d<-DotPlot(combined,assay="SCT",features=c("HLF", "ADGRG6", "CRHBP", "PCDH9","ICAM3","PROM1","BST2","CD164","CD37","CD44","KIT","MLLT3","EIF4A2","NPM1","RPL22"),group.by="seurat_clusters",cols=c("paleturquoise1", "violetred4"))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
print(d)
d<-DotPlot(combined,assay="SCT",features=GOI_supp,group.by="seurat_clusters",cols=c("paleturquoise1", "violetred4"))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
print(d)
d<-DotPlot(combined,assay="SCT",features=c("CSF3R","SMIM24","C1QTNF4","CSF1R","TOP2A","MKI67"),group.by="seurat_clusters",cols=c("paleturquoise1", "violetred4"))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
print(d)
d<-DotPlot(combined,assay="SCT",features=GOI_hayetal,group.by="seurat_clusters",cols=c("paleturquoise1", "violetred4"))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
print(d)
FeaturePlot(object = combined, features = c("AVP","SPINK2","NEGR1","CYGB","LGALS1","FCER1A"), cols=c("paleturquoise1","magenta4"),order=TRUE, ncol=1,reduction = "umap",pt.size = 2)&coord_equal()
dev.off()


#MLP=Multi Lineage Progenitor
#MEP= megakaryocyte/erythroid progenitors
#HSC= Hematopoetic Stem celldex
#MDP= monocyte/dendritic cell progenitors
#MKP= megakaryocyte progenitor
#GMP= granulocyte monocyte progenitors
#CLP= common lymphoid progenitors
#ERP= erythroid progenitors
#LMPP=lymphoid-primed multipotent progenitors
#PSC = patient specific cluster
Idents(combined)<-combined@meta.data$seurat_clusters
combined<-RenameIdents(object=combined,
            `0`="MLP",
            `1`="MEP",
            `2`="HSC1",
            `3`="PSC1",
            `4`="HSC2",
            `5`="HSC4",
            `6`="CLP_Ly2",
            `7`="MDP1",
            `8`="HSC3",
            `9`="MKP",
            `10`="GMP",
            `11`="HSC5",
            `12`="MDP2",
            `13`="Eo_B_Mast",
            `14`="PSC2",
            `15`="CLP_Ly2_B_PC",
            `16`="HSC6",
            `17`="LMPP_Ly1",
            `18`="ERP",
            `19`="HIST")

combined@meta.data$cell_type_specific<-Idents(combined)

saveRDS(combined,paste("all_integrated_merge_012022.rds",sep=""))


prop_plot1 <- proportion_plot_rj(
    seurat_obj = combined,
    feature1 = "cell_type_specific",
    feature2 = "sample"
)

prop_plot2 <- proportion_plot_rj(
    seurat_obj = combined,
    feature1 = "cell_type_specific",
    feature2 = "experiment"
)

pdf("score.pdf",height=12,width=10,useDingbats=FALSE)
print(DimPlot(combined, reduction = "umap",group.by="cell_type_specific",label=TRUE)&coord_equal())
print(prop_plot1)
print(prop_plot2)
dev.off()

################################
###SUBSET ONLY HSC POPULATIONS##
################################

HSC<-subset(x = combined, idents=c("HSC1","HSC2","HSC3","HSC4","HSC5","HSC6"))
HSC <- SCTransform(HSC, vars.to.regress = c("nCount_RNA","percent.mito"),return.only.var.genes = F)
HSC <- RunPCA(HSC, npcs = 30)
HSC <- RunUMAP(HSC, reduction = "pca", dims = 1:30)
HSC <- FindNeighbors(HSC, reduction = "pca", dims = 1:30,force.recalc=TRUE)
HSC <- FindClusters(HSC, resolution = 0.5)
HSC <- RunUMAP(HSC, dims = 1:30)
saveRDS(HSC,"HSC.rds")


##############################
##ADD ADDITIONAL ANNOTATIONS##
##############################

Idents(combined)<-combined@meta.data$experiment

combined<-RenameIdents(object=combined,
            `BL-8040`="Myeloma",
            `Healthy_BL-8040`="Healthy",
            `Healthy_Placebo`="Healthy",
            `Healthy_Plerixafor`="Healthy",
            `Placebo`="Myeloma",
            `Plerixafor`="Myeloma")

combined@meta.data$disease<-Idents(combined)


Idents(combined)<-combined@meta.data$experiment
combined<-RenameIdents(object=combined,
            `BL-8040`="BL-8040",
            `Healthy_BL-8040`="BL-8040",
            `Healthy_Placebo`="Placebo",
            `Healthy_Plerixafor`="Plerixafor",
            `Placebo`="Placebo",
            `Plerixafor`="Plerixafor")
combined@meta.data$treatment<-Idents(combined)

#############################################
#CREATE HEALTHY AND MYELOMA SPECIFIC OBJECTS#
#############################################

Idents(combined)<-combined@meta.data$disease
healthyonly<-subset(combined,idents=c("Healthy"))
healthyonly <- SCTransform(healthyonly, vars.to.regress = c("nCount_RNA","percent.mito"),return.only.var.genes = F)
healthyonly <- RunPCA(healthyonly, npcs = 30)
healthyonly <- RunUMAP(healthyonly, reduction = "pca", dims = 1:30)
healthyonly <- FindNeighbors(healthyonly, reduction = "pca", dims = 1:30,force.recalc=TRUE)
healthyonly <- FindClusters(healthyonly, resolution = 0.5)
healthyonly <- RunUMAP(healthyonly, dims = 1:30)
saveRDS(healthyonly,"healthy.rds")


myelomaonly<-subset(combined,idents=c("Myeloma"))
myelomaonly <- SCTransform(myelomaonly, vars.to.regress = c("nCount_RNA","percent.mito"),return.only.var.genes = F)
myelomaonly <- RunPCA(myelomaonly, npcs = 30)
myelomaonly <- RunUMAP(myelomaonly, reduction = "pca", dims = 1:30)
myelomaonly <- FindNeighbors(myelomaonly, reduction = "pca", dims = 1:30,force.recalc=TRUE)
myelomaonly <- FindClusters(myelomaonly, resolution = 0.5)
myelomaonly <- RunUMAP(myelomaonly, dims = 1:30)
saveRDS(myelomaonly,"myeloma.rds")

combined<-readRDS("~rjayasin/mmy_scratch/Dipersio/genesis/analysis/all_integrated_merge_012022.rds")
healthyonly<-readRDS("healthy.rds")

################################
#SELF RENEWAL AND QUIESCENCE####
################################

Idents(myeloma)<-myeloma@meta.data$cell_type_specific
myeloma_HSC<-subset(myeloma,idents=c("MLP","HSC1","HSC2","HSC3","HSC4","HSC5","HSC6"))
myeloma_HSC <- SCTransform(myeloma_HSC, vars.to.regress = c("nCount_RNA","percent.mito"),return.only.var.genes = F)
myeloma_HSC <- RunPCA(myeloma_HSC, npcs = 30)
myeloma_HSC <- RunUMAP(myeloma_HSC, reduction = "pca", dims = 1:30)
myeloma_HSC <- FindNeighbors(myeloma_HSC, reduction = "pca", dims = 1:30,force.recalc=TRUE)
myeloma_HSC <- FindClusters(myeloma_HSC, resolution = 0.5)
myeloma_HSC <- RunUMAP(myeloma_HSC, dims = 1:30)
saveRDS(myeloma_HSC,"myeloma_hsc.rds")


Idents(healthy)<-healthy@meta.data$cell_type_specific
healthy_HSC<-subset(healthy,idents=c("MLP","HSC1","HSC2","HSC3","HSC4","HSC5","HSC6"))
healthy_HSC <- SCTransform(healthy_HSC, vars.to.regress = c("nCount_RNA","percent.mito"),return.only.var.genes = F)
healthy_HSC <- RunPCA(healthy_HSC, npcs = 30)
healthy_HSC <- RunUMAP(healthy_HSC, reduction = "pca", dims = 1:30)
healthy_HSC <- FindNeighbors(healthy_HSC, reduction = "pca", dims = 1:30,force.recalc=TRUE)
healthy_HSC <- FindClusters(healthy_HSC, resolution = 0.5)
healthy_HSC <- RunUMAP(healthy_HSC, dims = 1:30)
saveRDS(healthy_HSC,"healthy_hsc.rds")

#https://pubmed.ncbi.nlm.nih.gov/34724565/
#selfreneweal and HSC quiescence genes
selfrenewal<-c("KAT7","MPL", "TEK", "GFI1B", "EGR1", "TAL1", "GATA2", "ERG", "PBX1", "MEIS1", "HOXA9","GATA1")

pdf(paste("plots.pdf", sep=""), width=10, height=8)
DotPlot(myeloma_hsc,assay="SCT",features=selfrenewal,split.by=c("treatment"), cols = "PiYG",group.by="cell_type_specific")+ggtitle("Myeloma HSC Only")+theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5))
#VlnPlot(myeloma_hsc,assay="SCT",features=selfrenewal,split.by=c("treatment"),group.by="cell_type_specific")+theme(legend.position="bottom")
#VlnPlot(myeloma_hsc,assay="SCT",features=selfrenewal,group.by="treatment")+theme(legend.position="bottom")
DotPlot(healthy_hsc,assay="SCT",features=selfrenewal,split.by=c("treatment"), cols = "PiYG",group.by="cell_type_specific")+ggtitle("Healthy HSC Only")+theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5))
#VlnPlot(healthy_hsc,assay="SCT",features=selfrenewal,split.by=c("treatment"),group.by="cell_type_specific")+theme(legend.position="bottom")
#VlnPlot(healthy_hsc,assay="SCT",features=selfrenewal,group.by="treatment")+theme(legend.position="bottom")
dev.off()

allHSC = merge(myeloma_hsc,y=c(healthy_hsc),project = "Merged")
allHSC <- SCTransform(allHSC, vars.to.regress = c("nCount_RNA","percent.mito"),return.only.var.genes = F)
allHSC <- RunPCA(allHSC, npcs = 30)
allHSC <- RunUMAP(allHSC, reduction = "pca", dims = 1:30)
allHSC <- FindNeighbors(allHSC, reduction = "pca", dims = 1:30,force.recalc=TRUE)
allHSC <- FindClusters(allHSC, resolution = 0.5)
allHSC <- RunUMAP(allHSC, dims = 1:30)
saveRDS(allHSC,"allHSC.rds")

selfrenewal<-c("KAT7","MPL", "TEK", "GFI1B", "EGR1", "TAL1", "GATA2", "ERG", "PBX1", "MEIS1", "HOXA9","GATA1")

papergenes<-c("EGR1","BTG2","JUNB","NR4A1","IER2","IL1B","TNFSF10","MYB")

pdf(paste("plots.pdf", sep=""), width=7, height=7)
print(DimPlot(allHSC,group.by=c("cell_type_specific"),split.by=c("treatment"))&coord_equal())
print(DimPlot(allHSC,group.by=c("cell_type_specific"),label=TRUE)&coord_equal())
print(DimPlot(allHSC,group.by=c("seurat_clusters"),label=TRUE)&coord_equal())
print(DimPlot(allHSC,group.by=c("cell_type_specific"),split.by=c("disease"))&coord_equal())
print(DotPlot(allHSC,assay="SCT",features=DEGs,split.by=c("treatment"), cols = "PiYG",group.by="cell_type_specific")+ggtitle("All HSC Only")+theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5)))
print(DotPlot(allHSC,assay="SCT",features=DEGs,split.by=c("treatment"), cols = "PiYG",group.by="disease")+ggtitle("All HSC Only")+theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5)))
print(DotPlot(allHSC,assay="SCT",features=selfrenewal,split.by=c("treatment"), cols = "PiYG",group.by="cell_type_specific")+ggtitle("All HSC Only")+theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5)))
print(DotPlot(allHSC,assay="SCT",features=selfrenewal,split.by=c("treatment"), cols = "PiYG",group.by="disease")+ggtitle("All HSC Only")+theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5)))
print(DotPlot(allHSC,assay="SCT",features=papergenes,split.by=c("treatment"), cols = "PiYG",group.by="cell_type_specific")+ggtitle("All HSC Only")+theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5)))
print(DotPlot(allHSC,assay="SCT",features=papergenes,split.by=c("treatment"), cols = "PiYG",group.by="disease")+ggtitle("All HSC Only")+theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5)))
print(DotPlot(allHSC,assay="RNA",features=papergenes,split.by=c("treatment2"), cols = "PiYG",group.by="disease")+ggtitle("All HSC Only")+theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5)))
print(DotPlot(allHSC,assay="SCT",features=c("TNF","NFKB1","NFKB2"),split.by=c("treatment"), cols = "PiYG",group.by="cell_type_specific")+ggtitle("All HSC Only")+theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5)))
print(DotPlot(allHSC,assay="SCT",features=c("TNF","NFKB1","NFKB2"),split.by=c("treatment"), cols = "PiYG",group.by="disease")+ggtitle("All HSC Only")+theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5)))
print(DotPlot(allHSC,assay="RNA",features=("TNF","NFKB1","NFKB2"),split.by=c("treatment2"), cols = "PiYG",group.by="disease")+ggtitle("All HSC Only")+theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5)))
dev.off()

###################
###Heatmap of DEGs#
###################

Idents(myeloma_HSC)<-myeloma_HSC@meta.data$treatment2
Idents(healthy_HSC)<-healthy_HSC@meta.data$treatment2

healthy_hsc<-readRDS("healthy_hsc.rds")
myeloma_hsc<-readRDS("myeloma_hsc.rds")


#Plerixafor(H) Motixafortide(H)         G-CSF(H)
datalist = list()

i=1
group1="Plerixafor(H)"
group2="G-CSF(H)"

i=2
group1="Motixafortide(H)"
group2="G-CSF(H)"

i=3
group1="Motixafortide(H)"
group2="Plerixafor(H)"

i=4
group1=c("Plerixafor(H)","Motixafortide(H)")
group2="G-CSF(H)"

df <- FindMarkers(object = healthy_hsc,only.pos=FALSE,min.pct=0.3,ident.1=c(paste0(group1)),ident.2=paste0(group2))
dfup = df[df$avg_log2FC  >= 1 & df$p_val_adj <= 0.05,]
dfdown = df[df$avg_log2FC  <= -1 & df$p_val_adj <= 0.05,]
dfup$ident1<-paste0(group1)
dfup$ident2<-paste0(group2)
dfdown$ident1<-paste0(group1)
dfdown$ident2<-paste0(group2)
dat<-rbind(dfup,dfdown)
dfdown$gene<-rownames(dfdown)
dfup$gene<-rownames(dfup)
dat<-rbind(dfup,dfdown)
rownames(dat)<-NULL
dat$i <- i 
datalist[[i]] <- dat

df <- FindMarkers(object = healthy_hsc,only.pos=FALSE,min.pct=0.3,ident.1=c(paste0(group1)),ident.2=paste0(group2))
dfup = df[df$avg_log2FC  >= 1 & df$p_val_adj <= 0.05,]
dfdown = df[df$avg_log2FC  <= -1 & df$p_val_adj <= 0.05,]
dfup$ident1<-paste0("CXCR4i")
dfup$ident2<-paste0(group2)
dfdown$ident1<-paste0("CXCR4i")
dfdown$ident2<-paste0(group2)
dfup$gene<-rownames(dfup)
dfdown$gene<-rownames(dfdown)
dat<-rbind(dfup,dfdown)
rownames(dat)<-NULL

dat$i <- i 
datalist[[i]] <- dat

big_data = do.call(rbind, datalist)
big_data$comparison<-paste0(big_data$ident1,"_",big_data$ident2)
#        p_val avg_log2FC pct.1 pct.2 p_val_adj        ident1   ident2 i
#H1FX        0   1.541079 0.990 0.963         0 Plerixafor(H) G-CSF(H) 1
#H1F0        0   1.255521 0.899 0.763         0 Plerixafor(H) G-CSF(H) 1
#TSC22D3     0   1.203081 0.889 0.805         0 Plerixafor(H) G-CSF(H) 1

p<-ggplot(big_data,aes(x=factor(comparison), 
  y=factor(gene), size=log10(p_val_adj+1), 
  fill=avg_log2FC,label=sprintf("%0.2f", round(avg_log2FC, digits = 2)))) +
        geom_point(alpha=0.8, shape=21) + geom_text(size=4)+
        scale_size(range = c(.1, 20), name="Proportion of Cells") +
        scale_fill_viridis(discrete=FALSE, option="plasma") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
        #ylab(paste0(feature2)) +
        #xlab(paste0(feature1))

pdf(paste("plots.pdf", sep=""), width=12, height=12)
print(p)
dev.off()

#########ORGANIZE GENE LIST FOR PLOTTING############

geneorder<-c("ZFP36","TXNIP","TSC22D3","NR4A1","EGR1","JUNB","FOSB","FOS","IER2","RPS26","H1FX","H1F0","S100A10","LGALS1","IFITM1","ID3","HSPA5","CRIP1","CD99","RNASET2")
big_data$geneorder <- as.character(big_data$gene)
big_data$geneorder <- factor(big_data$geneorder, levels=geneorder)

big_data$comparison <- as.character(big_data$comparison)
big_data$comparison <- factor(big_data$comparison, levels=c("Plerixafor(H)_G-CSF(H)","Motixafortide(H)_G-CSF(H)","CXCR4i_G-CSF(H)","Motixafortide(H)_Plerixafor(H)"))

###FINAL DEG LOG AVG FOLD CHANGE HEATMAP
pdf(paste("heatmap_DEG.pdf", sep=""), width=4, height=8)

#ggplot(big_data, aes(comparison, gene, fill= avg_log2FC,label=sprintf("%0.2f", round(avg_log2FC, digits = 2)))) + 
#  geom_tile() +geom_text(size=4)+
#  scale_fill_viridis(discrete=FALSE,option="plasma") +
#  theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

  ggplot(big_data, aes(comparison, geneorder, fill= avg_log2FC,label=sprintf("%0.2f", round(avg_log2FC, digits = 2)))) + 
  geom_tile() +geom_text(size=4)+
  scale_fill_gradient2(low="blue",mid = "white", high="red") +
  theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

dev.off()

highexpressors_Motixafortide<-c("ZFP36","TXNIP","TSC22D3","NR4A1","EGR1","JUNB","FOSB","FOS")
pdf(paste("plots.pdf", sep=""), width=6, height=6)
print(DotPlot(healthy_hsc,assay="SCT",features=highexpressors_Motixafortide,split.by=c("cell_type_specific"),group.by=("treatment2"), cols = "PiYG")+ggtitle("Healthy HSC Only")+theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5)))
dev.off()


gene_list<-unique(big_data$gene)

pdf(paste("plots.pdf", sep=""), width=8, height=6)
print(DotPlot(healthy_hsc,assay="SCT",features=gene_list,split.by=c("cell_type_specific"),group.by=("treatment2"), cols = "PiYG")+ggtitle("Healthy HSC Only")+theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5)))
print(DotPlot(myeloma_hsc,assay="RNA",features=gene_list,split.by=c("cell_type_specific"),group.by=("treatment2"), cols = "PiYG")+ggtitle("Myeloma HSC Only")+theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5)))
print(DotPlot(healthy_hsc,assay="SCT",features=gene_list,split.by=c("treatment2"), cols = "PiYG")+ggtitle("Healthy HSC Only")+theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5)))
print(DotPlot(myeloma_hsc,assay="RNA",features=gene_list,split.by=c("treatment2"), cols = "PiYG")+ggtitle("Myeloma HSC Only")+theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5)))
dev.off()


pdf(paste("plots.pdf", sep=""), width=12, height=12)
print(DimPlot(sobj,group.by=c("cell_type_specific"),cols=c("#73d055ff","#66cccc", "#1f968bff","#6633ff","#cc99ff","#996699","#660033",
 "#802582FF", "#952C80FF", "#AB337CFF", "#C03A76FF", "#D6456CFF",
"#E85362FF", "#F4685CFF", "#FA815FFF", "#FD9A6AFF", "#FEB37BFF", "#FECC8FFF",
"#FDE4A6FF", "#FCFDBFFF"),raster=TRUE,raster.dpi=c(1000, 1000),label=TRUE)&coord_equal())
print(DimPlot(allHSC,group.by=c("cell_type_specific"),split.by=c("treatment2"),cols=c("#73d055ff","#66cccc", "#1f968bff","#6633ff","#cc99ff","#996699","#660033",
"#481567ff"),ncol=3)&coord_equal()&scale_fill_viridis(option="magma"),raster=TRUE,raster.dpi=c(1000, 1000))
print(DimPlot(allHSC,group.by=c("cell_type_specific"),split.by=c("treatment2"),ncol=3,raster=TRUE,raster.dpi=c(1000, 1000))&coord_equal()&scale_fill_viridis(option="magma"))
dev.off()   

pdf(paste("plots.pdf", sep=""), width=12, height=12)
print(prop_plot1)
dev.off()


################################
####GENE LISTS STUFF############
################################


####Extract TNFA genes #filter out genes with lower than 10% of cells expressing gene of interest
#https://www.gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_TNFA_SIGNALING_VIA_NFKB.html
TNFA_genes<-c("ABCA1","ACKR3","AREG","ATF3","ATP2B1","B4GALT1","B4GALT5","BCL2A1","BCL3","BCL6","BHLHE40","BIRC2","BIRC3","BMP2","BTG1","BTG2","BTG3","CCL2","CCL20","CCL4","CCL5","CCN1","CCND1","CCNL1","CCRL2","CD44","CD69","CD80","CD83","CDKN1A","CEBPB","CEBPD","CFLAR","CLCF1","CSF1","CSF2","CXCL1","CXCL10","CXCL11","CXCL2","CXCL3","CXCL6","DDX58","DENND5A","DNAJB4","DRAM1","DUSP1","DUSP2","DUSP4","DUSP5","EDN1","EFNA1","EGR1","EGR2","EGR3","EHD1","EIF1","ETS2","F2RL1","F3","FJX1","FOS","FOSB","FOSL1","FOSL2","FUT4","G0S2","GADD45A","GADD45B","GCH1","GEM","GFPT2","GPR183","HBEGF","HES1","ICAM1","ICOSLG","ID2","IER2","IER3","IER5","IFIH1","IFIT2","IFNGR2","IL12B","IL15RA","IL18","IL1A","IL1B","IL23A","IL6","IL6ST","IL7R","INHBA","IRF1","IRS2","JAG1","JUN","JUNB","KDM6B","KLF10","KLF2","KLF4","KLF6","KLF9","KYNU","LAMB3","LDLR","LIF","LITAF","MAFF","MAP2K3","MAP3K8","MARCKS","MCL1","MSC","MXD1","MYC","NAMPT","NFAT5","NFE2L2","NFIL3","NFKB1","NFKB2","NFKBIA","NFKBIE","NINJ1","NR4A1","NR4A2","NR4A3","OLR1","PANX1","PDE4B","PDLIM5","PER1","PFKFB3","PHLDA1","PHLDA2","PLAU","PLAUR","PLEK","PLK2","PLPP3","PMEPA1","PNRC1","PPP1R15A","PTGER4","PTGS2","PTPRE","PTX3","RCAN1","REL","RELA","RELB","RHOB","RIPK2","RNF19B","SAT1","SDC4","SERPINB2","SERPINB8","SERPINE1","SGK1","SIK1","SLC16A6","SLC2A3","SLC2A6","SMAD3","SNN","SOCS3","SOD2","SPHK1","SPSB1","SQSTM1","STAT5A","TANK","TAP1","TGIF1","TIPARP","TLR2","TNC","TNF","TNFAIP2","TNFAIP3","TNFAIP6","TNFAIP8","TNFRSF9","TNFSF9","TNIP1","TNIP2","TRAF1","TRIB1","TRIP10","TSC22D1","TUBB2A","VEGFA","YRDC","ZBTB10","ZC3H12A","ZFP36")

a<-DotPlot(healthy_hsc,assay="SCT",features=TNFA_genes,split.by=c("treatment2"), cols = "PiYG")
a$data
top_TNFA_genes <- a$data %>% filter(pct.exp >= 10)
#do.call(paste, c(as.list(top_TNFA_genes$features.plot), sep = "','"))
GOI<-as.vector(unique(top_TNFA_genes$features.plot))
TNFAlist<-list(GOI)
healthy_hsc <- AddModuleScore(object = healthy_hsc,assay="SCT",features = TNFAlist,name = 'TNFA_NFKB')
pdf(paste("plots.pdf", sep=""), width=10, height=5)
#print(DotPlot(healthy_hsc,assay="SCT",features=GOI,split.by=c("treatment2"), cols = "PiYG")+theme(legend.position="bottom",axis.text.x = element_text(angle = 90, vjust = 0.5))+coord_flip())
VlnPlot(healthy_hsc,assay="SCT",features=c("TNFA_NFKB1"),group.by=c("treatment2"),split.by=c("cell_type_specific"))+theme(legend.position="bottom")
VlnPlot(healthy_hsc,assay="SCT",features=c("TNFA_NFKB1"),group.by=c("cell_type_specific"),split.by=c("treatment2"))+theme(legend.position="bottom")
VlnPlot(healthy_hsc,assay="SCT",features=c("TNFA_NFKB1"),group.by=c("treatment2"))+theme(legend.position="bottom")
VlnPlot(healthy_hsc,assay="SCT",features=c("TNFA_NFKB1"),group.by=c("cell_type_specific"))+theme(legend.position="bottom")
DotPlot(healthy_hsc,assay="SCT",features=c("TNFA_NFKB1"),group.by=c("treatment2"),split.by=c("cell_type_specific"), cols = "PiYG")
dev.off()

####Extract NAGASHIMA EGF UP genes #filter out genes with lower than 10% of cells expressing gene of interest
#https://www.gsea-msigdb.org/gsea/msigdb/cards/NAGASHIMA_EGF_SIGNALING_UP.html
#https://www.gsea-msigdb.org/gsea/msigdb/cards/AMIT_EGF_RESPONSE_40_HELA.html
EGF_genes<-c("AEN","ARC","AREG","ATF3","BCL10","BHLHE40","CCN1","CCN2","DLX2","DNAJB1","DUSP1","DUSP2","DUSP5","EDN1","EGR1","EGR3","EGR4","EPHA2","ETS2","F2RL1","FOS","FOSB","FOSL1","GEM","HBEGF","HES1","ID1","ID3","IER2","IER3","JUN","JUNB","KBTBD2","KDM6B","KLF10","KLF6","LIF","MAFF","MCL1","MIR22HG","MYC","NAB2","NAP1L1","NEDD9","NR4A1","NR4A2","NR4A3","PLEKHO2","RYBP","SIK1","SOWAHC","SPRED2","TIPARP","TNFRSF11B","TNFRSF12A","TRIB1","ZFP36","ATF3","ATL2","BCL3","BTG2","CAMKMT","CBX4","CCNA1","CCNL1","CDKN2AIP","CEBPB","DUSP1","DUSP2","DUSP5","EGR1","EGR3","HES1","ID1","IER2","IER3","IL6","JUN","JUNB","KCNJ12","KLF2","KLF6","LDLR","MBNL1","MBNL2","NFIB","NR4A1","PALMD","PIAS1","PIP5K1A","PPP1R15A","PTGS2","RNF141","SEC23A","SGK1","SLC2A3","UBE4B","UBR5","ZFP36")
EGF_uniqueM<-unique(EGF_genes)
a<-DotPlot(healthy_hsc,assay="SCT",features=EGF_uniqueM,split.by=c("treatment2"), cols = "PiYG")
top_EGF_genes <- a$data %>% filter(pct.exp >= 10)
#do.call(paste, c(as.list(top_TNFA_genes$features.plot), sep = "','"))
GOI<-as.vector(unique(top_EGF_genes$features.plot))
EGFlist<-list(GOI)
healthy_hsc <- AddModuleScore(object = healthy_hsc,assay="SCT",features = EGFlist,name = 'EGF_Score')
pdf(paste("plots.pdf", sep=""), width=10, height=5)
#print(DotPlot(healthy_hsc,assay="SCT",features=GOI,split.by=c("treatment2"), cols = "PiYG")+theme(legend.position="bottom",axis.text.x = element_text(angle = 90, vjust = 0.5))+coord_flip())
VlnPlot(healthy_hsc,assay="SCT",features=c("EGF_Score1"),group.by=c("treatment2"),split.by=c("cell_type_specific"))+theme(legend.position="bottom")
VlnPlot(healthy_hsc,assay="SCT",features=c("EGF_Score1"),group.by=c("cell_type_specific"),split.by=c("treatment2"))+theme(legend.position="bottom")
VlnPlot(healthy_hsc,assay="SCT",features=c("EGF_Score1"),group.by=c("treatment2"))+theme(legend.position="bottom")
VlnPlot(healthy_hsc,assay="SCT",features=c("EGF_Score1"),group.by=c("cell_type_specific"))+theme(legend.position="bottom")
DotPlot(healthy_hsc,assay="SCT",features=c("EGF_Score1"),group.by=c("treatment2"),split.by=c("cell_type_specific"), cols = "PiYG")
dev.off()

###Check Myeloma HSCs
myeloma_hsc <- AddModuleScore(object = myeloma_hsc,assay="SCT",features = TNFAlist,name = 'TNFA_NFKB')
myeloma_hsc <- AddModuleScore(object = myeloma_hsc,assay="SCT",features = EGFlist,name = 'EGF_Score')
allHSC <- AddModuleScore(object = allHSC,assay="SCT",features = TNFAlist,name = 'TNFA_NFKB')
allHSC <- AddModuleScore(object = allHSC,assay="SCT",features = EGFlist,name = 'EGF_Score')

pdf(paste("plots.pdf", sep=""), width=10, height=5)
#print(DotPlot(healthy_hsc,assay="SCT",features=GOI,split.by=c("treatment2"), cols = "PiYG")+theme(legend.position="bottom",axis.text.x = element_text(angle = 90, vjust = 0.5))+coord_flip())
VlnPlot(healthy_hsc,assay="SCT",features=c("TNFA_NFKB1"),group.by=c("cell_type_specific"),split.by=c("treatment2"),pt.size=0,cols=c("#fcfdbf","#fe9f6d","#de4968","#8c2981","#3b0f70","#000004"))+theme(legend.position="bottom")
VlnPlot(healthy_hsc,assay="SCT",features=c("EGF_Score1"),group.by=c("cell_type_specific"),split.by=c("treatment2"),pt.size=0,cols=c("#fcfdbf","#fe9f6d","#de4968","#8c2981","#3b0f70","#000004"))+theme(legend.position="bottom")
VlnPlot(myeloma_hsc,assay="SCT",features=c("TNFA_NFKB1"),group.by=c("cell_type_specific"),split.by=c("treatment2"),pt.size=0,cols=c("#fcfdbf","#fe9f6d","#de4968","#8c2981","#3b0f70","#000004"))+theme(legend.position="bottom")
VlnPlot(myeloma_hsc,assay="SCT",features=c("EGF_Score1"),group.by=c("cell_type_specific"),split.by=c("treatment2"),pt.size=0,cols=c("#fcfdbf","#fe9f6d","#de4968","#8c2981","#3b0f70","#000004"))+theme(legend.position="bottom")
VlnPlot(allHSC,assay="SCT",features=c("TNFA_NFKB1"),group.by=c("cell_type_specific"),split.by=c("treatment2"),pt.size=0,cols=c("#fcfdbf","#fe9f6d","#de4968","#8c2981","#3b0f70","#000004"))+theme(legend.position="bottom")
VlnPlot(allHSC,assay="SCT",features=c("EGF_Score1"),group.by=c("cell_type_specific"),split.by=c("treatment2"),pt.size=0,cols=c("#fcfdbf","#fe9f6d","#de4968","#8c2981","#3b0f70","#000004"))+theme(legend.position="bottom")
dev.off()

#Box/ViolinPlot
multi_dittoPlot(healthy_hsc, vars = c("HALLMARK_TNFA_SIGNALING_VIA_NFKB"), 
                group.by = "treatment2", plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Enrichment Scores", 
                theme = theme_classic())
multi_dittoPlot(healthy_hsc, vars = c("HALLMARK_TNFA_SIGNALING_VIA_NFKB"), 
                group.by = "cell_type_specific", plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Enrichment Scores", 
                theme = theme_classic())

s1<-read.table("overlap.txt",header=TRUE)

pdf(paste("GSetplots.pdf", sep=""), width=12, height=5)
ggplot(s1, aes(GeneSetName,-log10(FDRqvalue),size=Overlap,fill=FDRqvalue)) + 
  geom_point(alpha=0.8, shape=21)+
  scale_size(range = c(5,15))+
  scale_fill_viridis(discrete=FALSE,option="plasma") +
  theme_classic()+theme(legend.position="bottom",axis.text.x = element_text(angle = 90, vjust = 0.5))+coord_flip()
dev.off()   

saveRDS(sobj,"all_integrated_merge_012022.rds")

sobj<-readRDS("all_integrated_merge_012022.rds")
healthy_hsc<-readRDS("healthy_hsc.rds")
myeloma_hsc<-readRDS("myeloma_hsc.rds")


#conda activate monocle3_new
library(monocle3)
library(Seurat)
library(dplyr)
library(tibble)
library(stringr)
library(reshape2)
library(ggplot2)


#####################################
#######MONOCLE ANALYSIS##############
#####################################


#Try monocle plotting again but with less branching
#https://github.com/cole-trapnell-lab/monocle-release/issues/388
#####
#MAINTAIN ANCHORED SEURAT DATA!!!!!
#### Create a Monocle CDS Object
    # Project PC dimensions to whole data set

    sobj<-readRDS("all_integrated_merge_012022.rds")
    sobj <- ProjectDim(sobj, reduction = "pca")
  
    # Create an expression matrix
    expression_matrix <- sobj@assays$RNA@counts
    cell_metadata<-sobj@meta.data
    gene_annotation <- data.frame(gene_short_name=rownames(expression_matrix))
    rownames(gene_annotation) <- gene_annotation[,1]

    # Seurat-derived CDS
    my.cds <- new_cell_data_set(expression_matrix,
                                cell_metadata = cell_metadata,
                                gene_metadata = gene_annotation)
  
    # Transfer Seurat embeddings
    # Note that these may be calculated on the Integrated object, not the counts
    #   and thus will involve fewer genes
    reducedDim(my.cds, type = "PCA") <- sobj@reductions$pca@cell.embeddings 
    my.cds@preprocess_aux$prop_var_expl <- sobj@reductions$pca@stdev
    pdf(paste0("0_pcvariance.pdf"))
      plot_pc_variance_explained(my.cds)
    dev.off()
  
    # Transfer Seurat UMAP embeddings
    my.cds@int_colData@listData$reducedDims$UMAP <- sobj@reductions$umap@cell.embeddings
#    plot_cells(my.cds)
  
    # Copy cluster info from Seurat
    my.cds@clusters$UMAP_so$clusters <- sobj@meta.data$gt_tp_cell_type_integrated_.0.9

    #my.cds <- cluster_cells(my.cds, reduction_method = "UMAP", resolution = 1e-3)
    my.cds <- cluster_cells(my.cds, reduction_method = "UMAP")

    #DimPlot(sobj, reduction = "umap")
    #DimPlot(sobj, reduction = "umap")
    #plot_cells(my.cds, color_cells_by = "partition", group_label_size = 3.5)

  #Plot new dimensional reduction with cell_type calssifications
  pdf(paste0("1_Monocle_Cluster_seurat_v1.pdf"))
    a<-plot_cells(my.cds, label_groups_by_cluster=FALSE,  color_cells_by = "seurat_clusters")
    print(a)
  dev.off()
  pdf(paste0("1_Monocle_Cluster_origident_v1.pdf"))
    a<-plot_cells(my.cds, label_groups_by_cluster=FALSE,  color_cells_by = "experiment")
    print(a)
  dev.off()

  #Visualize genes of interest
#pdac_genes <- c("GRCh38-SELL","GRCh38-CCR7","GRCh38-LRRN3","GRCh38-GZMA","GRCh38-GZMB","GRCh38-PRF1")
#pdf(paste0("2_gene_expression.pdf"))
#  a<-plot_cells(my.cds,
#             genes=pdac_genes,
#             label_cell_groups=FALSE,
#             show_trajectory_graph=FALSE)
#  print(a)
#dev.off()


#Recluster cells
cds <- cluster_cells(my.cds,reduction_method = "UMAP")
pdf(paste0("3_cell_clusters_v1.pdf"))
  a<-plot_cells(my.cds, color_cells_by = "partition")
  print(a)
dev.off()

#Learn trajectory graph
my.cds <- learn_graph(my.cds, use_partition = TRUE)
pdf(paste0("4_trajectory_graph_clusters_v1.pdf"))
  a<-plot_cells(my.cds,
             color_cells_by = "cell_type_specific",
             label_groups_by_cluster=FALSE,
             label_leaves=FALSE,
             label_branch_points=FALSE)
  print(a)
dev.off()

pdf(paste0("4_trajectory_graph_experiment_v1.pdf"))
  a<-plot_cells(my.cds,
             color_cells_by = "treatment2",
             label_groups_by_cluster=FALSE,
             label_leaves=FALSE,
             label_branch_points=FALSE)
  print(a)
dev.off()


#ORder cells in pseudotime
pdf(paste0("5_pseudotime_nodes_v1.pdf"))
  a<-plot_cells(my.cds,
             color_cells_by = "cell_type_specific",
             label_cell_groups=FALSE,
             label_leaves=TRUE,
             label_branch_points=TRUE,
             graph_label_size=1.5)
  print(a)
dev.off()

#You need to specify the root node of the trajectory here

######################NODE FOR ADM 
get_earliest_principal_node <- function(my.cds, time_bin="HSC1"){
  cell_ids <- which(colData(my.cds)[, "cell_type_specific"] == time_bin)
  
  closest_vertex <-
  my.cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(my.cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(my.cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

get_earliest_principal_node(my.cds)
my.cds <- order_cells(my.cds, root_pr_nodes=get_earliest_principal_node(my.cds))

saveRDS(my.cds,"all_integrated_merge_042022_monocle3_v2.rds")

pdf(paste0("6_pseudotime_plot_v1.pdf"))
  a<-plot_cells(my.cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5, rasterize =TRUE)&coord_equal()
  print(a)
dev.off()

b<-DimPlot(sobj,group.by=c("cell_type_specific"),cols=c("#73d055ff","#66cccc", "#1f968bff","#6633ff","#cc99ff","#996699","#660033",
 "#802582FF", "#952C80FF", "#AB337CFF", "#C03A76FF", "#D6456CFF",
"#E85362FF", "#F4685CFF", "#FA815FFF", "#FD9A6AFF", "#FEB37BFF", "#FECC8FFF",
"#FDE4A6FF", "#FCFDBFFF"),raster=TRUE,label=TRUE)&coord_equal()&No_Legend()

pdf(paste0("6_pseudotime_plot_v1.pdf"),height=8,width=8,useDingbats=FALSE)
b + a
dev.off()


pdf(paste0("6_pseudotime_celltype_plot_v1.pdf"))
  a<-plot_cells(my.cds,
           color_cells_by = "cell_type_specific",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
  print(a)


  sobj@meta.data$cell_type_specific <- as.character(sobj@meta.data$cell_type_specific)
sobj@meta.data$cell_type_specific <- factor(sobj@meta.data$cell_type_specific, levels=c("HSC1","HSC2","HSC3","HSC4","HSC5","HSC6",
  "MLP","MKP","MEP","ERP","HIST","Eo_B_Mast",
  "GMP","MDP2","MDP1","CLP_Ly2_B_PC","CLP_Ly2","LMPP_Ly1",
  "PSC2","PSC1"))


###HEATMAP INSTEAD OF DOTPLOT FOR enriched populations
sobj<-readRDS("~rjayasin/mmy_scratch/Dipersio/genesis/analysis/all_integrated_merge_012022.rds")
#Make proportion table - extract proprotions and counts for plotting
data<-table(sobj@meta.data$sample,sobj@meta.data$cell_type)
#The value of each cell divided by the sum of the row cells - row is experiment in this case
data_prop<-prop.table(data,1)
data_adj<-data_prop*100
plot_data<-melt(data_adj)
plot_data$key<-paste0(plot_data$Var1,"_",plot_data$Var2)
colnames(plot_data)<-c("Sample","CellType","Proportion","key")
og_data<-melt(data)
og_data$key<-paste0(og_data$Var1,"_",og_data$Var2)
og_data$Var1<-NULL
og_data$Var2<-NULL
colnames(og_data)<-c("count","key")
merged_data<-merge(plot_data, og_data, by.x="key", by.y="key")


merged_data$Sample <- as.character(merged_data$Sample)
merged_data$Sample <- factor(merged_data$Sample, levels=c("S-005-G005","S-009-G009","D026-A026","DIV29-A029","BL010-BL010","BL011-BL011",
  "012-015-P15","012-037-P37","012-045-P45","012-055-P55","1096-A1096","1100-A1100","1102-A1102","1104-A1104","012-036-B36","012-049-B49","012-052-B012-052","012-059-B012-059"))


pdf(paste("heatmap_propplot.pdf", sep=""), width=12, height=8)
  ggplot(merged_data, aes(CellType, Sample, fill= Proportion,label=sprintf("%0.2f", round(Proportion, digits = 2)))) + 
  geom_tile() +geom_text(size=4)+
  scale_fill_gradient2(low="blue",mid = "white", high="red") +
  theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

dev.off()


###Gene list of interest for final plot
library(ggpubr)

healthy_hsc<-readRDS("healthy_hsc.rds")
myeloma_hsc<-readRDS("myeloma_hsc.rds")

Idents(healthy_hsc)<-healthy_hsc@meta.data$sample
healthy_hsc<-RenameIdents(object=healthy_hsc,
`012-037-P37`="G-CSF(M)",      `012-049-B49` ="Motixafortide+G-CSF(M)",     `1104-A1104` ="Plerixafor+G-CSF(M)",      `012-015-P15`="G-CSF(M)",
`012-045-P45`="G-CSF(M)",      `012-055-P55`="G-CSF(M)",      `012-036-B36`  ="Motixafortide+G-CSF(M)",    `1102-A1102`="Plerixafor+G-CSF(M)",
`1100-A1100` ="Plerixafor+G-CSF(M)",      `DIV29-A029`="Plerixafor(H)",       `012-052-B012-052`="Motixafortide+G-CSF(M)", `BL010-BL010`="Motixafortide(H)",
`S-005-G005` ="G-CSF(H)",      `012-059-B012-059`="Motixafortide+G-CSF(M)", `BL011-BL011`="Motixafortide(H)",      `S-009-G009`="G-CSF(H)",
`D026-A026`  ="Plerixafor(H)",      `1096-A1096`="Plerixafor+G-CSF(M)")
healthy_hsc@meta.data$treatment2<-Idents(healthy_hsc)

Idents(myeloma_hsc)<-myeloma_hsc@meta.data$sample
myeloma_hsc<-RenameIdents(object=myeloma_hsc,
`012-037-P37`="G-CSF(M)",      `012-049-B49` ="Motixafortide+G-CSF(M)",     `1104-A1104` ="Plerixafor+G-CSF(M)",      `012-015-P15`="G-CSF(M)",
`012-045-P45`="G-CSF(M)",      `012-055-P55`="G-CSF(M)",      `012-036-B36`  ="Motixafortide+G-CSF(M)",    `1102-A1102`="Plerixafor+G-CSF(M)",
`1100-A1100` ="Plerixafor+G-CSF(M)",      `DIV29-A029`="Plerixafor(H)",       `012-052-B012-052`="Motixafortide+G-CSF(M)", `BL010-BL010`="Motixafortide(H)",
`S-005-G005` ="G-CSF(H)",      `012-059-B012-059`="Motixafortide+G-CSF(M)", `BL011-BL011`="Motixafortide(H)",      `S-009-G009`="G-CSF(H)",
`D026-A026`  ="Plerixafor(H)",      `1096-A1096`="Plerixafor+G-CSF(M)")
myeloma_hsc@meta.data$treatment2<-Idents(myeloma_hsc)



myeloma_hsc@meta.data$treatment2 <- as.character(myeloma_hsc@meta.data$treatment2)
myeloma_hsc@meta.data$treatment2 <- factor(myeloma_hsc@meta.data$treatment2, levels=c("G-CSF(H)","Plerixafor(H)","Motixafortide(H)","G-CSF(M)","Plerixafor+G-CSF(M)","Motixafortide+G-CSF(M)"))

healthy_hsc@meta.data$treatment2 <- as.character(healthy_hsc@meta.data$treatment2)
healthy_hsc@meta.data$treatment2 <- factor(healthy_hsc@meta.data$treatment2, levels=c("G-CSF(H)","Plerixafor(H)","Motixafortide(H)","G-CSF(M)","Plerixafor+G-CSF(M)","Motixafortide+G-CSF(M)"))

myeloma_hsc@meta.data$cell_type_specific <- as.character(myeloma_hsc@meta.data$cell_type_specific)
myeloma_hsc@meta.data$cell_type_specific <- factor(myeloma_hsc@meta.data$cell_type_specific, levels=c("MLP","HSC6","HSC5","HSC4","HSC3","HSC2","HSC1"))

healthy_hsc@meta.data$cell_type_specific <- as.character(healthy_hsc@meta.data$cell_type_specific)
healthy_hsc@meta.data$cell_type_specific <- factor(healthy_hsc@meta.data$cell_type_specific, levels=c("MLP","HSC6","HSC5","HSC4","HSC3","HSC2","HSC1"))


healthy_hsc@meta.data$treatment3 <- paste0(healthy_hsc@meta.data$cell_type_specific,"_",healthy_hsc@meta.data$treatment2)
myeloma_hsc@meta.data$treatment3 <- paste0(myeloma_hsc@meta.data$cell_type_specific,"_",myeloma_hsc@meta.data$treatment2)

healthy_hsc@meta.data$treatment3 <- as.character(healthy_hsc@meta.data$treatment3)
healthy_hsc@meta.data$treatment3 <- factor(healthy_hsc@meta.data$treatment3, levels=c ("MLP_G-CSF(H)","MLP_Plerixafor(H)","MLP_Motixafortide(H)","HSC6_G-CSF(H)","HSC6_Plerixafor(H)","HSC6_Motixafortide(H)","HSC5_G-CSF(H)","HSC5_Plerixafor(H)","HSC5_Motixafortide(H)","HSC4_G-CSF(H)","HSC4_Plerixafor(H)","HSC4_Motixafortide(H)","HSC3_G-CSF(H)","HSC3_Plerixafor(H)","HSC3_Motixafortide(H)","HSC2_G-CSF(H)","HSC2_Plerixafor(H)","HSC2_Motixafortide(H)","HSC1_G-CSF(H)","HSC1_Plerixafor(H)","HSC1_Motixafortide(H)"))


myeloma_hsc@meta.data$treatment3 <- as.character(myeloma_hsc@meta.data$treatment3)
myeloma_hsc@meta.data$treatment3 <- factor(myeloma_hsc@meta.data$treatment3, levels=c("MLP_G-CSF(M)","MLP_Plerixafor+G-CSF(M)","MLP_Motixafortide+G-CSF(M)","HSC6_G-CSF(M)","HSC6_Plerixafor+G-CSF(M)","HSC6_Motixafortide+G-CSF(M)","HSC5_G-CSF(M)","HSC5_Plerixafor+G-CSF(M)","HSC5_Motixafortide+G-CSF(M)","HSC4_G-CSF(M)","HSC4_Plerixafor+G-CSF(M)","HSC4_Motixafortide+G-CSF(M)","HSC3_G-CSF(M)","HSC3_Plerixafor+G-CSF(M)","HSC3_Motixafortide+G-CSF(M)","HSC2_G-CSF(M)","HSC2_Plerixafor+G-CSF(M)","HSC2_Motixafortide+G-CSF(M)","HSC1_G-CSF(M)","HSC1_Plerixafor+G-CSF(M)","HSC1_Motixafortide+G-CSF(M)"))



#https://pubmed.ncbi.nlm.nih.gov/34724565/
#selfreneweal and HSC quiescence genes
selfrenewal<-c("EGR1","BTG2","JUNB","NR4A1","IER2","MYB", #https://stemcellres.biomedcentral.com/articles/10.1186/s13287-021-02498-0
  "EGR3",#https://www.nature.com/articles/s41556-020-0512-1
  #"FLT3", "CD34",#https://genome.cshlp.org/content/25/12/1860.long - decrease with age and increase with differentiation...
   "DNMT3B", #pluripotency related 25130491
   "TET2","BMI1", # https://ashpublications.org/blood/article/132/8/791/39470/Single-cell-approaches-identify-the-molecular - also mentions PBX1 and MEIS1
  "PBX1","MEIS1","KAT7","TAL1","GATA2","ERG","HOXA9","GATA1" #https://pubmed.ncbi.nlm.nih.gov/34724565/
  )
#https://pubmed.ncbi.nlm.nih.gov/31230859/#:~:text=We%20show%20that%20while%20inducing,functions%2C%20and%20poises%20HSCs%20for
#https://pubmed.ncbi.nlm.nih.gov/18397757/ - EGR1 specific 


pdf(paste("gene_dotplot_04192022.pdf", sep=""), width=15, height=5.5)
#a<-DotPlot(healthy_hsc,assay="SCT",features=selfrenewal,group.by=c("treatment3"), cols = "PiYG",scale.max = 100, scale.min = 0, col.min = -2.5,col.max = 2.5)+ggtitle("Healthy HSC Only")+theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5))&guides(x =  guide_axis(angle = 90))&theme_light()
#b<-DotPlot(myeloma_hsc,assay="SCT",features=selfrenewal,group.by=c("treatment3"), cols = "PiYG",scale.max = 100, scale.min = 0, col.min = -2.5,col.max = 2.5)+ggtitle("Myeloma HSC Only")+theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5))&guides(x =  guide_axis(angle = 90))&theme_light()
#a+b
#a<-DotPlot(healthy_hsc,assay="SCT",features=selfrenewal,group.by=c("treatment3"), cols = "PRGn",scale.max = 100, scale.min = 0, col.min = -2.5,col.max = 2.5)+ggtitle("Healthy HSC Only")+theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5))&guides(x =  guide_axis(angle = 90))&theme_light()
#b<-DotPlot(myeloma_hsc,assay="SCT",features=selfrenewal,group.by=c("treatment3"), cols = "PRGn",scale.max = 100, scale.min = 0, col.min = -2.5,col.max = 2.5)+ggtitle("Myeloma HSC Only")+theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5))&guides(x =  guide_axis(angle = 90))&theme_light()
#a+b
a<-DotPlot(healthy_hsc,assay="SCT",features=selfrenewal,group.by=c("treatment3"), cols = "RdGy",scale.max = 100, scale.min = 0, col.min = -2.5,col.max = 2.5)+ggtitle("Healthy HSC Only")+theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5))&guides(x =  guide_axis(angle = 90))&theme_light()
b<-DotPlot(myeloma_hsc,assay="SCT",features=selfrenewal,group.by=c("treatment3"), cols = "RdGy",scale.max = 100, scale.min = 0, col.min = -2.5,col.max = 2.5)+ggtitle("Myeloma HSC Only")+theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5))&guides(x =  guide_axis(angle = 90))&theme_light()
a+b
dev.off()


pdf(paste("SARM1.pdf", sep=""), width=4, height=5)
#a<-DotPlot(sobj,assay="SCT",features="SARM1",group.by=c("cell_type_specific"), cols = "RdGy",scale.max = 100, scale.min = 0, col.min = -2.5,col.max = 2.5)+ggtitle("Healthy HSC Only")+theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5))&guides(x =  guide_axis(angle = 90))&theme_light()
a<-DotPlot(sobj,assay="SCT",features="SARM1",group.by=c("cell_type_specific"), cols = "RdGy",scale.max = 100, scale.min = 0)+ggtitle("Healthy and Myeloma Cell Types")+theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5))&guides(x =  guide_axis(angle = 90))&theme_light()
#b<-DotPlot(myeloma_hsc,assay="SCT",features="SARM1",group.by=c("treatment3"), cols = "RdGy",scale.max = 100, scale.min = 0, col.min = -2.5,col.max = 2.5)+ggtitle("Myeloma HSC Only")+theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5))&guides(x =  guide_axis(angle = 90))&theme_light()
#a+b
print(a)
dev.off()

allHSC<-readRDS("allHSC.rds")

Idents(allHSC)<-allHSC@meta.data$sample
allHSC<-RenameIdents(object=allHSC,
`012-037-P37`="G-CSF(M)",      `012-049-B49` ="Motixafortide+G-CSF(M)",     `1104-A1104` ="Plerixafor+G-CSF(M)",      `012-015-P15`="G-CSF(M)",
`012-045-P45`="G-CSF(M)",      `012-055-P55`="G-CSF(M)",      `012-036-B36`  ="Motixafortide+G-CSF(M)",    `1102-A1102`="Plerixafor+G-CSF(M)",
`1100-A1100` ="Plerixafor+G-CSF(M)",      `DIV29-A029`="Plerixafor(H)",       `012-052-B012-052`="Motixafortide+G-CSF(M)", `BL010-BL010`="Motixafortide(H)",
`S-005-G005` ="G-CSF(H)",      `012-059-B012-059`="Motixafortide+G-CSF(M)", `BL011-BL011`="Motixafortide(H)",      `S-009-G009`="G-CSF(H)",
`D026-A026`  ="Plerixafor(H)",      `1096-A1096`="Plerixafor+G-CSF(M)")
allHSC@meta.data$treatment2<-Idents(allHSC)



####HSC Model DotPlot

 data<-table(allHSC@meta.data$cell_type_specific,allHSC@meta.data$treatment2)

    #The value of each cell divided by the sum of the row cells - row is feature2
    data_prop<-prop.table(data,2)
    data_adj<-data_prop*100
    plot_data<-melt(data_adj)


plot_data$Var2 <- as.character(plot_data$Var2)
plot_data$Var2 <- factor(plot_data$Var2, levels=c("G-CSF(H)","Plerixafor(H)","Motixafortide(H)","G-CSF(M)","Plerixafor+G-CSF(M)","Motixafortide+G-CSF(M)"))

plot_data$Var1 <- as.character(plot_data$Var1 )
plot_data$Var1  <- factor(plot_data$Var1,levels=c("MLP","HSC6","HSC5","HSC4","HSC3","HSC2","HSC1"))


 p<-ggplot(plot_data,aes(x=factor(Var2), y=factor(Var1), size=value, fill=factor(Var1),label=sprintf("%0.2f", round(value, digits = 2)))) +
        geom_point(alpha=0.8, shape=21, color="black") + geom_text(size=4)+
        scale_size(range = c(.1, 20), name="Proportion of Cells") +
        scale_fill_viridis(discrete=TRUE, option="A") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
        ylab("Treatment") +
        xlab("Cell Type")


pdf(paste("dotplot_cohort_model_04292022.pdf", sep=""), width=5.5, height=5.5)
print(p)
dev.off()


#############Heatmap enrichment plot
sobj<-readRDS("all_integrated_merge_012022.rds")

data<-table(sobj@meta.data$cell_type_specific,sobj@meta.data$sample)

#The value of each cell divided by the sum of the row cells - row is feature2
data_prop<-prop.table(data,2)
data_adj<-data_prop*100
    #plot_data<-melt(data_adj)

data_adj$Sample <- as.character(data_adj$Sample)
data_adj$Sample <- factor(data_adj$Sample, levels=c("S-005-G005","S-009-G009","D026-A026","DIV29-A029","BL010-BL010","BL011-BL011",
  "012-015-P15","012-037-P37","012-045-P45","012-055-P55","1096-A1096","1100-A1100","1102-A1102","1104-A1104","012-036-B36","012-049-B49","012-052-B012-052","012-059-B012-059"))

hclust_cols=c("S-005-G005","S-009-G009","D026-A026","DIV29-A029","BL010-BL010","BL011-BL011",
  "012-015-P15","012-037-P37","012-045-P45","012-055-P55","1096-A1096","1100-A1100","1102-A1102","1104-A1104","012-036-B36","012-049-B49","012-052-B012-052","012-059-B012-059")

hclust_rows=c("PSC1","PSC2","HIST","CLP_Ly2_B_PC","CLP_Ly2","LMPP_Ly1","MDP1","MDP2","GMP","MKP","Eo_B_Mast","ERP","MEP","MLP","HSC6","HSC5","HSC4","HSC3","HSC2","HSC1")

ordered_data<-data_adj[hclust_rows,hclust_cols]

pdf(paste("heatmap_sample_04292022.pdf", sep=""), width=8, height=5)
pheatmap(ordered_data, display_numbers = T,color = colorRampPalette(c("burlywood1", "darkred"))(50),cluster_rows = F, cluster_cols = F, fontsize_number = 9)
pheatmap(ordered_data, display_numbers = T,color = colorRampPalette(c("white", "darkred"))(50),cluster_rows = F, cluster_cols = F, fontsize_number = 9)
pheatmap(ordered_data, display_numbers = T,color = colorRampPalette(c("darkslategray1", "darkred"))(50),cluster_rows = F, cluster_cols = F, fontsize_number = 9)
pheatmap(ordered_data, display_numbers = T,color = colorRampPalette(c("azure", "darkslategray"))(50),cluster_rows = F, cluster_cols = F, fontsize_number = 9)
pheatmap(ordered_data, display_numbers = T,color = colorRampPalette(c("azure", "deeppink4"))(50),cluster_rows = F, cluster_cols = F, fontsize_number = 9)
dev.off()


###SCORE FIGURE

####Extract TNFA genes #filter out genes with lower than 10% of cells expressing gene of interest
#https://www.gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_TNFA_SIGNALING_VIA_NFKB.html
TNFA_genes<-c("ABCA1","ACKR3","AREG","ATF3","ATP2B1","B4GALT1","B4GALT5","BCL2A1","BCL3","BCL6","BHLHE40","BIRC2","BIRC3","BMP2","BTG1","BTG2","BTG3","CCL2","CCL20","CCL4","CCL5","CCN1","CCND1","CCNL1","CCRL2","CD44","CD69","CD80","CD83","CDKN1A","CEBPB","CEBPD","CFLAR","CLCF1","CSF1","CSF2","CXCL1","CXCL10","CXCL11","CXCL2","CXCL3","CXCL6","DDX58","DENND5A","DNAJB4","DRAM1","DUSP1","DUSP2","DUSP4","DUSP5","EDN1","EFNA1","EGR1","EGR2","EGR3","EHD1","EIF1","ETS2","F2RL1","F3","FJX1","FOS","FOSB","FOSL1","FOSL2","FUT4","G0S2","GADD45A","GADD45B","GCH1","GEM","GFPT2","GPR183","HBEGF","HES1","ICAM1","ICOSLG","ID2","IER2","IER3","IER5","IFIH1","IFIT2","IFNGR2","IL12B","IL15RA","IL18","IL1A","IL1B","IL23A","IL6","IL6ST","IL7R","INHBA","IRF1","IRS2","JAG1","JUN","JUNB","KDM6B","KLF10","KLF2","KLF4","KLF6","KLF9","KYNU","LAMB3","LDLR","LIF","LITAF","MAFF","MAP2K3","MAP3K8","MARCKS","MCL1","MSC","MXD1","MYC","NAMPT","NFAT5","NFE2L2","NFIL3","NFKB1","NFKB2","NFKBIA","NFKBIE","NINJ1","NR4A1","NR4A2","NR4A3","OLR1","PANX1","PDE4B","PDLIM5","PER1","PFKFB3","PHLDA1","PHLDA2","PLAU","PLAUR","PLEK","PLK2","PLPP3","PMEPA1","PNRC1","PPP1R15A","PTGER4","PTGS2","PTPRE","PTX3","RCAN1","REL","RELA","RELB","RHOB","RIPK2","RNF19B","SAT1","SDC4","SERPINB2","SERPINB8","SERPINE1","SGK1","SIK1","SLC16A6","SLC2A3","SLC2A6","SMAD3","SNN","SOCS3","SOD2","SPHK1","SPSB1","SQSTM1","STAT5A","TANK","TAP1","TGIF1","TIPARP","TLR2","TNC","TNF","TNFAIP2","TNFAIP3","TNFAIP6","TNFAIP8","TNFRSF9","TNFSF9","TNIP1","TNIP2","TRAF1","TRIB1","TRIP10","TSC22D1","TUBB2A","VEGFA","YRDC","ZBTB10","ZC3H12A","ZFP36")
a<-DotPlot(sobj,assay="SCT",features=TNFA_genes,split.by=c("treatment2"), cols = "PiYG")
top_TNFA_genes <- a$data %>% filter(pct.exp >= 10)
GOI<-as.vector(unique(top_TNFA_genes$features.plot))
TNFAlist<-list(GOI)
####Extract NAGASHIMA EGF UP genes #filter out genes with lower than 10% of cells expressing gene of interest
#https://www.gsea-msigdb.org/gsea/msigdb/cards/NAGASHIMA_EGF_SIGNALING_UP.html
#https://www.gsea-msigdb.org/gsea/msigdb/cards/AMIT_EGF_RESPONSE_40_HELA.html
EGF_genes<-c("AEN","ARC","AREG","ATF3","BCL10","BHLHE40","CCN1","CCN2","DLX2","DNAJB1","DUSP1","DUSP2","DUSP5","EDN1","EGR1","EGR3","EGR4","EPHA2","ETS2","F2RL1","FOS","FOSB","FOSL1","GEM","HBEGF","HES1","ID1","ID3","IER2","IER3","JUN","JUNB","KBTBD2","KDM6B","KLF10","KLF6","LIF","MAFF","MCL1","MIR22HG","MYC","NAB2","NAP1L1","NEDD9","NR4A1","NR4A2","NR4A3","PLEKHO2","RYBP","SIK1","SOWAHC","SPRED2","TIPARP","TNFRSF11B","TNFRSF12A","TRIB1","ZFP36","ATF3","ATL2","BCL3","BTG2","CAMKMT","CBX4","CCNA1","CCNL1","CDKN2AIP","CEBPB","DUSP1","DUSP2","DUSP5","EGR1","EGR3","HES1","ID1","IER2","IER3","IL6","JUN","JUNB","KCNJ12","KLF2","KLF6","LDLR","MBNL1","MBNL2","NFIB","NR4A1","PALMD","PIAS1","PIP5K1A","PPP1R15A","PTGS2","RNF141","SEC23A","SGK1","SLC2A3","UBE4B","UBR5","ZFP36")
EGF_uniqueM<-unique(EGF_genes)
a<-DotPlot(sobj,assay="SCT",features=EGF_uniqueM,split.by=c("treatment2"), cols = "PiYG")
top_EGF_genes <- a$data %>% filter(pct.exp >= 10)
#do.call(paste, c(as.list(top_TNFA_genes$features.plot), sep = "','"))
GOI<-as.vector(unique(top_EGF_genes$features.plot))
EGFlist<-list(GOI)
###Check Myeloma HSCs
sobj <- AddModuleScore(object = sobj,assay="SCT",features = TNFAlist,name = 'TNFA_NFKB')
sobj <- AddModuleScore(object = sobj,assay="SCT",features = EGFlist,name = 'EGF_Score')

##HBO1/KAT7 Score
#https://pubmed.ncbi.nlm.nih.gov/34724565/
#selfreneweal and HSC quiescence genes
selfrenewal<-c("KAT7","MPL", "TEK", "GFI1B", "EGR1", "TAL1", "GATA2", "ERG", "PBX1", "MEIS1", "HOXA9","GATA1")
a<-DotPlot(sobj,assay="SCT",features=selfrenewal,split.by=c("treatment2"), cols = "PiYG")
top_HBO1_genes <- a$data %>% filter(pct.exp >= 10)
GOI<-as.vector(unique(top_HBO1_genes$features.plot))
HBO1list<-list(GOI)
sobj <- AddModuleScore(object = sobj,assay="SCT",features = HBO1list,name = 'HBOScore')

Idents(sobj)<-sobj@meta.data$cell_type_specific
HSConly<-subset(sobj,idents=c("HSC1","HSC2","HSC3","HSC4","HSC5","HSC6","MLP"))


HSConly@meta.data$treatment2 <- as.character(HSConly@meta.data$treatment2)
HSConly@meta.data$treatment2 <- factor(HSConly@meta.data$treatment2, levels=c("G-CSF(H)","Plerixafor(H)","Motixafortide(H)","G-CSF(M)","Plerixafor+G-CSF(M)","Motixafortide+G-CSF(M)"))

HSConly@meta.data$cell_type_specific <- as.character(HSConly@meta.data$cell_type_specific)
HSConly@meta.data$cell_type_specific <- factor(HSConly@meta.data$cell_type_specific,levels=c("HSC1","HSC2","HSC3","HSC4","HSC5","HSC6","MLP"))

a<-VlnPlot(HSConly,assay="SCT",features=c("TNFA_NFKB1"),group.by=c("cell_type_specific"),split.by=c("treatment2"),pt.size=0,cols=c("#fcfdbf","#fe9f6d","#de4968","#8c2981","#3b0f70","#000004"))+theme(legend.position="right")
b<-VlnPlot(HSConly,assay="SCT",features=c("EGF_Score1"),group.by=c("cell_type_specific"),split.by=c("treatment2"),pt.size=0,cols=c("#fcfdbf","#fe9f6d","#de4968","#8c2981","#3b0f70","#000004"))+theme(legend.position="right")
c<-VlnPlot(HSConly,assay="SCT",features=c("HBOScore1"),group.by=c("cell_type_specific"),split.by=c("treatment2"),pt.size=0,cols=c("#fcfdbf","#fe9f6d","#de4968","#8c2981","#3b0f70","#000004"))+theme(legend.position="right")

pdf(paste("Gene_Scores.pdf", sep=""), width=10, height=10)
ggarrange(a,b,c,
          ncol = 1, nrow = 3)
dev.off()


### Add pvalue to plots with gene score
sobj<-readRDS("~rjayasin/mmy_scratch/Dipersio/genesis/analysis/all_integrated_merge_012022.rds")

DEGlist<-c(
"CRHBP", "HLF", "DUSP1","ADGRG6","PCDH9","AVP", #HSC
"C1QTNF4","SPINK2","CSF3R","SMIM24", #MLP
"ITGA2B", "PF4", "VWF","PLEK","TOP2A","MKI67","MYOZ3","CMTM5", #MKP
"ELANE", "MPO", "LYZ", "CSF1R", "CTSG", "PRTN3", "AZU1" ,"MS4A6A", "ANXA2", #GMP
"FCER1A","TIMP3",#MEP
"SLC40A1","PKIG","APOC1","BLVRB", "CA1", "HBB", "KLF1", "TFR2", #ERP
"TCL1A","IRF4","SCT","CCR2", "IRF8", "MPEG1", "SPIB", "IGKC", #MDP1
"SAMHD1","LGALS1","TMEM154", #MDP2
"RGS1", "NPTX2", "DDIT4", "ID2","SH2D1A",#CLP_Ly2
"DNTT","VPREB3","CYGB","ADA","LTB","HMBG1","LEF1","MS4A1","RAG1", "RAG2" ,"CD79A", "VPREB1","TGFB1", #CLP_Ly2_B_PC/CLP_Ly2
"GZMA","ACY3","NEGR1",#LMPP_Ly1
"CLC", "CPA3", "HDC") #Eo_B_Mast

#cstem<-list(c("LINC00511","CNTN5","KALRN","PDE4D","SLC26A7","CCSER1","ANKH","GALNT13","SNTB1","NCALD","DMD","SOX6","MAML2","FGFR2","POLR2F","UNC5C","SLC5A8","MEIS2","AGMO","TGFBR3","TRIM2","PDZRN3","PHLDA1","GLS","FBXW7","MAP2","RIMS2","PODXL","DOCK4","CNTN1","NKAIN3","SEMA3D","DSC2","ZFPM2","PHLPP1","ST3GAL1","C10orf90","CEP85L","ITPR2","NCMAP","PREP","ITGA6","MELTF","PDE4B","ROCR","SHROOM3","NIPAL3","ESRRG","PADI2","SLC26A2","PYGB","KCND2","LAMA4","FAM189A2","CEMIP2","TMEM108","ST3GAL6","DNAH11","LSAMP","PTPRG","NPAS2","ZNF462","KHDRBS3","CST3","ZNF521","STOX2","AC078923.1","DCLK2","MTSS1","ATL2","ARHGEF28","MUC7","PRICKLE2","EHBP1","SFMBT2","NFIB","CEACAM1","UNC80","CCDC14","KANK4","ADCY2","SFRP1","MID1","UGP2","PLCB4","AC002066.1","C4orf19","LRRC7","PDGFC","ART3","PPARGC1A","AC024230.1","GNE","GAB1","HIBCH","CMTM7","UBE2E2","PIK3C2G","C5orf17","GABRP","GK5","MBNL1","PPP2R2B","TOX","FGF13","HNF4G","ERC2","DOCK7","ITGBL1","CDK6","FNDC3B","MGLL","LMO4","TGFB2","CELF2","PTAR1","BACH2","BACE2","CNNM2","ZBTB20","IL12RB2","HOMER2","KIF1B","AMD1","SNX1","LMCD1","LMCD1-AS1","DAPK1","NALCN","DEPTOR","RRBP1","CARMIL1","AQP5","KCNK2","SLC35F3","SCHLAP1","ADAM12","FMNL2","KIF26B","CLIC4","PABPC1","ARHGEF9","JAZF1","WWTR1","WWP2","EPHA7","PAG1","ATP6V1C2","BTBD3","TMEM87A","ABCA13","LINC01122","SNCA","ARFGAP3","AC009478.1","CDK19","GXYLT2","AL117329.1","OGFRL1","SPATA6","PPP2R3A","EXT1","AL592295.3","APBA1","FOXN2","MT-ND2","GPR156","CNOT4","MBNL2","PRICKLE1","CREB3L2","IFRD1","GPM6B","FAM13A","MEF2A","MT-ND3","RAB7A","MTMR12","FAT1","ARRDC3","CAV1","MGAM2","CLCN5","BCL11A","ITGB4","SLC24A3","CPNE4","CD44","DSCAM","DIAPH2","AC025470.2","EPHB1","SH3YL1","FGF14","OSBPL8","CEP350","ATXN7L1","USO1","MKLN1","ITIH6","ADAMTS17","GNAI1","SEC31A","LRP6","MRAS","PLCH1","TPTEP2-CSNK1E","GRHL1","CRISPLD1","PPM1L","KRT15","RANBP17","KAT2B","ANKRD44","RNF13","DYNC1I1","DIP2A","ARL17B","RNF152","AHCYL2","CRADD","RNF145","TENM2","MAP4K4","NR3C2","RNF217","GRIP1","ADCY8","PTPRJ","C1GALT1","ANKRD17","ARHGEF38","TNRC6C","SH3PXD2B","TNRC6B","LINC02473","ROPN1B","PTPRK","UTP20","PLCE1","UGT8","CCDC50","NR2F2-AS1","NRG1","FRMD3","RNPC3","DDX17","SFT2D2","ITGA9","NRG2","GARNL3","GRIA2","ARFIP1","ATR","PCDH7","SPTSSB","FNBP1L","MFGE8","ELF5","SOX10","HNRNPD","DTNB","DAPK2","ELMO1","CTTNBP2","PIK3CB","LDLRAD4","KANSL1L","UNC13B","CNGA1","11-Sep","SLC39A8","CTNNB1","TRIQK","USP54","PDE1C","SYN2","AC136616.1","P4HB","PRMT2","VPS13D","TCF20","FAM160A1","CPNE8","VEGFA","EVA1A","PDSS2","CYP7B1","DGKI","MCTP2","OXR1","CLSTN1","NUP205","LINC00964","AC083843.3","EP300","PPP3CA","MTUS1","LAMA1","SMARCA2","PCAT1","ACSL1","PPARA","SAT1","CLINT1","SNX22","PCDH11X","MAP7D3","ETV5"))
#combined <- AddModuleScore(object = combined,assay="SCT",features = cstem,name = 'cstem_score')

pdf("cellmarkers.pdf",height=5,width=13,useDingbats=FALSE)
a<-DotPlot(sobj,assay="SCT",features=DEGlist,group.by=c("cell_type_specific"), cols = "RdGy",scale.max = 100, scale.min = 0)+theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5))&guides(x =  guide_axis(angle = 90))&theme_light()
print(a)
dev.off()


mycolors <- colorRampPalette((brewer.pal(9, "YlOrRd")))(1024)

a<-DotPlot(HSConly,assay="SCT",features=c("TNFA_NFKB1"),group.by=c("cell_type_specific"),split.by=c("treatment2"),cols=c("#fcfdbf","#fe9f6d","#de4968","#8c2981","#3b0f70","#000004"))+theme(legend.position="right")
b<-DotPlot(HSConly,assay="SCT",features=c("EGF_Score1"),group.by=c("cell_type_specific"),split.by=c("treatment2"),cols=c("#fcfdbf","#fe9f6d","#de4968","#8c2981","#3b0f70","#000004"))+theme(legend.position="right")
c<-DotPlot(HSConly,assay="SCT",features=c("HBOScore1"),group.by=c("cell_type_specific"),split.by=c("treatment2"),cols=c("#fcfdbf","#fe9f6d","#de4968","#8c2981","#3b0f70","#000004"))+theme(legend.position="right")


a_plot<- a$data %>% separate(id, c("celltype", "group"), "_", extra = "merge") %>% 
  mutate(Group2 = substr(group, nchar(group)-3+1, nchar(group)))%>%
      ggplot(aes(x = celltype, y = group,fill=avg.exp,size=pct.exp,label=sprintf("%0.2f", round(avg.exp, digits = 2)))) + 
      geom_point( shape=21) +
      facet_grid(Group2~., scales = 'free', space = 'free') +
      #coord_equal(ratio = 10/8) +
      #Creating color range
      scale_fill_gradientn(colors=mycolors, guide="colorbar", na.value = 'white' ) +
      scale_size(range = c(0,10))+
      geom_text(size=3)+
      #Rotating labels
      theme_classic()+
      #ggtitle(prot) +
      theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
            panel.grid.major = element_line(color = 'grey80', size = 0.1),
            #axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text.y = element_text(size = 10, colour = 'grey20'),
            axis.title = element_blank(),  
            plot.margin=unit(c(0,0,0,0), "cm"))+ggtitle("TNFA_NFKB Score")

      b_plot<- b$data %>% separate(id, c("celltype", "group"), "_", extra = "merge") %>% 
  mutate(Group2 = substr(group, nchar(group)-3+1, nchar(group)))%>%
      ggplot(aes(x = celltype, y = group,fill=avg.exp,size=pct.exp,label=sprintf("%0.2f", round(avg.exp, digits = 2)))) + 
      geom_point( shape=21) +
      facet_grid(Group2~., scales = 'free', space = 'free') +
      #coord_equal(ratio = 10/8) +
      #Creating color range
      scale_fill_gradientn(colors=mycolors, guide="colorbar", na.value = 'white' ) +
      scale_size(range = c(0,10))+
      geom_text(size=3)+
      #Rotating labels
      theme_classic()+
      #ggtitle(prot) +
      theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
            panel.grid.major = element_line(color = 'grey80', size = 0.1),
            #axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text.y = element_text(size = 10, colour = 'grey20'),
            axis.title = element_blank(),  
            plot.margin=unit(c(0,0,0,0), "cm"))+ggtitle("EGF Score")

c_plot<- c$data %>% separate(id, c("celltype", "group"), "_", extra = "merge") %>% 
  mutate(Group2 = substr(group, nchar(group)-3+1, nchar(group)))%>%
      ggplot(aes(x = celltype, y = group,fill=avg.exp,size=pct.exp,label=sprintf("%0.2f", round(avg.exp, digits = 2)))) + 
      geom_point( shape=21) +
      facet_grid(Group2~., scales = 'free', space = 'free') +
      #coord_equal(ratio = 10/8) +
      #Creating color range
      scale_fill_gradientn(colors=mycolors, guide="colorbar", na.value = 'white' ) +
      scale_size(range = c(0,10))+
      geom_text(size=3)+
      #Rotating labels
      theme_classic()+
      #ggtitle(prot) +
      theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
            panel.grid.major = element_line(color = 'grey80', size = 0.1),
            #axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text.y = element_text(size = 10, colour = 'grey20'),
            axis.title = element_blank(),  
            plot.margin=unit(c(0,0,0,0), "cm"))+ggtitle("HBO Score")


pdf(paste("Gene_Scores_Version2.pdf", sep=""), width=15, height=4)
ggarrange(a_plot, b_plot,c_plot,
          ncol = 3, nrow = 1)
dev.off()

####################################
###########UPLOAD DATA to GEO#######
####################################

###OUTPUT RNA Counts for publication
library(Seurat)
library(Matrix)

fullobj<-readRDS("~rjayasin/mmy_scratch/Dipersio/genesis/analysis/all_integrated_merge_012022.rds")

pdf(paste("check.pdf", sep=""), width=10, height=10)
  b<-DimPlot(fullobj,group.by=c("cell_type_specific"),cols=c("#73d055ff","#66cccc", "#1f968bff","#6633ff","#cc99ff","#996699","#660033",
 "#802582FF", "#952C80FF", "#AB337CFF", "#C03A76FF", "#D6456CFF",
"#E85362FF", "#F4685CFF", "#FA815FFF", "#FD9A6AFF", "#FEB37BFF", "#FECC8FFF",
"#FDE4A6FF", "#FCFDBFFF"),raster=TRUE,label=TRUE)&coord_equal()
  print(b)
dev.off()

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
ifnb.list <- SplitObject(fullobj, split.by = "orig.ident")

# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  print(x)
  sample=x@meta.data$orig.ident[1]
  dir.create(paste0(sample))
  writeMM(x@assays$RNA@counts, file = paste0(sample,'/matrix.mtx'))
  # save gene and cell names
  #Convert gene names to upper case if mouse
  write(x = toupper(rownames(x@assays$RNA@counts)), file = paste0(sample,"/features.tsv"))
  write(x = colnames(x@assays$RNA@counts), file = paste0(sample,"/barcodes.tsv"))
  ###Extract METADATA
  x@meta.data$barcodes = rownames(x@meta.data)
  df = x@meta.data[, c('barcodes', 'cell_type_specific')]
  write.table(df, file =paste0(sample,'/',sample,'_meta_data.tsv'), sep = '\t', quote = F, row.names = F)
})



