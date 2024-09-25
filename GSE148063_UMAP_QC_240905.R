# Documents > mouse_retina > GSE148063

library(dplyr)
library(Seurat)
set.seed(0)

whole.data <- ReadMtx(mtx = "GSE148063_matrix.mtx.gz",
                      features = "GSE148063_genes.tsv.gz",
                      cells = "GSE148063_barcodes.tsv.gz")

whole <- CreateSeuratObject(counts = whole.data, project = "GSE148063",
                            min.cells = 3, min.features = 200)

whole[["percent.mt"]] <- PercentageFeatureSet(whole, pattern = "^mt-")

VlnPlot(whole, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)

whole <- subset(whole, subset = nFeature_RNA > 2500 &
                  nFeature_RNA < 7500 & percent.mt < 5)

whole <- NormalizeData(whole, normalization.method = "LogNormalize",
                       scale.factor = 10000)

whole <- FindVariableFeatures(whole, selection.method = "vst",
                              nfeatures = 2000)

all.genes <- rownames(whole)
whole <- ScaleData(whole, features = all.genes)

whole <- RunPCA(whole, features = VariableFeatures(object = whole))
DimPlot(whole, reduction = "pca")

whole <- JackStraw(whole, num.replicate = 100)
whole <- ScoreJackStraw(whole, dims = 1:20)
ElbowPlot(whole)

whole <- FindNeighbors(whole, dims = 1:15)
whole <- FindClusters(whole, resolution = 0.7)
whole <- RunUMAP(whole, dims = 1:7)
DimPlot(whole, reduction = "umap", label = TRUE)

stage <- read.csv("stage.csv", header = TRUE)
stage <- stage[,8]
CellsMeta <- whole@meta.data
CellsMeta["stage"] <- stage
whole <- AddMetaData(whole, CellsMeta)

Notch <- c("Notch1", "Notch2", "Notch3", "Notch4")
DotPlot(object = whole, features = Notch)

ForDis <- c("Hes1", "Hes3", "Hes5", "Hes6", "Hes7", "Hey1", "Hey2", "Heyl")
DotPlot(object = whole, features = ForDis)

markers <- c("Pax6", "Six3", "Rax", "Lhx2",
             "Twist1", "Col3a1", "Alx3", "Prrx2", "Acta2",
             "Krt18", "Krt8", "Dlx5",
             "Hbb-bs", "Pecam1", "Cd34", "Lyz2", "Ncf1")
DotPlot(object = whole, features = markers)

new.cluster.ids <- c("OV1", "OV2", "OV3", "MES1","OV4", 
                     "MES2", "MES3", "OV5", "MES4",
                     "ECTO", "VAS", "OV6", "MES5", "OV7",
                     "ERY", "MES6", "LEU", "MES7")

names(new.cluster.ids) <- levels(whole)
whole <- RenameIdents(whole, new.cluster.ids)
DimPlot(whole, reduction = "umap", label = TRUE)
DotPlot(object = whole, features = markers)

cluster.markers <- FindMarkers(whole, ident.1 = "LEU")
head(cluster.markers, n = 15)

DimPlot(whole, split.by = "stage")

saveRDS(whole, file = "GSE148063_UMAP_QC.rds")