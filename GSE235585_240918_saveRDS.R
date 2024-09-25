# GSE235585

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

D020 <- read.csv("GSM7505844_D20_counts.csv.gz", row.names = 1)
D020 <- CreateSeuratObject(counts = D020, project = "D020",
                          min.cells = 3, min.features = 1000)

D033 <- read.csv("GSM7505846_D33_counts.csv.gz", row.names = 1)
D033 <- CreateSeuratObject(counts = D033, project = "D033",
                                    min.cells = 3, min.features = 1000)

D090 <- read.csv("GSM7505849_D90_counts.csv.gz", row.names = 1)
D090 <- CreateSeuratObject(counts = D090, project = "D090",
                          min.cells = 3, min.features = 1000)

D210 <- read.csv("GSM7505845_D210_counts.csv.gz", row.names = 1)
D210 <- CreateSeuratObject(counts = D210, project = "D210",
                          min.cells = 3, min.features = 1000)

organoids <- merge(D020, y = c(D033, D090, D210),
                   add.cell.ids = c("D020", "D033", "D090", "D210"),
                   project = "oranoids")

organoids[["percent.mt"]] <- PercentageFeatureSet(organoids, pattern = "^MT-")

VlnPlot(organoids, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

organoids <- NormalizeData(organoids, normalization.method = "LogNormalize", scale.factor = 10000)

organoids <- FindVariableFeatures(organoids, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(organoids)
organoids <- ScaleData(organoids, features = all.genes)

organoids <- RunPCA(organoids, features = VariableFeatures(object = organoids))

DimPlot(organoids, reduction = "pca") 

ElbowPlot(organoids)

organoids <- FindNeighbors(organoids, dims = 1:15)
organoids <- FindClusters(organoids, resolution = 0.4)

organoids <- RunUMAP(organoids, dims = 1:10)
DimPlot(organoids, reduction = "umap")

DimPlot(organoids, split.by = "orig.ident")

markers <- c("PAX6", "SIX3", "RAX", "LHX2",
             "TWIST1", "COL3A1", "ALX3", "PRRX2", "ACTA2",
             "KRT18", "KRT8", "DLX5",
             "HBB-BS", "PECAM1", "CD34", "LYZ2", "NCF1")
DotPlot(object = organoids, features = markers)

CellTypeMarkers <- c("RAX", "VSX2", "MKI67", "TOP2A", "RLBP1", "PAX2",
                     "GFAP", "CRYBB1", "SOX1", "FOXG1", "ELAVL4", "DLX6",
                     "GDF7", "LMX1A", "HES6", "ATOH7", "POU4F2", "ISL1",
                     "PRDM13", "TFAP2A", "PROX1", "ONECUT1", "CRX",
                     "OTX2", "ARR3", "OPN1SW", "NRL", "VSX1", "GRM6", "GRIK1",
                     "MAP2")
DotPlot(object = organoids, features = CellTypeMarkers)

notch <- c("NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4")
DotPlot(object = organoids, features = notch)

target <- c("HES1", "HES3", "HES5", "HES6", "HES7", "HEY1", "HEY2", "HEYL")
DotPlot(object = organoids, features = target)

integrated_markers <- c("LHX2", "PAX6", "RAX", "SIX3",   
                        "ELAVL3", "RCVRN",
                        "VSX2", "MAP2", "SOX2", "NES", "MKI67",
                        "TOP2A", "RLBP1", "PAX2",
                        "GFAP", "CRYBB1", "SOX1", "FOXG1", "ELAVL4", "DLX6",
                        "GDF7", "LMX1A", "HES6", "ATOH7", "POU4F2", "ISL1",
                        "PRDM13", "TFAP2A", "PROX1", "ONECUT1", "CRX",
                        "OTX2", "ARR3", "OPN1SW", "NRL", "VSX1", "GRM6", "GRIK1")
DotPlot(object = organoids, features = integrated_markers) + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20)) +
  RotatedAxis()

dorothea <- c("ESR1", "KLF5", "MYC", "SNAI2", "ZBTB7A", "ZNF263")
DotPlot(object = organoids, features = dorothea)

new.cluster.ids <- c("CM1", "Late RPC", "RGC1", "ROD", "MULLER",
                     "CM2", "BIPOLAR Pre", "Early RPC", "INTER", "NRPC",
                     "CONE", "Transient NRPC", "PHOTO Pre", "RGC2")

names(new.cluster.ids) <- levels(organoids)
organoids <- RenameIdents(organoids, new.cluster.ids)
DimPlot(organoids, reduction = "umap", label = TRUE)
DimPlot(organoids, split.by = "orig.ident", label = TRUE)

DotPlot(object = organoids, features = markers)
DotPlot(object = organoids, features = CellTypeMarkers)
DotPlot(object = organoids, features = notch)
DotPlot(object = organoids, features = target)

saveRDS(organoids, file = "GSE235585_UMAP_QC_240918.rds")