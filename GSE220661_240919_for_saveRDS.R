# GSE220661

library(dplyr)
library(Seurat)
library(patchwork)

D028.count <- ReadMtx(mtx = "GSM6808628_D28_matrix.mtx.gz",
                     cells = "GSM6808628_D28_barcodes.tsv.gz",
                     features = "GSM6808628_D28_genes.tsv.gz")

D028 <- CreateSeuratObject(counts = D028.count, project = "D028",
                          min.cells = 3, min.features = 200)

D078.count <- ReadMtx(mtx = "GSM6808629_D78_matrix.mtx.gz",
                     cells = "GSM6808629_D78_barcodes.tsv.gz",
                     features = "GSM6808629_D78_genes.tsv.gz")

D078 <- CreateSeuratObject(counts = D078.count, project = "D078",
                          min.cells = 3, min.features = 200)

D185.count <- ReadMtx(mtx = "GSM6808630_D185_matrix.mtx.gz",
                     cells = "GSM6808630_D185_barcodes.tsv.gz",
                     features = "GSM6808630_D185_genes.tsv.gz")

D185 <- CreateSeuratObject(counts = D185.count, project = "D185",
                          min.cells = 3, min.features = 200)

organoids <- merge(D028, y = c(D078, D185),
                   add.cell.ids = c("D028", "D078", "D185"),
                   project = "oranoids")

organoids[["percent.mt"]] <- PercentageFeatureSet(organoids, pattern = "^MT-")

VlnPlot(organoids, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

organoids <- subset(organoids, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15)

organoids <- NormalizeData(organoids, normalization.method = "LogNormalize", scale.factor = 10000)

organoids <- FindVariableFeatures(organoids, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(organoids)
organoids <- ScaleData(organoids, features = all.genes)

organoids <- RunPCA(organoids, features = VariableFeatures(object = organoids))

DimPlot(organoids, reduction = "pca") 

ElbowPlot(organoids)

organoids <- FindNeighbors(organoids, dims = 1:15)
organoids <- FindClusters(organoids, resolution = 0.7)

organoids <- RunUMAP(organoids, dims = 1:10)
DimPlot(organoids, reduction = "umap", label = TRUE)

DimPlot(organoids, split.by = "orig.ident", label = TRUE)

markers <- c("LHX2", "PAX6", "RAX", "SIX3",   
             "POU4F2", "ELAVL3", "OTX2", "RCVRN",
             "VSX2", "MAP2")
DotPlot(object = organoids, features = markers)

CellTypeMarkers <- c("RAX", "VSX2", "MKI67", "TOP2A", "RLBP1", "PAX2",
                    "GFAP", "CRYBB1", "SOX1", "FOXG1", "ELAVL4", "DLX6",
                    "GDF7", "LMX1A", "HES6", "ATOH7", "POU4F2", "ISL1",
                    "PRDM13", "TFAP2A", "PROX1", "ONECUT1", "CRX",
                    "OTX2", "ARR3", "OPN1SW", "NRL", "VSX1", "GRM6", "GRIK1")
DotPlot(object = organoids, features = CellTypeMarkers)

notch <- c("NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4")
DotPlot(object = organoids, features = notch)

target <- c("HES1", "HES3", "HES5", "HES6", "HES7", "HEY1", "HEY2", "HEYL")
DotPlot(object = organoids, features = target)

new.cluster.ids <- c("RPC1", "Npre", "RPC2", "CONE", "RPC3",
                     "RPC4", "ROD", "RPC5", "LENS", "RPC6", "MULLER",
                     "ECTO1", "RPC7", "BIPOLAR", "HOR", "ECTO2")

names(new.cluster.ids) <- levels(organoids)
organoids <- RenameIdents(organoids, new.cluster.ids)
DimPlot(organoids, reduction = "umap", label = TRUE)
DimPlot(organoids, split.by = "orig.ident", label = TRUE)

#organoids_temp <- JoinLayers(organoids_temp)
#cluster.markers <- FindMarkers(organoids_temp, ident.1 = "ASTR")
#head(cluster.markers, n = 50)

#DotPlot(object = organoids, features = markers)
#DotPlot(object = organoids, features = CellTypeMarkers)
#DotPlot(object = organoids, features = notch)
#DotPlot(object = organoids, features = target)

#library(dorothea)
#library(dplyr)

#dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))
#write.csv(dorothea_regulon_human, "regulon_human.csv")

# dorothea
# line49320 ESR1 A NOTCH2
# line124906 KLF5 C NOTCH1
# line147328 MYC A NOTCH4
# line220080 SNAI2 C NOTCH1
# line287736 ZBTB7A C NOTCH1
# line333389 ZNF263 B NOTCH1

#dorothea <- c("ESR1", "KLF5", "MYC", "SNAI2", "ZBTB7A", "ZNF263")
#DotPlot(object = organoids, features = dorothea)

saveRDS(organoids, file = "GSE220661_UMAP_QC_240919.rds")