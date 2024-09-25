# human_retina_organoids > GSE235585

library(dplyr)
library(Seurat)
library(ggplot2)

whole <- readRDS("GSE235585_UMAP_QC_240918.rds")

integrated_markers <- c("LHX2", "PAX6", "RAX", "SIX3",   
                        "ELAVL3", "RCVRN",
                        "VSX2", "MAP2", "SOX2", "NES", "MKI67",
                        "TOP2A", "RLBP1", "PAX2",
                        "GFAP", "CRYBB1", "SOX1", "FOXG1", "ELAVL4", "DLX6",
                        "GDF7", "LMX1A", "HES6", "ATOH7", "POU4F2", "ISL1",
                        "PRDM13", "TFAP2A", "PROX1", "ONECUT1", "CRX",
                        "OTX2", "ARR3", "OPN1SW", "NRL", "VSX1", "GRM6", "GRIK1")
DotPlot(object = whole, features = integrated_markers) + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20)) +
  RotatedAxis() +
  scale_y_discrete(
    limits = c("Late RPC", "Early RPC", "Transient NRPC", "NRPC",
               "RGC1", "RGC2", "CM1", "CM2", "MULLER", "INTER", "PHOTO Pre",
               "ROD", "CONE", "BIPOLAR Pre")
  )

new.cluster.ids <- c("CM1", "Late RPC", "RGC1", "ROD", "MULLER",
                     "CM2", "BIPOLAR Pre", "Early RPC", "INTER", "NRPC",
                     "CONE", "Transient NRPC", "PHOTO Pre", "RGC2")
names(new.cluster.ids) <- levels(whole)
whole <- RenameIdents(whole, new.cluster.ids)

DimPlot(whole, reduction = "umap", label = TRUE, label.size = 10, repel = TRUE) + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20))

RPCs <- subset(whole, idents = c("Late RPC", "Early RPC"))

DimPlot(whole, split.by = "orig.ident", label = TRUE, label.size = 10, repel = TRUE) + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20))

notch <- c("NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4")
DotPlot(object = whole, features = notch) + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20)) +
  scale_y_discrete(
    limits = c("Late RPC", "Early RPC", "Transient NRPC", "NRPC",
               "RGC1", "RGC2", "CM1", "CM2", "MULLER", "INTER", "PHOTO Pre",
               "ROD", "CONE", "BIPOLAR Pre")
  )

target <- c("HES1", "HES3", "HES5", "HES6", "HES7", "HEY1", "HEY2", "HEYL")
DotPlot(object = whole, features = target) + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20)) +
  RotatedAxis() +
  scale_y_discrete(
    limits = c("Late RPC", "Early RPC", "Transient NRPC", "NRPC",
               "RGC1", "RGC2", "CM1", "CM2", "MULLER", "INTER", "PHOTO Pre",
               "ROD", "CONE", "BIPOLAR Pre")
  )

library(CSCORE)
stage = RPCs
mean_exp = rowMeans(stage@assays$RNA$counts/stage$nCount_RNA)
genes_selected = names(sort.int(mean_exp, decreasing = T))[1:5000]
CSCORE_result <- CSCORE(stage, genes = genes_selected)

CSCORE_coexp <- CSCORE_result$est
CSCORE_p <- CSCORE_result$p_value
p_matrix_BH = matrix(0, length(genes_selected), length(genes_selected))
p_matrix_BH[upper.tri(p_matrix_BH)] = p.adjust(CSCORE_p[upper.tri(CSCORE_p)], method = "BH")
p_matrix_BH <- p_matrix_BH + t(p_matrix_BH)
CSCORE_coexp[p_matrix_BH > 0.05] <- 0

library(WGCNA)
adj = WGCNA::adjacency.fromSimilarity(abs(CSCORE_coexp), power = 1)
TOM = WGCNA::TOMsimilarity(adj)
dissTOM = 1-TOM
rownames(dissTOM) <- colnames(dissTOM) <- genes_selected
hclust_dist = hclust(as.dist(dissTOM), method = "average")
memb = dynamicTreeCut::cutreeDynamic(dendro = hclust_dist,
                                     distM = dissTOM,
                                     deepSplit = 2,
                                     pamRespectsDendro = FALSE,
                                     minClusterSize = 10)
names(memb) = genes_selected
memb_tab <- table(memb)
module_list = lapply(sort(unique(memb)), function(i_k) names(which(memb == i_k)))

list_to_write <- sapply(module_list, function(x) paste(unlist(x), collapse = ", "))
writeLines(list_to_write, "module_list_GSE235585_RPCs.txt")

library(SeuratWrappers)
library(monocle3)

whole <- JoinLayers(whole)
cds <- as.cell_data_set(whole)
cds <- cluster_cells(cds, resolution=1e-3)

p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
plot(p1)
plot(p2)

integrated.sub <- subset(as.Seurat(cds, assay = NULL))
cds <- as.cell_data_set(integrated.sub)

cds <- learn_graph(cds, use_partition = FALSE, verbose = FALSE)

plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

#DimPlot(whole, split.by = "stage")
#somite_12 <- subset(whole, subset = stage == "somite_12")
#DimPlot(somite_12)

#FeaturePlot(integrated.sub, features = "PAX6")
#FeaturePlot(integrated.sub, features = "MKI67")

cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 21]))
plot_cells(cds,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "grey60")

integrated.sub <- as.Seurat(cds, assay = NULL)
FeaturePlot(integrated.sub, "monocle3_pseudotime")

integrated.sub <- SetIdent(integrated.sub, 
                           value = integrated.sub@meta.data$seurat_clusters)

df2 <- FetchData(object = integrated.sub, vars = c("monocle3_pseudotime"), layer = "count")
df3 <- FetchData(object = whole, vars = c("NOTCH1", "NOTCH3", "HES1", "HES5", "HES6", "JAG1",
                                         "JAG2", "DLL1", "DLL3", "DLL4", "ident"), layer = "scale.data")
df4 <- cbind(df2, df3)
df5 <- subset(df4, df4$ident == "Late RPC" | df4$ident == "Early RPC")

# ggplot(df4) +
#   geom_point(aes(monocle3_pseudotime, NOTCH1, color = "NOTCH1")) +
#   geom_point(aes(monocle3_pseudotime, NOTCH3, color = "NOTCH3")) +
#   scale_color_manual(values = c("NOTCH1" = "red", "NOTCH3" = "blue")) +
#   guides(color = guide_legend(title = NULL)) +
#   FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
#   theme(legend.text = element_text(size = 20), axis.title.y = element_blank(),
#         panel.background = element_rect(fill = "transparent", color = NA),
#         panel.grid.minor = element_line(color = "grey"))

ggplot(df5) +
  geom_point(aes(monocle3_pseudotime, NOTCH1, color = "NOTCH1")) +
  geom_point(aes(monocle3_pseudotime, NOTCH3, color = "NOTCH3")) +
  scale_color_manual(values = c("NOTCH1" = "red", "NOTCH3" = "blue")) +
  guides(color = guide_legend(title = NULL)) +
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20), axis.title.y = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.minor = element_line(color = "grey"))