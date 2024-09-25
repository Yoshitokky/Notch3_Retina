# mouse_retina > GSE148063

library(dplyr)
library(Seurat)
library(ggrepel)

set.seed(0)
whole <- readRDS("GSE148063_UMAP_QC.rds")

markers <- c("Pax6", "Six3", "Rax", "Lhx2", "Mki67",
             "Twist1", "Col3a1", "Alx3", "Prrx2", "Acta2",
             "Krt18", "Krt8", "Dlx5",
             "Hbb-bs", "Pecam1", "Cd34", "Lyz2", "Ncf1")

DotPlot(object = whole, features = markers) + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20)) +
  RotatedAxis() +
  scale_y_discrete(
    limits = c("OV1", "OV2", "OV3", "OV4", "OV5", "OV6", "OV7",
               "MES1", "MES2", "MES3", "MES4", "MES5", "MES6", "MES7",
               "ECTO", "ERY", "VAS", "LEU")
  )

integrated_markers <- c("Lhx2", "Pax6", "Rax", "Six3",   
                        "Elavl3", "Rcvrn",
                        "Vsx2", "Map2", "Sox2", "Nes", "Mki67",
                        "Top2a", "Rlbp1", "Pax2",
                        "Gfap", "Crybb1", "Sox1", "Foxg1", "Elavl4", "Dlx6",
                        "Gdf7", "Lmx1a", "Hes6", "Atoh7", "Pou4f2", "Isl1",
                        "Prdm13", "Tfap2a", "Prox1", "Onecut1", "Crx",
                        "Otx2", "Arr3", "Opn1sw", "Nrl", "Vsx1", "Grm6", "Grik1")
DotPlot(object = whole, features = integrated_markers) + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20)) +
  RotatedAxis()

new.cluster.ids <- c("OV1", "OV2", "OV3", "MES1","OV4", 
                     "MES2", "MES3", "OV5", "MES4",
                     "ECTO", "VAS", "OV6", "MES5", "OV7",
                     "ERY", "MES6", "LEU", "MES7")

names(new.cluster.ids) <- levels(whole)
whole <- RenameIdents(whole, new.cluster.ids)
DimPlot(whole, reduction = "umap", label = TRUE, label.size = 10, repel = TRUE) + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20))
OV <- subset(whole, idents = c("OV1", "OV2", "OV3", "OV4",
                               "OV5", "OV6", "OV7"))

library(CSCORE)
stage = OV
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
writeLines(list_to_write, "module_list_GSE148063_OV.txt")

library(SeuratWrappers)
library(monocle3)

cds <- as.cell_data_set(OV)
cds <- cluster_cells(cds, resolution=1e-3)

p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
plot(p1, label = TRUE)
plot(p2)

integrated.sub <- subset(as.Seurat(cds, assay = NULL), monocle3_partitions == 1)
cds <- as.cell_data_set(integrated.sub)

cds <- learn_graph(cds, use_partition = TRUE, verbose = FALSE)

plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

DimPlot(whole, split.by = "stage", label = TRUE, label.size = 6, repel = TRUE) + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20))

cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 2]))

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
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20))

#head(integrated.sub[[]])
#head(whole[[]])
#colData(cds)
#levels(integrated.sub)

integrated.sub <- SetIdent(integrated.sub, 
                           value = integrated.sub@meta.data$seurat_clusters)

notch <- c("Notch1", "Notch2", "Notch3", "Notch4")
target <- c("Hes1", "Hes2", "Hes3", "Hes5", "Hes6", "Hes7",
           "Hey1", "Hey2", "Heyl")
delta <- c("Dll1", "Dll3", "Dll4", "Jag1", "Jag2")

new.cluster.ids_2 <- c("OV1", "OV2", "OV3", "OV4", 
                       "OV5", "OV6", "OV7")

names(new.cluster.ids_2) <- levels(integrated.sub)
integrated.sub <- RenameIdents(integrated.sub, new.cluster.ids_2)

DotPlot(integrated.sub, features = notch) + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20))
DotPlot(integrated.sub, features = target) + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20))

df2 <- FetchData(object = integrated.sub, vars = c("monocle3_pseudotime"), layer = "count")
df3 <- FetchData(object = OV, vars = c("Notch1", "Notch3", "Hes1", "Hes5", "Hes6", "Zbtb7a", "Zfp263"), layer = "scale.data")
df4 <- cbind(df2, df3)

library(ggplot2)
# ggplot(df4, aes(monocle3_pseudotime, Notch1), ) + geom_point(color = "red") +
#   guides(color = guide_legend(title = NULL)) +
#   FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
#   theme(legend.text = element_text(size = 20), axis.title.y = element_blank())
# ggplot(df4, aes(monocle3_pseudotime, Notch3)) + geom_point()
# ggplot(df4, aes(monocle3_pseudotime, Zbtb7a)) + geom_point()
# ggplot(df4, aes(monocle3_pseudotime, Zfp263)) + geom_point()
# ggplot(df4, aes(monocle3_pseudotime, Hes6)) + geom_point()
ggplot(df4) +
  geom_point(aes(monocle3_pseudotime, Notch1, color = "Notch1")) +
  geom_point(aes(monocle3_pseudotime, Notch3, color = "Notch3")) +
  scale_color_manual(values = c("Notch1" = "red", "Notch3" = "blue")) +
  guides(color = guide_legend(title = NULL)) +
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20), axis.title.y = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.minor = element_line(color = "grey"))

# correlation <- cor(df4$Notch1, df4$Notch3) # 0.117
# ggplot(df4, aes(Notch1, Notch3)) + geom_point()
# correlation <- cor(df4$Notch1, df4$Zbtb7a) # 0.045
# correlation <- cor(df4$Notch1, df4$Zfp263) # 0.015
# 
# markers <- c("Pax6", "Six3", "Rax", "Lhx2",
#              "Twist1", "Col3a1", "Alx3", "Prrx2", "Acta2",
#              "Krt18", "Krt8", "Dlx5",
#              "Hbb-bs", "Pecam1", "Cd34", "Lyz2", "Ncf1")
# DotPlot(object = whole, features = markers)
# 
# # dorothea
# # line 2734 Esr1 A Notch2
# # line 4997 Klf5 C Notch1
# # line 9287 Snai2 C Notch1
# # line 12143 Zbtb7a C Notch1
# # line 12509 Zfp263 B Notch1
# dorothea <- c("Esr1", "Klf5", "Snai2", "Zbtb7a", "Zfp263")
# DotPlot(integrated.sub, dorothea)
# 
# df <- FetchData(object = integrated.sub, vars = c("Notch1", "Zbtb7a", "Zfp263", "monocle3_pseudotime"), layer = "counts")
# df_Notch1 <- subset(df, df$Notch1 > 0)
# df_Notch1_Zbtb7a_Zfp263 <- subset(df_Notch1, df_Notch1$Zbtb7a > 0 | Zfp263 > 0)
# df_Notch1_Zbtb7a_Zfp263_0 <- subset(df_Notch1, df_Notch1$Zbtb7a == 0 & Zfp263 == 0)
# df_Notch1_0 <- subset(df, df$Notch1 == 0)
# df_Notch1_0_Zbtb7a_Zfp263_0 <- subset(df_Notch1_0, df_Notch1_0$Zbtb7a == 0 & Zfp263 == 0)
# df_Notch1_0_Zbtb7a_Zfp263 <- subset(df_Notch1_0, df_Notch1_0$Zbtb7a > 0 | Zfp263 > 0)
# df_Zbtb7a_Zfp263 <- subset(df, df$Zbtb7a > 0 & Zfp263 > 0)
# df_Zbtb7a <- subset(df, df$Zbtb7a > 0)
# 
# lengths(df_Notch1_0_Zbtb7a_Zfp263_0) # 1622 / 4042 = 0.401
# lengths(df_Notch1_0_Zbtb7a_Zfp263) # 1079 / 4042 = 0.266, 1079 / 2420 = 0.445
# lengths(df_Notch1_Zbtb7a_Zfp263_0) # 651 / 4042 = 0.161, 651 / 2420 = 0.269
# lengths(df_Notch1_Zbtb7a_Zfp263) # 690 / 4042 = 0.170, 690 / 2420 = 0.285
# vx = matrix(c(1662, 1079, 651, 690), nrow = 2, byrow = TRUE)
# chisq.test(vx) # X-squared = 53.099, df = 1, p-value = 3.172e-13
# 
# lengths(df_Zbtb7a_Zfp263)
# lengths(df_Zbtb7a)
# 
# library(ggVennDiagram)
# x <- list(Notch1 = 651:1341, TF = 1079:1769)
# ggVennDiagram(x)
# 
# library(beeswarm)
# 
# jointlist <- append(df_Notch1_Zbtb7a_Zfp263_0, df_Notch1_0_Zbtb7a_Zfp263)
# 
# par(mfrow = c(1, 2))
# beeswarm(df_Notch1_Zbtb7a_Zfp263_0$monocle3_pseudotime)
# beeswarm(df_Notch1_0_Zbtb7a_Zfp263$monocle3_pseudotime)
# 
# library(ggplot2)
# 
# boxplot(df_Notch1_Zbtb7a_Zfp263_0$monocle3_pseudotime,
#         df_Notch1_0_Zbtb7a_Zfp263$monocle3_pseudotime)
# 
# wilcox.test(df_Notch1_Zbtb7a_Zfp263_0$monocle3_pseudotime,
#             df_Notch1_0_Zbtb7a_Zfp263$monocle3_pseudotime)
# 
# FeaturePlot(object = OV, features = c("Zbtb7a", "Zfp263", "Notch1", "Hes1"))