# human_retina_organoids > GSE220661

library(dplyr)
library(Seurat)
library(ggplot2)

whole <- readRDS("GSE220661_UMAP_QC_240919.rds")

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
    limits = c("RPC1", "RPC2", "RPC3", "RPC4", "RPC5", "RPC6", "RPC7",
               "MULLER", "LENS", "ECTO1", "ECTO2", "Npre", "HOR", "CONE",
               "ROD", "BIPOLAR")
  )

new.cluster.ids <- c("RPC1", "Npre", "RPC2", "CONE", "RPC3",
                     "RPC4", "ROD", "RPC5", "LENS", "RPC6", "MULLER",
                     "ECTO1", "RPC7", "BIPOLAR", "HOR", "ECTO2")
names(new.cluster.ids) <- levels(whole)
whole <- RenameIdents(whole, new.cluster.ids)

DimPlot(whole, reduction = "umap", label = TRUE, label.size = 10, repel = TRUE) + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20))

RPCs <- subset(whole, idents = c("RPC1", "RPC2", "RPC3", "RPC4",
                               "RPC5", "RPC6", "RPC7"))

DimPlot(whole, split.by = "orig.ident", label = TRUE, label.size = 10, repel = TRUE) + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20))

notch <- c("NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4")
DotPlot(object = whole, features = notch) + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20)) +
  scale_y_discrete(
    limits = c("RPC1", "RPC2", "RPC3", "RPC4", "RPC5", "RPC6", "RPC7",
               "MULLER", "LENS", "ECTO1", "ECTO2", "Npre", "HOR", "CONE",
               "ROD", "BIPOLAR")
  )

target <- c("HES1", "HES3", "HES5", "HES6", "HES7", "HEY1", "HEY2", "HEYL")
DotPlot(object = whole, features = target) + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20)) +
  RotatedAxis() +
  scale_y_discrete(
    limits = c("RPC1", "RPC2", "RPC3", "RPC4", "RPC5", "RPC6", "RPC7",
               "MULLER", "LENS", "ECTO1", "ECTO2", "Npre", "HOR", "CONE",
               "ROD", "BIPOLAR")
  )

#delta <- c("JAG1", "JAG2", "DLL1", "DLL3", "DLL4")
#DotPlot(object = whole, features = delta) + 
#  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
#  theme(legend.text = element_text(size = 20))

#cluster.markers <- FindMarkers(whole, ident.1 = "LEU")
#head(cluster.markers, n = 15)

library(CSCORE)
RPCs_joined <- JoinLayers(RPCs)
stage = RPCs_joined
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
writeLines(list_to_write, "module_list_GSE220661_RPCs.txt")

library(SeuratWrappers)
library(monocle3)

whole_joined <- JoinLayers(whole)
cds <- as.cell_data_set(whole_joined)
cds <- cluster_cells(cds, resolution=1e-3)

p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
plot(p1, label = TRUE)
plot(p2)

integrated.sub <- subset(as.Seurat(cds, assay = NULL), monocle3_partitions == 1 | 2)
cds <- as.cell_data_set(integrated.sub)

cds <- learn_graph(cds, use_partition = FALSE, verbose = FALSE)

plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 3]))

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
FeaturePlot(integrated.sub, "monocle3_pseudotime") + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20))

#head(integrated.sub[[]])
#head(whole[[]])
#colData(cds)
#levels(integrated.sub)

integrated.sub <- SetIdent(integrated.sub, 
                           value = integrated.sub@meta.data$seurat_clusters)

df2 <- FetchData(object = integrated.sub, vars = c("monocle3_pseudotime"), layer = "count")
df3 <- FetchData(object = whole, vars = c("NOTCH1", "NOTCH3", "HES1", "HES5", "HES6", "JAG1",
                                          "JAG2", "DLL1", "DLL3", "DLL4", "ident"), layer = "scale.data")
df4 <- cbind(df2, df3)
df5 <- subset(df4, df4$ident == "RPC1" | df4$ident == "RPC2" | df4$ident == "RPC3"
              | df4$ident == "RPC4" | df4$ident == "RPC5" | df4$ident == "RPC6"
              | df4$ident == "RPC7")

library(ggplot2)
# ggplot(df4, aes(monocle3_pseudotime, NOTCH1)) + geom_point() + geom_smooth()
# ggplot(df4, aes(monocle3_pseudotime, NOTCH3)) + geom_point()
# ggplot(df4, aes(monocle3_pseudotime, HES6)) + geom_point() 
# ggplot(df4, aes(monocle3_pseudotime, JAG1)) + geom_point() 
# ggplot(df4, aes(monocle3_pseudotime, JAG2)) + geom_point() 
# ggplot(df4, aes(monocle3_pseudotime, DLL1)) + geom_point() 
# ggplot(df4, aes(monocle3_pseudotime, DLL3)) + geom_point() 
# ggplot(df4, aes(monocle3_pseudotime, DLL4)) + geom_point() 
# ggplot(df4) +
#   geom_point(aes(monocle3_pseudotime, NOTCH1), color = "red") +
#   geom_point(aes(monocle3_pseudotime, NOTCH3), color = "blue")+
#   geom_point(aes(monocle3_pseudotime, HES1), color = "green")
# ggplot(df4) +
#   geom_point(aes(monocle3_pseudotime, NOTCH1), color = "red") +
#   geom_point(aes(monocle3_pseudotime, JAG2), color = "blue")
# ggplot(df4) +
#   geom_point(aes(monocle3_pseudotime, DLL1), color = "red") +
#   geom_point(aes(monocle3_pseudotime, DLL3), color = "blue")+
#   geom_point(aes(monocle3_pseudotime, DLL4), color = "green")
# correlation <- cor(df4$NOTCH1, df4$NOTCH3) # 0.117
# ggplot(df4, aes(NOTCH1, NOTCH3)) + geom_point()
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

#correlation <- cor(df4$Notch1, df4$Zbtb7a)
#correlation <- cor(df4$Notch1, df4$Zfp263) 

#export as tsv
# whole_joined <- JoinLayers(whole)
# mat <- GetAssayData(object = whole_joined, assay = "RNA", layer = "data")
# write.csv(mat, "GSE220661_240913.csv")
# 
# plot_cells(integrated.sub,
#            genes=delta,
#            label_cell_groups=FALSE,
#            show_trajectory_graph=FALSE)

# dorothea
# line 2734 Esr1 A Notch2
# line 4997 Klf5 C Notch1
# line 9287 Snai2 C Notch1
# line 12143 Zbtb7a C Notch1
# line 12509 Zfp263 B Notch1
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

# > sessionInfo()
# R version 4.4.1 (2024-06-14)
# Platform: x86_64-pc-linux-gnu
# Running under: Ubuntu 22.04.4 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0 
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=ja_JP.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=ja_JP.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=ja_JP.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=ja_JP.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Asia/Tokyo
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] beeswarm_0.4.0              ggVennDiagram_1.5.2         ggplot2_3.5.1               monocle3_1.3.7             
# [5] SingleCellExperiment_1.26.0 SummarizedExperiment_1.34.0 GenomicRanges_1.56.1        GenomeInfoDb_1.40.1        
# [9] IRanges_2.38.1              S4Vectors_0.42.1            MatrixGenerics_1.16.0       matrixStats_1.4.1          
# [13] Biobase_2.64.0              BiocGenerics_0.50.0         SeuratWrappers_0.3.5        WGCNA_1.72-5               
# [17] fastcluster_1.2.6           dynamicTreeCut_1.63-1       CSCORE_0.0.0.9000           Seurat_5.1.0               
# [21] SeuratObject_5.0.2          sp_2.1-4                    dplyr_1.1.4                
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.22        splines_4.4.1           later_1.3.2             tibble_3.2.1           
# [5] R.oo_1.26.0             polyclip_1.10-7         preprocessCore_1.66.0   rpart_4.1.23           
# [9] fastDummies_1.7.4       lifecycle_1.0.4         doParallel_1.0.17       globals_0.16.3         
# [13] lattice_0.22-5          MASS_7.3-61             backports_1.5.0         magrittr_2.0.3         
# [17] Hmisc_5.1-3             plotly_4.10.4           rmarkdown_2.28          remotes_2.5.0          
# [21] httpuv_1.6.15           sctransform_0.4.1       spam_2.10-0             spatstat.sparse_3.1-0  
# [25] reticulate_1.39.0       minqa_1.2.8             cowplot_1.1.3           pbapply_1.7-2          
# [29] DBI_1.2.3               RColorBrewer_1.1-3      abind_1.4-8             zlibbioc_1.50.0        
# [33] Rtsne_0.17              purrr_1.0.2             R.utils_2.12.3          nnet_7.3-19            
# [37] GenomeInfoDbData_1.2.12 ggrepel_0.9.6           irlba_2.3.5.1           listenv_0.9.1          
# [41] spatstat.utils_3.1-0    goftest_1.2-3           RSpectra_0.16-2         spatstat.random_3.3-1  
# [45] fitdistrplus_1.2-1      parallelly_1.38.0       DelayedArray_0.30.1     leiden_0.4.3.1         
# [49] codetools_0.2-19        tidyselect_1.2.1        UCSC.utils_1.0.0        lme4_1.1-35.5          
# [53] base64enc_0.1-3         spatstat.explore_3.3-2  jsonlite_1.8.8          progressr_0.14.0       
# [57] Formula_1.2-5           ggridges_0.5.6          survival_3.7-0          iterators_1.0.14       
# [61] foreach_1.5.2           tools_4.4.1             ica_1.0-3               Rcpp_1.0.13            
# [65] glue_1.7.0              SparseArray_1.4.8       gridExtra_2.3           xfun_0.47              
# [69] withr_3.0.1             BiocManager_1.30.25     fastmap_1.2.0           boot_1.3-30            
# [73] fansi_1.0.6             rsvd_1.0.5              digest_0.6.37           R6_2.5.1               
# [77] mime_0.12               colorspace_2.1-1        scattermore_1.2         GO.db_3.19.1           
# [81] tensor_1.5              spatstat.data_3.1-2     RSQLite_2.3.7           R.methodsS3_1.8.2      
# [85] utf8_1.2.4              tidyr_1.3.1             generics_0.1.3          data.table_1.16.0      
# [89] S4Arrays_1.4.1          httr_1.4.7              htmlwidgets_1.6.4       uwot_0.2.2             
# [93] pkgconfig_2.0.3         gtable_0.3.5            blob_1.2.4              lmtest_0.9-40          
# [97] impute_1.78.0           XVector_0.44.0          htmltools_0.5.8.1       dotCall64_1.1-1        
# [101] scales_1.3.0            png_0.1-8               spatstat.univar_3.0-1   knitr_1.48             
# [105] rstudioapi_0.16.0       reshape2_1.4.4          nloptr_2.1.1            checkmate_2.3.2        
# [109] nlme_3.1-165            cachem_1.1.0            zoo_1.8-12              stringr_1.5.1          
# [113] KernSmooth_2.23-24      parallel_4.4.1          miniUI_0.1.1.1          foreign_0.8-86         
# [117] AnnotationDbi_1.66.0    pillar_1.9.0            grid_4.4.1              vctrs_0.6.5            
# [121] RANN_2.6.2              promises_1.3.0          xtable_1.8-4            cluster_2.1.6          
# [125] htmlTable_2.4.3         evaluate_0.24.0         cli_3.6.3               compiler_4.4.1         
# [129] rlang_1.1.4             crayon_1.5.3            future.apply_1.11.2     plyr_1.8.9             
# [133] stringi_1.8.4           viridisLite_0.4.2       deldir_2.0-4            munsell_0.5.1          
# [137] Biostrings_2.72.1       lazyeval_0.2.2          spatstat.geom_3.3-2     Matrix_1.6-5           
# [141] RcppHNSW_0.6.0          patchwork_1.2.0         bit64_4.0.5             future_1.34.0          
# [145] KEGGREST_1.44.1         shiny_1.9.1             ROCR_1.0-11             igraph_2.0.3           
# [149] memoise_2.0.1           bit_4.0.5 