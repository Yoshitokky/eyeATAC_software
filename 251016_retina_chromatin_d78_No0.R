# GSE183684

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

D78 <- ReadMtx(
  mtx = "GSM5567528_d78_matrix.mtx.gz",
  features = "GSM5567528_d78_features.tsv.gz",
  cells = "GSM5567528_d78_barcodes.tsv.gz",
  feature.column = 2
)

RNA_d78 <- CreateSeuratObject(counts = D78, min.cells = 3, min.features = 200)

RNA_d78[["percent.mt"]] <- PercentageFeatureSet(RNA_d78, pattern = "^MT-")

VlnPlot(RNA_d78, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, layer = "counts")

RNA_d78 <- subset(RNA_d78, subset = nFeature_RNA > 1250 & nFeature_RNA < 5000 & percent.mt < 5)

RNA_d78 <- NormalizeData(RNA_d78, normalization.method = "LogNormalize", scale.factor = 10000)

RNA_d78 <- FindVariableFeatures(RNA_d78, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(RNA_d78)
RNA_d78 <- ScaleData(RNA_d78, features = all.genes)

RNA_d78 <- RunPCA(RNA_d78, features = VariableFeatures(object = RNA_d78))

ElbowPlot(RNA_d78)

RNA_d78 <- FindNeighbors(RNA_d78, dims = 1:15)

# RNA_d78_05_UMAP <- FindClusters(RNA_d78, resolution = 0.5)
# RNA_d78_05_UMAP <- RunUMAP(RNA_d78_05_UMAP, dims = 1:15)
# DimPlot(RNA_d78_05_UMAP, reduction = "umap")
 
# RNA_d78_06_UMAP <- FindClusters(RNA_d78, resolution = 0.6)
# RNA_d78_06_UMAP <- RunUMAP(RNA_d78_06_UMAP, dims = 1:15)
# DimPlot(RNA_d78_06_UMAP, reduction = "umap", label = TRUE)

RNA_d78_065_UMAP <- FindClusters(RNA_d78, resolution = 0.65)
RNA_d78_065_UMAP <- RunUMAP(RNA_d78_065_UMAP, dims = 1:15)
DimPlot(RNA_d78_065_UMAP, reduction = "umap", label = TRUE)

# RNA_d78_07_UMAP <- FindClusters(RNA_d78, resolution = 0.7)
# RNA_d78_07_UMAP <- RunUMAP(RNA_d78_07_UMAP, dims = 1:15)
# DimPlot(RNA_d78_07_UMAP, reduction = "umap", label =TRUE)

# saveRDS(RNA_d78_065_UMAP, file = "241108_GSE183684_D78.rds")

RNA_d78_065_UMAP <- readRDS("241108_GSE183684_D78.rds")

anno_markers <- c("LHX2", "PAX6", "RAX", "VSX2", "SIX3", "SOX2",
                  "MAP2","ELAVL2", "ELAVL4", "ONECUT1","ONECUT2", "MKI67",
                  "PAX2", "EEF1A1", "MEIS2", "LHX1", "PTF1A")
# Annotation strategy
# LHX2+ PAX6+ RAX+ VSX2+ --> RPC
# SIX3 SOX2 --> RPC, cited from PMID: 33723039
# MAP2+ --> Neural cells
# ELAVL2,4 --> terminally differentiate into AMA
# ONECUT1 --> terminally differentiate into RGC, HOR
# MKI67 --> Prolifiliration
# EEF1A1 --> Muller glia, rod photoreceptor, immature amacrine and ganglion

RNA_d78_065_anno <- readRDS("241108_GSE183684_D78.rds")

FeaturePlot(RNA_d78_065_anno, features = anno_markers)

DotPlot(RNA_d78_065_UMAP, features = anno_markers)  + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20))

new.cluster.ids <- c("Early markers+ RPC1_d78", "Early markers+ RPC2_d78",
                     "Early markers+ RPC3_d78", "ELAVL2+/4+ RPC2_d78", "ELAVL2+/4+ RPC1_d78",
                     "ONECUT1+/2+/LHX1+ RPC_d78", "ONECUT1+/2+/PTF1A+ RPC1_d78",
                     "ELAVL2+/4+ RPC3_d78", "ONECUT1+/MEIS2+ RPC3_d78",
                     "ONECUT1+/2+/PTF1A+ RPC2_d78", "ONECUT1+/MEIS2+ RPC1_d78", "ONECUT1+/MEIS2+ RPC2_d78",
                     "ONECUT1+/MEIS2+ RPC4_d78", "PAX6+/RAX+ RPC_d78",
                     "MKI67+ RPC2_d78", "MKI67+ RPC1_d78",
                     "VSX2++ RPC_d78", "EEF1A1+ RPC_d78")
names(new.cluster.ids) <- levels(RNA_d78_065_UMAP)
RNA_d78_065_anno <- RenameIdents(RNA_d78_065_UMAP, new.cluster.ids)

DotPlot(RNA_d78_065_anno, features = anno_markers) + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20)) +
  RotatedAxis() +
  scale_y_discrete(
    limits = c("Early markers+ RPC1_d78", "Early markers+ RPC2_d78", "Early markers+ RPC3_d78", "MKI67+ RPC1_d78", "MKI67+ RPC2_d78", 
               "PAX6+/RAX+ RPC_d78", "ONECUT1+/2+/PTF1A+ RPC1_d78",
               "ONECUT1+/2+/PTF1A+ RPC2_d78", "ONECUT1+/2+/LHX1+ RPC_d78",
               "ONECUT1+/MEIS2+ RPC1_d78", "ONECUT1+/MEIS2+ RPC2_d78",
               "ONECUT1+/MEIS2+ RPC3_d78", "ONECUT1+/MEIS2+ RPC4_d78",
               "VSX2++ RPC_d78", "ELAVL2+/4+ RPC1_d78","ELAVL2+/4+ RPC2_d78",
               "ELAVL2+/4+ RPC3_d78", "EEF1A1+ RPC_d78")
  )

# FeaturePlot(RNA_d78_065_anno, features = "MEIS2")
# PMID 30018341, RGC marker
 
# FeaturePlot(RNA_d78_065_anno, features = "LHX1")
# FeaturePlot(RNA_d78_065_anno, features = "PTF1A")
# PMID 27886389, horizontal cell marker

# FeaturePlot(RNA_d78_065_anno, features = anno_markers)

DimPlot(RNA_d78_065_anno, reduction = "umap", label =TRUE, repel = TRUE)  + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20))

FeaturePlot(RNA_d78_065_anno, features = "NOTCH1") + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20))
FeaturePlot(RNA_d78_065_anno, features = "NOTCH2") + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20))
FeaturePlot(RNA_d78_065_anno, features = "NOTCH3") + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20))
FeaturePlot(RNA_d78_065_anno, features = "NOTCH4") + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20))

notch <- c("NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4")
DotPlot(RNA_d78_065_anno, features = notch) + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20)) +
  RotatedAxis() +
  scale_y_discrete(
    limits = c("Early markers+ RPC1_d78", "Early markers+ RPC2_d78", "Early markers+ RPC3_d78", "MKI67+ RPC1_d78", "MKI67+ RPC2_d78", 
               "PAX6+/RAX+ RPC_d78", "ONECUT1+/2+/PTF1A+ RPC1_d78",
               "ONECUT1+/2+/PTF1A+ RPC2_d78", "ONECUT1+/2+/LHX1+ RPC_d78",
               "ONECUT1+/MEIS2+ RPC1_d78", "ONECUT1+/MEIS2+ RPC2_d78",
               "ONECUT1+/MEIS2+ RPC3_d78", "ONECUT1+/MEIS2+ RPC4_d78",
               "VSX2++ RPC_d78", "ELAVL2+/4+ RPC1_d78","ELAVL2+/4+ RPC2_d78",
               "ELAVL2+/4+ RPC3_d78", "EEF1A1+ RPC_d78")
  )

library(SeuratWrappers)
library(monocle3)

# the next line is necessary if the data are bound together
# whole <- JoinLayers(whole)
cds <- as.cell_data_set(RNA_d78_065_anno)

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

cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 1]))
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

library(viridis)

FeaturePlot(integrated.sub, "monocle3_pseudotime") + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20)) +
  scale_color_viridis(option = "C") 

# sessionInfo()
# R version 4.4.1 (2024-06-14)
# Platform: x86_64-pc-linux-gnu
# Running under: Ubuntu 22.04.4 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0 
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=ja_JP.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=ja_JP.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=ja_JP.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=ja_JP.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Asia/Tokyo
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] viridis_0.6.5               viridisLite_0.4.2          
# [3] monocle3_1.3.7              SingleCellExperiment_1.26.0
# [5] SummarizedExperiment_1.34.0 GenomicRanges_1.56.2       
# [7] GenomeInfoDb_1.40.1         IRanges_2.38.1             
# [9] S4Vectors_0.42.1            MatrixGenerics_1.16.0      
# [11] matrixStats_1.4.1           Biobase_2.64.0             
# [13] BiocGenerics_0.50.0         SeuratWrappers_0.3.5       
# [15] ggplot2_3.5.1               patchwork_1.3.0            
# [17] Seurat_5.1.0                SeuratObject_5.0.2         
# [19] sp_2.1-4                    dplyr_1.1.4                
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3      rstudioapi_0.17.1       jsonlite_1.8.9         
# [4] magrittr_2.0.3          spatstat.utils_3.1-0    ggbeeswarm_0.7.2       
# [7] nloptr_2.1.1            farver_2.1.2            zlibbioc_1.50.0        
# [10] vctrs_0.6.5             ROCR_1.0-11             minqa_1.2.8            
# [13] spatstat.explore_3.3-3  S4Arrays_1.4.1          htmltools_0.5.8.1      
# [16] SparseArray_1.4.8       sctransform_0.4.1       parallelly_1.38.0      
# [19] KernSmooth_2.23-24      htmlwidgets_1.6.4       ica_1.0-3              
# [22] plyr_1.8.9              plotly_4.10.4           zoo_1.8-12             
# [25] igraph_2.1.1            mime_0.12               lifecycle_1.0.4        
# [28] pkgconfig_2.0.3         rsvd_1.0.5              Matrix_1.6-5           
# [31] R6_2.5.1                fastmap_1.2.0           GenomeInfoDbData_1.2.12
# [34] fitdistrplus_1.2-1      future_1.34.0           shiny_1.9.1            
# [37] digest_0.6.37           colorspace_2.1-1        tensor_1.5             
# [40] RSpectra_0.16-2         irlba_2.3.5.1           labeling_0.4.3         
# [43] progressr_0.15.0        fansi_1.0.6             spatstat.sparse_3.1-0  
# [46] httr_1.4.7              polyclip_1.10-7         abind_1.4-8            
# [49] compiler_4.4.1          proxy_0.4-27            remotes_2.5.0          
# [52] withr_3.0.2             fastDummies_1.7.4       R.utils_2.12.3         
# [55] MASS_7.3-61             DelayedArray_0.30.1     tools_4.4.1            
# [58] vipor_0.4.7             lmtest_0.9-40           beeswarm_0.4.0         
# [61] httpuv_1.6.15           future.apply_1.11.3     goftest_1.2-3          
# [64] R.oo_1.26.0             glue_1.8.0              nlme_3.1-165           
# [67] promises_1.3.0          grid_4.4.1              Rtsne_0.17             
# [70] cluster_2.1.6           reshape2_1.4.4          generics_0.1.3         
# [73] gtable_0.3.6            spatstat.data_3.1-2     R.methodsS3_1.8.2      
# [76] tidyr_1.3.1             data.table_1.16.2       XVector_0.44.0         
# [79] utf8_1.2.4              spatstat.geom_3.3-3     RcppAnnoy_0.0.22       
# [82] ggrepel_0.9.6           RANN_2.6.2              pillar_1.9.0           
# [85] stringr_1.5.1           spam_2.11-0             RcppHNSW_0.6.0         
# [88] later_1.3.2             splines_4.4.1           lattice_0.22-5         
# [91] survival_3.7-0          deldir_2.0-4            tidyselect_1.2.1       
# [94] miniUI_0.1.1.1          pbapply_1.7-2           gridExtra_2.3          
# [97] scattermore_1.2         leidenbase_0.1.31       UCSC.utils_1.0.0       
# [100] stringi_1.8.4           boot_1.3-30             lazyeval_0.2.2         
# [103] codetools_0.2-19        tibble_3.2.1            BiocManager_1.30.25    
# [106] cli_3.6.3               uwot_0.2.2              xtable_1.8-4           
# [109] reticulate_1.39.0       munsell_0.5.1           Rcpp_1.0.13            
# [112] globals_0.16.3          spatstat.random_3.3-2   png_0.1-8              
# [115] ggrastr_1.0.2           spatstat.univar_3.0-1   parallel_4.4.1         
# [118] assertthat_0.2.1        dotCall64_1.2           lme4_1.1-35.5          
# [121] listenv_0.9.1           scales_1.3.0            ggridges_0.5.6         
# [124] crayon_1.5.3            leiden_0.4.3.1          purrr_1.0.2            
# [127] rlang_1.1.4             cowplot_1.1.3