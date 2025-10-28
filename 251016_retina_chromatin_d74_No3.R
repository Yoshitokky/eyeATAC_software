# GSE183684

library(Signac)
library(Seurat)
library(GenomicRanges)
library(ggplot2)
library(patchwork)

set.seed(1)

multiome_d74 <- readRDS("241108_retina_chromatin_d74_No2.rds")

multiome_d74 <- RunTFIDF(multiome_d74)
multiome_d74 <- FindTopFeatures(multiome_d74, min.cutoff = 'q0')
multiome_d74 <- RunSVD(multiome_d74)

DepthCor(multiome_d74)

multiome_d74 <- RunUMAP(object = multiome_d74, reduction = 'lsi', dims = 2:30)
multiome_d74 <- FindNeighbors(object = multiome_d74, reduction = 'lsi', dims = 2:30)
multiome_d74 <- FindClusters(object = multiome_d74, verbose = FALSE, algorithm = 3)
DimPlot(object = multiome_d74, label = TRUE) + NoLegend()

gene.activities <- GeneActivity(multiome_d74)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
multiome_d74[['RNA']] <- CreateAssayObject(counts = gene.activities)
multiome_d74 <- NormalizeData(
  object = multiome_d74,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(multiome_d74$nCount_RNA)
)

DefaultAssay(multiome_d74) <- 'RNA'

FeaturePlot(
  object = multiome_d74,
  features = c('NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)

# RNAseq

RNA_d74 <- readRDS("241106_GSE183684_D74.rds")

new.cluster.ids <- c("Early markers+ RPC1_d74", "ELAVL2+/4+ RPC2_d74", "Early markers+ RPC2_d74", 
                     "Early markers+ RPC3_d74", "ELAVL2+/4+ RPC3_d74",
                     "ONECUT1+/2+/LHX1+ RPC_d74", 
                     "ELAVL2+/4+ RPC1_d74", "ONECUT1+/MEIS2+ RPC2_d74",
                     "Early markers+ RPC4_d74", "ONECUT1+/2+/PTF1A+ RPC_d74", "ONECUT1+/MEIS2+ RPC1_d74", 
                     "MKI67+ RPC_d74", "ONECUT1+/MEIS2+ RPC3_d74", "VSX2++ RPC_d74", "EEF1A1+ RPC_d74")
names(new.cluster.ids) <- levels(RNA_d74)
RNA_d74 <- RenameIdents(RNA_d74, new.cluster.ids)
RNA_d74$celltype <- Idents(RNA_d74)

DimPlot(RNA_d74, reduction = "umap", label = TRUE, repel = TRUE)

cell_counts <- table(RNA_d74$celltype)
print(cell_counts)
# Early markers+ RPC1_d74        ELAVL2+/4+ RPC2_d74    Early markers+ RPC2_d74 
# 1446                       1232                        989 
# Early markers+ RPC3_d74        ELAVL2+/4+ RPC3_d74  ONECUT1+/2+/LHX1+ RPC_d74 
# 957                        854                        690 
# ELAVL2+/4+ RPC1_d74   ONECUT1+/MEIS2+ RPC2_d74    Early markers+ RPC4_d74 
# 562                        455                        436 
# ONECUT1+/2+/PTF1A+ RPC_d74   ONECUT1+/MEIS2+ RPC1_d74             MKI67+ RPC_d74 
# 314                        306                        267 
# ONECUT1+/MEIS2+ RPC3_d74             VSX2++ RPC_d74            EEF1A1+ RPC_d74 
# 261                        158                         63 

RNA_d74 <- UpdateSeuratObject(RNA_d74)

transfer.anchors <- FindTransferAnchors(
  reference = RNA_d74,
  query = multiome_d74,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = RNA_d74$celltype,
  weight.reduction = multiome_d74[['lsi']],
  dims = 2:30
)

multiome_d74 <- AddMetaData(object = multiome_d74, metadata = predicted.labels)

# RNA_d74 <- readRDS("241108_retina_chromatin_d74_RNAseq_No3.rds")
# multiome_d74 <- readRDS("241108_retina_chromatin_d74_ATACseq_No3.rds")

celltypes <- levels(RNA_d74)
library(scales)
color_vector <- hue_pal()(length(celltypes))
rna_colors <- setNames(color_vector, celltypes)
show_col(rna_colors)
multiome_d74$celltype <- factor(multiome_d74$predicted.id, levels = names(rna_colors))

plot1 <- DimPlot(
  object = RNA_d74,
  group.by = 'celltype',
  label = TRUE,
  repel = TRUE,
  cols = rna_colors) + NoLegend() + ggtitle('scRNA-seq') + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20))

plot2 <- DimPlot(
  object = multiome_d74,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE,
  cols = rna_colors) + NoLegend() + ggtitle('scATAC-seq') + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20))

plot1
plot2

saveRDS(object = RNA_d74, file = "250829_retina_chromatin_d74_RNAseq_No3.rds")
saveRDS(object = multiome_d74, file = "250829_retina_chromatin_d74_ATACseq_No3.rds")

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
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=ja_JP.UTF-8       
# [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=ja_JP.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=ja_JP.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=ja_JP.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Asia/Tokyo
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] scales_1.4.0         future_1.58.0        patchwork_1.3.1      ggplot2_3.5.2       
# [5] GenomicRanges_1.56.2 GenomeInfoDb_1.40.1  IRanges_2.38.1       S4Vectors_0.42.1    
# [9] BiocGenerics_0.50.0  Seurat_5.3.0         SeuratObject_5.1.0   sp_2.2-0            
# [13] Signac_1.14.0       
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3      rstudioapi_0.17.1       jsonlite_2.0.0          magrittr_2.0.3         
# [5] spatstat.utils_3.1-4    farver_2.1.2            zlibbioc_1.50.0         vctrs_0.6.5            
# [9] ROCR_1.0-11             spatstat.explore_3.4-3  Rsamtools_2.20.0        RcppRoll_0.3.1         
# [13] htmltools_0.5.8.1       sctransform_0.4.2       parallelly_1.45.0       KernSmooth_2.23-24     
# [17] htmlwidgets_1.6.4       ica_1.0-3               plyr_1.8.9              plotly_4.11.0          
# [21] zoo_1.8-14              igraph_2.1.4            mime_0.13               lifecycle_1.0.4        
# [25] pkgconfig_2.0.3         Matrix_1.6-5            R6_2.6.1                fastmap_1.2.0          
# [29] GenomeInfoDbData_1.2.12 fitdistrplus_1.2-4      shiny_1.11.1            digest_0.6.37          
# [33] colorspace_2.1-1        tensor_1.5.1            RSpectra_0.16-2         irlba_2.3.5.1          
# [37] labeling_0.4.3          progressr_0.15.1        spatstat.sparse_3.1-0   httr_1.4.7             
# [41] polyclip_1.10-7         abind_1.4-8             compiler_4.4.1          withr_3.0.2            
# [45] BiocParallel_1.38.0     fastDummies_1.7.5       MASS_7.3-61             tools_4.4.1            
# [49] lmtest_0.9-40           httpuv_1.6.16           future.apply_1.20.0     goftest_1.2-3          
# [53] glue_1.8.0              nlme_3.1-165            promises_1.3.3          grid_4.4.1             
# [57] Rtsne_0.17              cluster_2.1.6           reshape2_1.4.4          generics_0.1.4         
# [61] gtable_0.3.6            spatstat.data_3.1-6     tidyr_1.3.1             data.table_1.17.8      
# [65] XVector_0.44.0          spatstat.geom_3.4-1     RcppAnnoy_0.0.22        ggrepel_0.9.6          
# [69] RANN_2.6.2              pillar_1.11.0           stringr_1.5.1           spam_2.11-1            
# [73] RcppHNSW_0.6.0          later_1.4.2             splines_4.4.1           dplyr_1.1.4            
# [77] lattice_0.22-5          survival_3.7-0          deldir_2.0-4            tidyselect_1.2.1       
# [81] Biostrings_2.72.1       miniUI_0.1.2            pbapply_1.7-2           gridExtra_2.3          
# [85] scattermore_1.2         matrixStats_1.5.0       stringi_1.8.7           UCSC.utils_1.0.0       
# [89] lazyeval_0.2.2          codetools_0.2-19        tibble_3.3.0            cli_3.6.5              
# [93] uwot_0.2.3              xtable_1.8-4            reticulate_1.42.0       dichromat_2.0-0.1      
# [97] Rcpp_1.1.0              globals_0.18.0          spatstat.random_3.4-1   png_0.1-8              
# [101] spatstat.univar_3.1-4   parallel_4.4.1          dotCall64_1.2           bitops_1.0-9           
# [105] listenv_0.9.1           viridisLite_0.4.2       ggridges_0.5.6          purrr_1.1.0            
# [109] crayon_1.5.3            rlang_1.1.6             cowplot_1.2.0           fastmatch_1.1-6