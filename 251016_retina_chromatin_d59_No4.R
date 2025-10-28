# GSE183684

library(Signac)
library(Seurat)
library(GenomicRanges)
library(ggplot2)
library(patchwork)
library(EnsDb.Hsapiens.v86)
library(GenomicRanges)
library(GenomeInfoDb)

RNA_d59 <- readRDS("250827_retina_chromatin_d59_RNAseq_No3.rds")
multiome_d59 <- readRDS("250827_retina_chromatin_d59_ATACseq_No3.rds")

# change back to working with peaks instead of gene activities
DefaultAssay(multiome_d59) <- 'peaks'

Idents(object = multiome_d59) <- "predicted.id"

multiome_d59 <- SortIdents(multiome_d59)

peaks <- granges(multiome_d59[["peaks"]])

gene_anno <- Annotation(multiome_d59)

# NOTCH1
region_NOTCH1 <- gene_anno[gene_anno$gene_name == "NOTCH1"]
region_NOTCH1 <- Extend(region_NOTCH1, upstream = 100000, downstream = 100000)
peaks_in_region_NOTCH1 <- subsetByOverlaps(peaks, region_NOTCH1)
peaks_in_region_NOTCH1
# write.csv(peaks_in_region_NOTCH1, "ATAC_peaks_NOTCH1_d59.csv")
result_NOTCH1 <- read.csv("ATAC_peaks_NOTCH1_d59.csv")
result_NOTCH1

target_region_NOTCH1 <- GRanges(seqnames = "chr9",
                            ranges = IRanges(start = 136500000, end = 136610000))
peaks_in_target_region_NOTCH1 <- subsetByOverlaps(peaks, target_region_NOTCH1)
peaks_in_target_region_NOTCH1
# write.csv(peaks_in_target_region_NOTCH1, "ATAC_peaks_NOTCH1_roi_d59.csv")
result_NOTCH1_roi <- read.csv("ATAC_peaks_NOTCH1_roi_d59.csv")
result_NOTCH1_roi

# Idents(multiome_d59) <- factor(Idents(multiome_d59),
#                                levels = c("Early markers+ RPC1_d59",
#                                           "Early markers+ RPC2_d59",
#                                           "Early markers+ RPC3_d59",
#                                           "MKI67+ RPC1_d59", 
#                                           "MKI67+ RPC2_d59", 
#                                           "ONECUT1+/2+/LHX1+/PTF1A+ RPC_d59",
#                                           "PAX6+/RAX+ RPC_d59",
#                                           "ELAVL2+/4+ RPC1_d59",
#                                           "ELAVL2+/4+ RPC2_d59",
#                                           "ELAVL2+/4+ RPC3_d59",
#                                           "ELAVL2+/4+ RPC4_d59",
#                                           "EEF1A1+ RPC_d59",
#                                           "GLIA_d59"))

CoveragePlot(
  object = multiome_d59,
  region = "NOTCH1",
  extend.upstream = 100000,
  extend.downstream = 100000
)

CoveragePlot(
  object = multiome_d59,
  region = "NOTCH1",
  extend.upstream = 1000,
  extend.downstream = 1000
)

# Notch1 arrow1
# library(JASPAR2020)
# library(TFBSTools)
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(dplyr)
# 
# pfm <- getMatrixSet(
#   x = JASPAR2020,
#   opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
# )
# 
# test <- AddMotifs(
#   object = multiome_d59,
#   genome = N1A1,
#   pfm = pfm,
#   assay = "peaks"
# )

# earlyRPC <- subset(multiome_d59,
#                    multiome_d59@meta.data$predicted.id == c("Early markers+ RPC1_d59",
#                                                             "Early markers+ RPC2_d59",
#                                                             "Early markers+ RPC3_d59",
#                                                             "MKI67+ RPC1_d59",
#                                                             "MKI67+ RPC2_d59"))
# 
# N1A1_LHX2 <- GRanges(seqnames = "chr9",
#                 ranges = IRanges(start = 136542391, end = 136542396))
# accessible_counts_N1A1_LHX2 <- FeatureMatrix(
#   fragments = Fragments(earlyRPC),
#   features = N1A1_LHX2,
#   cells = colnames(earlyRPC)
# )
# head(accessible_counts_N1A1_LHX2)
# df_N1A1_LHX2 <- data.frame(cell = colnames(accessible_counts_N1A1_LHX2),
#                  count = as.numeric(accessible_counts_N1A1_LHX2[1, ]),
#                  cluster = earlyRPC@meta.data$predicted.id)
# head(df_N1A1_LHX2)
# sum(df_N1A1_LHX2$count > 0)
# # [1] 53
# sum(df_N1A1_LHX2$count == 0)
# # [1] 412
# nrow(df_N1A1_LHX2)
# # [1] 465
# sum(df_N1A1_LHX2$count != 0)

# NOTCH2
CoveragePlot(
  object = multiome_d59,
  region = "NOTCH2",
  extend.upstream = 100000,
  extend.downstream = 100000
)

CoveragePlot(
  object = multiome_d59,
  region = "NOTCH2",
  extend.upstream = 1000,
  extend.downstream = 1000
)

# NOTCH3
region_NOTCH3 <- gene_anno[gene_anno$gene_name == "NOTCH3"]
region_NOTCH3 <- Extend(region_NOTCH3, upstream = 100000, downstream = 100000)
peaks_in_region_NOTCH3 <- subsetByOverlaps(peaks, region_NOTCH3)
peaks_in_region_NOTCH3
# write.csv(peaks_in_region_NOTCH3, "ATAC_peaks_NOTCH3_d59.csv")
result_NOTCH3 <- read.csv("ATAC_peaks_NOTCH3_d59.csv")
result_NOTCH3

target_region_NOTCH3 <- GRanges(seqnames = "chr19",
                                ranges = IRanges(start = 15100000, end = 15120000))
peaks_in_target_region_NOTCH3 <- subsetByOverlaps(peaks, target_region_NOTCH3)
peaks_in_target_region_NOTCH3
# write.csv(peaks_in_target_region_NOTCH3, "ATAC_peaks_NOTCH3_roi_d59.csv")
result_NOTCH3_roi <- read.csv("ATAC_peaks_NOTCH3_roi_d59.csv")
result_NOTCH3_roi

CoveragePlot(
  object = multiome_d59,
  region = "NOTCH3",
  extend.upstream = 100000,
  extend.downstream = 100000
)

CoveragePlot(
  object = multiome_d59,
  region = "NOTCH3",
  extend.upstream = 1000,
  extend.downstream = 1000
)

# NOTCH4
CoveragePlot(
  object = multiome_d59,
  region = "NOTCH4",
  extend.upstream = 100000,
  extend.downstream = 100000
)

CoveragePlot(
  object = multiome_d59,
  region = "NOTCH4",
  extend.upstream = 1000,
  extend.downstream = 1000
)

region_NOTCH4 <- gene_anno[gene_anno$gene_name == "NOTCH4"]
region_NOTCH4 <- Extend(region_NOTCH4, upstream = 100000, downstream = 100000)
peaks_in_region_NOTCH4 <- subsetByOverlaps(peaks, region_NOTCH4)
peaks_in_region_NOTCH4
# write.csv(peaks_in_region_NOTCH4, "ATAC_peaks_NOTCH4_d59.csv")
result_NOTCH4 <- read.csv("ATAC_peaks_NOTCH4_d59.csv")
result_NOTCH4

target_region_NOTCH4 <- GRanges(seqnames = "chr6",
                                ranges = IRanges(start = 32110000, end = 32200000))
peaks_in_target_region_NOTCH4 <- subsetByOverlaps(peaks, target_region_NOTCH4)
peaks_in_target_region_NOTCH4
# write.csv(peaks_in_target_region_NOTCH4, "ATAC_peaks_NOTCH4_roi_d59.csv")
result_NOTCH4_roi <- read.csv("ATAC_peaks_NOTCH4_roi_d59.csv")
result_NOTCH4_roi

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
#   [1] EnsDb.Hsapiens.v86_2.99.0 ensembldb_2.28.1          AnnotationFilter_1.28.0   GenomicFeatures_1.56.0   
# [5] AnnotationDbi_1.66.0      Biobase_2.64.0            patchwork_1.3.1           ggplot2_3.5.2            
# [9] GenomicRanges_1.56.2      GenomeInfoDb_1.40.1       IRanges_2.38.1            S4Vectors_0.42.1         
# [13] BiocGenerics_0.50.0       Seurat_5.3.0              SeuratObject_5.1.0        sp_2.2-0                 
# [17] Signac_1.14.0            
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.22            splines_4.4.1               later_1.4.2                
# [4] BiocIO_1.14.0               bitops_1.0-9                tibble_3.3.0               
# [7] polyclip_1.10-7             XML_3.99-0.18               fastDummies_1.7.5          
# [10] lifecycle_1.0.4             globals_0.18.0              lattice_0.22-5             
# [13] MASS_7.3-61                 magrittr_2.0.3              plotly_4.11.0              
# [16] yaml_2.3.10                 httpuv_1.6.16               sctransform_0.4.2          
# [19] spam_2.11-1                 spatstat.sparse_3.1-0       reticulate_1.42.0          
# [22] cowplot_1.2.0               pbapply_1.7-2               DBI_1.2.3                  
# [25] RColorBrewer_1.1-3          abind_1.4-8                 zlibbioc_1.50.0            
# [28] Rtsne_0.17                  purrr_1.1.0                 RCurl_1.98-1.17            
# [31] GenomeInfoDbData_1.2.12     ggrepel_0.9.6               irlba_2.3.5.1              
# [34] listenv_0.9.1               spatstat.utils_3.1-4        goftest_1.2-3              
# [37] RSpectra_0.16-2             spatstat.random_3.4-1       fitdistrplus_1.2-4         
# [40] parallelly_1.45.0           codetools_0.2-19            DelayedArray_0.30.1        
# [43] RcppRoll_0.3.1              tidyselect_1.2.1            UCSC.utils_1.0.0           
# [46] farver_2.1.2                matrixStats_1.5.0           spatstat.explore_3.4-3     
# [49] GenomicAlignments_1.40.0    jsonlite_2.0.0              progressr_0.15.1           
# [52] ggridges_0.5.6              survival_3.7-0              tools_4.4.1                
# [55] ica_1.0-3                   Rcpp_1.1.0                  glue_1.8.0                 
# [58] gridExtra_2.3               SparseArray_1.4.8           MatrixGenerics_1.16.0      
# [61] dplyr_1.1.4                 withr_3.0.2                 fastmap_1.2.0              
# [64] digest_0.6.37               R6_2.6.1                    mime_0.13                  
# [67] colorspace_2.1-1            scattermore_1.2             tensor_1.5.1               
# [70] dichromat_2.0-0.1           spatstat.data_3.1-6         RSQLite_2.4.1              
# [73] tidyr_1.3.1                 generics_0.1.4              data.table_1.17.8          
# [76] rtracklayer_1.64.0          httr_1.4.7                  htmlwidgets_1.6.4          
# [79] S4Arrays_1.4.1              uwot_0.2.3                  pkgconfig_2.0.3            
# [82] gtable_0.3.6                blob_1.2.4                  lmtest_0.9-40              
# [85] XVector_0.44.0              htmltools_0.5.8.1           dotCall64_1.2              
# [88] ProtGenerics_1.36.0         scales_1.4.0                png_0.1-8                  
# [91] spatstat.univar_3.1-4       rstudioapi_0.17.1           reshape2_1.4.4             
# [94] rjson_0.2.23                nlme_3.1-165                curl_6.4.0                 
# [97] zoo_1.8-14                  cachem_1.1.0                stringr_1.5.1              
# [100] KernSmooth_2.23-24          parallel_4.4.1              miniUI_0.1.2               
# [103] restfulr_0.0.16             pillar_1.11.0               grid_4.4.1                 
# [106] vctrs_0.6.5                 RANN_2.6.2                  promises_1.3.3             
# [109] xtable_1.8-4                cluster_2.1.6               cli_3.6.5                  
# [112] compiler_4.4.1              Rsamtools_2.20.0            rlang_1.1.6                
# [115] crayon_1.5.3                future.apply_1.20.0         plyr_1.8.9                 
# [118] stringi_1.8.7               viridisLite_0.4.2           deldir_2.0-4               
# [121] BiocParallel_1.38.0         Biostrings_2.72.1           lazyeval_0.2.2             
# [124] spatstat.geom_3.4-1         Matrix_1.6-5                RcppHNSW_0.6.0             
# [127] bit64_4.6.0-1               future_1.58.0               KEGGREST_1.44.1            
# [130] shiny_1.11.1                SummarizedExperiment_1.34.0 ROCR_1.0-11                
# [133] igraph_2.1.4                memoise_2.0.1               fastmatch_1.1-6            
# [136] bit_4.6.0 