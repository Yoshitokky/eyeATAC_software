library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

# in terminal

# gzip -d GSM5567518_d59_fragments.tsv.gz
# bgzip GSM5567518_d59_fragments.tsv
# tabix -p bed GSM5567518_d59_fragments.tsv.gz

RNA_d59 <- ReadMtx(mtx = "GSM5567526_d59_matrix.mtx.gz",
                   cells = "GSM5567526_d59_barcodes.tsv.gz",
                   features = "GSM5567526_d59_features.tsv.gz")

data_d59 <- CreateSeuratObject(counts = RNA_d59, assay = "RNA")

fragpath <- 'GSM5567518_d59_fragments.tsv.gz'

total_counts <- CountFragments(fragpath)

cutoff <- 1000 # 1000 in tutorial

barcodes <- total_counts[total_counts$frequency_count > cutoff, ]$CB

frags <- CreateFragmentObject(path = fragpath, cells = barcodes)

peaks <- CallPeaks(frags, macs2.path = "/home/notch/anaconda3/bin/macs2")

counts <- FeatureMatrix(fragments = frags, features = peaks, cells = barcodes)

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

data_d59 <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks"
)

saveRDS(object = data_d59, file = "241108_retina_chromatin_d59_No1.rds")

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
#   [1] BSgenome.Hsapiens.UCSC.hg38_1.4.5 BSgenome_1.72.0                  
# [3] rtracklayer_1.64.0                BiocIO_1.14.0                    
# [5] Biostrings_2.72.1                 XVector_0.44.0                   
# [7] EnsDb.Hsapiens.v86_2.99.0         ensembldb_2.28.1                 
# [9] AnnotationFilter_1.28.0           GenomicFeatures_1.56.0           
# [11] AnnotationDbi_1.66.0              Biobase_2.64.0                   
# [13] GenomicRanges_1.56.2              GenomeInfoDb_1.40.1              
# [15] IRanges_2.38.1                    S4Vectors_0.42.1                 
# [17] BiocGenerics_0.50.0               Signac_1.14.0                    
# [19] Seurat_5.1.0                      SeuratObject_5.0.2               
# [21] sp_2.1-4                         
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.22            splines_4.4.1               later_1.3.2                
# [4] bitops_1.0-9                tibble_3.2.1                polyclip_1.10-7            
# [7] rpart_4.1.23                XML_3.99-0.17               fastDummies_1.7.4          
# [10] lifecycle_1.0.4             globals_0.16.3              lattice_0.22-5             
# [13] MASS_7.3-61                 backports_1.5.0             magrittr_2.0.3             
# [16] rmarkdown_2.28              Hmisc_5.2-0                 plotly_4.10.4              
# [19] yaml_2.3.10                 httpuv_1.6.15               sctransform_0.4.1          
# [22] spam_2.11-0                 spatstat.sparse_3.1-0       reticulate_1.39.0          
# [25] cowplot_1.1.3               pbapply_1.7-2               DBI_1.2.3                  
# [28] RColorBrewer_1.1-3          abind_1.4-8                 zlibbioc_1.50.0            
# [31] Rtsne_0.17                  purrr_1.0.2                 biovizBase_1.52.0          
# [34] RCurl_1.98-1.16             nnet_7.3-19                 VariantAnnotation_1.50.0   
# [37] GenomeInfoDbData_1.2.12     ggrepel_0.9.6               irlba_2.3.5.1              
# [40] listenv_0.9.1               spatstat.utils_3.1-0        goftest_1.2-3              
# [43] RSpectra_0.16-2             spatstat.random_3.3-2       fitdistrplus_1.2-1         
# [46] parallelly_1.38.0           leiden_0.4.3.1              codetools_0.2-19           
# [49] DelayedArray_0.30.1         RcppRoll_0.3.1              tidyselect_1.2.1           
# [52] UCSC.utils_1.0.0            farver_2.1.2                base64enc_0.1-3            
# [55] matrixStats_1.4.1           spatstat.explore_3.3-3      GenomicAlignments_1.40.0   
# [58] jsonlite_1.8.9              Formula_1.2-5               progressr_0.15.0           
# [61] ggridges_0.5.6              survival_3.7-0              tools_4.4.1                
# [64] ica_1.0-3                   Rcpp_1.0.13                 glue_1.8.0                 
# [67] gridExtra_2.3               SparseArray_1.4.8           xfun_0.48                  
# [70] MatrixGenerics_1.16.0       dplyr_1.1.4                 withr_3.0.2                
# [73] fastmap_1.2.0               fansi_1.0.6                 digest_0.6.37              
# [76] R6_2.5.1                    mime_0.12                   colorspace_2.1-1           
# [79] scattermore_1.2             tensor_1.5                  dichromat_2.0-0.1          
# [82] spatstat.data_3.1-2         RSQLite_2.3.7               utf8_1.2.4                 
# [85] tidyr_1.3.1                 generics_0.1.3              data.table_1.16.2          
# [88] httr_1.4.7                  htmlwidgets_1.6.4           S4Arrays_1.4.1             
# [91] uwot_0.2.2                  pkgconfig_2.0.3             gtable_0.3.6               
# [94] blob_1.2.4                  lmtest_0.9-40               htmltools_0.5.8.1          
# [97] dotCall64_1.2               ProtGenerics_1.36.0         scales_1.3.0               
# [100] png_0.1-8                   spatstat.univar_3.0-1       knitr_1.48                 
# [103] rstudioapi_0.17.1           reshape2_1.4.4              rjson_0.2.23               
# [106] checkmate_2.3.2             nlme_3.1-165                curl_5.2.3                 
# [109] zoo_1.8-12                  cachem_1.1.0                stringr_1.5.1              
# [112] KernSmooth_2.23-24          parallel_4.4.1              miniUI_0.1.1.1             
# [115] foreign_0.8-86              restfulr_0.0.15             pillar_1.9.0               
# [118] grid_4.4.1                  vctrs_0.6.5                 RANN_2.6.2                 
# [121] promises_1.3.0              xtable_1.8-4                cluster_2.1.6              
# [124] htmlTable_2.4.3             evaluate_1.0.1              cli_3.6.3                  
# [127] compiler_4.4.1              Rsamtools_2.20.0            rlang_1.1.4                
# [130] crayon_1.5.3                future.apply_1.11.3         plyr_1.8.9                 
# [133] stringi_1.8.4               viridisLite_0.4.2           deldir_2.0-4               
# [136] BiocParallel_1.38.0         munsell_0.5.1               lazyeval_0.2.2             
# [139] spatstat.geom_3.3-3         Matrix_1.6-5                RcppHNSW_0.6.0             
# [142] patchwork_1.3.0             bit64_4.5.2                 future_1.34.0              
# [145] ggplot2_3.5.1               KEGGREST_1.44.1             shiny_1.9.1                
# [148] SummarizedExperiment_1.34.0 ROCR_1.0-11                 igraph_2.1.1               
# [151] memoise_2.0.1               fastmatch_1.1-4             bit_4.5.0