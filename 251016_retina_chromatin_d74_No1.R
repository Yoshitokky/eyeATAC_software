library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

# in terminal

# gzip -d GSM5567519_d74_fragments.tsv.gz
# bgzip GSM5567519_d74_fragments.tsv
# tabix -p bed GSM5567519_d74_fragments.tsv.gz

RNA_d74 <- ReadMtx(mtx = "GSM5567527_d74_matrix.mtx.gz",
                   cells = "GSM5567527_d74_barcodes.tsv.gz",
                   features = "GSM5567527_d74_features.tsv.gz")

data_d74 <- CreateSeuratObject(counts = RNA_d74, assay = "RNA")

fragpath <- 'GSM5567519_d74_fragments.tsv.gz'

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

data_d74 <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks"
)

saveRDS(object = data_d74, file = "241108_retina_chromatin_d74_No1.rds")

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
#   [1] RcppAnnoy_0.0.22            splines_4.4.1              
# [3] later_1.3.2                 bitops_1.0-9               
# [5] tibble_3.2.1                polyclip_1.10-7            
# [7] rpart_4.1.23                XML_3.99-0.17              
# [9] fastDummies_1.7.4           lifecycle_1.0.4            
# [11] globals_0.16.3              lattice_0.22-5             
# [13] MASS_7.3-61                 backports_1.5.0            
# [15] magrittr_2.0.3              rmarkdown_2.28             
# [17] Hmisc_5.2-0                 plotly_4.10.4              
# [19] yaml_2.3.10                 httpuv_1.6.15              
# [21] sctransform_0.4.1           spam_2.11-0                
# [23] spatstat.sparse_3.1-0       reticulate_1.39.0          
# [25] cowplot_1.1.3               pbapply_1.7-2              
# [27] DBI_1.2.3                   RColorBrewer_1.1-3         
# [29] abind_1.4-8                 zlibbioc_1.50.0            
# [31] Rtsne_0.17                  purrr_1.0.2                
# [33] biovizBase_1.52.0           RCurl_1.98-1.16            
# [35] nnet_7.3-19                 VariantAnnotation_1.50.0   
# [37] GenomeInfoDbData_1.2.12     ggrepel_0.9.6              
# [39] irlba_2.3.5.1               listenv_0.9.1              
# [41] spatstat.utils_3.1-0        goftest_1.2-3              
# [43] RSpectra_0.16-2             spatstat.random_3.3-2      
# [45] fitdistrplus_1.2-1          parallelly_1.38.0          
# [47] leiden_0.4.3.1              codetools_0.2-19           
# [49] DelayedArray_0.30.1         RcppRoll_0.3.1             
# [51] tidyselect_1.2.1            UCSC.utils_1.0.0           
# [53] farver_2.1.2                base64enc_0.1-3            
# [55] matrixStats_1.4.1           spatstat.explore_3.3-3     
# [57] GenomicAlignments_1.40.0    jsonlite_1.8.9             
# [59] Formula_1.2-5               progressr_0.15.0           
# [61] ggridges_0.5.6              survival_3.7-0             
# [63] tools_4.4.1                 ica_1.0-3                  
# [65] Rcpp_1.0.13                 glue_1.8.0                 
# [67] gridExtra_2.3               SparseArray_1.4.8          
# [69] xfun_0.48                   MatrixGenerics_1.16.0      
# [71] dplyr_1.1.4                 withr_3.0.2                
# [73] fastmap_1.2.0               fansi_1.0.6                
# [75] digest_0.6.37               R6_2.5.1                   
# [77] mime_0.12                   colorspace_2.1-1           
# [79] scattermore_1.2             tensor_1.5                 
# [81] dichromat_2.0-0.1           spatstat.data_3.1-2        
# [83] RSQLite_2.3.7               utf8_1.2.4                 
# [85] tidyr_1.3.1                 generics_0.1.3             
# [87] data.table_1.16.2           httr_1.4.7                 
# [89] htmlwidgets_1.6.4           S4Arrays_1.4.1             
# [91] uwot_0.2.2                  pkgconfig_2.0.3            
# [93] gtable_0.3.6                blob_1.2.4                 
# [95] lmtest_0.9-40               htmltools_0.5.8.1          
# [97] dotCall64_1.2               ProtGenerics_1.36.0        
# [99] scales_1.3.0                png_0.1-8                  
# [101] spatstat.univar_3.0-1       knitr_1.48                 
# [103] rstudioapi_0.17.1           reshape2_1.4.4             
# [105] rjson_0.2.23                checkmate_2.3.2            
# [107] nlme_3.1-165                curl_5.2.3                 
# [109] zoo_1.8-12                  cachem_1.1.0               
# [111] stringr_1.5.1               KernSmooth_2.23-24         
# [113] parallel_4.4.1              miniUI_0.1.1.1             
# [115] foreign_0.8-86              restfulr_0.0.15            
# [117] pillar_1.9.0                grid_4.4.1                 
# [119] vctrs_0.6.5                 RANN_2.6.2                 
# [121] promises_1.3.0              xtable_1.8-4               
# [123] cluster_2.1.6               htmlTable_2.4.3            
# [125] evaluate_1.0.1              cli_3.6.3                  
# [127] compiler_4.4.1              Rsamtools_2.20.0           
# [129] rlang_1.1.4                 crayon_1.5.3               
# [131] future.apply_1.11.3         plyr_1.8.9                 
# [133] stringi_1.8.4               viridisLite_0.4.2          
# [135] deldir_2.0-4                BiocParallel_1.38.0        
# [137] munsell_0.5.1               lazyeval_0.2.2             
# [139] spatstat.geom_3.3-3         Matrix_1.6-5               
# [141] RcppHNSW_0.6.0              patchwork_1.3.0            
# [143] bit64_4.5.2                 future_1.34.0              
# [145] ggplot2_3.5.1               KEGGREST_1.44.1            
# [147] shiny_1.9.1                 SummarizedExperiment_1.34.0
# [149] ROCR_1.0-11                 igraph_2.1.1               
# [151] memoise_2.0.1               fastmatch_1.1-4            
# [153] bit_4.5.0