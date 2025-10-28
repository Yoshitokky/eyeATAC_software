# GSE183684

library(Signac)
library(Seurat)
library(GenomicRanges)
library(ggplot2)
library(patchwork)

RNAseq_ATACseq_merged <- readRDS("241108_retina_chromatin_d59_No1.rds")
# This data is derived from the following code;
# 241108_retina_chromatin_d59_No1.R

str(RNAseq_ATACseq_merged)
head(RNAseq_ATACseq_merged)

RNAseq_ATACseq_merged[['peaks']]

grange <- granges(RNAseq_ATACseq_merged)

peaks.keep <- seqnames(granges(RNAseq_ATACseq_merged)) %in% standardChromosomes(granges(RNAseq_ATACseq_merged))
RNAseq_ATACseq_merged <- RNAseq_ATACseq_merged[as.vector(peaks.keep), ]

# compute nucleosome signal score per cell
RNAseq_ATACseq_merged <- NucleosomeSignal(object = RNAseq_ATACseq_merged)

# compute TSS enrichment score per cell
RNAseq_ATACseq_merged <- TSSEnrichment(object = RNAseq_ATACseq_merged)

# not run
# peak_ranges should be a set of genomic ranges spanning the set of peaks to be quantified per cell
peak_matrix <- FeatureMatrix(
  fragments = Fragments(RNAseq_ATACseq_merged),
  features = grange
)

# not run
total_fragments <- CountFragments('GSM5567518_d59_fragments.tsv.gz')
rownames(total_fragments) <- total_fragments$CB
RNAseq_ATACseq_merged$fragments <- total_fragments[colnames(RNAseq_ATACseq_merged), "frequency_count"]

RNAseq_ATACseq_merged <- FRiP(
  object = RNAseq_ATACseq_merged,
  assay = 'peaks',
  total.fragments = 'fragments'
)

# add blacklist ratio
RNAseq_ATACseq_merged$blacklist_ratio <- FractionCountsInRegion(
  object = RNAseq_ATACseq_merged, 
  assay = 'peaks',
  regions = blacklist_hg38_unified
)

DensityScatter(RNAseq_ATACseq_merged, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

RNAseq_ATACseq_merged$nucleosome_group <- ifelse(RNAseq_ATACseq_merged$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = RNAseq_ATACseq_merged, group.by = 'nucleosome_group')

VlnPlot(
  object = RNAseq_ATACseq_merged,
  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'FRiP'),
  pt.size = 0.1,
  ncol = 5
)

multiome_d59 <- subset(
  x = RNAseq_ATACseq_merged,
  subset = nCount_peaks > 2000 &
    nCount_peaks < 30000 &
    FRiP > 0.2 &
    blacklist_ratio < 0.01 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2.5
)
multiome_d59

saveRDS(object = multiome_d59, file = "241108_retina_chromatin_d59_No2.rds")

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
#   [1] patchwork_1.3.0      ggplot2_3.5.1        GenomicRanges_1.56.2 GenomeInfoDb_1.40.1 
# [5] IRanges_2.38.1       S4Vectors_0.42.1     BiocGenerics_0.50.0  Seurat_5.1.0        
# [9] SeuratObject_5.0.2   sp_2.1-4             Signac_1.14.0       
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3      rstudioapi_0.17.1       jsonlite_1.8.9         
# [4] magrittr_2.0.3          ggbeeswarm_0.7.2        spatstat.utils_3.1-0   
# [7] farver_2.1.2            zlibbioc_1.50.0         vctrs_0.6.5            
# [10] ROCR_1.0-11             spatstat.explore_3.3-3  Rsamtools_2.20.0       
# [13] RcppRoll_0.3.1          htmltools_0.5.8.1       sctransform_0.4.1      
# [16] parallelly_1.38.0       KernSmooth_2.23-24      htmlwidgets_1.6.4      
# [19] ica_1.0-3               plyr_1.8.9              plotly_4.10.4          
# [22] zoo_1.8-12              igraph_2.1.1            mime_0.12              
# [25] lifecycle_1.0.4         pkgconfig_2.0.3         Matrix_1.6-5           
# [28] R6_2.5.1                fastmap_1.2.0           GenomeInfoDbData_1.2.12
# [31] fitdistrplus_1.2-1      future_1.34.0           shiny_1.9.1            
# [34] digest_0.6.37           colorspace_2.1-1        tensor_1.5             
# [37] RSpectra_0.16-2         irlba_2.3.5.1           labeling_0.4.3         
# [40] progressr_0.15.0        fansi_1.0.6             spatstat.sparse_3.1-0  
# [43] httr_1.4.7              polyclip_1.10-7         abind_1.4-8            
# [46] compiler_4.4.1          withr_3.0.2             BiocParallel_1.38.0    
# [49] fastDummies_1.7.4       MASS_7.3-61             tools_4.4.1            
# [52] vipor_0.4.7             lmtest_0.9-40           beeswarm_0.4.0         
# [55] httpuv_1.6.15           future.apply_1.11.3     goftest_1.2-3          
# [58] glue_1.8.0              nlme_3.1-165            promises_1.3.0         
# [61] grid_4.4.1              Rtsne_0.17              cluster_2.1.6          
# [64] reshape2_1.4.4          generics_0.1.3          gtable_0.3.6           
# [67] spatstat.data_3.1-2     tidyr_1.3.1             data.table_1.16.2      
# [70] utf8_1.2.4              XVector_0.44.0          spatstat.geom_3.3-3    
# [73] RcppAnnoy_0.0.22        ggrepel_0.9.6           RANN_2.6.2             
# [76] pillar_1.9.0            stringr_1.5.1           spam_2.11-0            
# [79] RcppHNSW_0.6.0          later_1.3.2             splines_4.4.1          
# [82] dplyr_1.1.4             lattice_0.22-5          survival_3.7-0         
# [85] deldir_2.0-4            tidyselect_1.2.1        Biostrings_2.72.1      
# [88] miniUI_0.1.1.1          pbapply_1.7-2           gridExtra_2.3          
# [91] scattermore_1.2         matrixStats_1.4.1       stringi_1.8.4          
# [94] UCSC.utils_1.0.0        lazyeval_0.2.2          codetools_0.2-19       
# [97] tibble_3.2.1            cli_3.6.3               uwot_0.2.2             
# [100] xtable_1.8-4            reticulate_1.39.0       munsell_0.5.1          
# [103] Rcpp_1.0.13             globals_0.16.3          spatstat.random_3.3-2  
# [106] png_0.1-8               ggrastr_1.0.2           spatstat.univar_3.0-1  
# [109] parallel_4.4.1          dotCall64_1.2           bitops_1.0-9           
# [112] listenv_0.9.1           viridisLite_0.4.2       scales_1.3.0           
# [115] ggridges_0.5.6          leiden_0.4.3.1          purrr_1.0.2            
# [118] crayon_1.5.3            rlang_1.1.4             cowplot_1.1.3          
# [121] fastmatch_1.1-4