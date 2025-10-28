# GSE183684

library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
# library(BSgenome.Mmusculus.UCSC.mm10)
# library(patchwork)
library(ggplot2)

multiome_d78 <- readRDS("251007_retina_chromatin_d78_ATACseq_No3.rds")
multiome_d78
# An object of class Seurat 
# 189768 features across 11468 samples within 2 assays 
# Active assay: RNA (19607 features, 0 variable features)
# 2 layers present: counts, data
# 1 other assay present: peaks
# 2 dimensional reductions calculated: lsi, umap

p1 <- DimPlot(multiome_d78, label = TRUE, pt.size = 0.1) + NoLegend()
p1

p2 <- DimPlot(
  object = multiome_d78,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq') + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20))
p2

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
multiome_d78 <- AddMotifs(
  object = multiome_d78,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm,
  assay = "peaks"
)

DefaultAssay(multiome_d78) <- 'peaks'

Idents(object = multiome_d78) <- "predicted.id"

multiome_d78 <- SortIdents(multiome_d78)

# saveRDS(multiome_d78, "251007_retina_chromatin_d78_No5.rds")
# multiome_d78 <- readRDS("251007_retina_chromatin_d78_No5.rds")

da_peaks <- FindMarkers(
  object = multiome_d78,
  ident.1 = c("Early markers+ RPC1_d78", "Early markers+ RPC2_d78",
              "Early markers+ RPC3_d78", "MKI67+ RPC1_d78", "MKI67+ RPC2_d78"),
  assay = "peaks",
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

# saveRDS(da_peaks, "251007_retina_chromatin_d78_No5_da_peaks.rds")
# da_peaks <- readRDS("251007_retina_chromatin_d78_No5_da_peaks.rds")

# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005 & da_peaks$pct.1 > 0.2, ])

# test enrichment
enriched.motifs <- FindMotifs(
  object = multiome_d78,
  features = top.da.peak,
  assay = "peaks"
)

# Selecting background regions to match input sequence characteristics
# Matching GC.percent distribution
# Testing motif enrichment in 1132 regions

head(enriched.motifs)
motifs_of_interest <- c("MA0700.2", "MA0069.1", "MA0718.1", "MA0726.1")
# MA0700.2 = LHX2
# MA0069.1 = PAX6
# MA0718.1 = RAX
# MA0726.1 = VSX2

enrichment_moi <- enriched.motifs[motifs_of_interest, , drop = FALSE]
print(enrichment_moi)
# motif observed background percent.observed percent.background fold.enrichment
# MA0700.2 MA0700.2      150       2147         13.25088             5.3675        2.468725
# MA0069.1 MA0069.1       77       1626          6.80212             4.0650        1.673338
# MA0718.1 MA0718.1      130       1994         11.48410             4.9850        2.303731
# MA0726.1 MA0726.1      131       1876         11.57244             4.6900        2.467471
# pvalue motif.name     p.adjust
# MA0700.2 1.209630e-24       LHX2 3.418120e-24
# MA0069.1 8.705281e-06       PAX6 9.237752e-06
# MA0718.1 7.315325e-19        RAX 1.421154e-18
# MA0726.1 1.597430e-21       VSX2 3.724009e-21

# write.csv(enrichment_moi, "251007_d78_motif_enrichment.csv")

DefaultAssay(multiome_d78) <- "peaks"

MotifPlot(
  object = multiome_d78,
  motifs = head(rownames(enriched.motifs))
)

MotifPlot(
  object = multiome_d78,
  motifs = c("LHX2", "PAX6", "RAX", "VSX2")
)

# footprinting
library(Signac)
library(Seurat)

multiome <- readRDS("251007_retina_chromatin_d78_ATACseq_No3.rds")
multiome
# An object of class Seurat 
# 189768 features across 11468 samples within 2 assays 
# Active assay: RNA (19607 features, 0 variable features)
# 2 layers present: counts, data
# 1 other assay present: peaks
# 2 dimensional reductions calculated: lsi, umap

DefaultAssay(multiome) <- 'peaks'
Idents(object = multiome) <- "predicted.id"
multiome <- SortIdents(multiome)

DimPlot(multiome, label = TRUE)

library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(EnsDb.Hsapiens.v86)

# extract position frequency matrices for the motifs
pwm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

# add motif information
multiome <- AddMotifs(multiome, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pwm)

# gather the footprinting information for sets of motifs
multiome <- Footprint(
  object = multiome,
  motif.name = c("MA0700.2", "MA0069.1", "MA0718.1", "MA0726.1"),
  genome = BSgenome.Hsapiens.UCSC.hg38
)
# MA0700.2 = LHX2
# MA0069.1 = PAX6
# MA0718.1 = RAX
# MA0726.1 = VSX2

# plot the footprint data for each group of cells
p2 <- PlotFootprint(multiome, features = c("MA0700.2", "MA0069.1", "MA0718.1", "MA0726.1"))
p2 + patchwork::plot_layout(ncol = 1)

PlotFootprint(multiome, features = "MA0700.2", label.top = 5)
# Warning messages:
#   1: Removed 1 row containing missing values or values outside the scale range (`geom_line()`). 
# 2: Removed 8687 rows containing missing values or values outside the scale range
# (`geom_label_repel()`). 
PlotFootprint(multiome, features = "MA0069.1", label.top = 5)
# Warning message:
#   Removed 8738 rows containing missing values or values outside the scale range
# (`geom_label_repel()`). 
PlotFootprint(multiome, features = "MA0718.1", label.top = 5)
# Warning message:
#   Removed 8670 rows containing missing values or values outside the scale range
# (`geom_label_repel()`). 
PlotFootprint(multiome, features = "MA0726.1", label.top = 5)
# Warning message:
#   Removed 8636 rows containing missing values or values outside the scale range
# (`geom_label_repel()`). 

# saveRDS(multiome, "251010_retina_chromatin_d78_No5_FP.rds")
multiome <- readRDS("251010_retina_chromatin_d78_No5_FP.rds")

# motifs <- Motifs(multiome)
# motifs
# # A Motif object containing 633 motifs in 170161 regions
# 
# motifs_pos <- motifs@positions
# motifs_pos
# 
# # MA0700.2 = LHX2
# # MA0069.1 = PAX6
# # MA0718.1 = RAX
# # MA0726.1 = VSX2
# 
# # LHX2
# LHX2_pos <- motifs_pos[["MA0700.2"]]
# LHX2_pos
# # LHX2_NOTCH1
# LHX2_chr9 <- LHX2_pos[seqnames(LHX2_pos) == "chr9"]
# LHX2_chr9
# # write.csv(LHX2_chr9, "FootPrint_LHX2_NOTCH1_d78.csv")
# read.csv("FootPrint_LHX2_NOTCH1_d78.csv")
# # LHX2_NOTCH3
# LHX2_chr19 <- LHX2_pos[seqnames(LHX2_pos) == "chr19"]
# LHX2_chr19
# # write.csv(LHX2_chr19, "FootPrint_LHX2_NOTCH3_d78.csv")
# read.csv("FootPrint_LHX2_NOTCH3_d78.csv")
# # LHX2_NOTCH4
# LHX2_chr6 <- LHX2_pos[seqnames(LHX2_pos) == "chr6"]
# LHX2_chr6
# # write.csv(LHX2_chr6, "FootPrint_LHX2_NOTCH4_d78.csv")
# read.csv("FootPrint_LHX2_NOTCH4_d78.csv")
# # PAX6
# PAX6_pos <- motifs_pos[["MA0069.1"]]
# PAX6_pos
# # PAX6_NOTCH1
# PAX6_chr9 <- PAX6_pos[seqnames(PAX6_pos) == "chr9"]
# PAX6_chr9
# # write.csv(PAX6_chr9, "FootPrint_PAX6_NOTCH1_d78.csv")
# read.csv("FootPrint_PAX6_NOTCH1_d78.csv")
# # PAX6_NOTCH3
# PAX6_chr19 <- PAX6_pos[seqnames(PAX6_pos) == "chr19"]
# PAX6_chr19
# # write.csv(PAX6_chr19, "FootPrint_PAX6_NOTCH3_d78.csv")
# read.csv("FootPrint_PAX6_NOTCH3_d78.csv")
# # PAX6_NOTCH4
# PAX6_chr6 <- PAX6_pos[seqnames(PAX6_pos) == "chr6"]
# PAX6_chr6
# # write.csv(PAX6_chr6, "FootPrint_PAX6_NOTCH4_d78.csv")
# read.csv("FootPrint_PAX6_NOTCH4_d78.csv")
# 
# # RAX
# RAX_pos <- motifs_pos[["MA0718.1"]]
# RAX_pos
# # RAX_NOTCH1
# RAX_chr9 <- RAX_pos[seqnames(RAX_pos) == "chr9"]
# RAX_chr9
# # write.csv(RAX_chr9, "FootPrint_RAX_NOTCH1_d78.csv")
# read.csv("FootPrint_RAX_NOTCH1_d78.csv")
# # RAX_NOTCH3
# RAX_chr19 <- RAX_pos[seqnames(RAX_pos) == "chr19"]
# RAX_chr19
# # write.csv(RAX_chr19, "FootPrint_RAX_NOTCH3_d78.csv")
# read.csv("FootPrint_RAX_NOTCH3_d78.csv")
# # RAX_NOTCH4
# RAX_chr6 <- RAX_pos[seqnames(RAX_pos) == "chr6"]
# RAX_chr6
# # write.csv(RAX_chr6, "FootPrint_RAX_NOTCH4_d78.csv")
# read.csv("FootPrint_RAX_NOTCH4_d78.csv")
# 
# # VSX2
# VSX2_pos <- motifs_pos[["MA0726.1"]]
# VSX2_pos
# # RAX_NOTCH1
# VSX2_chr9 <- VSX2_pos[seqnames(VSX2_pos) == "chr9"]
# VSX2_chr9
# # write.csv(VSX2_chr9, "FootPrint_VSX2_NOTCH1_d78.csv")
# read.csv("FootPrint_VSX2_NOTCH1_d78.csv")
# # VSX2_NOTCH3
# VSX2_chr19 <- VSX2_pos[seqnames(VSX2_pos) == "chr19"]
# VSX2_chr19
# # write.csv(VSX2_chr19, "FootPrint_VSX2_NOTCH3_d78.csv")
# read.csv("FootPrint_VSX2_NOTCH3_d78.csv")
# # VSX2_NOTCH4
# VSX2_chr6 <- VSX2_pos[seqnames(VSX2_pos) == "chr6"]
# VSX2_chr6
# # write.csv(VSX2_chr6, "FootPrint_VSX2_NOTCH4_d78.csv")
# read.csv("FootPrint_VSX2_NOTCH4_d78.csv")

sessionInfo()
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
#   [1] EnsDb.Hsapiens.v86_2.99.0         ensembldb_2.28.1                 
# [3] AnnotationFilter_1.28.0           GenomicFeatures_1.56.0           
# [5] AnnotationDbi_1.66.0              Biobase_2.64.0                   
# [7] ggplot2_4.0.0                     BSgenome.Hsapiens.UCSC.hg38_1.4.5
# [9] BSgenome_1.72.0                   rtracklayer_1.64.0               
# [11] BiocIO_1.14.0                     Biostrings_2.72.1                
# [13] XVector_0.44.0                    GenomicRanges_1.56.2             
# [15] GenomeInfoDb_1.40.1               IRanges_2.38.1                   
# [17] S4Vectors_0.42.1                  BiocGenerics_0.50.0              
# [19] TFBSTools_1.42.0                  JASPAR2020_0.99.10               
# [21] motifmatchr_1.26.0                Seurat_5.3.0                     
# [23] SeuratObject_5.2.0                sp_2.2-0                         
# [25] Signac_1.15.0                    
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.22            splines_4.4.1               later_1.4.4                
# [4] bitops_1.0-9                R.oo_1.27.1                 tibble_3.3.0               
# [7] polyclip_1.10-7             DirichletMultinomial_1.46.0 XML_3.99-0.19              
# [10] fastDummies_1.7.5           lifecycle_1.0.4             pwalign_1.0.0              
# [13] globals_0.18.0              lattice_0.22-5              MASS_7.3-61                
# [16] magrittr_2.0.4              plotly_4.11.0               yaml_2.3.10                
# [19] httpuv_1.6.16               sctransform_0.4.2           spam_2.11-1                
# [22] spatstat.sparse_3.1-0       reticulate_1.43.0           CNEr_1.40.0                
# [25] cowplot_1.2.0               pbapply_1.7-4               DBI_1.2.3                  
# [28] RColorBrewer_1.1-3          abind_1.4-8                 zlibbioc_1.50.0            
# [31] Rtsne_0.17                  R.utils_2.13.0              purrr_1.1.0                
# [34] RCurl_1.98-1.17             GenomeInfoDbData_1.2.12     ggrepel_0.9.6              
# [37] irlba_2.3.5.1               listenv_0.9.1               spatstat.utils_3.2-0       
# [40] seqLogo_1.70.0              goftest_1.2-3               RSpectra_0.16-2            
# [43] annotate_1.82.0             spatstat.random_3.4-2       fitdistrplus_1.2-4         
# [46] parallelly_1.45.1           codetools_0.2-19            DelayedArray_0.30.1        
# [49] RcppRoll_0.3.1              tidyselect_1.2.1            UCSC.utils_1.0.0           
# [52] farver_2.1.2                matrixStats_1.5.0           spatstat.explore_3.5-3     
# [55] GenomicAlignments_1.40.0    jsonlite_2.0.0              progressr_0.16.0           
# [58] ggridges_0.5.7              survival_3.7-0              tools_4.4.1                
# [61] TFMPvalue_0.0.9             ica_1.0-3                   Rcpp_1.1.0                 
# [64] glue_1.8.0                  gridExtra_2.3               SparseArray_1.4.8          
# [67] MatrixGenerics_1.16.0       dplyr_1.1.4                 withr_3.0.2                
# [70] fastmap_1.2.0               caTools_1.18.3              digest_0.6.37              
# [73] R6_2.6.1                    mime_0.13                   colorspace_2.1-2           
# [76] GO.db_3.19.1                scattermore_1.2             poweRlaw_1.0.0             
# [79] gtools_3.9.5                tensor_1.5.1                dichromat_2.0-0.1          
# [82] spatstat.data_3.1-8         RSQLite_2.4.3               R.methodsS3_1.8.2          
# [85] tidyr_1.3.1                 generics_0.1.4              data.table_1.17.8          
# [88] httr_1.4.7                  htmlwidgets_1.6.4           S4Arrays_1.4.1             
# [91] uwot_0.2.3                  pkgconfig_2.0.3             gtable_0.3.6               
# [94] blob_1.2.4                  lmtest_0.9-40               S7_0.2.0                   
# [97] htmltools_0.5.8.1           dotCall64_1.2               ProtGenerics_1.36.0        
# [100] scales_1.4.0                png_0.1-8                   spatstat.univar_3.1-4      
# [103] rstudioapi_0.17.1           tzdb_0.5.0                  reshape2_1.4.4             
# [106] rjson_0.2.23                nlme_3.1-165                curl_7.0.0                 
# [109] zoo_1.8-14                  cachem_1.1.0                stringr_1.5.2              
# [112] KernSmooth_2.23-24          parallel_4.4.1              miniUI_0.1.2               
# [115] restfulr_0.0.16             pillar_1.11.1               grid_4.4.1                 
# [118] vctrs_0.6.5                 RANN_2.6.2                  promises_1.3.3             
# [121] xtable_1.8-4                cluster_2.1.6               readr_2.1.5                
# [124] cli_3.6.5                   compiler_4.4.1              Rsamtools_2.20.0           
# [127] rlang_1.1.6                 crayon_1.5.3                future.apply_1.20.0        
# [130] labeling_0.4.3              plyr_1.8.9                  stringi_1.8.7              
# [133] viridisLite_0.4.2           deldir_2.0-4                BiocParallel_1.38.0        
# [136] lazyeval_0.2.2              spatstat.geom_3.6-0         Matrix_1.6-5               
# [139] RcppHNSW_0.6.0              hms_1.1.3                   patchwork_1.3.2            
# [142] bit64_4.6.0-1               future_1.67.0               KEGGREST_1.44.1            
# [145] shiny_1.11.1                SummarizedExperiment_1.34.0 ROCR_1.0-11                
# [148] igraph_2.1.4                memoise_2.0.1               fastmatch_1.1-6            
# [151] bit_4.6.0