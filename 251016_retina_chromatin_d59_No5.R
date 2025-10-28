# GSE183684

library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
# library(BSgenome.Mmusculus.UCSC.mm10)
# library(patchwork)
library(ggplot2)

multiome_d59 <- readRDS("250827_retina_chromatin_d59_ATACseq_No3.rds")
multiome_d59
# An object of class Seurat
# 206672 features across 6925 samples within 2 assays
# Active assay: RNA (19607 features, 0 variable features)
# 2 layers present: counts, data
# 1 other assay present: peaks
# 2 dimensional reductions calculated: lsi, umap

p1 <- DimPlot(multiome_d59, label = TRUE, pt.size = 0.1) + NoLegend()
p1

p2 <- DimPlot(
  object = multiome_d59,
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
multiome_d59 <- AddMotifs(
  object = multiome_d59,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm,
  assay = "peaks"
)

DefaultAssay(multiome_d59) <- 'peaks'

Idents(object = multiome_d59) <- "predicted.id"

multiome_d59 <- SortIdents(multiome_d59)

# saveRDS(multiome_d59, "251003_retina_chromatin_No5.rds")
# multiome_d59 <- readRDS("251003_retina_chromatin_No5.rds")

da_peaks <- FindMarkers(
  object = multiome_d59,
  ident.1 = c("Early markers+ RPC1_d59", "Early markers+ RPC2_d59",
              "Early markers+ RPC3_d59", "MKI67+ RPC1_d59", "MKI67+ RPC2_d59"),
  assay = "peaks",
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

# saveRDS(da_peaks, "251003_retina_chromatin_No5_da_peaks.rds")
# da_peaks <- readRDS("251003_retina_chromatin_No5_da_peaks.rds")

# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005 & da_peaks$pct.1 > 0.2, ])

# test enrichment
enriched.motifs <- FindMotifs(
  object = multiome_d59,
  features = top.da.peak,
  assay = "peaks"
)

# Selecting background regions to match input sequence characteristics
# Matching GC.percent distribution
# Testing motif enrichment in 1273 regions

head(enriched.motifs)
motifs_of_interest <- c("MA0700.2", "MA0069.1", "MA0718.1", "MA0726.1")
# MA0700.2 = LHX2
# MA0069.1 = PAX6
# MA0718.1 = RAX
# MA0726.1 = VSX2

enrichment_moi <- enriched.motifs[motifs_of_interest, , drop = FALSE]
print(enrichment_moi)
# motif observed background percent.observed percent.background fold.enrichment
# MA0700.2 MA0700.2      583       2339        12.209424             5.8475        2.087973
# MA0069.1 MA0069.1      326       1689         6.827225             4.2225        1.616868
# MA0718.1 MA0718.1      518       1883        10.848168             4.7075        2.304443
# MA0726.1 MA0726.1      527       1956        11.036649             4.8900        2.256983
# pvalue motif.name     p.adjust
# MA0700.2 4.296954e-72       LHX2 8.782269e-72
# MA0069.1 3.717354e-19       PAX6 3.950351e-19
# MA0718.1 2.028112e-79        RAX 4.641017e-79
# MA0726.1 2.325237e-77       VSX2 5.162580e-77

# write.csv(enrichment_moi, "251003_d59_motif_enrichment.csv")

DefaultAssay(multiome_d59) <- "peaks"

MotifPlot(
  object = multiome_d59,
  motifs = head(rownames(enriched.motifs))
)

MotifPlot(
  object = multiome_d59,
  motifs = c("LHX2", "PAX6", "RAX", "VSX2")
)

# footprinting
library(Signac)
library(Seurat)

multiome <- readRDS("250827_retina_chromatin_d59_ATACseq_No3.rds")
multiome
# An object of class Seurat
# 206672 features across 6925 samples within 2 assays
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
# Warning message:
#   Removed 7154 rows containing missing values or values outside the scale range (`geom_label_repel()`). 
PlotFootprint(multiome, features = "MA0069.1", label.top = 5)
# Warning message:
#   Removed 7196 rows containing missing values or values outside the scale range (`geom_label_repel()`).
PlotFootprint(multiome, features = "MA0718.1", label.top = 5)
# Warning message:
#   Removed 7140 rows containing missing values or values outside the scale range (`geom_label_repel()`).
PlotFootprint(multiome, features = "MA0726.1", label.top = 5)
# Warning message:
#   Removed 7112 rows containing missing values or values outside the scale range (`geom_label_repel()`). 

# saveRDS(multiome, "251010_retina_chromatin_d59_No5_FP.rds")
multiome <- readRDS("251010_retina_chromatin_d59_No5_FP.rds")

# motifs <- Motifs(multiome)
# motifs
# # A Motif object containing 633 motifs in 187065 regions
# 
# motifs_pos <- motifs@positions
# motifs_pos
# # GRangesList object of length 633:
# #   $MA0030.1
# # GRanges object with 12732 ranges and 1 metadata column:
# #   seqnames            ranges strand |     score
# # <Rle>         <IRanges>  <Rle> | <numeric>
# #   [1]     chr1     778470-778483      + |   15.1064
# # [2]     chr1     944553-944566      + |   11.2755
# # [3]     chr1     975823-975836      + |   12.5552
# # [4]     chr1   1237592-1237605      + |   13.0563
# # [5]     chr1   1540204-1540217      + |   12.4298
# # ...      ...               ...    ... .       ...
# # [12728]     chrY   6897272-6897285      - |   12.0851
# # [12729]     chrY 12422476-12422489      + |   13.2440
# # [12730]     chrY 14247154-14247167      - |   11.2326
# # [12731]     chrY 17331643-17331656      - |   12.1101
# # [12732]     chrY 18952173-18952186      - |   14.9138
# # -------
# #   seqinfo: 24 sequences from an unspecified genome; no seqlengths
# # 
# # ...
# # <632 more elements>
# 
# # MA0700.2 = LHX2
# # MA0069.1 = PAX6
# # MA0718.1 = RAX
# # MA0726.1 = VSX2
# 
# # LHX2
# LHX2_pos <- motifs_pos[["MA0700.2"]]
# LHX2_pos
# # GRanges object with 15769 ranges and 1 metadata column:
# #   seqnames            ranges strand |     score
# # <Rle>         <IRanges>  <Rle> | <numeric>
# #   [1]     chr1     854346-854356      + |   11.8913
# # [2]     chr1     940196-940206      - |   12.3744
# # [3]     chr1     941963-941973      - |   12.2462
# # [4]     chr1   1039503-1039513      + |   12.9397
# # [5]     chr1   1039494-1039504      - |   12.2545
# # ...      ...               ...    ... .       ...
# # [15765]     chrY 13802996-13803006      - |   12.1932
# # [15766]     chrY 14748356-14748366      + |   11.8214
# # [15767]     chrY 15299042-15299052      - |   11.8913
# # [15768]     chrY 15419302-15419312      + |   11.5475
# # [15769]     chrY 19018560-19018570      + |   11.5659
# # -------
# #   seqinfo: 24 sequences from an unspecified genome; no seqlengths
# 
# # LHX2_NOTCH1
# LHX2_chr9 <- LHX2_pos[seqnames(LHX2_pos) == "chr9"]
# LHX2_chr9
# # GRanges object with 631 ranges and 1 metadata column:
# #   seqnames              ranges strand |     score
# # <Rle>           <IRanges>  <Rle> | <numeric>
# #   [1]     chr9       304448-304458      - |   13.7504
# # [2]     chr9       346565-346575      + |   13.0553
# # [3]     chr9       590044-590054      - |   11.5844
# # [4]     chr9       686589-686599      - |   11.6755
# # [5]     chr9       899508-899518      - |   12.3744
# # ...      ...                 ...    ... .       ...
# # [627]     chr9 135363937-135363947      - |   11.7794
# # [628]     chr9 135456919-135456929      + |   12.1034
# # [629]     chr9 135802884-135802894      + |   12.5259
# # [630]     chr9 136542636-136542646      - |   12.5490
# # [631]     chr9 136886184-136886194      + |   12.0710
# # -------
# #   seqinfo: 24 sequences from an unspecified genome; no seqlengths
# # write.csv(LHX2_chr9, "FootPrint_LHX2_NOTCH1_d59.csv")
# read.csv("FootPrint_LHX2_NOTCH1_d59.csv")
# 
# # peaks <- granges(multiome[["peaks"]])
# # gene_anno <- Annotation(multiome)
# # region_NOTCH1 <- gene_anno[gene_anno$gene_name == "NOTCH1"]
# # region_NOTCH1 <- Extend(region_NOTCH1, upstream = 100000, downstream = 100000)
# # 
# # LHX2_NOTCH1 <- subsetByOverlaps(LHX2_pos, region_NOTCH1)
# # 
# # fp <- Footprint(
# #   object = multiome,
# #   regions = LHX2_NOTCH1,
# #   genome = BSgenome.Hsapiens.UCSC.hg38,
# #   assay = "peaks",
# #   key = "LHX2_NOTCH1_"
# # )
# # fp
# # # An object of class Seurat 
# # # 206672 features across 6925 samples within 2 assays 
# # # Active assay: peaks (187065 features, 187065 variable features)
# # # 2 layers present: counts, data
# # # 1 other assay present: RNA
# # # 2 dimensional reductions calculated: lsi, umap
# 
# # LHX2_NOTCH3
# LHX2_chr19 <- LHX2_pos[seqnames(LHX2_pos) == "chr19"]
# LHX2_chr19
# # write.csv(LHX2_chr19, "FootPrint_LHX2_NOTCH3_d59.csv")
# read.csv("FootPrint_LHX2_NOTCH3_d59.csv")
# 
# # LHX2_NOTCH4
# LHX2_chr6 <- LHX2_pos[seqnames(LHX2_pos) == "chr6"]
# LHX2_chr6
# # write.csv(LHX2_chr6, "FootPrint_LHX2_NOTCH4_d59.csv")
# read.csv("FootPrint_LHX2_NOTCH4_d59.csv")
# 
# # PAX6
# PAX6_pos <- motifs_pos[["MA0069.1"]]
# PAX6_pos
# # PAX6_NOTCH1
# PAX6_chr9 <- PAX6_pos[seqnames(PAX6_pos) == "chr9"]
# PAX6_chr9
# # write.csv(PAX6_chr9, "FootPrint_PAX6_NOTCH1_d59.csv")
# read.csv("FootPrint_PAX6_NOTCH1_d59.csv")
# # PAX6_NOTCH3
# PAX6_chr19 <- PAX6_pos[seqnames(PAX6_pos) == "chr19"]
# PAX6_chr19
# # write.csv(PAX6_chr19, "FootPrint_PAX6_NOTCH3_d59.csv")
# read.csv("FootPrint_PAX6_NOTCH3_d59.csv")
# # PAX6_NOTCH4
# PAX6_chr6 <- PAX6_pos[seqnames(PAX6_pos) == "chr6"]
# PAX6_chr6
# # write.csv(PAX6_chr6, "FootPrint_PAX6_NOTCH4_d59.csv")
# read.csv("FootPrint_PAX6_NOTCH4_d59.csv")
# 
# # RAX
# RAX_pos <- motifs_pos[["MA0718.1"]]
# RAX_pos
# # RAX_NOTCH1
# RAX_chr9 <- RAX_pos[seqnames(RAX_pos) == "chr9"]
# RAX_chr9
# # write.csv(RAX_chr9, "FootPrint_RAX_NOTCH1_d59.csv")
# read.csv("FootPrint_RAX_NOTCH1_d59.csv")
# # RAX_NOTCH3
# RAX_chr19 <- RAX_pos[seqnames(RAX_pos) == "chr19"]
# RAX_chr19
# # write.csv(RAX_chr19, "FootPrint_RAX_NOTCH3_d59.csv")
# read.csv("FootPrint_RAX_NOTCH3_d59.csv")
# # RAX_NOTCH4
# RAX_chr6 <- RAX_pos[seqnames(RAX_pos) == "chr6"]
# RAX_chr6
# # write.csv(RAX_chr6, "FootPrint_RAX_NOTCH4_d59.csv")
# read.csv("FootPrint_RAX_NOTCH4_d59.csv")
# 
# # VSX2
# VSX2_pos <- motifs_pos[["MA0726.1"]]
# VSX2_pos
# # RAX_NOTCH1
# VSX2_chr9 <- VSX2_pos[seqnames(VSX2_pos) == "chr9"]
# VSX2_chr9
# # write.csv(VSX2_chr9, "FootPrint_VSX2_NOTCH1_d59.csv")
# read.csv("FootPrint_VSX2_NOTCH1_d59.csv")
# # VSX2_NOTCH3
# VSX2_chr19 <- VSX2_pos[seqnames(VSX2_pos) == "chr19"]
# VSX2_chr19
# # write.csv(VSX2_chr19, "FootPrint_VSX2_NOTCH3_d59.csv")
# read.csv("FootPrint_VSX2_NOTCH3_d59.csv")
# # VSX2_NOTCH4
# VSX2_chr6 <- VSX2_pos[seqnames(VSX2_pos) == "chr6"]
# VSX2_chr6
# # write.csv(VSX2_chr6, "FootPrint_VSX2_NOTCH4_d59.csv")
# read.csv("FootPrint_VSX2_NOTCH4_d59.csv")

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
#   [1] EnsDb.Hsapiens.v86_2.99.0         ensembldb_2.28.1                  AnnotationFilter_1.28.0          
# [4] GenomicFeatures_1.56.0            AnnotationDbi_1.66.0              Biobase_2.64.0                   
# [7] ggplot2_4.0.0                     BSgenome.Hsapiens.UCSC.hg38_1.4.5 BSgenome_1.72.0                  
# [10] rtracklayer_1.64.0                BiocIO_1.14.0                     Biostrings_2.72.1                
# [13] XVector_0.44.0                    GenomicRanges_1.56.2              GenomeInfoDb_1.40.1              
# [16] IRanges_2.38.1                    S4Vectors_0.42.1                  BiocGenerics_0.50.0              
# [19] TFBSTools_1.42.0                  JASPAR2020_0.99.10                motifmatchr_1.26.0               
# [22] Seurat_5.3.0                      SeuratObject_5.2.0                sp_2.2-0                         
# [25] Signac_1.15.0                    
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.22            splines_4.4.1               later_1.4.4                 bitops_1.0-9               
# [5] R.oo_1.27.1                 tibble_3.3.0                polyclip_1.10-7             XML_3.99-0.19              
# [9] DirichletMultinomial_1.46.0 fastDummies_1.7.5           lifecycle_1.0.4             pwalign_1.0.0              
# [13] globals_0.18.0              lattice_0.22-5              MASS_7.3-61                 magrittr_2.0.4             
# [17] plotly_4.11.0               yaml_2.3.10                 httpuv_1.6.16               sctransform_0.4.2          
# [21] spam_2.11-1                 spatstat.sparse_3.1-0       reticulate_1.43.0           cowplot_1.2.0              
# [25] pbapply_1.7-4               DBI_1.2.3                   CNEr_1.40.0                 RColorBrewer_1.1-3         
# [29] abind_1.4-8                 zlibbioc_1.50.0             Rtsne_0.17                  R.utils_2.13.0             
# [33] purrr_1.1.0                 RCurl_1.98-1.17             GenomeInfoDbData_1.2.12     ggrepel_0.9.6              
# [37] irlba_2.3.5.1               listenv_0.9.1               spatstat.utils_3.2-0        seqLogo_1.70.0             
# [41] goftest_1.2-3               RSpectra_0.16-2             annotate_1.82.0             spatstat.random_3.4-2      
# [45] fitdistrplus_1.2-4          parallelly_1.45.1           codetools_0.2-19            DelayedArray_0.30.1        
# [49] RcppRoll_0.3.1              tidyselect_1.2.1            UCSC.utils_1.0.0            farver_2.1.2               
# [53] matrixStats_1.5.0           spatstat.explore_3.5-3      GenomicAlignments_1.40.0    jsonlite_2.0.0             
# [57] progressr_0.16.0            ggridges_0.5.7              survival_3.7-0              tools_4.4.1                
# [61] TFMPvalue_0.0.9             ica_1.0-3                   Rcpp_1.1.0                  glue_1.8.0                 
# [65] gridExtra_2.3               SparseArray_1.4.8           MatrixGenerics_1.16.0       dplyr_1.1.4                
# [69] withr_3.0.2                 fastmap_1.2.0               caTools_1.18.3              digest_0.6.37              
# [73] R6_2.6.1                    mime_0.13                   colorspace_2.1-2            GO.db_3.19.1               
# [77] scattermore_1.2             poweRlaw_1.0.0              gtools_3.9.5                tensor_1.5.1               
# [81] dichromat_2.0-0.1           spatstat.data_3.1-8         RSQLite_2.4.3               R.methodsS3_1.8.2          
# [85] tidyr_1.3.1                 generics_0.1.4              data.table_1.17.8           httr_1.4.7                 
# [89] htmlwidgets_1.6.4           S4Arrays_1.4.1              uwot_0.2.3                  pkgconfig_2.0.3            
# [93] gtable_0.3.6                blob_1.2.4                  lmtest_0.9-40               S7_0.2.0                   
# [97] htmltools_0.5.8.1           dotCall64_1.2               ProtGenerics_1.36.0         scales_1.4.0               
# [101] png_0.1-8                   spatstat.univar_3.1-4       rstudioapi_0.17.1           tzdb_0.5.0                 
# [105] reshape2_1.4.4              rjson_0.2.23                nlme_3.1-165                curl_7.0.0                 
# [109] zoo_1.8-14                  cachem_1.1.0                stringr_1.5.2               KernSmooth_2.23-24         
# [113] parallel_4.4.1              miniUI_0.1.2                restfulr_0.0.16             pillar_1.11.1              
# [117] grid_4.4.1                  vctrs_0.6.5                 RANN_2.6.2                  promises_1.3.3             
# [121] xtable_1.8-4                cluster_2.1.6               readr_2.1.5                 cli_3.6.5                  
# [125] compiler_4.4.1              Rsamtools_2.20.0            rlang_1.1.6                 crayon_1.5.3               
# [129] future.apply_1.20.0         labeling_0.4.3              plyr_1.8.9                  stringi_1.8.7              
# [133] viridisLite_0.4.2           deldir_2.0-4                BiocParallel_1.38.0         lazyeval_0.2.2             
# [137] spatstat.geom_3.6-0         Matrix_1.6-5                RcppHNSW_0.6.0              hms_1.1.3                  
# [141] patchwork_1.3.2             bit64_4.6.0-1               future_1.67.0               KEGGREST_1.44.1            
# [145] shiny_1.11.1                SummarizedExperiment_1.34.0 ROCR_1.0-11                 igraph_2.1.4               
# [149] memoise_2.0.1               fastmatch_1.1-6             bit_4.6.0