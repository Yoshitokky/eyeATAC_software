# GSE183684

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

D59 <- ReadMtx(
 mtx = "GSM5567526_d59_matrix.mtx.gz",
 features = "GSM5567526_d59_features.tsv.gz",
 cells = "GSM5567526_d59_barcodes.tsv.gz",
 feature.column = 2
)

RNA_d59 <- CreateSeuratObject(counts = D59, min.cells = 3, min.features = 200)

RNA_d59[["percent.mt"]] <- PercentageFeatureSet(RNA_d59, pattern = "^MT-")

VlnPlot(RNA_d59, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, layer = "counts")

RNA_d59 <- subset(RNA_d59, subset = percent.mt < 20)

RNA_d59 <- NormalizeData(RNA_d59, normalization.method = "LogNormalize", scale.factor = 10000)

RNA_d59 <- FindVariableFeatures(RNA_d59, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(RNA_d59)
RNA_d59 <- ScaleData(RNA_d59, features = all.genes)

RNA_d59 <- RunPCA(RNA_d59, features = VariableFeatures(object = RNA_d59))

ElbowPlot(RNA_d59)

RNA_d59 <- FindNeighbors(RNA_d59, dims = 1:15)

# RNA_d59_04 <- FindClusters(RNA_d59, resolution = 0.4)
# RNA_d59_04 <- RunTSNE(RNA_d59_04, dims = 1:10)
# DimPlot(RNA_d59_04, reduction = "tsne")

# RNA_d59_05 <- FindClusters(RNA_d59, resolution = 0.5)
# RNA_d59_05 <- RunTSNE(RNA_d59_05, dims = 1:10)
# DimPlot(RNA_d59_05, reduction = "tsne", label = TRUE)

# RNA_d59_07 <- FindClusters(RNA_d59, resolution = 0.7)
# RNA_d59_07 <- RunTSNE(RNA_d59_07, dims = 1:10)
# DimPlot(RNA_d59_07, reduction = "tsne")

# RNA_d59_09 <- FindClusters(RNA_d59, resolution = 0.9)
# RNA_d59_09 <- RunTSNE(RNA_d59_09, dims = 1:10)
# DimPlot(RNA_d59_09, reduction = "tsne")

RNA_d59_05_UMAP <- FindClusters(RNA_d59, resolution = 0.5)
RNA_d59_05_UMAP <- RunUMAP(RNA_d59_05_UMAP, dims = 1:15)
DimPlot(RNA_d59_05_UMAP, reduction = "umap")

# RNA_d59_07_UMAP <- FindClusters(RNA_d59, resolution = 0.7)
# RNA_d59_07_UMAP <- RunUMAP(RNA_d59_07_UMAP, dims = 1:15)
# DimPlot(RNA_d59_07_UMAP, reduction = "umap")

saveRDS(RNA_d59_05_UMAP, file = "241108_GSE183684_D59_No0.rds")
# RNA_d59_05_UMAP <- readRDS("241108_GSE183684_D59_No0.rds")

# # finding cluster markers 
# cluster0.markers <- FindMarkers(RNA_d59_05_UMAP, ident.1 = 0, logfc.threshold = 1.00, test.use = "roc", only.pos = TRUE)
# cluster0.markers

# DotPlot(RNA_d59_05, features = "EEF1A1")
# # eukaryotic translation elongation factor 1 alpha 1
# # according to NIH data
# FeaturePlot(RNA_d59_05, features = "EEF1A1")

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
# cluster13.markers <- FindMarkers(RNA_d59_05, ident.1 = 13, logfc.threshold = 1.00, test.use = "roc", only.pos = TRUE)
# cluster13.markers
# PAX2, Glial cells, PMID: 20437530

DotPlot(RNA_d59_05_UMAP, features = anno_markers)
# FeaturePlot(RNA_d59_05, features = anno_markers)

new.cluster.ids <- c("EEF1A1+ RPC_d59", "Early markers+ RPC2_d59", "ONECUT1+/2+/LHX1+/PTF1A+ RPC_d59", 
                     "ELAVL2+/4+ RPC1_d59", "MKI67+ RPC1_d59", "ELAVL2+/4+ RPC4_d59", 
                     "ELAVL2+/4+ RPC2_d59", "ONECUT1+/MEIS2+ RPC_d59",
                     "PAX6+/RAX+ RPC_d59", "ELAVL2+/4+ RPC3_d59", "Early markers+ RPC1_d59", "Early markers+ RPC3_d59", 
                     "MKI67+ RPC2_d59", "GLIA_d59")
names(new.cluster.ids) <- levels(RNA_d59_05_UMAP)
RNA_d59_05_anno <- RenameIdents(RNA_d59_05_UMAP, new.cluster.ids)
DimPlot(RNA_d59_05_anno, reduction = "umap", label = TRUE, repel = TRUE) + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20))

FeaturePlot(RNA_d59_05_anno, features = anno_markers)

DotPlot(RNA_d59_05_anno, features = anno_markers) + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20)) +
  RotatedAxis() +
  scale_y_discrete(
    limits = c("Early markers+ RPC1_d59", "Early markers+ RPC2_d59", "Early markers+ RPC3_d59",  "MKI67+ RPC1_d59",
               "MKI67+ RPC2_d59", "ONECUT1+/2+/LHX1+/PTF1A+ RPC_d59", "ONECUT1+/MEIS2+ RPC_d59", "PAX6+/RAX+ RPC_d59",
               "ELAVL2+/4+ RPC1_d59", "ELAVL2+/4+ RPC2_d59", "ELAVL2+/4+ RPC3_d59",
               "ELAVL2+/4+ RPC4_d59", "EEF1A1+ RPC_d59", "GLIA_d59")
  )

FeaturePlot(RNA_d59_05_anno, features = "NOTCH1") + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20))
FeaturePlot(RNA_d59_05_anno, features = "NOTCH2") + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20))
FeaturePlot(RNA_d59_05_anno, features = "NOTCH3") + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20))
FeaturePlot(RNA_d59_05_anno, features = "NOTCH4") + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20))

notch <- c("NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4")
DotPlot(RNA_d59_05_anno, features = notch) + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20)) +
  RotatedAxis() +
  scale_y_discrete(
    limits = c("Early markers+ RPC1_d59", "Early markers+ RPC2_d59", "Early markers+ RPC3_d59",  "MKI67+ RPC1_d59",
               "MKI67+ RPC2_d59", "ONECUT1+/2+/LHX1+/PTF1A+ RPC_d59", "ONECUT1+/MEIS2+ RPC_d59", "PAX6+/RAX+ RPC_d59",
               "ELAVL2+/4+ RPC1_d59", "ELAVL2+/4+ RPC2_d59", "ELAVL2+/4+ RPC3_d59",
               "ELAVL2+/4+ RPC4_d59", "EEF1A1+ RPC_d59", "GLIA_d59")
  )

library(SeuratWrappers)
library(monocle3)

# the next line is necessary if the data are bound together
# whole <- JoinLayers(whole)
cds <- as.cell_data_set(RNA_d59_05_anno)

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

library(viridis)

FeaturePlot(integrated.sub, "monocle3_pseudotime") + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20)) +
  scale_color_viridis(option = "C") 