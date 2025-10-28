# GSE183684

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

D74 <- ReadMtx(
  mtx = "GSM5567527_d74_matrix.mtx.gz",
  features = "GSM5567527_d74_features.tsv.gz",
  cells = "GSM5567527_d74_barcodes.tsv.gz",
  feature.column = 2
)

RNA_d74 <- CreateSeuratObject(counts = D74, min.cells = 3, min.features = 200)

RNA_d74[["percent.mt"]] <- PercentageFeatureSet(RNA_d74, pattern = "^MT-")

VlnPlot(RNA_d74, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, layer = "counts")

RNA_d74 <- subset(RNA_d74, subset = nFeature_RNA > 1250 & nFeature_RNA < 6000 & percent.mt < 10)

RNA_d74 <- NormalizeData(RNA_d74, normalization.method = "LogNormalize", scale.factor = 10000)

RNA_d74 <- FindVariableFeatures(RNA_d74, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(RNA_d74)
RNA_d74 <- ScaleData(RNA_d74, features = all.genes)

RNA_d74 <- RunPCA(RNA_d74, features = VariableFeatures(object = RNA_d74))

ElbowPlot(RNA_d74)

RNA_d74 <- FindNeighbors(RNA_d74, dims = 1:15)

# RNA_d74_05_UMAP <- FindClusters(RNA_d74, resolution = 0.5)
# RNA_d74_05_UMAP <- RunUMAP(RNA_d74_05_UMAP, dims = 1:15)
# DimPlot(RNA_d74_05_UMAP, reduction = "umap")

RNA_d74_07_UMAP <- FindClusters(RNA_d74, resolution = 0.7)
RNA_d74_07_UMAP <- RunUMAP(RNA_d74_07_UMAP, dims = 1:15)
DimPlot(RNA_d74_07_UMAP, reduction = "umap", label =TRUE)

# RNA_d74_08_UMAP <- FindClusters(RNA_d74, resolution = 0.8)
# RNA_d74_08_UMAP <- RunUMAP(RNA_d74_08_UMAP, dims = 1:15)
# DimPlot(RNA_d74_08_UMAP, reduction = "umap", label =TRUE)

# RNA_d74_09_UMAP <- FindClusters(RNA_d74, resolution = 0.9)
# RNA_d74_09_UMAP <- RunUMAP(RNA_d74_09_UMAP, dims = 1:15)
# DimPlot(RNA_d74_09_UMAP, reduction = "umap", label =TRUE)

# saveRDS(RNA_d74_07_UMAP, file = "241108_GSE183684_D74.rds")
RNA_d74_07_UMAP <- readRDS("241108_GSE183684_D74.rds")

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

DotPlot(RNA_d74_07_UMAP, features = anno_markers)

new.cluster.ids <- c("Early markers+ RPC1_d74", "ELAVL2+/4+ RPC2_d74", "Early markers+ RPC2_d74", 
                     "Early markers+ RPC3_d74", "ELAVL2+/4+ RPC3_d74",
                     "ONECUT1+/2+/LHX1+ RPC_d74", 
                     "ELAVL2+/4+ RPC1_d74", "ONECUT1+/MEIS2+ RPC2_d74",
                     "Early markers+ RPC4_d74", "ONECUT1+/2+/PTF1A+ RPC_d74", "ONECUT1+/MEIS2+ RPC1_d74", 
                     "MKI67+ RPC_d74", "ONECUT1+/MEIS2+ RPC3_d74", "VSX2++ RPC_d74", "EEF1A1+ RPC_d74")
names(new.cluster.ids) <- levels(RNA_d74_07_UMAP)
RNA_d74_07_anno <- RenameIdents(RNA_d74_07_UMAP, new.cluster.ids)
DotPlot(RNA_d74_07_anno, features = anno_markers) + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20)) +
  RotatedAxis() +
  scale_y_discrete(
    limits = c("Early markers+ RPC1_d74", "Early markers+ RPC2_d74",
               "Early markers+ RPC3_d74",  "Early markers+ RPC4_d74",
               "MKI67+ RPC_d74", "ONECUT1+/2+/PTF1A+ RPC_d74",
               "ONECUT1+/2+/LHX1+ RPC_d74", "ONECUT1+/MEIS2+ RPC1_d74",
               "ONECUT1+/MEIS2+ RPC2_d74", "ONECUT1+/MEIS2+ RPC3_d74", "VSX2++ RPC_d74",
               "ELAVL2+/4+ RPC1_d74", "ELAVL2+/4+ RPC2_d74",
               "ELAVL2+/4+ RPC3_d74", "EEF1A1+ RPC_d74")
  )

# FeaturePlot(RNA_d74_07_anno, features = "MEIS2")
# PMID 30018341, RGC marker

# FeaturePlot(RNA_d74_07_anno, features = "LHX1")
# FeaturePlot(RNA_d74_07_anno, features = "PTF1A")
# PMID 27486389, horizontal cell marker

FeaturePlot(RNA_d74_07_anno, features = anno_markers)

DimPlot(RNA_d74_07_anno, reduction = "umap", label =TRUE, repel = TRUE)  + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20))

FeaturePlot(RNA_d74_07_anno, features = "NOTCH1") + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20))
FeaturePlot(RNA_d74_07_anno, features = "NOTCH2") + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20))
FeaturePlot(RNA_d74_07_anno, features = "NOTCH3") + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20))
FeaturePlot(RNA_d74_07_anno, features = "NOTCH4") + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20))

notch <- c("NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4")
DotPlot(RNA_d74_07_anno, features = notch) + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20)) +
  RotatedAxis() +
  scale_y_discrete(
    limits = c("Early markers+ RPC1_d74", "Early markers+ RPC2_d74",
               "Early markers+ RPC3_d74",  "Early markers+ RPC4_d74",
               "MKI67+ RPC_d74", "ONECUT1+/2+/PTF1A+ RPC_d74",
               "ONECUT1+/2+/LHX1+ RPC_d74", "ONECUT1+/MEIS2+ RPC1_d74",
               "ONECUT1+/MEIS2+ RPC2_d74", "ONECUT1+/MEIS2+ RPC3_d74", "VSX2++ RPC_d74",
               "ELAVL2+/4+ RPC1_d74", "ELAVL2+/4+ RPC2_d74",
               "ELAVL2+/4+ RPC3_d74", "EEF1A1+ RPC_d74")
  )

library(SeuratWrappers)
library(monocle3)

cds <- as.cell_data_set(RNA_d74_07_anno)

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