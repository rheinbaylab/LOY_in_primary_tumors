Yfeatures <- c('DDX3Y', 'EIF1AY', 'KDM5D', 'RPS4Y1', 'USP9Y', 'UTY', 'ZFY')
gene_lis <- read.csv('./genelist.housekeeping.txt')
chr_Y <- read.csv('./chrY.noPARgenes.txt', header = F)
stable_lis <- read.csv("./stable_expressed.csv", header = F)
Yfeatures <- c('DDX3Y', 'EIF1AY', 'KDM5D', 'RPS4Y1', 'USP9Y', 'UTY', 'ZFY')

# libraries
library(ggplot2)
library(Seurat)
library(dplyr)
library(stringr)
library(ggthemes)
library(matrixStats)
library(pheatmap)
library(Rtsne)


###############################################################
# 1. load data
# data are downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139829
folders <- list.files('./GSE139829/')

sceList <- lapply(folders,function(folder){
  CreateSeuratObject(counts = Read10X(paste('./GSE139829/', folder, sep = "")),
                     min.features = 100,
                     project = folder )
})


uvm <- merge(sceList[[1]],
             y = c(sceList[[2]],sceList[[3]],sceList[[4]],
                   sceList[[5]],sceList[[6]],
                   sceList[[7]],sceList[[8]],sceList[[9]],
                   sceList[[10]],sceList[[11]]),
             add.cell.ids = folders)


###############################################################
# 2. Quality control
# filter mt and cell
uvm[["percent.mt"]] <- PercentageFeatureSet(uvm, pattern = "^MT-")

uvm$log10GenesPerUMI <- log10(uvm$nFeature_RNA) / log10(uvm$nCount_RNA)

uvm_filter <- subset(uvm, subset = nFeature_RNA > 250 & nFeature_RNA < 8000 & nCount_RNA > 500 & percent.mt < 10 & log10GenesPerUMI > 0.8)

clinical <- read.csv("./clinical.csv")

uvm_filter@meta.data$cell <- row.names(uvm_filter@meta.data)
meta <- merge(uvm_filter@meta.data, clinical, by.x = "orig.ident", by.y = "Sample")
row.names(meta) <- meta$cell

uvm_filter@meta.data <- meta

uvm_filter <- NormalizeData(uvm_filter)
uvm_filter <- FindVariableFeatures(uvm_filter, selection.method = "vst", nfeatures = 2000)
uvm_filter <- ScaleData(uvm_filter)
uvm_filter <- RunPCA(uvm_filter)



###############################################################
# 3. annotate cell types
# Determine the K-nearest neighbor graph
uvm_filter <- FindNeighbors(object = uvm_filter, dims = 1:20)
uvm_filter <- FindClusters(uvm_filter, resolution = 1.5)
uvm_filter <- RunUMAP(uvm_filter, dims = 1:20)
uvm_filter <- RunTSNE(uvm_filter, dims = 1:20)

DimPlot(uvm_filter, reduction = "tsne", label = T) # 46 clusters

# B cells - CD19, CD79A, MS4A1
VlnPlot(object = uvm_filter,
        features = c("CD19", "CD79A", "MS4A1"))
FeaturePlot(uvm_filter,
            reduction = "umap",
            features = c("CD19", "CD79A", "MS4A1"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)


# T cells - CD3D, CD3E, CD8A
VlnPlot(object = uvm_filter,
        features = c("CD3D", "CD3E", "CD8A"))
FeaturePlot(uvm_filter,
            reduction = "umap",
            features = c("CD3D", "CD3E", "CD8A"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

# Plasma cells - IGHG1, MZB1, SDC1, CD79A
VlnPlot(object = uvm_filter,
        features = c("IGHG1", "MZB1", "SDC1", "CD79A"))
FeaturePlot(uvm_filter,
            reduction = "umap",
            features = c("IGHG1", "MZB1", "SDC1", "CD79A"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

# Monocytes and macrophages - CD68, CD163, CD14
VlnPlot(object = uvm_filter,
        features = c("CD68", "CD163", "CD14"))
FeaturePlot(uvm_filter,
            reduction = "umap",
            features = c("CD68", "CD163", "CD14"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

# NK Cells - FGFBP2, FCG3RA, CX3CR1
VlnPlot(object = uvm_filter,
        features = c("FGFBP2", "FCG3RA", "CX3CR1"))
FeaturePlot(uvm_filter,
            reduction = "umap",
            features = c("FGFBP2", "FCG3RA", "CX3CR1"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

# tumor cells - "MLANA", "MITF", "DCT"
VlnPlot(object = uvm_filter,
        features = c("MLANA", "MITF", "DCT"))
FeaturePlot(uvm_filter,
            reduction = "umap",
            features = c("MLANA", "MITF", "DCT"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

VlnPlot(object = uvm_filter,
        features = c("ROBO1", "EIF1AX", "SF3B1", "HTR2B", "ECM1", "CDH1"))

# Retinal pigment epithelium - RPE65
VlnPlot(object = uvm_filter,
        features = c("RPE65"))
FeaturePlot(uvm_filter,
            reduction = "umap",
            features = c("RPE65"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

# Photoreceptor cells - RCVRN
VlnPlot(object = uvm_filter,
        features = c("RCVRN"))
FeaturePlot(uvm_filter,
            reduction = "umap",
            features = c("RCVRN"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)


# Endothelial cells - PECAM1, VWF
VlnPlot(object = uvm_filter,
        features = c("PECAM1", "VWF"))
FeaturePlot(uvm_filter,
            reduction = "umap",
            features = c("PECAM1", "VWF"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

uvm_filter <- RenameIdents(object = uvm_filter,
                           "0" = "Tumor cells",
                           "1" = "Tumor cells",
                           "2" = "Tumor cells",
                           "3" = "Tumor cells",
                           "4" = "Tumor cells",
                           "5" = "Tumor cells",
                           "6" = "T cells",
                           "7" = "Tumor cells",
                           "8" = "Tumor cells",
                           "9" = "Monocytes & Macrophages",
                           "10" = "Monocytes & Macrophages",
                           "11" = "Tumor cells",
                           "12" = "T cells",
                           "13" = "T cells",
                           "14" = "Tumor cells",
                           "15" = "Tumor cells",
                           "16" = "Tumor cells",
                           "17" = "Tumor cells",
                           "18" = "Tumor cells",
                           "19" = "Tumor cells",
                           "20" = "Tumor cells",
                           "21" = "Tumor cells",
                           "22" = "T cells",
                           "23" = "T cells",
                           "24" = "Plasma cells",
                           "25" = "T cells",
                           "26" = "T cells",
                           "27" = "Tumor cells",
                           "28" = "Monocytes & Macrophages",
                           "29" = "T cells",
                           "30" = "NK cells",
                           "31" = "Tumor cells",
                           "32" = "T cells",
                           "33" = "T cells",
                           "34" = "B cells",
                           "35" = "Endothelial cells",
                           "36" = "Tumor cells",
                           "37" = "Photoreceptor cells",
                           "38" = "Tumor cells",
                           "39" = "Tumor cells",
                           "40" = "T cells",
                           "41" = "Monocytes & Macrophages",
                           "42" = "Tumor cells",
                           "43" = "Retinal Pigment Epithelium cells",
                           "44" = "Tumor cells",
                           "45" = "Tumor cells")

uvm_filter@meta.data$cell_type <- Idents(uvm_filter)

# Group cells: LOY vs WT
uvm_Y <- subset(uvm_filter, features = Yfeatures)
uvm_filter$Y_avg <- rowMeans(t(as.matrix(uvm_Y@assays$RNA@data)))
thre <- mean(uvm_filter@meta.data[uvm_filter@meta.data$Sex=="F",]$Y_avg) + 3*sd(uvm_filter@meta.data[uvm_filter@meta.data$Sex=="F",]$Y_avg)
uvm_filter$loy_pred <- ifelse(uvm_filter$Y_avg>thre, "Not_LOY", "LOY")

# subset data
uvm_Y <- subset(uvm_filter, features = Yfeatures)

uvm_Y$nCount_RNA_total <- uvm_filter@meta.data$nCount_RNA
uvm_Y$nFeature_RNA_total <- uvm_filter@meta.data$nFeature_RNA

# housekeeping genes
housekeeping <- subset(uvm_filter, features = gene_lis$housekeeping)
uvm_Y@meta.data$housekeeping_avg <- colMeans(as.matrix(housekeeping@assays$RNA@data))

# stable gene
stable <- subset(uvm_filter, features = stable_lis$V1)
uvm_Y@meta.data$stable_avg <- colMeans(as.matrix(stable@assays$RNA@data))

uvm_exp <- as.data.frame(t(as.matrix(uvm_Y@assays$RNA@data)))
uvm_meta <- uvm_Y@meta.data

# save the data
write.csv(uvm_exp, "uvm_exp.csv")
write.csv(uvm_meta, "uvm_meta.csv")
saveRDS(uvm_filter, "uvm_processed.rds")

