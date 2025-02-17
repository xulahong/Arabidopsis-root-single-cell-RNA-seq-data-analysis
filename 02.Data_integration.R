library(Seurat)
PTI.combined <-  readRDS( "./PTI.combined.rds")
DefaultAssay(PTI.combined) <- "RNA"

list <- SplitObject(PTI.combined, split.by = "orig.ident")
list <- lapply(X = list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = list)
list <- lapply(X = list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
rm(PTI.combined)

# anchors <- FindIntegrationAnchors(object.list = list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = list, reference = c(1, 2), reduction = "rpca", dims = 1:50)

# this command creates an 'integrated' data assay
PTI.combined <- IntegrateData(anchorset = anchors)
# PTI.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)

# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(PTI.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
PTI.combined <- ScaleData(PTI.combined, verbose = FALSE)
PTI.combined <- RunPCA(PTI.combined, npcs = 50, verbose = FALSE)
PTI.combined <- RunUMAP(PTI.combined, reduction = "pca", dims = 1:50)
DimPlot(PTI.combined)
head(PTI.combined@meta.data)
PTI.combined <- FindNeighbors(PTI.combined, reduction = "pca", dims = 1:50)
head(PTI.combined@meta.data)

PTI.combined <- FindClusters(PTI.combined, resolution = 0.1)

saveRDS(PTI.combined, file = "./PTI.combined.rds")
