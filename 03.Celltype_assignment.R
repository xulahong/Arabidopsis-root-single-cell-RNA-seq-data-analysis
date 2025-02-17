library(Seurat)
PTI.combined <-  readRDS( "./PTI.combined.rds")
DefaultAssay(PTI.combined) <- "integrated"
PTI.combined <- FindClusters(PTI.combined, resolution = 0.1)
Idents(PTI.combined) <- PTI.combined$integrated_snn_res.0.1
DimPlot(PTI.combined, reduction = "umap",label = T)
PTI.combined <- RenameIdents(object = PTI.combined, 
                             `0` = "Epidermis", `1` = "LRC", `2` = "Epidermis", `3`='Stele', `4`='Cortex',`5`='Columella',
                             `6`='Endodermis', `7`='Columella', `8`='Stele')

DimPlot(PTI.combined, reduction = "umap",label = T)

###### ###### ###### ------------------- order the cell type ------------------ ###### ###### ###### 
PTI.combined$celltype6 <- Idents(PTI.combined)
PTI.combined$celltype6 <- factor(PTI.combined$celltype6, level=c("Columella", "LRC", "Epidermis", "Cortex", "Endodermis", "Stele"))
table(PTI.combined$celltype6, PTI.combined$orig.ident)
colnames(PTI.combined@meta.data)

DimPlot(PTI.combined, reduction = "umap",label = F, split.by = "orig.ident", group.by = "celltype6")
saveRDS(PTI.combined, file = "PTI.combined.rds")
