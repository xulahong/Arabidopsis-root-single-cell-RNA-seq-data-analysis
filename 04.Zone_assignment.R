library(Seurat)
PTI.combined <- readRDS("./PTI.combined.rds")
PTI_Bargmann2013 <- readRDS("PTI_zone_metadata_Bargmann2013.rds")

PTI.combined$zone_bargman <- PTI_Bargmann2013$PTI_class_data
PTI.combined$zone4 <- as.character(PTI.combined$zone_bargman)
table(PTI.combined$zone4)

PTI.combined$zone4[PTI.combined$celltype6 == 'Columella' | PTI.combined$celltype6 == 'LRC'] <- "Root cap"
PTI.combined$zone4 <- gsub("meristematic zone", "Meristem", PTI.combined$zone4)
PTI.combined$zone4 <- gsub("elongation zone", "Elongation", PTI.combined$zone4)
PTI.combined$zone4 <- gsub("maturation zone", "Maturation", PTI.combined$zone4)
table(PTI.combined$zone4)
PTI.combined$zone4 <- factor(PTI.combined$zone4, level=c("Root cap", "Meristem", "Elongation", "Maturation"))
DimPlot(PTI.combined, reduction = "umap", group.by = "zone4")
table(PTI.combined$zone4 )
colors = c("#bee6ae", "#c4c3be", "#c8cbde", "#575756")
DimPlot(PTI.combined, reduction = "umap", group.by = "zone4", cols = colors)

DimPlot(PTI.combined, reduction = "umap", group.by = "zone4", cols = colors, split.by = "orig.ident")

saveRDS(PTI.combined, file = "PTI.combined.rds")

