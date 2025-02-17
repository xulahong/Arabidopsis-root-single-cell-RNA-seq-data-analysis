library(Seurat)
library(SingleCellExperiment)
library(condiments)
library(slingshot)
# For data manipulation
library(dplyr)
library(tidyr)
# For visualization
library(ggplot2)
library(RColorBrewer)
library(viridis)
set.seed(2071)
library(tradeSeq)
library(RColorBrewer)
library(cowplot)
library(scales)
library(pheatmap)
############################################### 
PTI.combined <- readRDS("./PTI.combined.rds")
DimPlot(PTI.combined, reduction = "umap",label = T)
Idents(PTI.combined) <- "celltype6"
table(PTI.combined$celltype6)
############################################### flg22
DDD.data <- list()

for (i in c('Epidermis', 'Cortex', 'Endodermis', 'Stele')) {
  PTI <- subset(PTI.combined, idents = i)
  Idents(PTI) <- "orig.ident"
  table(PTI$orig.ident)
  
  PTI <- subset(PTI, idents = c("Ctrl1", "flg22"))
  
  PTI$conditions <- "Ctrl1"
  PTI$conditions[PTI$orig.ident == "flg22"] <- "flg22"
  table(PTI$conditions)
  
  table(PTI$zone4)
  
  PTI$spatial_id <- "Meristem"
  PTI$spatial_id[PTI$zone4 == "Elongation"] <- "Elongation"
  PTI$spatial_id[PTI$zone4 == "Maturation"] <- "Maturation"
  
  table(PTI$spatial_id)
  table(PTI$orig.ident)
  table(PTI$celltype6)
  table(PTI$conditions)
  
  
  #### convert back to singleCellExperiment ####
  tgfb <- as.SingleCellExperiment(PTI, assay = "RNA")
  
  
  #### run slingshot
  tgfb <- slingshot(tgfb, reducedDim = 'UMAP', clusterLabels = colData(tgfb)$spatial_id, 
                    start.clus = 'Meristem', approx_points = 150)
  
  
  # Differential progression test
  prog_res <- progressionTest(SlingshotDataSet(tgfb), conditions = tgfb$conditions)
  prog_res
  ks.test(a$Pseudotime[a$orig.ident == "Ctrl1"], a$Pseudotime[a$orig.ident == "flg22"])
  
  ks.test(slingPseudotime(tgfb)[colData(tgfb)$conditions == "Ctrl1",1],
          slingPseudotime(tgfb)[colData(tgfb)$conditions == "flg22",1])
  
  df <- bind_cols(
    as.data.frame(reducedDims(tgfb)$UMAP),
    as.data.frame(colData(tgfb)@listData[-44])
  ) %>%
    sample_frac(1)
  
  curve <- slingCurves(tgfb)[[1]]
  
  DDD.data[[i]] <- list(coldata=df, prog_res=prog_res, curve=curve)
}

save(DDD.data, file = "DDD.data.flg22.Rdata")
############################################### pep1

DDD.data <- list()

for (i in c('Epidermis', 'Cortex', 'Endodermis', 'Stele')) {
  PTI <- subset(PTI.combined, idents = i)
  Idents(PTI) <- "orig.ident"
  table(PTI$orig.ident)
  
  PTI <- subset(PTI, idents = c("Ctrl1", "Pep1"))
  
  PTI$conditions <- "Ctrl1"
  PTI$conditions[PTI$orig.ident == "Pep1"] <- "Pep1"
  table(PTI$conditions)
  
  table(PTI$zone4)
  
  PTI$spatial_id <- "Meristem"
  PTI$spatial_id[PTI$zone4 == "Elongation"] <- "Elongation"
  PTI$spatial_id[PTI$zone4 == "Maturation"] <- "Maturation"
  
  table(PTI$spatial_id)
  table(PTI$orig.ident)
  table(PTI$celltype6)
  table(PTI$conditions)
  
  
  #### convert back to singleCellExperiment ####
  tgfb <- as.SingleCellExperiment(PTI, assay = "RNA")
  
  
  #### run slingshot
  tgfb <- slingshot(tgfb, reducedDim = 'UMAP', clusterLabels = colData(tgfb)$spatial_id, 
                    start.clus = 'Meristem', approx_points = 150)
  
  
  # Differential progression test
  prog_res <- progressionTest(SlingshotDataSet(tgfb), conditions = tgfb$conditions)
  
  df <- bind_cols(
    as.data.frame(reducedDims(tgfb)$UMAP),
    as.data.frame(colData(tgfb)@listData[-44])
  ) %>%
    sample_frac(1)
  
  curve <- slingCurves(tgfb)[[1]]
  
  DDD.data[[i]] <- list(coldata=df, prog_res=prog_res, curve=curve)
}

save(DDD.data, file = "DDD.data.Pep1.Rdata")

############################################### 

############################################### Chitin
DDD.data <- list()

for (i in c('Epidermis', 'Cortex', 'Endodermis', 'Stele')) {
  PTI <- subset(PTI.combined, idents = i)
  Idents(PTI) <- "orig.ident"
  table(PTI$orig.ident)
  
  PTI <- subset(PTI, idents = c("Ctrl2", "Chitin"))
  
  PTI$conditions <- "Ctrl2"
  PTI$conditions[PTI$orig.ident == "Chitin"] <- "Chitin"
  table(PTI$conditions)
  
  table(PTI$zone4)
  
  PTI$spatial_id <- "Meristem"
  PTI$spatial_id[PTI$zone4 == "Elongation"] <- "Elongation"
  PTI$spatial_id[PTI$zone4 == "Maturation"] <- "Maturation"
  
  table(PTI$spatial_id)
  table(PTI$orig.ident)
  table(PTI$celltype6)
  table(PTI$conditions)
  
  
  #### convert back to singleCellExperiment ####
  tgfb <- as.SingleCellExperiment(PTI, assay = "RNA")
  
  
  #### run slingshot
  tgfb <- slingshot(tgfb, reducedDim = 'UMAP', clusterLabels = colData(tgfb)$spatial_id, 
                    start.clus = 'Meristem', approx_points = 150)
  
  
  # Differential progression test
  prog_res <- progressionTest(SlingshotDataSet(tgfb), conditions = tgfb$conditions)
  
  df <- bind_cols(
    as.data.frame(reducedDims(tgfb)$UMAP),
    as.data.frame(colData(tgfb)@listData[-44])
  ) %>%
    sample_frac(1)
  
  curve <- slingCurves(tgfb)[[1]]
  
  DDD.data[[i]] <- list(coldata=df, prog_res=prog_res, curve=curve)
}

save(DDD.data, file = "DDD.data.Chitin.Rdata")

############################################### SCOOP10
DDD.data <- list()

for (i in c('Epidermis', 'Cortex', 'Endodermis', 'Stele')) {
  PTI <- subset(PTI.combined, idents = i)
  Idents(PTI) <- "orig.ident"
  table(PTI$orig.ident)
  
  PTI <- subset(PTI, idents = c("Ctrl2", "SCOOP10"))
  
  PTI$conditions <- "Ctrl2"
  PTI$conditions[PTI$orig.ident == "SCOOP10"] <- "SCOOP10"
  table(PTI$conditions)
  
  PTI$spatial_id <- "Meristem"
  PTI$spatial_id[PTI$zone4 == "Elongation"] <- "Elongation"
  PTI$spatial_id[PTI$zone4 == "Maturation"] <- "Maturation"
  
  
  #### convert back to singleCellExperiment ####
  tgfb <- as.SingleCellExperiment(PTI, assay = "RNA")
  
  
  #### run slingshot
  tgfb <- slingshot(tgfb, reducedDim = 'UMAP', clusterLabels = colData(tgfb)$spatial_id, 
                    start.clus = 'Meristem', approx_points = 150)
  
  
  # Differential progression test
  prog_res <- progressionTest(SlingshotDataSet(tgfb), conditions = tgfb$conditions)
  
  df <- bind_cols(
    as.data.frame(reducedDims(tgfb)$UMAP),
    as.data.frame(colData(tgfb)@listData[-44])
  ) %>%
    sample_frac(1)
  
  curve <- slingCurves(tgfb)[[1]]
  
  DDD.data[[i]] <- list(coldata=df, prog_res=prog_res, curve=curve)
}

save(DDD.data, file = "DDD.data.SCOOP10.Rdata")

############################################### SCREW1
DDD.data <- list()

for (i in c('Epidermis', 'Cortex', 'Endodermis', 'Stele')) {
  PTI <- subset(PTI.combined, idents = i)
  Idents(PTI) <- "orig.ident"
  table(PTI$orig.ident)
  
  PTI <- subset(PTI, idents = c("Ctrl2", "SCREW1"))
  
  PTI$conditions <- "Ctrl2"
  PTI$conditions[PTI$orig.ident == "SCREW1"] <- "SCREW1"
  table(PTI$conditions)
  
  PTI$spatial_id <- "Meristem"
  PTI$spatial_id[PTI$zone4 == "Elongation"] <- "Elongation"
  PTI$spatial_id[PTI$zone4 == "Maturation"] <- "Maturation"
  
  
  #### convert back to singleCellExperiment ####
  tgfb <- as.SingleCellExperiment(PTI, assay = "RNA")
  
  
  #### run slingshot
  tgfb <- slingshot(tgfb, reducedDim = 'UMAP', clusterLabels = colData(tgfb)$spatial_id, 
                    start.clus = 'Meristem', approx_points = 150)
  
  
  # Differential progression test
  prog_res <- progressionTest(SlingshotDataSet(tgfb), conditions = tgfb$conditions)
  
  df <- bind_cols(
    as.data.frame(reducedDims(tgfb)$UMAP),
    as.data.frame(colData(tgfb)@listData[-44])
  ) %>%
    sample_frac(1)
  
  curve <- slingCurves(tgfb)[[1]]
  
  DDD.data[[i]] <- list(coldata=df, prog_res=prog_res, curve=curve)
}

save(DDD.data, file = "DDD.data.SCREW1.Rdata")
