library(Seurat)
ctrl1.data <- Read10X("/scratch/user/xulahong/scRNA_seq_flg22_pep1/Root_CTL/Ctrl_count_04132021/outs/filtered_feature_bc_matrix")
flg22.data <- Read10X("/scratch/user/xulahong/scRNA_seq_flg22_pep1/Root_Flgn1/Flg22_count_04132021/outs/filtered_feature_bc_matrix")
pep1.data <- Read10X("/scratch/user/xulahong/scRNA_seq_flg22_pep1/Root_Pep1/Pep1_count_04132021/outs/filtered_feature_bc_matrix")

ctrl2.data <- Read10X("/scratch/user/xulahong/scRNA_seq_chitin_scoop_screw1/Arabadopsis_sc_02_2021/Ctrl_count_04122021/outs/filtered_feature_bc_matrix")
chitin.data <- Read10X("/scratch/user/xulahong/scRNA_seq_chitin_scoop_screw1/Arabadopsis_sc_02_2021/Chitin_count_04122021/outs/filtered_feature_bc_matrix")
scoop10.data <- Read10X("/scratch/user/xulahong/scRNA_seq_chitin_scoop_screw1/Arabadopsis_sc_02_2021/SCOOP10_count_04122021/outs/filtered_feature_bc_matrix")
screw1.data <- Read10X("/scratch/user/xulahong/scRNA_seq_chitin_scoop_screw1/Arabadopsis_sc_02_2021/SCREW1_count_04122021/outs/filtered_feature_bc_matrix")

ctrl1 <- CreateSeuratObject(counts = ctrl1.data, project = "ctrl1")
flg22 <- CreateSeuratObject(counts = flg22.data, project = "flg22")
pep1 <- CreateSeuratObject(counts = pep1.data, project = "pep1")

ctrl2 <- CreateSeuratObject(counts = ctrl2.data, project = "ctrl2")
chitin <- CreateSeuratObject(counts = chitin.data, project = "chitin")
scoop10 <- CreateSeuratObject(counts = scoop10.data, project = "scoop10")
screw1 <- CreateSeuratObject(counts = screw1.data, project = "screw1")

# Merging three seurat object
PTI.combined <- merge(ctrl1, y = c(flg22, pep1, ctrl2, chitin, screw1, scoop10), add.cell.ids = c("ctrl1", "flg22", "pep1", "ctrl2", "chitin","screw1", "scoop10"), project = "PTI")
table(PTI.combined$orig.ident)

#save file
PTI.combined.raw <- PTI.combined
saveRDS(PTI.combined.raw, file = "./PTI.combined.raw.rds")

# Ploting
Idents(PTI) <- "orig.ident"
DefaultAssay(PTI) <- "RNA"
PTI[["percent.mt"]] <- PercentageFeatureSet(PTI, pattern = "^ATMG")
VlnPlot(PTI, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(PTI, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(PTI, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# Filtering
PTI <- subset(PTI, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & nCount_RNA < 150000 & percent.mt < 15)
# PTI.combined <- subset(PTI.combined, subset = nFeature_RNA > 200 & nFeature_RNA < 115000 & nCount_RNA < 200000 & percent.mt < 40)
VlnPlot(PTI, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
table(PTI$orig.ident)


# Ploting
Idents(PTI2) <- "orig.ident"
DefaultAssay(PTI2) <- "RNA"
PTI2[["percent.mt"]] <- PercentageFeatureSet(PTI2, pattern = "^ATMG")
VlnPlot(PTI2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(PTI2, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(PTI2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# Filtering
PTI2 <- subset(PTI2, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & nCount_RNA < 75000 & percent.mt < 30)
VlnPlot(PTI2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
table(PTI2$orig.ident)

table(PTI$orig.ident)
table(PTI2$orig.ident)
head(PTI@meta.data)
head(PTI2@meta.data)

PTI.combined <- merge(PTI, PTI2)
head(PTI.combined@meta.data)
table(PTI.combined$orig.ident)

## remove some genes
counts <- GetAssayData(PTI.combined, assay = "RNA")
mito.genes <- grep(pattern = "^ATMG", x = rownames(x = counts), value = TRUE)
chlo.genes <- grep(pattern = "^ATCG", x = rownames(x = counts), value = TRUE)
proto.genes <- read.csv("~/lahong/scRNA-seq/scRNA_seq_flg22_pep1/cell_type/cell type/356genes_protoplast.csv")
proto.genes <- proto.genes$GENEID
counts <- counts[-(which(rownames(counts) %in% c(mito.genes, chlo.genes, proto.genes))),]
PTI.combined <- subset(PTI.combined, features = rownames(counts))
dim(PTI.combined)

saveRDS(PTI.combined, file = "./PTI.combined.rds")
