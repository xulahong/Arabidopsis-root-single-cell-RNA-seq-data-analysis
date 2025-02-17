library(Seurat)
PTI.combined <- readRDS("./PTI.combined.rds")
DimPlot(PTI.combined, reduction = "umap", group.by = "zone4")
DimPlot(PTI.combined, reduction = "umap", group.by = "celltype6")
PTI.combined$zone.celltype.stim <- paste(PTI.combined$zone4, PTI.combined$celltype6, PTI.combined$orig.ident, sep = "_")
Idents(PTI.combined) <- "zone.celltype.stim"
table(PTI.combined@active.ident)

####-------------------------------------- flg22 VS control-----------------------
Rootcap_LRC_flg22 <- FindMarkers(PTI.combined, ident.1 = "Root cap_LRC_flg22", ident.2 = "Root cap_LRC_Ctrl1", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Rootcap_Columella_flg22 <- FindMarkers(PTI.combined, ident.1 = "Root cap_Columella_flg22", ident.2 = "Root cap_Columella_Ctrl1", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)


Meristem_Epidermis_flg22 <- FindMarkers(PTI.combined, ident.1 = "Meristem_Epidermis_flg22", ident.2 = "Meristem_Epidermis_Ctrl1", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Meristem_Cortex_flg22 <- FindMarkers(PTI.combined, ident.1 = "Meristem_Cortex_flg22", ident.2 = "Meristem_Cortex_Ctrl1", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Meristem_Endodermis_flg22 <- FindMarkers(PTI.combined, ident.1 = "Meristem_Endodermis_flg22", ident.2 = "Meristem_Endodermis_Ctrl1", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Meristem_Stele_flg22 <- FindMarkers(PTI.combined, ident.1 = "Meristem_Stele_flg22", ident.2 = "Meristem_Stele_Ctrl1", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)

Elongation_Epidermis_flg22 <- FindMarkers(PTI.combined, ident.1 = "Elongation_Epidermis_flg22", ident.2 = "Elongation_Epidermis_Ctrl1", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Elongation_Cortex_flg22 <- FindMarkers(PTI.combined, ident.1 = "Elongation_Cortex_flg22", ident.2 = "Elongation_Cortex_Ctrl1", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Elongation_Endodermis_flg22 <- FindMarkers(PTI.combined, ident.1 = "Elongation_Endodermis_flg22", ident.2 = "Elongation_Endodermis_Ctrl1", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Elongation_Stele_flg22 <- FindMarkers(PTI.combined, ident.1 = "Elongation_Stele_flg22", ident.2 = "Elongation_Stele_Ctrl1", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)


Maturation_Epidermis_flg22 <- FindMarkers(PTI.combined, ident.1 = "Maturation_Epidermis_flg22", ident.2 = "Maturation_Epidermis_Ctrl1", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Maturation_Cortex_flg22 <- FindMarkers(PTI.combined, ident.1 = "Maturation_Cortex_flg22", ident.2 = "Maturation_Cortex_Ctrl1", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Maturation_Endodermis_flg22 <- FindMarkers(PTI.combined, ident.1 = "Maturation_Endodermis_flg22", ident.2 = "Maturation_Endodermis_Ctrl1", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Maturation_Stele_flg22 <- FindMarkers(PTI.combined, ident.1 = "Maturation_Stele_flg22", ident.2 = "Maturation_Stele_Ctrl1", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)


DE.flg22 <- list(Rootcap_LRC_flg22=Rootcap_LRC_flg22, Rootcap_Columella_flg22=Rootcap_Columella_flg22,
                 Meristem_Epidermis_flg22=Meristem_Epidermis_flg22, Meristem_Cortex_flg22=Meristem_Cortex_flg22, Meristem_Endodermis_flg22=Meristem_Endodermis_flg22, Meristem_Stele_flg22=Meristem_Stele_flg22,
                 Elongation_Epidermis_flg22=Elongation_Epidermis_flg22, Elongation_Cortex_flg22=Elongation_Cortex_flg22, Elongation_Endodermis_flg22=Elongation_Endodermis_flg22, Elongation_Stele_flg22=Elongation_Stele_flg22,
                 Maturation_Epidermis_flg22=Maturation_Epidermis_flg22, Maturation_Cortex_flg22=Maturation_Cortex_flg22, Maturation_Endodermis_flg22=Maturation_Endodermis_flg22, Maturation_Stele_flg22=Maturation_Stele_flg22)


####-------------------------------------- Pep1 VS control-----------------------
Rootcap_LRC_Pep1 <- FindMarkers(PTI.combined, ident.1 = "Root cap_LRC_Pep1", ident.2 = "Root cap_LRC_Ctrl1", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Rootcap_Columella_Pep1 <- FindMarkers(PTI.combined, ident.1 = "Root cap_Columella_Pep1", ident.2 = "Root cap_Columella_Ctrl1", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)


Meristem_Epidermis_Pep1 <- FindMarkers(PTI.combined, ident.1 = "Meristem_Epidermis_Pep1", ident.2 = "Meristem_Epidermis_Ctrl1", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Meristem_Cortex_Pep1 <- FindMarkers(PTI.combined, ident.1 = "Meristem_Cortex_Pep1", ident.2 = "Meristem_Cortex_Ctrl1", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Meristem_Endodermis_Pep1 <- FindMarkers(PTI.combined, ident.1 = "Meristem_Endodermis_Pep1", ident.2 = "Meristem_Endodermis_Ctrl1", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Meristem_Stele_Pep1 <- FindMarkers(PTI.combined, ident.1 = "Meristem_Stele_Pep1", ident.2 = "Meristem_Stele_Ctrl1", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)

Elongation_Epidermis_Pep1 <- FindMarkers(PTI.combined, ident.1 = "Elongation_Epidermis_Pep1", ident.2 = "Elongation_Epidermis_Ctrl1", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Elongation_Cortex_Pep1 <- FindMarkers(PTI.combined, ident.1 = "Elongation_Cortex_Pep1", ident.2 = "Elongation_Cortex_Ctrl1", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Elongation_Endodermis_Pep1 <- FindMarkers(PTI.combined, ident.1 = "Elongation_Endodermis_Pep1", ident.2 = "Elongation_Endodermis_Ctrl1", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Elongation_Stele_Pep1 <- FindMarkers(PTI.combined, ident.1 = "Elongation_Stele_Pep1", ident.2 = "Elongation_Stele_Ctrl1", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)


Maturation_Epidermis_Pep1 <- FindMarkers(PTI.combined, ident.1 = "Maturation_Epidermis_Pep1", ident.2 = "Maturation_Epidermis_Ctrl1", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Maturation_Cortex_Pep1 <- FindMarkers(PTI.combined, ident.1 = "Maturation_Cortex_Pep1", ident.2 = "Maturation_Cortex_Ctrl1", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Maturation_Endodermis_Pep1 <- FindMarkers(PTI.combined, ident.1 = "Maturation_Endodermis_Pep1", ident.2 = "Maturation_Endodermis_Ctrl1", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Maturation_Stele_Pep1 <- FindMarkers(PTI.combined, ident.1 = "Maturation_Stele_Pep1", ident.2 = "Maturation_Stele_Ctrl1", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)

DE.Pep1 <- list(Rootcap_LRC_Pep1=Rootcap_LRC_Pep1, Rootcap_Columella_Pep1=Rootcap_Columella_Pep1, 
                Meristem_Epidermis_Pep1=Meristem_Epidermis_Pep1, Meristem_Cortex_Pep1=Meristem_Cortex_Pep1, Meristem_Endodermis_Pep1=Meristem_Endodermis_Pep1, Meristem_Stele_Pep1=Meristem_Stele_Pep1,
                Elongation_Epidermis_Pep1=Elongation_Epidermis_Pep1, Elongation_Cortex_Pep1=Elongation_Cortex_Pep1, Elongation_Endodermis_Pep1=Elongation_Endodermis_Pep1, Elongation_Stele_Pep1=Elongation_Stele_Pep1,
                Maturation_Epidermis_Pep1=Maturation_Epidermis_Pep1, Maturation_Cortex_Pep1=Maturation_Cortex_Pep1, Maturation_Endodermis_Pep1=Maturation_Endodermis_Pep1, Maturation_Stele_Pep1=Maturation_Stele_Pep1)

####-------------------------------------- Chitin VS control-----------------------
Rootcap_LRC_Chitin <- FindMarkers(PTI.combined, ident.1 = "Root cap_LRC_Chitin", ident.2 = "Root cap_LRC_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Rootcap_Columella_Chitin <- FindMarkers(PTI.combined, ident.1 = "Root cap_Columella_Chitin", ident.2 = "Root cap_Columella_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)


Meristem_Epidermis_Chitin <- FindMarkers(PTI.combined, ident.1 = "Meristem_Epidermis_Chitin", ident.2 = "Meristem_Epidermis_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Meristem_Cortex_Chitin <- FindMarkers(PTI.combined, ident.1 = "Meristem_Cortex_Chitin", ident.2 = "Meristem_Cortex_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Meristem_Endodermis_Chitin <- FindMarkers(PTI.combined, ident.1 = "Meristem_Endodermis_Chitin", ident.2 = "Meristem_Endodermis_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Meristem_Stele_Chitin <- FindMarkers(PTI.combined, ident.1 = "Meristem_Stele_Chitin", ident.2 = "Meristem_Stele_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)

Elongation_Epidermis_Chitin <- FindMarkers(PTI.combined, ident.1 = "Elongation_Epidermis_Chitin", ident.2 = "Elongation_Epidermis_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Elongation_Cortex_Chitin <- FindMarkers(PTI.combined, ident.1 = "Elongation_Cortex_Chitin", ident.2 = "Elongation_Cortex_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Elongation_Endodermis_Chitin <- FindMarkers(PTI.combined, ident.1 = "Elongation_Endodermis_Chitin", ident.2 = "Elongation_Endodermis_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Elongation_Stele_Chitin <- FindMarkers(PTI.combined, ident.1 = "Elongation_Stele_Chitin", ident.2 = "Elongation_Stele_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)


Maturation_Epidermis_Chitin <- FindMarkers(PTI.combined, ident.1 = "Maturation_Epidermis_Chitin", ident.2 = "Maturation_Epidermis_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Maturation_Cortex_Chitin <- FindMarkers(PTI.combined, ident.1 = "Maturation_Cortex_Chitin", ident.2 = "Maturation_Cortex_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Maturation_Endodermis_Chitin <- FindMarkers(PTI.combined, ident.1 = "Maturation_Endodermis_Chitin", ident.2 = "Maturation_Endodermis_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Maturation_Stele_Chitin <- FindMarkers(PTI.combined, ident.1 = "Maturation_Stele_Chitin", ident.2 = "Maturation_Stele_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)

DE.Chitin <- list(Rootcap_LRC_Chitin=Rootcap_LRC_Chitin, Rootcap_Columella_Chitin=Rootcap_Columella_Chitin,
                  Meristem_Epidermis_Chitin=Meristem_Epidermis_Chitin, Meristem_Cortex_Chitin=Meristem_Cortex_Chitin, Meristem_Endodermis_Chitin=Meristem_Endodermis_Chitin, Meristem_Stele_Chitin=Meristem_Stele_Chitin,
                  Elongation_Epidermis_Chitin=Elongation_Epidermis_Chitin, Elongation_Cortex_Chitin=Elongation_Cortex_Chitin, Elongation_Endodermis_Chitin=Elongation_Endodermis_Chitin, Elongation_Stele_Chitin=Elongation_Stele_Chitin,
                  Maturation_Epidermis_Chitin=Maturation_Epidermis_Chitin, Maturation_Cortex_Chitin=Maturation_Cortex_Chitin, Maturation_Endodermis_Chitin=Maturation_Endodermis_Chitin, Maturation_Stele_Chitin=Maturation_Stele_Chitin)

####-------------------------------------- SCOOP10 VS control-----------------------
Rootcap_LRC_SCOOP10 <- FindMarkers(PTI.combined, ident.1 = "Root cap_LRC_SCOOP10", ident.2 = "Root cap_LRC_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Rootcap_Columella_SCOOP10 <- FindMarkers(PTI.combined, ident.1 = "Root cap_Columella_SCOOP10", ident.2 = "Root cap_Columella_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)


Meristem_Epidermis_SCOOP10 <- FindMarkers(PTI.combined, ident.1 = "Meristem_Epidermis_SCOOP10", ident.2 = "Meristem_Epidermis_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Meristem_Cortex_SCOOP10 <- FindMarkers(PTI.combined, ident.1 = "Meristem_Cortex_SCOOP10", ident.2 = "Meristem_Cortex_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Meristem_Endodermis_SCOOP10 <- FindMarkers(PTI.combined, ident.1 = "Meristem_Endodermis_SCOOP10", ident.2 = "Meristem_Endodermis_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Meristem_Stele_SCOOP10 <- FindMarkers(PTI.combined, ident.1 = "Meristem_Stele_SCOOP10", ident.2 = "Meristem_Stele_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)


Elongation_Epidermis_SCOOP10 <- FindMarkers(PTI.combined, ident.1 = "Elongation_Epidermis_SCOOP10", ident.2 = "Elongation_Epidermis_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Elongation_Cortex_SCOOP10 <- FindMarkers(PTI.combined, ident.1 = "Elongation_Cortex_SCOOP10", ident.2 = "Elongation_Cortex_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Elongation_Endodermis_SCOOP10 <- FindMarkers(PTI.combined, ident.1 = "Elongation_Endodermis_SCOOP10", ident.2 = "Elongation_Endodermis_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Elongation_Stele_SCOOP10 <- FindMarkers(PTI.combined, ident.1 = "Elongation_Stele_SCOOP10", ident.2 = "Elongation_Stele_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)


Maturation_Epidermis_SCOOP10 <- FindMarkers(PTI.combined, ident.1 = "Maturation_Epidermis_SCOOP10", ident.2 = "Maturation_Epidermis_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Maturation_Cortex_SCOOP10 <- FindMarkers(PTI.combined, ident.1 = "Maturation_Cortex_SCOOP10", ident.2 = "Maturation_Cortex_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Maturation_Endodermis_SCOOP10 <- FindMarkers(PTI.combined, ident.1 = "Maturation_Endodermis_SCOOP10", ident.2 = "Maturation_Endodermis_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Maturation_Stele_SCOOP10 <- FindMarkers(PTI.combined, ident.1 = "Maturation_Stele_SCOOP10", ident.2 = "Maturation_Stele_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)


DE.SCOOP10 <- list(Rootcap_LRC_SCOOP10=Rootcap_LRC_SCOOP10, Rootcap_Columella_SCOOP10=Rootcap_Columella_SCOOP10,
                   Meristem_Epidermis_SCOOP10=Meristem_Epidermis_SCOOP10, Meristem_Cortex_SCOOP10=Meristem_Cortex_SCOOP10, Meristem_Endodermis_SCOOP10=Meristem_Endodermis_SCOOP10, Meristem_Stele_SCOOP10=Meristem_Stele_SCOOP10,
                   Elongation_Epidermis_SCOOP10=Elongation_Epidermis_SCOOP10, Elongation_Cortex_SCOOP10=Elongation_Cortex_SCOOP10, Elongation_Endodermis_SCOOP10=Elongation_Endodermis_SCOOP10, Elongation_Stele_SCOOP10=Elongation_Stele_SCOOP10,
                   Maturation_Epidermis_SCOOP10=Maturation_Epidermis_SCOOP10, Maturation_Cortex_SCOOP10=Maturation_Cortex_SCOOP10, Maturation_Endodermis_SCOOP10=Maturation_Endodermis_SCOOP10, Maturation_Stele_SCOOP10=Maturation_Stele_SCOOP10)


####-------------------------------------- SCREW1 VS control-----------------------
Rootcap_LRC_SCREW1 <- FindMarkers(PTI.combined, ident.1 = "Root cap_LRC_SCREW1", ident.2 = "Root cap_LRC_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Rootcap_Columella_SCREW1 <- FindMarkers(PTI.combined, ident.1 = "Root cap_Columella_SCREW1", ident.2 = "Root cap_Columella_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)

Meristem_Epidermis_SCREW1 <- FindMarkers(PTI.combined, ident.1 = "Meristem_Epidermis_SCREW1", ident.2 = "Meristem_Epidermis_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Meristem_Cortex_SCREW1 <- FindMarkers(PTI.combined, ident.1 = "Meristem_Cortex_SCREW1", ident.2 = "Meristem_Cortex_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Meristem_Endodermis_SCREW1 <- FindMarkers(PTI.combined, ident.1 = "Meristem_Endodermis_SCREW1", ident.2 = "Meristem_Endodermis_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Meristem_Stele_SCREW1 <- FindMarkers(PTI.combined, ident.1 = "Meristem_Stele_SCREW1", ident.2 = "Meristem_Stele_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)


Elongation_Epidermis_SCREW1 <- FindMarkers(PTI.combined, ident.1 = "Elongation_Epidermis_SCREW1", ident.2 = "Elongation_Epidermis_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Elongation_Cortex_SCREW1 <- FindMarkers(PTI.combined, ident.1 = "Elongation_Cortex_SCREW1", ident.2 = "Elongation_Cortex_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Elongation_Endodermis_SCREW1 <- FindMarkers(PTI.combined, ident.1 = "Elongation_Endodermis_SCREW1", ident.2 = "Elongation_Endodermis_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Elongation_Stele_SCREW1 <- FindMarkers(PTI.combined, ident.1 = "Elongation_Stele_SCREW1", ident.2 = "Elongation_Stele_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)


Maturation_Epidermis_SCREW1 <- FindMarkers(PTI.combined, ident.1 = "Maturation_Epidermis_SCREW1", ident.2 = "Maturation_Epidermis_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Maturation_Cortex_SCREW1 <- FindMarkers(PTI.combined, ident.1 = "Maturation_Cortex_SCREW1", ident.2 = "Maturation_Cortex_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Maturation_Endodermis_SCREW1 <- FindMarkers(PTI.combined, ident.1 = "Maturation_Endodermis_SCREW1", ident.2 = "Maturation_Endodermis_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)
Maturation_Stele_SCREW1 <- FindMarkers(PTI.combined, ident.1 = "Maturation_Stele_SCREW1", ident.2 = "Maturation_Stele_Ctrl2", verbose = FALSE, test.use = "MAST", logfc.threshold = 0)


DE.SCREW1 <- list(Rootcap_LRC_SCREW1=Rootcap_LRC_SCREW1, Rootcap_Columella_SCREW1=Rootcap_Columella_SCREW1,
                  Meristem_Epidermis_SCREW1=Meristem_Epidermis_SCREW1, Meristem_Cortex_SCREW1=Meristem_Cortex_SCREW1, Meristem_Endodermis_SCREW1=Meristem_Endodermis_SCREW1, Meristem_Stele_SCREW1=Meristem_Stele_SCREW1,
                  Elongation_Epidermis_SCREW1=Elongation_Epidermis_SCREW1, Elongation_Cortex_SCREW1=Elongation_Cortex_SCREW1, Elongation_Endodermis_SCREW1=Elongation_Endodermis_SCREW1, Elongation_Stele_SCREW1=Elongation_Stele_SCREW1,
                  Maturation_Epidermis_SCREW1=Maturation_Epidermis_SCREW1, Maturation_Cortex_SCREW1=Maturation_Cortex_SCREW1, Maturation_Endodermis_SCREW1=Maturation_Endodermis_SCREW1, Maturation_Stele_SCREW1=Maturation_Stele_SCREW1)

####------------------------------------------------------------------
DE_results <- list()
DE_results$DE.flg22 <- DE.flg22
DE_results$DE.Chitin <- DE.Chitin
DE_results$DE.Pep1 <- DE.Pep1
DE_results$DE.SCOOP10 <- DE.SCOOP10
DE_results$DE.SCREW1 <- DE.SCREW1

saveRDS(DE_results, file = "DE.rootzone4.celltype6.results.rds")


################################################################----------------------------------------------------------------------
deg <- DE.flg22

sigdeg <- list()
for (i in names(deg)) {
  a <- deg[[i]]$p_val_adj < 0.05 & abs(deg[[i]]$avg_logFC) >= 0.25
  b <- deg[[i]][a,]
  sigdeg[[i]] <- row.names(b)
}

sigdeg.flg22 <- sigdeg

# saveRDS(sigdeg.flg22, file = "sigDE.id.celltype7.flg22.rds")


################################################################----------------------------------------------------------------------
deg <- DE.Pep1

sigdeg <- list()
for (i in names(deg)) {
  a <- deg[[i]]$p_val_adj < 0.05 & abs(deg[[i]]$avg_logFC) >= 0.25
  b <- deg[[i]][a,]
  sigdeg[[i]] <- row.names(b)
}

sigdeg.pep1 <- sigdeg

# saveRDS(sigdeg.pep1, file = "sigDE.id.celltype7.pep1.rds")

################################################################----------------------------------------------------------------------
deg <- DE.Chitin

sigdeg <- list()
for (i in names(deg)) {
  a <- deg[[i]]$p_val_adj < 0.05 & abs(deg[[i]]$avg_logFC) >= 0.25
  b <- deg[[i]][a,]
  sigdeg[[i]] <- row.names(b)
}

sigdeg.chitin <- sigdeg

# saveRDS(sigdeg.chitin, file = "sigDE.id.celltype7.chitin.rds")

################################################################----------------------------------------------------------------------
deg <- DE.SCREW1

sigdeg <- list()
for (i in names(deg)) {
  a <- deg[[i]]$p_val_adj < 0.05 & abs(deg[[i]]$avg_logFC) >= 0.25
  b <- deg[[i]][a,]
  sigdeg[[i]] <- row.names(b)
}

sigdeg.edp1 <- sigdeg

# saveRDS(sigdeg.edp1, file = "sigDE.id.celltype7.edp1.rds")

################################################################----------------------------------------------------------------------
deg <- DE.SCOOP10

sigdeg <- list()
for (i in names(deg)) {
  a <- deg[[i]]$p_val_adj < 0.05 & abs(deg[[i]]$avg_logFC) >= 0.25
  b <- deg[[i]][a,]
  sigdeg[[i]] <- row.names(b)
}

sigdeg.scoop10 <- sigdeg

# saveRDS(sigdeg.scoop10, file = "sigDE.id.celltype7.scoop10.rds")


DE_results$sig.DE.flg22 <- sigdeg.flg22
DE_results$sig.DE.Chitin <- sigdeg.chitin
DE_results$sig.DE.Pep1 <- sigdeg.pep1
DE_results$sig.DE.SCOOP10 <- sigdeg.scoop10
DE_results$sig.DE.SCREW1 <- sigdeg.edp1

saveRDS(DE_results, file = "DE.rootzone4.celltype6.results.rds")

#########################################################################
############################################################ up and down
#########################################################################

setwd("~/lahong/scRNA-seq/scRNA_seq_all_3/08.root_zone4_celltype6_DE_data/")
####-------------------- flg22 sigDE up and down gene list -------------------------------
# DE.flg22 <- readRDS("DE.celltype7.flg22.rds")

degUp  <- list()
degDown  <- list()
for (i in names(DE.flg22)){
  up <- DE.flg22[[i]][DE.flg22[[i]]$p_val_adj < 0.05 & DE.flg22[[i]]$avg_logFC >= 0.25,]
  down <- DE.flg22[[i]][DE.flg22[[i]]$p_val_adj < 0.05 & DE.flg22[[i]]$avg_logFC <= (-0.25),]
  degUp[[i]] <- row.names(up)
  degDown[[i]] <- row.names(down)
}
degUp
degDown

flg22.up.gene <- degUp
flg22.down.gene <- degDown

# saveRDS(flg22.up.gene, file = "sigDE.id.celltype7.flg22.up.rds")
# saveRDS(flg22.down.gene, file = "sigDE.id.celltype7.flg22.down.rds")

####------------------- pep1 sigDE up and down gene list -------------------------------
# DE.Pep1 <- readRDS("DE.celltype7.pep1.rds")

degUp  <- list()
degDown  <- list()
for (i in names(DE.Pep1)){
  up <- DE.Pep1[[i]][DE.Pep1[[i]]$p_val_adj < 0.05 & DE.Pep1[[i]]$avg_logFC >= 0.25,]
  down <- DE.Pep1[[i]][DE.Pep1[[i]]$p_val_adj < 0.05 & DE.Pep1[[i]]$avg_logFC <= (-0.25),]
  degUp[[i]] <- row.names(up)
  degDown[[i]] <- row.names(down)
}
degUp
degDown

pep1.up.gene <- degUp
pep1.down.gene <- degDown

# saveRDS(pep1.up.gene, file = "sigDE.id.celltype7.pep1.up.rds")
# saveRDS(pep1.down.gene, file = "sigDE.id.celltype7.pep1.down.rds")

####-------------------- chitin sigDE up and down gene list -------------------------------
# DE.Chitin <- readRDS("DE.celltype7.chitin.rds")

degUp  <- list()
degDown  <- list()
for (i in names(DE.Chitin)){
  up <- DE.Chitin[[i]][DE.Chitin[[i]]$p_val_adj < 0.05 & DE.Chitin[[i]]$avg_logFC >= 0.25,]
  down <- DE.Chitin[[i]][DE.Chitin[[i]]$p_val_adj < 0.05 & DE.Chitin[[i]]$avg_logFC <= (-0.25),]
  degUp[[i]] <- row.names(up)
  degDown[[i]] <- row.names(down)
}
degUp
degDown

chitin.up.gene <- degUp
chitin.down.gene <- degDown

# saveRDS(chitin.up.gene, file = "sigDE.id.celltype7.chitin.up.rds")
# saveRDS(chitin.down.gene, file = "sigDE.id.celltype7.chitin.down.rds")

####-------------------- edp1 sigDE up and down gene list -------------------------------
# DE.SCREW1 <- readRDS("DE.celltype7.edp1.rds")

degUp  <- list()
degDown  <- list()
for (i in names(DE.SCREW1)){
  up <- DE.SCREW1[[i]][DE.SCREW1[[i]]$p_val_adj < 0.05 & DE.SCREW1[[i]]$avg_logFC >= 0.25,]
  down <- DE.SCREW1[[i]][DE.SCREW1[[i]]$p_val_adj < 0.05 & DE.SCREW1[[i]]$avg_logFC <= (-0.25),]
  degUp[[i]] <- row.names(up)
  degDown[[i]] <- row.names(down)
}
degUp
degDown

edp1.up.gene <- degUp
edp1.down.gene <- degDown

# saveRDS(edp1.up.gene, file = "sigDE.id.celltype7.edp1.up.rds")
# saveRDS(edp1.down.gene, file = "sigDE.id.celltype7.edp1.down.rds")

####-------------------- scoop10 sigDE up and down gene list -------------------------------
# DE.SCOOP10 <- readRDS("DE.celltype7.scoop10.rds")

degUp  <- list()
degDown  <- list()
for (i in names(DE.SCOOP10)){
  up <- DE.SCOOP10[[i]][DE.SCOOP10[[i]]$p_val_adj < 0.05 & DE.SCOOP10[[i]]$avg_logFC >= 0.25,]
  down <- DE.SCOOP10[[i]][DE.SCOOP10[[i]]$p_val_adj < 0.05 & DE.SCOOP10[[i]]$avg_logFC <= (-0.25),]
  degUp[[i]] <- row.names(up)
  degDown[[i]] <- row.names(down)
}
degUp
degDown

scoop10.up.gene <- degUp
scoop10.down.gene <- degDown

# saveRDS(scoop10.up.gene, file = "sigDE.id.celltype7.scoop10.up.rds")
# saveRDS(scoop10.down.gene, file = "sigDE.id.celltype7.scoop10.down.rds")
##############################################################################
DE_results$sig.DE.flg22.up <- flg22.up.gene
DE_results$sig.DE.flg22.down <- flg22.down.gene

DE_results$sig.DE.Chitin.up <- chitin.up.gene
DE_results$sig.DE.Chitin.down <- chitin.down.gene

DE_results$sig.DE.Pep1.up <- pep1.up.gene
DE_results$sig.DE.Pep1.down <- pep1.down.gene

DE_results$sig.DE.SCOOP10.up <- scoop10.up.gene
DE_results$sig.DE.SCOOP10.down <- scoop10.down.gene

DE_results$sig.DE.SCREW1.up <- edp1.up.gene
DE_results$sig.DE.SCREW1.down <- edp1.down.gene

saveRDS(DE_results, file = "DE.rootzone4.celltype6.results.rds")
