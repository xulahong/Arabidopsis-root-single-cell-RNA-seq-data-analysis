library(clusterProfiler)
library(RColorBrewer)
DE_result <- readRDS("./DE.rootzone4.celltype6.results.rds")

sigdeg.flg22 <- DE_result[["sig.DE.flg22"]]
sigdeg.chitin <- DE_result[["sig.DE.Chitin"]]
sigdeg.pep1 <- DE_result[["sig.DE.Pep1"]]
sigdeg.scoop10 <- DE_result[["sig.DE.SCOOP10"]]
sigdeg.edp1 <- DE_result[["sig.DE.SCREW1"]]

deg <- c(sigdeg.flg22, sigdeg.chitin, sigdeg.pep1, sigdeg.scoop10, sigdeg.edp1)

go_result <- list()
for (i in names(deg)) {
  gene <- deg[[i]]
  ego <- enrichGO(gene, OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")
  go_result[[i]] <- ego@result[ego@result$pvalue < 0.01, ]
}


go_term <- list()
for (i in names(go_result)) {
  pvalue  <- go_result[[i]][, c(1,5)]
  colnames(pvalue)[2] <- paste0(i, "_", colnames(pvalue)[2])
  go_term[[i]] <- pvalue
}

df <- Reduce(function(...) merge(..., all=TRUE), go_term)

go_result_all <- do.call("rbind", go_result)
go_result_all <- go_result_all[, 1:2]
b <- go_result_all
b= b[!duplicated(b[,1]), ]
b = as.data.frame(b[,-c(3:6)])
head(b)

library(dplyr)
df.merge <- left_join(df, b, by="ID")
df.merge[is.na(df.merge)] <- 1
head(df.merge)
df.merge[, 2:64] <- -log10(df.merge[, 2:64])
head(df.merge)

df.merge <- df.merge[rowSums(df.merge[, 2:64])>15,]

library(tidyr)
df.gather <- gather(df.merge, celltype, pvalue, 2:64)
head(df.gather)

### top 10
top10 <- df.gather %>% group_by(celltype) %>% top_n(n = 3, wt = pvalue)
head(top10)
length(unique(top10$ID))

df <- df.merge[which(df.merge$ID %in% unique(top10$ID)), -1]
head(df)
row.names(df) <- df$Description
df <- df[, -64]

# color = colorRampPalette(rev(brewer.pal(10, "PiYG")))(10)

library(pheatmap)
palette = colorRampPalette(brewer.pal(n = 4, name = "RdPu"))(100)
color = c("#d4d4d2", palette[1:6], palette[60:100])
pheatmap(df, color = color, border_color = NA, clustering_method = "ward.D2")

# col_group2 <- substr(colnames(pep1.FC), 1, 3)
col_group <- strsplit(colnames(df), "[_]")
col_group <- sapply(col_group,function(x) x[1])
col_group1 <- strsplit(colnames(df), "[_]")
col_group1 <- sapply(col_group1,function(x) x[2])
col_group2 <- strsplit(colnames(df), "[_]")
col_group2 <- sapply(col_group2,function(x) x[3])

# Data frame with column annotations.
mat_col <- data.frame(Type= col_group1, Zone=col_group, Treatment= col_group2)
rownames(mat_col) <- colnames(df)
mat_col$Zone <- gsub("Rootcap", "Root cap", mat_col$Zone)

table(mat_col$Type)
mat_col$Type <- factor(mat_col$Type, level=c("Epidermis", "Cortex", "Endodermis", "Stele", "LRC", "Columella"))
table(mat_col$Treatment)
mat_col$Treatment <- factor(mat_col$Treatment, level=c("flg22", "Chitin", "Pep1", "SCOOP10", "SCREW1"))
table(mat_col$Zone)
mat_col$Zone <- factor(mat_col$Zone, level=c("Root cap", "Meristem", "Elongation", "Maturation"))

# List with colors for each annotation.
type_colors = c("#377EB8", "#E41A1C", "#984EA3", "#FF7F00", "#4DAF4A", "#A65628")
zone_colors = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")
mat_colors <- list(Type = type_colors, Zone = zone_colors,  Treatment= brewer.pal(5, "Dark2"))
                                     
# mat_colors <- list(Type = brewer.pal(7, "Set1"), Treatment= c(rep("#FFFFFF", times = 5)))

names(mat_colors$Zone) <- unique(mat_col$Zone)
names(mat_colors$Type) <- unique(mat_col$Type)
names(mat_colors$Treatment) <- unique(mat_col$Treatment)
mat_colors

p <- pheatmap(
  mat               = as.matrix(df),
  gaps_col          = c(12,22,36,49),
  border_color      = NA,
  color             = color,
  show_colnames     = FALSE,
  show_rownames     = T,
  cluster_cols      = FALSE, 
  cluster_rows      = TRUE, 
  scale             = "none",
  clustering_method = "ward.D2",
  annotation_col    = mat_col,
  annotation_colors = mat_colors,
  fontsize          = 8,
  main              = " ",
  legend = T,
  annotation_names_col = F,
  annotation_legend = T
)


p


### save to pdf
pdf(file = 'heatmap_GO.pdf',  width = 12, height = 25)
print(p)
dev.off()

