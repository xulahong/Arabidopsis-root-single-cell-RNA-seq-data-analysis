setwd("~/lahong/scRNA-seq/scRNA_seq_all_2/12.GRNs/")
library(dplyr)
library(tidyr)
library(tidyverse)
library(plyr)
###############################################
if(T){
  load("./TF.FunTFBS.RData")
  
  pep1.grn <- read.csv("pep1/pep1_regulons.csv")
  colnames(pep1.grn) <- c('TF1', 'targets', 'cluster')
  head(pep1.grn)
  
  ###########
  pep1.TF.FunTFBS <- TF.FunTFBS[["pep1.TF.FunTFBS"]]
  pep1.TF.FunTFBS <- pep1.TF.FunTFBS[, c( "TF", "Gene", "Cell_type")]
  head(pep1.TF.FunTFBS)
  
  a <- left_join(pep1.grn, pep1.TF.FunTFBS, by= c("TF1" = "TF", "targets" = "Gene"))
  a <- na.omit(a)
  head(a)
  keep <- a[which(a$cluster == a$Cell_type),]
  length(unique(keep$TF))
  keep <- keep[, -4]
  deduped.data <- unique( keep[ , 1:3 ] )
  pep1.grn <- deduped.data
  
  ###############################################
  load("./TF.FunTFBS.RData")
  
  scoop10.grn <- read.csv("scoop10/scoop10_regulons.csv")
  colnames(scoop10.grn) <- c('TF1', 'targets', 'cluster')
  head(scoop10.grn)
  
  ###########
  scoop10.TF.FunTFBS <- TF.FunTFBS[["scoop10.TF.FunTFBS"]]
  scoop10.TF.FunTFBS <- scoop10.TF.FunTFBS[, c( "TF", "Gene", "Cell_type")]
  head(scoop10.TF.FunTFBS)
  
  a <- left_join(scoop10.grn, scoop10.TF.FunTFBS, by= c("TF1" = "TF", "targets" = "Gene"))
  a <- na.omit(a)
  head(a)
  keep <- a[which(a$cluster == a$Cell_type),]
  length(unique(keep$TF))
  keep <- keep[, -4]
  deduped.data <- unique( keep[ , 1:3 ] )
  scoop10.grn <- deduped.data
  
  ###############################################
  load("./TF.FunTFBS.RData")
  
  screw1.grn <- read.csv("screw1/screw1_regulons.csv")
  colnames(screw1.grn) <- c('TF1', 'targets', 'cluster')
  head(screw1.grn)
  
  ###########
  screw1.TF.FunTFBS <- TF.FunTFBS[["screw1.TF.FunTFBS"]]
  screw1.TF.FunTFBS <- screw1.TF.FunTFBS[, c( "TF", "Gene", "Cell_type")]
  head(screw1.TF.FunTFBS)
  
  a <- left_join(screw1.grn, screw1.TF.FunTFBS, by= c("TF1" = "TF", "targets" = "Gene"))
  a <- na.omit(a)
  head(a)
  keep <- a[which(a$cluster == a$Cell_type),]
  length(unique(keep$TF))
  keep <- keep[, -4]
  deduped.data <- unique( keep[ , 1:3 ] )
  screw1.grn <- deduped.data
}

#############################################################################################################################################
if(T){
  pep1.grn.list <- list()
  pep1.grn.list.2 <- list()
  
  for (i in pep1.grn$cluster){
    x <- pep1.grn[grep(i, pep1.grn$cluster), ]
    
    for (j in x$targets) {
      a <- x[grep(j, x$targets), 1]
      pep1.grn.list[[i]][[j]] <- length(a)
    }
    
    y <- t(as.data.frame(pep1.grn.list[[i]]))
    colnames(y) <- "TF"
    y <- as.data.frame(y)
    y$targets <- row.names(y)
    pep1.grn.list.2[[i]] <- y
  }
  
  df <- ldply (pep1.grn.list.2, data.frame)
  colnames(df) <- c("Celltype", "TF", "Target")
  df$Treatment <- "pep1"
  head(df)
  df.pep1 <- df
}

if(T){
  scoop10.grn.list <- list()
  scoop10.grn.list.2 <- list()
  
  for (i in scoop10.grn$cluster){
    x <- scoop10.grn[grep(i, scoop10.grn$cluster), ]
    
    for (j in x$targets) {
      a <- x[grep(j, x$targets), 1]
      scoop10.grn.list[[i]][[j]] <- length(a)
    }
    
    y <- t(as.data.frame(scoop10.grn.list[[i]]))
    colnames(y) <- "TF"
    y <- as.data.frame(y)
    y$targets <- row.names(y)
    scoop10.grn.list.2[[i]] <- y
  }
  
  df <- ldply (scoop10.grn.list.2, data.frame)
  colnames(df) <- c("Celltype", "TF", "Target")
  df$Treatment <- "scoop10"
  head(df)
  df.scoop10 <- df
}

if(T){
  screw1.grn.list <- list()
  screw1.grn.list.2 <- list()
  
  for (i in screw1.grn$cluster){
    x <- screw1.grn[grep(i, screw1.grn$cluster), ]
    
    for (j in x$targets) {
      a <- x[grep(j, x$targets), 1]
      screw1.grn.list[[i]][[j]] <- length(a)
    }
    
    y <- t(as.data.frame(screw1.grn.list[[i]]))
    colnames(y) <- "TF"
    y <- as.data.frame(y)
    y$targets <- row.names(y)
    screw1.grn.list.2[[i]] <- y
  }
  
  df <- ldply (screw1.grn.list.2, data.frame)
  colnames(df) <- c("Celltype", "TF", "Target")
  df$Treatment <- "screw1"
  head(df)
  df.screw1 <- df
}

#######################################
head(df.pep1)
head(df.scoop10)
head(df.screw1)

df <- rbind(df.pep1, df.scoop10, df.screw1)
df$celltype <- paste0(df$Treatment, "_", df$Celltype)
head(df)
df <- df[, c(2,3,5)]
df <- spread(df, celltype, TF)

df[is.na(df)] <- 0
head(df)

save(df, file = "DEGs.Number.of.TFs.RData")










