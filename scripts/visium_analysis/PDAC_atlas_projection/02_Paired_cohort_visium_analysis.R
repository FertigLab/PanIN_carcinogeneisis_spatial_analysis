title: "Paired_cohort_visium_analysis"
author: "Alexander Bell"
date: "2/3/2024"

library(pracma)
library(plyr)
library(dplyr)
library(spatstat)
library(Seurat)
library(tidyverse)
library(qvalue)
library(spatialEco)
library(UpSetR)
library(harmony)
library(ComplexHeatmap)
library(BBmisc)
library(ggpubr)
library(ggplot2)
library(CoGAPS)
library(forcats)
library(monocle3)
library(Seurat)
library(projectR)
library(dittoSeq)
library(rstatix)
source("..Custom_Functions.R")
data_merged <- readRDS("..objects_and_metadata/full_unprepped_object.rds")
prepped <- readRDS(file = "..objects_and_metadata/prepped_full_object.rds")
sessionInfo()
# BAR PLOT OF CODA ANNOTATIONS BY LOUVAIN CLUSTER (Supplemental 10B) ####
my_list <- list()
for(x in 1:length(unique(data_merged$post_harmony_clusters))){
  data_merged@meta.data %>% filter(post_harmony_clusters == x-1) %>% select("islets", "normal epithelium", "smooth muscle", "fat", "acini", "collagen", "panin") -> clust
  colSums(clust) -> clustcolsums
  my_list[[x]] <- clustcolsums/sum(clustcolsums)
}

tmplist <- list()
for(x in 1:length(my_list)){
  as.data.frame((my_list[[x]])) -> tmp
  tmp$celltype <- rownames(tmp)
  colnames(tmp) <- c("percent", "celltype")
  #tmp$percent <- tmp$percent*100
  tmp$cluster <- as.character(x-1)
  tmplist[[x]] <- tmp
}

final_table <- do.call(rbind, tmplist[0:15])

final_table$cluster <- factor(final_table$cluster, levels=c('0','1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14'))

pdf("..PreProcessing/Louvain_CODA_Barplot.pdf")
ggplot(final_table, aes(fill=celltype, y=percent, x=cluster)) + 
  geom_bar(position="fill", stat="identity") + scale_fill_manual(values = c("#9523B8", "#FFC2F5", "#FFFF00", "#79F8FC","#0000FF", "#FF0000", "#50ED50"))
dev.off()

# BAR PLOTS OF CODA ANNOTATION BEFORE AND AFTER CODA FILTRATION (Supplemental 10C) ####
my_list <- list()
for(x in 1:length(unique(prepped$post_harmony_clusters))){
  clust <- prepped@meta.data %>% filter(post_harmony_clusters == x-1) %>% select("islets", "normal epithelium", "smooth muscle", "fat", "acini", "collagen", "panin")
  clustcolsums <- colSums(clust)
  my_list[[x]] <- clustcolsums/sum(clustcolsums)
}

tmplist <- list()
for(x in 1:length(my_list)){
  tmp <- as.data.frame((my_list[[x]]))
  tmp$celltype <- rownames(tmp)
  colnames(tmp) <- c("percent", "celltype")
  #tmp$percent <- tmp$percent*100
  tmp$cluster <- as.character(x-1)
  tmplist[[x]] <- tmp
}


final_table <- do.call(rbind, tmplist[0:15])
final_table$cluster <- factor(final_table$cluster, levels=c('0','1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14'))
final_table_trimmed <- final_table %>% filter(cluster == 6 | cluster == 8 | cluster == 11)


my_list <- list()
for(x in 1:length(unique(prepped$trueidents))){
  tmp <- unique(prepped$trueidents)[x]
  clust <- prepped@meta.data %>% filter(trueidents == tmp) %>% select("islets", "normal epithelium", "smooth muscle", "fat", "acini", "collagen", "panin")
  clustcolsums <- colSums(clust)
  my_list[[x]] <- clustcolsums/sum(clustcolsums)
}

tmplist <- list()
for(x in 1:length(my_list)){
  tmp <- as.data.frame((my_list[[x]]))
  tmp$celltype <- rownames(tmp)
  colnames(tmp) <- c("percent", "celltype")
  #tmp$percent <- tmp$percent*100
  tmp$cluster <- unique(prepped$trueidents)[x]
  tmplist[[x]] <- tmp
}

final_table_CODA <- do.call(rbind, tmplist[1:4])
final_table_CODA$cluster <- factor(final_table_CODA$cluster, levels=c("nocells", "trueneoplasia", "trueduct", "trueislets"))
final_table_CODA_trimmed <- final_table_CODA %>% filter(cluster == "trueduct" | cluster == "trueneoplasia" | cluster == "trueislets")

pdf("..Louvain_CODA_Barplot_before_CODA.pdf")
ggplot(final_table_trimmed, aes(fill=celltype, y=percent, x=factor(cluster, levels = c(8,11,6)))) + 
  geom_bar(position="fill", stat="identity") + scale_fill_manual(values = c("#9523B8", "#FFC2F5", "#FFFF00", "#79F8FC","#0000FF", "#FF0000", "#50ED50")) +
  labs(fill = "CODA Annotation", x = "", y = "Percentage") + 
  scale_x_discrete(labels=c("8" = "Cluster 8", "11" = "Cluster 11", "6" = "Cluster 6"))
dev.off()

pdf("..Louvain_CODA_Barplot_after_CODA.pdf")
ggplot(final_table_CODA_trimmed, aes(fill=celltype, y=percent, x=factor(cluster, levels = c("trueduct", "trueneoplasia", "trueislets")))) + 
  geom_bar(position="fill", stat="identity") + scale_fill_manual(values = c("#9523B8", "#FFC2F5", "#FFFF00", "#79F8FC","#0000FF", "#FF0000", "#50ED50")) +
  labs(fill = "CODA Annotation", x = "", y = "Percentage") +scale_x_discrete(labels=c("trueduct" = "Cluster 8", "trueneoplasia" = "Cluster 11", "trueislets" = "Cluster 6"))
dev.off()

# BOX PLOTS OF GROUP MARKER GENES BEFORE AND AFTER CODA FILTRATION (Supplemental 10D) ####

genelist <- c("INS", "PRSS1", "CTSE")
trueidentsbarcodes <- rownames(prepped@meta.data %>% filter(trueidents == "trueislets" | trueidents == "trueduct"| trueidents == "trueneoplasia"))
seuratclusterbarcodes <- rownames(prepped@meta.data %>% filter(post_harmony_clusters == 6 | post_harmony_clusters == 8 | post_harmony_clusters == 11))
trueidentslist <- prepped@meta.data %>% filter(trueidents == "trueislets" | trueidents == "trueduct"| trueidents == "trueneoplasia") %>% dplyr::select(trueidents) %>% as.vector() %>% unlist() %>% as.vector()
seuratclusterlist <- prepped@meta.data %>% filter(post_harmony_clusters == 6 | post_harmony_clusters == 8 | post_harmony_clusters == 11) %>% dplyr::select(post_harmony_clusters) %>% as.vector() %>% unlist() %>% as.vector()

# making a vector of labels generated by intersecting CODA and Louvain clustering
trueidentslist <- str_replace(trueidentslist, "trueduct", "Duct")
trueidentslist <- str_replace(trueidentslist, "trueneoplasia", "Neoplasia")
trueidentslist <- str_replace(trueidentslist, "trueislets", "Islet")

# making a vector of labels bgenerated by Louvain clustering alone (before CODA)
seuratclusterlist <- str_replace(seuratclusterlist, "6", "Islet")
seuratclusterlist <- str_replace(seuratclusterlist, "8", "Duct")
seuratclusterlist <- str_replace(seuratclusterlist, "11", "Neoplasia")

# making a dataframe of the SCT data
dfprepped <- prepped@assays$SCT@data %>% as.data.frame()

# making tables for the genes of interest and the spots of interest
trueidentstable <- dfprepped[genelist, trueidentsbarcodes]
seuratclustertable <- dfprepped[genelist, seuratclusterbarcodes]

for(loopedgene in genelist){
  the_group <- c(rep("Before CODA", ncol(seuratclustertable)), rep("After CODA", ncol(trueidentstable)))
  the_identities <- c(seuratclusterlist, trueidentslist)
  the_gene <- as.numeric(unlist(c(seuratclustertable[loopedgene,], trueidentstable[loopedgene,])))
  the_table <- as.data.frame(cbind(the_group, the_gene, the_identities))
  colnames(the_table) <- c("Group", "Gene", "Identity")
  the_table$Gene <- as.numeric(the_table$Gene)
  the_table$Identity <- factor(the_table$Identity, levels = c("Duct", "Neoplasia", "Islet"))
  the_table$Group <- factor(the_table$Group, levels = c("Before CODA", "After CODA"))
  
  stat.test <- the_table %>%
    group_by(Identity) %>%
    wilcox_test(Gene ~ Group)
  
  stat.test <- stat.test %>%
    add_xy_position(x = "Identity", dodge = 0.8)
  pdf(paste0("..geneplots/", loopedgene,".pdf"))
  print(ggplot(the_table, aes(x=factor(Identity, levels = c("Duct", "Neoplasia", "Islet")), y=Gene, color = factor(Group, levels = c("Before CODA", "After CODA")))) + 
          geom_boxplot() +
          ggtitle(label = paste(loopedgene)) +
          theme(plot.title = element_text(hjust = 0.5)) +
          labs(color = "Spot filtration", x = "", y = "Normalized count") + 
          stat_pvalue_manual(
            stat.test,  label = "p", tip.length = 0))
  dev.off()
}

# SPATIAL PLOTS SHOWING IDENTITIES BEFORE AND AFTER CODA FILTRATION (Supplemental 10A) ####
prepped <- SetIdent(prepped, value = "trueidents")  
cleanedpanin <- CellsByIdentities(prepped, idents = c("trueneoplasia"))
cleanedduct <- CellsByIdentities(prepped, idents = c("trueduct"))

prepped@meta.data %>% filter(trueidents == "trueneoplasia" & paningrade == "high_grade") %>% nrow()
prepped@meta.data %>% filter(trueidents == "trueneoplasia" & paningrade == "low_grade") %>% nrow()
pdf("..filtered_panin.pdf")
SpatialDimPlot(prepped, images = "slice1.6", cells.highlight = cleanedpanin, cols.highlight = c("red", "#363838"),
               pt.size.factor = 2.2)
dev.off()

pdf("..filtered_duct.pdf")
SpatialDimPlot(prepped, images = "slice1.6", cells.highlight = cleanedduct, cols.highlight = c("red", "#363838"),
               pt.size.factor = 2.2)
dev.off()

prepped <- SetIdent(prepped, value = "post_harmony_clusters")  
unfilteredpanin <- CellsByIdentities(prepped, idents = c(11))
unfilteredduct <- CellsByIdentities(prepped, idents = c(8))

pdf("..unfiltered_panin.pdf")
SpatialDimPlot(prepped, images = "slice1.6", cells.highlight = unfilteredpanin, cols.highlight = c("red", "#363838"),
               pt.size.factor = 2.2)
dev.off()

pdf("..unfiltered_duct.pdf")
SpatialDimPlot(prepped, images = "slice1.6", cells.highlight = unfilteredduct, cols.highlight = c("red", "#363838"),
               pt.size.factor = 2.2)
dev.off()

# LOUVAIN CLUSTER SPATIAL PLOTS (Supplemental 3-9) ####
pdf("..Spatial_Louvain.pdf")
SpatialDimPlot(data_merged, group.by = "post_harmony_clusters", images = "slice1" )
SpatialDimPlot(data_merged, group.by = "post_harmony_clusters", images = "slice1.1" )
SpatialDimPlot(data_merged, group.by = "post_harmony_clusters", images = "slice1.2" )
SpatialDimPlot(data_merged, group.by = "post_harmony_clusters", images = "slice1.3" )
SpatialDimPlot(data_merged, group.by = "post_harmony_clusters", images = "slice1.4" )
SpatialDimPlot(data_merged, group.by = "post_harmony_clusters", images = "slice1.5" )
SpatialDimPlot(data_merged, group.by = "post_harmony_clusters", images = "slice1.6" )
dev.off()

data_merged <- SetIdent(data_merged, value = "post_harmony_clusters")

# HEATMAP OF LOUVAIN CLUSTER MARKER GENES (Supplemental 2) ####
# Get markers for each post-harmony louvain cluster
markers <- FindAllMarkers(prepped, test.use = 'MAST', only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
# Get top 10 marker genes per cluster
markers_df <- as.data.frame(markers %>% group_by(cluster) %>% slice_max(n = 10, order_by = avg_log2FC))
# Make heatmap
hm <- DoHeatmap(prepped,features = markers_df$gene) + theme(text = element_text(size = 4))
pdf(file= "../heatmaps/Louvain_Clusters.pdf", width = 14)
hm
dev.off()

# panCAF MODULE SCORE SPATIAL PLOTS (Main figure 2B) ####

# gene list
pancaf_genes <- list(c('LUM','DCN','COL1A1','VIM','CD39','FSP1','FAP','ACTA2','PDPN'))

# make module score
prepped <- AddModuleScore(prepped, features = pancaf_genes, name = 'panCAF')

# make spatial plots
pdf("..panCAF_spatial.pdf")
SpatialFeaturePlot(prepped, features = 'panCAF1', images = 'slice1', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'panCAF1', images = 'slice1.1', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'panCAF1', images = 'slice1.2', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'panCAF1', images = 'slice1.3', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'panCAF1', images = 'slice1.4', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'panCAF1', images = 'slice1.5', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'panCAF1', images = 'slice1.6', pt.size.factor = 2.3)
dev.off()

#iCAF MODULE SCORE SPATIAL PLOTS (Main figure 2D) ####

# gene list
iCAF_genes <- list(c('CXCL1','CXCL2','CCL2','CXCL12','PDGFRA','CFD','LMNA','DPT','HAS1','HAS2'))

# make module score
prepped <- AddModuleScore(prepped, features = iCAF_genes, name = 'iCAF')

# make spatial plots
pdf("..iCAF_spatial.pdf")
SpatialFeaturePlot(prepped, features = 'iCAF1', images = 'slice1', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'iCAF1', images = 'slice1.1', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'iCAF1', images = 'slice1.2', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'iCAF1', images = 'slice1.3', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'iCAF1', images = 'slice1.4', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'iCAF1', images = 'slice1.5', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'iCAF1', images = 'slice1.6', pt.size.factor = 2.3)
dev.off()

# myCAF MODULE SCORE SPATIAL PLOTS (Main figure 2C) ####

# gene list
myCAF_genes <- list(c('TAGLN','MYL9','TPM2','MMP11','POSTN','HOPX','TWIST1','SOX4'))

# make module score
prepped <- AddModuleScore(prepped, features = myCAF_genes, name = 'myCAF')

# make spatial plots
pdf("..myCAF_spatial.pdf")
SpatialFeaturePlot(prepped, features = 'myCAF1', images = 'slice1', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'myCAF1', images = 'slice1.1', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'myCAF1', images = 'slice1.2', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'myCAF1', images = 'slice1.3', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'myCAF1', images = 'slice1.4', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'myCAF1', images = 'slice1.5', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'myCAF1', images = 'slice1.6', pt.size.factor = 2.3)
dev.off()

# apCAF MODULE SCORE SPATIAL PLOTS (Main figure 2E) ####

# gene list
apCAF_genes <- list(c('HLA-DRA','HLA-DPA1','CD74','HLA-DQ','SLPI'))

# make module score
prepped <- AddModuleScore(prepped, features = apCAF_genes, name = 'apCAF')

# make spatial plots
pdf("..apCAF_spatial.pdf")
SpatialFeaturePlot(prepped, features = 'apCAF1', images = 'slice1', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'apCAF1', images = 'slice1.1', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'apCAF1', images = 'slice1.2', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'apCAF1', images = 'slice1.3', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'apCAF1', images = 'slice1.4', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'apCAF1', images = 'slice1.5', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'apCAF1', images = 'slice1.6', pt.size.factor = 2.3)
dev.off()

# PTPRC (CD45) SPATIAL PLOTS (Main figure 2F) ####

# make spatial plots
pdf("..PTPRC_spatial.pdf")
SpatialFeaturePlot(prepped, features = 'PTPRC', images = 'slice1', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'PTPRC', images = 'slice1.1', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'PTPRC', images = 'slice1.2', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'PTPRC', images = 'slice1.3', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'PTPRC', images = 'slice1.4', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'PTPRC', images = 'slice1.5', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'PTPRC', images = 'slice1.6', pt.size.factor = 2.3)
dev.off()

# CLASSICAL PDAC MODULE SCORE SPATIAL PLOTS (Main figure 4A) ####

# classical subtype gene list
classical_genes <- list(c('BTNL8','FAM3D','ATAD4','AGR3','CTSE','LOC400573','LYZ','TFF2',
                          'TFF1','ANXA10','LGALS4','PLA2G10','CEACAM6','VSIG2','TSPAN8',
                          'ST6GALNAC1','AGR2','TFF3','CYP3A7','MYO1A','CLRN3','KRT20',
                          'CDH17','SPINK4','REG4'))

# make module score
prepped <- AddModuleScore(prepped, features = classical_genes, name = 'classical')

# make spatial plots
pdf("..classical_spatial.pdf")
SpatialFeaturePlot(prepped, features = 'classical1', images = 'slice1', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'classical1', images = 'slice1.1', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'classical1', images = 'slice1.2', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'classical1', images = 'slice1.3', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'classical1', images = 'slice1.4', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'classical1', images = 'slice1.5', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'classical1', images = 'slice1.6', pt.size.factor = 2.3)
dev.off()

# BASAL PDAC MODULE SCORE SPATIAL PLOTS (Main figure 4B) ####

# Basal subtype gene list
basal_genes <- list(c('VGLL','UCA1','S100A2','LY6D','SPRR3','SPRR1B','LEMD1','KRT15',
                      'CTSL2','DHRS9','AREG','CST6','SERPINB3','KRT6C','KRT6A','SERPINB4',
                      'FAM83A','SCEL','FGFBP1','KRT7','KRT17','GPR87','TNS4','SLC2A1',
                      'ANXA8L2'))

# make module score
prepped <- AddModuleScore(prepped, features = basal_genes, name = 'basal')

# make spatial plots
pdf("..basal_spatial.pdf")
SpatialFeaturePlot(prepped, features = 'basal1', images = 'slice1', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'basal1', images = 'slice1.1', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'basal1', images = 'slice1.2', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'basal1', images = 'slice1.3', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'basal1', images = 'slice1.4', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'basal1', images = 'slice1.5', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'basal1', images = 'slice1.6', pt.size.factor = 2.3)
dev.off()

# CANCER STEM CELL MODULE SCORE SPATIAL PLOTs (Main figure 4C) ####

# cancer stem cell gene list
csc_genes <- list(c('ABCG2','ALDH1A1','CD24','CD44','EPCAM','PROM1','CXCR4','NES','DCLK1','SOX9','NANOG'))

# make module score
prepped <- AddModuleScore(prepped, features = csc_genes, name = 'CSC')

# make spatial plots
pdf("CSC_spatial.pdf")
SpatialFeaturePlot(prepped, features = 'CSC1', images = 'slice1', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'CSC1', images = 'slice1.1', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'CSC1', images = 'slice1.2', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'CSC1', images = 'slice1.3', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'CSC1', images = 'slice1.4', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'CSC1', images = 'slice1.5', pt.size.factor = 2.3)
SpatialFeaturePlot(prepped, features = 'CSC1', images = 'slice1.6', pt.size.factor = 2.3)
dev.off()

# NEOPLASIA VS DUCTAL DIFFERENTIAL EXPRESSION (Main figure 4D) ####

# perform differential gene expression
trueneoduct <- FindMarkers(prepped, group.by = "trueidents", ident.1 = "trueneoplasia", ident.2 = "trueduct", test.use = "MAST", latent.vars = "Slide_ID")

# make volcano plot (volcanoPlot from Custom_Functions)
volcanoPlot(trueneoduct, title_1 = "Neoplasia", title_2 = "Duct", batch = "CODA_Louvain", directory = "../neo_vs_duct/")

# HEATMAP OF NEOPLASIA VS DUCTAL GENES (Main figure 4E) ####

# this file is generated by running the volcanoPlot function (figure 3D)
neo_duct_degs <- read.csv("../neo_vs_duct/CODA_Louvain_NeoplasiavsDuct_DEGS.csv")

# top 20 DEGs
top20 <- head(neo_duct_degs, n = 20)

# bottom 20 DEGs
bottom20 <- tail(neo_duct_degs, n = 20)

# combine DEGs
neo_duct_degs <- rbind(top20, bottom20)

# subset seurat object for only neoplasia and ducts
prepped_subs <- subset(prepped, subset = trueidents %in% c("trueneoplasia", "trueduct"))

# data slot for the subset seurat object
the_table <- prepped_subs@assays$SCT@data

# making table of significant degs
the_table[rownames(the_table) %in% neo_duct_degs$Genes,] -> significant_degs_table

# making data for heatmap
datScale <- t(apply(significant_degs_table,1,scale))

### making metadata for the heatmap

# columns
celltype <- prepped_subs@meta.data$trueidents
slidenumber <- prepped_subs@meta.data$Slide_Number
qc <- as.data.frame(cbind(celltype, slidenumber))

# rows
neoductgenes <- neo_duct_degs %>% select("Genes", "Group")
neoductgenes$Group[neoductgenes$Group == "Uppers"] <- "Neoplasia_DEGs"
neoductgenes$Group[neoductgenes$Group == "Lowers"] <- "Duct_DEGs"

# Making legend for heatmap
qc$slidenumber <- factor(qc$slidenumber, levels = c("slide_113_1", "slide_113_2", "slide_113_3", "slide_113_4", "slide_114_1", "slide_114_2", "slide_114_4"))
qc$tissue[qc$celltype == "trueneoplasia"] <- "Neoplasia"
qc$tissue[qc$celltype == "trueduct"] <- "Duct"
column_order = rev(order(qc$tissue))
row_order <- match(neoductgenes$Genes, rownames(datScale))
datScale <- datScale[row_order,]
row_order <- 1:length(rownames(datScale))

# Generating and saving heatmap
pdf("..neoplasia_duct_heatmap_top20.pdf")
Heatmap(datScale, show_row_names = T, show_column_names = F,
        clustering_distance_rows = 'pearson',
        row_names_gp = grid::gpar(fontsize = 6),
        column_order = column_order,
        row_order = row_order,
        top_annotation = columnAnnotation(df = qc[,c("tissue",'slidenumber')], col = list(tissue = c("Duct" = "#0000FF", "Neoplasia" = "#FF0000"))),
        left_annotation = rowAnnotation(Genes = neoductgenes[,c("Group")], col = list(Genes = c("Duct_DEGs" = "#0000FF", "Neoplasia_DEGs" = "#FF0000"))))
dev.off()

# HIGH GRADE VS LOW GRADE DIFFERENTIAL GENE EXPRESSION ANALYSIS VOLCANO PLOT (Supplemental figure 16) ####

# differential gene expression
degtable <- FindMarkers(prepped, group.by = "paningrade", ident.1 = "high_grade", ident.2 = "low_grade", test.use = "MAST", latent.vars = "Slide_ID")

# make volcano plots
volcanoPlot(degtable, title_1 = "High Grade", title_2 = "Low Grade", batch = "All_Slides_Louvain", directory = "..hg_vs_lg/")

# HIGH/LOW GRADE PANIN SPATIAL PLOT (Main figure 6E) ####

# paningrade is the ident that shows HG and LG spots
data_merged <- SetIdent(data_merged, value = "paningrade")

# plot HG and LG
pdf("..paningrade.pdf")
SpatialDimPlot(data_merged, cells.highlight = CellsByIdentities(data_merged, idents = c("high_grade", "low_grade", "normal")), images = "slice1.6", cols.highlight = c("#FF0000", "#FFAF2D", "#FEEFA7", "#D4D4D4"))
dev.off()

# MUCL3 SPATIAL PLOT (Main figure 6G) ####
pdf("..MUCL3.pdf")
SpatialFeaturePlot(data_merged, features = "MUCL3", images = "slice1.6")
dev.off()
# TSPAN1 SPATIAL PLOT (Main figure 6F) ####
pdf("..TSPAN1.pdf")
SpatialFeaturePlot(data_merged, features = "TSPAN1", images = "slice1.6")
dev.off()
# TFF1 SPATIAL PLOT (Supplemental figure 15, Main figure 6H) ####
pdf("..TFF1.pdf")
SpatialFeaturePlot(data_merged, features = "TFF1", images = "slice1")
SpatialFeaturePlot(data_merged, features = "TFF1", images = "slice1.1")
SpatialFeaturePlot(data_merged, features = "TFF1", images = "slice1.2")
SpatialFeaturePlot(data_merged, features = "TFF1", images = "slice1.3")
SpatialFeaturePlot(data_merged, features = "TFF1", images = "slice1.4")
SpatialFeaturePlot(data_merged, features = "TFF1", images = "slice1.5")
SpatialFeaturePlot(data_merged, features = "TFF1", images = "slice1.6")
dev.off()

# GENE SET ENRICHMENT ANALYSIS OF DIFFERENTIAL EXPRESSION RESULTS (Supplemental 14) ####

### PREPPING PRE-RANKED GENE LISTS FOR GSEA 

# optional importing of data to this step 
trueneoduct <- FindMarkers(prepped, group.by = "trueidents", ident.1 = "trueneoplasia", ident.2 = "trueduct", test.use = "MAST", latent.vars = "Slide_ID", logfc.threshold = 0)

# ranking by statistic and direction
statrank <- trueneoduct %>% mutate(neg.log.padj = -1*log10(p_val_adj)) %>% mutate(statistic = avg_log2FC/abs(avg_log2FC) * neg.log.padj)
statrank <- statrank[order(statrank$statistic, decreasing = T), ] %>% select(statistic)

# ranking by LFC
lfcrank <- trueneoduct[order(trueneoduct$avg_log2FC, decreasing = T),] %>% select(avg_log2FC)

# this input is used for GSEA GUI with the following settings: pre-ranked by LFC, weighted, collapsed.
#write.csv(lfcrank, "..neo_duct_lfc_ranked.csv")

### Plotting significant pathways
# these are the outputs from GSEA
positives <- read.table("..gsea_report_for_na_pos_1651694308838.tsv", sep = '\t', header = T)
negatives <- read.table("..gsea_report_for_na_neg_1651694308838.tsv", sep = '\t', header = T)

positives <- positives %>% mutate(group = "positives")
negatives <- negatives %>% mutate(group = "negatives")
positives <- positives[nrow(positives):1, ] 

sig_pos <- positives %>% filter(NOM.p.val < 0.05) %>% select("NAME")
sig_neg <- negatives %>% filter(NOM.p.val < 0.05) %>% select("NAME")

# Saving the significant pathways
write(sig_pos$NAME, "..Significant_Positives.txt")
write(sig_neg$NAME, "..Significant_Negatives.txt")

thedata <- rbind(negatives, positives)
thedata <- thedata %>% select("NAME", "NES", "group")
thedata$NAME <- factor(thedata$NAME, levels = thedata$NAME)

# Make lollipop graph

pdf("..GSEA_NEO_DUCT.pdf")
ggplot(thedata, aes(x=NAME, y=NES)) +
  geom_segment(aes(x=NAME, xend=NAME, y=0, yend=NES), size=0.5, alpha=0.5) +
  geom_point(size = 3, aes(fill= group), alpha=0.7, shape=21, stroke=1.4) +
  scale_fill_manual(values = c("blue", "red")) +
  theme_light() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
  ) +
  xlab("") +
  ylab("Normalized Enrichment Score") +
  coord_flip() +
  theme(axis.text.y = element_text(size = 6.8)) +
  scale_colour_manual("red", "blue")
dev.off()









