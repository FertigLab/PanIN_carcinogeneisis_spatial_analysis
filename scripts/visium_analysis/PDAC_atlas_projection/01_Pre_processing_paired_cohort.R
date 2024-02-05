title: "Pre_processing_paired_cohort"
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

source("..Custom_Functions.R")

sessionInfo()

# PREPROCESSING 1 - LOAD DATA ####
file_list <- c("1152735",
               "1153847",
               "1153848",
               "Batch_2/JHU_P2_114_D1",
               "Batch_2/JHU_P2_113_C1",
               "Batch_2/JHU_P2_114_A1",
               "Batch_2/JHU_P2_113_D1")

# importing all the .h5 files
n=0
namvec <- c()
for(file in file_list) {
  n = n+1
  data_dir <- paste0('..10X_Files/',file)
  nam = paste0("data", n)
  namvec <- append(namvec, nam)
  tmp <- Load10X_Spatial(data_dir, filename = "filtered_feature_bc_matrix.h5")
  assign(nam, tmp)
}

# PREPROCESSING 2 - PRE-PROCESS, MERGE, AND ANNOTATE SEURAT OBJECTS ####

# this loop SCTransforms and then merges all the Seurat objects
for(x in 1:length(namvec)) {
  data <- get(namvec[x])
  data <- SCTransform(data, assay = "Spatial", verbose = FALSE, return.only.var.genes = FALSE)
  if(x == 1){
    vf = VariableFeatures(data)
    nam = paste0("data",x)
    assign(nam, data)
    next
  }
  if(x == 2) {
    vf = append(vf, VariableFeatures(data))
    data_merged <- merge(data, data1)
    nam = paste0("data",x)
    assign(nam, "1")
    data1 <- "1"
    next
  }
  if(x > 2){
    vf = append(vf, VariableFeatures(data))
    data_merged <- merge(data_merged, data)
    nam = paste0("data",x)
    assign(nam, "1")
  }
}

# making slide ID  metadata column
data_merged@meta.data <- cbind(data_merged@meta.data, substr(rownames(data_merged@meta.data), 18, 25))
colnames(data_merged@meta.data)[length(colnames(data_merged@meta.data))] <- "Slide_ID"

# Making "Slide_Name" metadata  column
testlist <- data_merged@meta.data$Slide_ID
slides <- unique(data_merged@meta.data$Slide_ID) 
for(x in 1:length(slides)) {
  tmp <- grep(paste0("^",slides[x], "$"), testlist)
  testlist <- replace(testlist, tmp, paste0("Slide", " ", x))
  print(tmp)
}

data_merged@meta.data <- cbind(data_merged@meta.data, testlist)
colnames(data_merged@meta.data)[length(colnames(data_merged@meta.data))] <- "Slide_Name"

# making SCT the default assay
DefaultAssay(data_merged) <- "SCT"


# PREPROCESSING 3 - PCA, CLUSTERS, UMAP ####

# Variable features from merged data
VariableFeatures(data_merged) <- vf

# Run PCA
data_merged <- RunPCA(data_merged, verbose = FALSE)

# Elbow plot
ElbowPlot(data_merged)

# Run UMAP
data_merged <- RunUMAP(data_merged, dims = 1:20)

# Find Neighbors
data_merged <- FindNeighbors(data_merged, dims = 1:20)

# Louvain Cluster
data_merged <- FindClusters(data_merged, verbose = FALSE, resolution = 0.5)

# Renaming the metadata column for the clusters
colnames(data_merged@meta.data)[which(colnames(data_merged@meta.data) == "SCT_snn_res.0.5")] <- "pre_harmony_clusters"

#Saving UMAP prior to batch adjustment with harmony, grouped by Louvain Clusters
pdf("..UMAP_Clusters_Prior_To_Harmony_Correction.pdf")
DimPlot(data_merged, reduction = "umap", group.by = c("pre_harmony_clusters"))
dev.off()



# Saving UMAP prior to batch adjustment with harmony, grouped by slides
pdf("..PreProcessing//UMAP_Slide_Names_Prior_To_Harmony_Correction.pdf")
DimPlot(data_merged, reduction = "umap", group.by = c("Slide_Name"))
dev.off()




# PREPROCESSING 4 - HARMONY REDUCTION ####

## adjusting with harmony--OPTIONAL--can use original harmony reduction by importing (see blow)
#data_merged <- RunHarmony(data_merged, group.by.vars = "Slide_ID", assay.use = "SCT", plot_convergence = T)
#harmony_embeddings <- Embeddings(data_merged, 'harmony')

# using original harmony reduction
data_merged@reductions$harmony <- readRDS("..objects_and_metadata/harmony_reduction.rds")

# Plotting Harmony reduction
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = data_merged, reduction = "harmony", pt.size = .1, group.by = "Slide_Name")
p2 <- VlnPlot(object = data_merged, features = "harmony_1", group.by = "Slide_Name", pt.size = .1)
pdf("..PreProcessing/Plotting_Harmony_Reduction.pdf", height = 5, width = 10)
plot_grid(p1,p2)
dev.off()


# Performing UMAP and clustering on Harmony reduction
data_merged <- data_merged %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

# using the original clusters because they are going to be different every time
data_merged$SCT_snn_res.0.5 <- readRDS("..objects_and_metadata/harmony_louvain.rds")

# renaming the metadata column for the clusters
colnames(data_merged@meta.data)[which(colnames(data_merged@meta.data) == "SCT_snn_res.0.5")] <- "post_harmony_clusters"

# Plotting UMAP with ident of harmony clusters 
pdf("..PreProcessing/UMAP_of_Louvain_Clusters_After_Harmony.pdf")
DimPlot(data_merged,  group.by = "post_harmony_clusters")
dev.off()

# Plotting UMAP with ident of slide name
pdf("..PreProcessing/UMAP_Slide_Name_After_Harmony.pdf")
DimPlot(data_merged,  group.by = "Slide_Name")
dev.off()

# Save Seurat Object
#saveRDS(data_merged, "..objects_and_metadata/save_1.rds")
# optional, start script here by reading in this object
data_merged <- readRDS("..objects_and_metadata/save_1.rds")

# Save meta data
#saveRDS(data_merged@meta.data, "..objects_and_metadata/metadata_1.rds")

# PREPROCESSING 5 - MERGING VISIUM WITH CODA ####
data_merged1 <- data_merged

# Making the "Slide_Number" identity, which will be useful for importing the CODA spot composition annotations
new_vec <- data_merged@meta.data$Slide_Name
new_vec <- replace(new_vec, which(new_vec == "Slide 1"), "slide_113_1")
new_vec <- replace(new_vec, which(new_vec == "Slide 2"), "slide_114_2")
new_vec <- replace(new_vec, which(new_vec == "Slide 3"), "slide_113_2")
new_vec <- replace(new_vec, which(new_vec == "Slide 4"), "slide_114_4")
new_vec <- replace(new_vec, which(new_vec == "Slide 5"), "slide_113_3")
new_vec <- replace(new_vec, which(new_vec == "Slide 6"), "slide_114_1")
new_vec <- replace(new_vec, which(new_vec == "Slide 7"), "slide_113_4")
data_merged@meta.data$Slide_Number <- new_vec

# Importing the CODA annotations
ashleyslidelist <- c("113_1", "113_2", "113_3", "113_4", "114_1", "114_2", "114_4")
tablelist <- list()
for(x in 1:length(ashleyslidelist)) {
  tmp <- readxl::read_xlsx(paste0("..April12_Ashley_CODA_With_Fat/", ashleyslidelist[x], "_tissue_positions_list_cellular_compositions_fat.xlsx"))
  tmp <- as.data.frame(tmp)
  assign(paste0("slide_",ashleyslidelist[x],"_CODA_compositions"), tmp)
  tablelist[[x]] <- tmp
  colnames(tablelist[[x]])[1:6] <- c("barcode", "slide", "spot_row", "spot_col", "pixel_row", "pixel_col") #check later to make sure this is right
}

names(tablelist) <- c("slide_113_1", "slide_113_2", "slide_113_3", "slide_113_4", "slide_114_1", "slide_114_2", "slide_114_4")

# ADDING CODA CELLULAR COMPOSITION TO METADATA
for(x in 1:7){
  data_merged@meta.data %>% filter(Slide_Number == names(tablelist[x])) -> tmpdata
  tmpdata %>%
    rownames() %>%
    substring(first = 0, last = 18) -> tmp1
  tmpdata$barcode <- tmp1
  as.data.frame(tablelist[[x]]) -> thetable
  leftjoined <- left_join(tmpdata, thetable, by = "barcode")
  rownames(leftjoined) <- rownames(tmpdata)
  leftjoined <- leftjoined %>%
    select("orig.ident", "nCount_Spatial", "nFeature_Spatial", "nCount_SCT", "nFeature_SCT", "pre_harmony_clusters", "seurat_clusters", "Slide_ID", "post_harmony_clusters", "Slide_Name", "Slide_Number", "slide", "spot_row", "spot_col", "pixel_row", "pixel_col", "islets", "normal epithelium", "smooth muscle", "fat", "acini", "collagen", "nontissue", "panin")
  colnames(leftjoined) <- c("orig.ident", "nCount_Spatial", "nFeature_Spatial", "nCount_SCT", "nFeature_SCT", "pre_harmony_clusters", "seurat_clusters", "Slide_ID", "post_harmony_clusters", "Slide_Name", "Slide_Number", "slide", "spot_row", "spot_col", "pixel_row", "pixel_col", "islets", "normal epithelium", "smooth muscle", "fat", "acini", "collagen", "nontissue", "panin")
  assign(paste0(names(tablelist[x]),"_meta_data"), leftjoined)
}


# Making new metadata sheet from the individual slides
new_metadata <- rbind(slide_113_1_meta_data, slide_114_2_meta_data, slide_113_2_meta_data, slide_114_4_meta_data, slide_113_3_meta_data, slide_114_1_meta_data, slide_113_4_meta_data)
data_merged@meta.data <- new_metadata

# Generating cutoffs for CODA spot composition
panin_index <- which(data_merged@meta.data$panin > 70)
acini_index <- which(data_merged@meta.data$acini > 70)
islets_index <- which(data_merged@meta.data$islets > 70)
epithelium_index <- which(data_merged@meta.data$`normal epithelium` > 70)
muscle_index <- which(data_merged@meta.data$`smooth muscle` > 70) 
collagen_index <- which(data_merged@meta.data$collagen > 70) 
fat_index <- which(data_merged@meta.data$fat > 70)

# Labelling CODA defined cell types
data_merged@meta.data$coda70 <- "nocells"
data_merged@meta.data$coda70[panin_index] <- "panin"
data_merged@meta.data$coda70[acini_index] <- "acini"
data_merged@meta.data$coda70[islets_index] <- "islets"
data_merged@meta.data$coda70[epithelium_index] <- "normal epithelium"
data_merged@meta.data$coda70[muscle_index] <- "smooth muscle"
data_merged@meta.data$coda70[collagen_index] <- "collagen"
data_merged@meta.data$coda70[fat_index] <- "fat"

# TAKING THE OVERLAP OF CODA AND LOUVAIN CLUSTERS. "TRUEIDENTS" identity contains the groups that will be used for comparing neoplasia to normal duct
# annotating the gene subsets that overlap the CODA cell types and the Seurat cell types
data_merged@meta.data$trueidents <- "nocells"
data_merged@meta.data %>% filter(post_harmony_clusters == 6 & coda70 == "islets") %>% rownames() -> truecells
data_merged$trueidents[rownames(data_merged@meta.data) %in% truecells] <- "trueislets"
data_merged@meta.data %>% filter(post_harmony_clusters == 8 & coda70 == "normal epithelium") %>% rownames() -> truecells
data_merged$trueidents[rownames(data_merged@meta.data) %in% truecells] <- "trueduct"
data_merged@meta.data %>% filter(post_harmony_clusters == 11 & coda70 == "panin") %>% rownames() -> truecells
data_merged$trueidents[rownames(data_merged@meta.data) %in% truecells] <- "trueneoplasia"

# PREPROCESSING 6 - ANNOTATING HIGH AND LOW GRADE NEOPLASIA ####
# LG/HG annotation with SpatialDimChoose (import previously annotated data)
# FUNCTION: SpatialDimChoose()
# object is data
# image is the slide (i.e. "slice1.1), 
# group.by is the ident that is displayed on the dimplot
# new_identity is the new identity that the subset will go into. can be an identity that already exists
# group_name is the name of the subset for that round of clicking spots (i.e. low_grade

# data_merged <- SpatialDimChoose(data_merged, group.by = "post_harmony_clusters", new_identity = "paningrade", group_name = "low_grade", image = "slice1")
data_merged$paningrade <- readRDS("..objects_and_metadata/paningrade.rds")

# slide specific LG/HG with SpatialDimChoose (import previously annotated data)
# data_merged <- SpatialDimChoose(data_merged, group.by = "post_harmony_clusters", new_identity = "paningradeslide", group_name = "LG_113_1", image = "slice1")
data_merged$paningradeslide <- readRDS("..objects_and_metadata/paningradeslide.rds")

# saving non-SCT prepped seurat object
saveRDS(data_merged, "..objects_and_metadata/full_unprepped_object.rds")

# PREPROCESSING 7 - SCT PREPPING DATA TO BE USED WITH CROSS-SLIDE COMPARISONS ####
# Generate new seurat object that has been SCT Prepped
prepped <- PrepSCTFindMarkers(data_merged, assay = "SCT")
prepped$paningrade <- data_merged$paningrade
prepped$paningradeslide <- data_merged$paningradeslide

# saving full SCT prepped object
saveRDS(prepped, file = " ..objects_and_metadata/prepped_full_object.rds")

