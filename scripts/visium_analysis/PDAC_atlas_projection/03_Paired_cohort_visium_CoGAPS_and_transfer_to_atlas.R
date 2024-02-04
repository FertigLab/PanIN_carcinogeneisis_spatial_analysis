title: "Paired_cohort_visium_CoGAPS_and_transfer_to_atlas"
author: "Alexander Bell"
date: "2/3/2024"

library(tidyverse)
library(CoGAPS)
library(Seurat)
library(projectR)
library(viridis)
library(ggplot2)
library(ggpubr)
library(fgsea)
library(msigdbr)
library(biomaRt)
library(forcats)

prepped <- readRDS(file = "../objects_and_metadata/prepped_full_object.rds")
sessionInfo()
# SUBSETTING EPITHELIAL SPOT POPULATION WITH CODA ####

# getting epithelial spots which are both >70% CODA and belong to the dominant Louvain cluster
epi_index <- which(prepped@meta.data$trueidents == "trueneoplasia" | prepped@meta.data$trueidents == "trueduct" )
epi_mat <- prepped@assays$SCT@data[,epi_index]
epi_dat <- as.matrix(as.data.frame(epi_mat))

# saving matrix of epithelial spots
saveRDS(epi_dat, "..epithelial_matrix.rds")


# PERFOMING COGAPS ON VISIUM (CODA-DEFINED) EPITHELIAL SPOTS ####
epi_dat <- readRDS("..epithelial_matrix.rds")

params <- new("CogapsParams")
params <- CogapsParams(
  sparseOptimization = TRUE,
  nPatterns = 2,
  seed = 367,
  geneNames = rownames(epi_dat),
  sampleNames = colnames(epi_dat),
  nIterations = 50000,
  distributed = NULL
)
# Update the sets for the object
#params <- setDistributedParams(params, nSets = 3, cut = 4, minNS = 1, maxNS = 5)
epi_cogaps <- CoGAPS(data = epi_dat, params = params)
saveRDS(epi_cogaps, "../epi_cogaps_n2")

params <- new("CogapsParams")
params <- CogapsParams(
  sparseOptimization = TRUE,
  nPatterns = 4,
  seed = 367,
  geneNames = rownames(epi_dat),
  sampleNames = colnames(epi_dat),
  nIterations = 50000,
  distributed = NULL
)
# Update the sets for the object
#params <- setDistributedParams(params, nSets = 3, cut = 4, minNS = 1, maxNS = 5)
epi_cogaps <- CoGAPS(data = epi_dat, params = params)
saveRDS(epi_cogaps, "../epi_cogaps_n4")

params <- new("CogapsParams")
params <- CogapsParams(
  sparseOptimization = TRUE,
  nPatterns = 6,
  seed = 367,
  geneNames = rownames(epi_dat),
  sampleNames = colnames(epi_dat),
  nIterations = 50000,
  distributed = NULL
)
# Update the sets for the object
#params <- setDistributedParams(params, nSets = 3, cut = 4, minNS = 1, maxNS = 5)
epi_cogaps <- CoGAPS(data = epi_dat, params = params)
saveRDS(epi_cogaps, "../epi_cogaps_n6")

params <- new("CogapsParams")
params <- CogapsParams(
  sparseOptimization = TRUE,
  nPatterns = 8,
  seed = 367,
  geneNames = rownames(epi_dat),
  sampleNames = colnames(epi_dat),
  nIterations = 50000,
  distributed = NULL
)
# Update the sets for the object
#params <- setDistributedParams(params, nSets = 3, cut = 4, minNS = 1, maxNS = 5)
epi_cogaps <- CoGAPS(data = epi_dat, params = params)
saveRDS(epi_cogaps, "../epi_cogaps_n8")

params <- new("CogapsParams")
params <- CogapsParams(
  sparseOptimization = TRUE,
  nPatterns = 10,
  seed = 367,
  geneNames = rownames(epi_dat),
  sampleNames = colnames(epi_dat),
  nIterations = 50000,
  distributed = NULL
)
# Update the sets for the object
#params <- setDistributedParams(params, nSets = 3, cut = 4, minNS = 1, maxNS = 5)
epi_cogaps <- CoGAPS(data = epi_dat, params = params)
saveRDS(epi_cogaps, "../epi_cogaps_n10")

# PREPPING PLOTS OF ATLAS AND ST-PATTERNS BY PANIN GRADE ####

# subsetting the seurat object for the projection of the atlas patterns onto coda-defined epithelial spots
alex_panin_epithelial_atlas <- subset(prepped, subset = (trueidents == "trueneoplasia" | trueidents == "trueduct"))

# isolating metadata because it contains LG, HG classification and module scores
panin_meta <- alex_panin_epithelial_atlas@meta.data
panin_meta$Grade <- 'other'
panin_meta$Grade[panin_meta$paningrade == 'high_grade'] <- 'HG'
panin_meta$Grade[panin_meta$paningrade == 'low_grade'] <- 'LG'
panin_meta$Grade[panin_meta$paningrade == 'normal'] <- 'N'
alex_panin_epithelial_atlas@meta.data <- panin_meta

# filtering for hg/lg/normal spots
alex_panin_epithelial_atlas <- subset(alex_panin_epithelial_atlas, Grade != 'other')
alex_panin_epithelial_atlas$Grade <- factor(alex_panin_epithelial_atlas$Grade,
                                      levels = c("N", "LG", "HG"))
# Palette for grade colors
# yellow, orange, red
Grade_pal <- c("#ebf4a5", "#f2a93b", "#ea3f25")

# importing CoGAPS results from PDAC scRNA-Seq atlas
atlas_cogaps <- readRDS(".../epiMat-e74511f6-8a1d-4929-ab42-1cb67c5a4fb6-result-8pattern.rds")

# projecitng atlas patterns onto CODA-defined epithelial cells 
projected_atlas_patterns <- projectR(data = alex_panin_epithelial_atlas@assays$SCT@scale.data, loadings = atlas_cogaps)

# add projected patterns to seurat object
atlas_pattern_weights <- as.data.frame(t(projected_atlas_patterns))

# adding atlas patterns to seurat metadata
for(pattern in colnames(atlas_pattern_weights)) {
  alex_panin_epithelial_atlas <- AddMetaData(object = alex_panin_epithelial_atlas,
                                       metadata = atlas_pattern_weights[[pattern]],
                                       col.name = pattern)
}

# making second seurat object to analyze the ST-CoGAPS data
alex_panin_epithelial_st <- subset(prepped, subset = (trueidents == "trueneoplasia" | trueidents == "trueduct"))

# importing ST CoGAPS restuls for nPatterns = 6
st_cogaps <- readRDS("..epi_cogaps_n6.rds")

# appending st_cogaps results to the seurat metadata
alex_panin_epithelial_st@meta.data <- cbind(alex_panin_epithelial_st@meta.data, st_cogaps@sampleFactors)

# list of Panin Images
panin_images <- Images(alex_panin_epithelial_atlas)

# isolating metadata because it contains LG, HG classification and module scores
panin_meta <- alex_panin_epithelial_st@meta.data
panin_meta$Grade <- 'other'
panin_meta$Grade[panin_meta$paningrade == 'high_grade'] <- 'HG'
panin_meta$Grade[panin_meta$paningrade == 'low_grade'] <- 'LG'
panin_meta$Grade[panin_meta$paningrade == 'normal'] <- 'N'
alex_panin_epithelial_st@meta.data <- panin_meta

# filter for hg/lg/normal stops
alex_panin_epithelial_st <- subset(alex_panin_epithelial_st, Grade != 'other')
alex_panin_epithelial_st$Grade <- factor(alex_panin_epithelial_st$Grade,
                                      levels = c("N", "LG", "HG"))

# Function for generating violin plots by panin grade
plotViolin <- function(seurat, feature, title, show = TRUE){
  violin <- VlnPlot(seurat, features = feature, group.by = "Grade",
                    cols = Grade_pal)
  # save minimum and maximum y axis limits for feature to ensure that space is left
  # for the ggpubr significance marks
  
  feat_range <- layer_scales(violin)$y$get_limits()
  feat_dif <- feat_range[2] - feat_range[1]
  y_min <- feat_range[1]
  y_max <- feat_range[2] + feat_dif/2.5
  
  plot <- violin +
    ylim(c(y_min, y_max)) +
    stat_compare_means(comparisons = list(
      c("N", "LG"),
      c("N", "HG"),
      c("LG", "HG"))) +
    ggtitle(title) + 
    ylab(paste0("projected ", feature, " weight")) +
    xlab("epithelial lesion grade")
  if(show){
    print(plot)
  }
  return(plot)
}

# function for generating Spatial Plots of CoGAPS patterns
plotPatternSpatial <- function(seurat, feature, image, show = TRUE){
  plot <- SpatialFeaturePlot(seurat, features = feature, images = image,
                             pt.size.factor = 2.5, crop = T) + scale_color_viridis(option = "magma")
  if(show){
    print(plot)
  }
  return(plot)
}

Vlnpath = "..outputs/violin/"

# violin plots of the weights of the ST-Patterns by panin grade (Supplemental 18A) ####
patterns <- colnames(st_cogaps@sampleFactors)
for(p in patterns){
  vln_plot <- plotViolin(alex_panin_epithelial_st, feature = p, title = NULL)
  ggsave(filename = paste0(Vlnpath, "/", p, "violin_ST_CoGAPS.pdf"),
         plot = vln_plot,
         width = unit(6, "in"),
         height = unit(6, "in"))
}

# violin plots of projected weights of the atlas CoGAPS patterns by panin grade (Main Figure 6) ####
patterns <- colnames(atlas_cogaps@sampleFactors)
for(p in patterns){
  vln_plot <- plotViolin(alex_panin_epithelial_atlas, feature = p, title = NULL)
  ggsave(filename = paste0(Vlnpath, "/", p, "violin_atlas_CoGAPS.pdf"),
         plot = vln_plot,
         width = unit(6, "in"),
         height = unit(6, "in"))
}

# spatial plots of projected atlas patterns onto ST (Main Figure 6) ####
PatternSpatialpath = "..outputs/spatial/"
patterns <- colnames(atlas_cogaps@sampleFactors)
for(i in panin_images){
  i_name <- gsub("\\.", "_", i)
  for(p in patterns){
    pat_plot <- plotPatternSpatial(alex_panin_epithelial_atlas, feature = p, image = i) +
      theme(legend.position="right") + scale_fill_viridis(option = "magma")
    
    pdf(paste0(PatternSpatialpath, "/", p, "_", i_name, "spatial_atlas_CoGAPS.pdf"))
    print(pat_plot)
    dev.off()
  }
}







# scatter plots of projected ST patterns onto atlas (Supplemental 18C) ####
library(projectR)
# reading in pre-processed PDAC atlas matrix prepared for CoGAPS
epimat <- Matrix::readMM('../epiMat.mtx')

# gene and sample names for the matrix
geneNames <- readRDS("../geneNames.rds")
sampleNames <- readRDS("../sampleNames.rds")

# reading in ST-CoGAPS results
st_cogaps <- readRDS("..epi_cogaps_n6.rds")

# projecting ST CoGAPS onto the PDAC atlas
project_embed <- projectR(data = epimat, loadings = st_cogaps)

# making scatterplots of the projected weights of each ST pattern on the PDAC atlas
for(x in 1:6){
  thepat <- paste0("Pattern_",x)
  zeplot <- FeaturePlot(thedata, features = thepat) + scale_color_viridis(option = "magma")
  pdf(paste0("..atlas_projection_pattern_",x,".pdf"))
  print(zeplot)
  dev.off()
}

# MSIGDBR OVERREPRESENTATION ANALYSIS OF ST-COGAPS PATTERN MARKERS (Supplemental 18B) ####
# List of MsigDB hallmarks and the genes in each set
hallmark_df <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- hallmark_df %>% split(x = .$gene_symbol, f = .$gs_name)

# obtain human hgnc gene symbols from ensembl
mart <- useEnsembl('genes', dataset = "hsapiens_gene_ensembl") #GRCh38
Hs_genes <- getBM("hgnc_symbol", mart = mart)

# prepar objects for overrepresentation analysis
theobject <- st_cogaps
thename <- "cogaps_n6"
patternMarkerResults <- patternMarkers(theobject, threshold = "cut", lp = NA, axis = 1)
names(patternMarkerResults$PatternMarkers) <- colnames(patternMarkerResults$PatternMarkerRanks)

# List of each gene as a pattern marker
PMlist <- patternMarkerResults$PatternMarkers

patternmarks <- data.frame(lapply(patternMarkerResults$PatternMarkers, "length<-",
                                  max(lengths(patternMarkerResults$PatternMarkers))))

#### looping through making bar plots of the overrepresented pattern markers within the ST-Patterns 
for(pattern in names(PMlist)) {
    
  # Over-representation analysis of pattern markers 
  suppressWarnings(
    result <- fora(pathways = hallmark_list,
                   genes = PMlist[[pattern]],
                   universe = Hs_genes$hgnc_symbol,
                   maxSize=2038)
  )
  
  # Append ratio of overlap/# of genes in hallmark (k/K)
  result$"k/K" <- result$overlap/result$size
  
  # Append log-transformed HB-adjusted q value (-10 * log(adjusted p value))
  result$"neg.log.q" <- (-10) * log10(result$padj)
  
  # Append shortened version of hallmark's name without "HALLMARK_"
  result$MsigDB_Hallmark <- substr(result$pathway, 10, nchar(result$pathway))
  
  # Reorder dataframe by ascending adjusted p value
  result <- mutate(result, MsigDB_Hallmark=fct_reorder(MsigDB_Hallmark, - padj))
  
  # rearrange columns to move overlap genes to the end and convert to vector for
  # compatibility with saving as csv
  result <- relocate(result, "overlapGenes", .after = "MsigDB_Hallmark")
  result <- relocate(result, "neg.log.q", .after = "padj")
  result <- result %>% mutate(overlapGenes = sapply(overlapGenes, toString))
  
  # save pattern marker overlaps as a csv file
  write.csv(result,
            file = paste0("../",thename,"/ORA_", pattern,"_MsigDB.csv"))
  
  # for plotting, limit the results to significant over-representation
  result <- result[result$padj < 0.05,]
  
  # cut off the pathways to the top 10 most significant
  result <- head(result, n = 10)
  thepaths <- as.vector(result$pathway)
  #plot and save the waterfall plot of ORA p-values
  plot <- ggplot(result, aes_string(y = "neg.log.q", x = "MsigDB_Hallmark", fill = "MsigDB_Hallmark")) +
    ## Specifies barplot
    geom_col() +
    ## Rename y axis
    ylab("-10*log10(FDR q-value)") + 
    ## Flips the coordinates
    coord_flip() +
    ## Makes the background white
    theme_minimal() +
    ## Add title
    ggtitle(paste0(pattern, ": Overrepresented MsigDB Hallmarks")) +
    ## This creates the dotted line at .05 value 
    geom_hline(yintercept=c(13.0103), linetype="dotted") + # Add veritcle line to show significances
    ## Adds the q values
    #geom_text(aes(label=format(signif(padj, 4))), hjust = -.04) +
    ## Removes legend
    theme(legend.position = "none", axis.text.y.left = element_text(color = "black", size=11, hjust=0.5)) +
    ## specifies limits 
    ylim(0, ceiling(max(result$"neg.log.q")) + (max(result$"neg.log.q")/4))
  ggsave(paste0("../",thename,"/ORA_", pattern,"_MsigDB.pdf"),
         plot = plot,
         width = unit(10, "in"),
         height = unit(6, "in"),
         device = "pdf")

}


