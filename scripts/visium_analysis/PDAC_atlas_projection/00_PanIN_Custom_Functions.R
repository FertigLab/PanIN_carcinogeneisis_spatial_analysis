title: "PanIN Custom Functions"
author: "Alexander Bell"
date: "2/3/2024"

library(Seurat)
library(miniUI)
library(grid)
library(ggplot2)
library(shiny)
library(cowplot)
library(patchwork)
library(hdf5r)
library(dplyr)
library(pheatmap)
library(data.table)
library(rqdatatable)
library(gdata)
library(ggrepel)
library(extrafont)

sessionInfo()
loadfonts()

# Pattern Weight UMAP Plots for the PDAC Atlas ####
plotUMAP <- function(cds, color_cells_by, title, show = TRUE) {
  plot <- plot_cells(cds, color_cells_by = color_cells_by,
                     label_cell_groups = FALSE, show_trajectory_graph = FALSE) +
    theme(legend.position = "bottom") +
    ggtitle(title)
  if(show){
    print(plot)
  }
  return(plot)
}

# Violin Plots: Projected Pattern Weight in Each Tumor Grade ####
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

# Spatial Plots: Projected Pattern Weight in Spots on Slides ####
plotPatternSpatial <- function(seurat, feature, image, show = TRUE){
  plot <- SpatialFeaturePlot(seurat, features = feature, images = image,
                             pt.size.factor = 3)
  if(show){
    print(plot)
  }
  return(plot)
}

# Spatial Plots: Tumor Grade in Spots on Slides ####
plotGradeSpatial <- function(seurat, group, image, show = TRUE){
  plot <- SpatialDimPlot(seurat, group.by = group, images = image,
                         pt.size.factor = 3, cols = Grade_pal)
  if(show){
    print(plot)
  }
  return(plot)
}

# FUNCTION: MAKE VOLCANO PLOTS ####
volcanoPlot <- function(input, title_1, title_2, batch, directory) {
  # Making uppers/mids/downers
  data <- input %>%
    mutate(neg.log.padj = -1*log10(p_val_adj)) %>% 
    filter(!is.na(neg.log.padj)) %>%
    mutate(topstat = neg.log.padj * abs(avg_log2FC * 2))
  #data$rn <- rownames(data)
  #data
  Uppers <- data %>% 
    filter((avg_log2FC >= 0.58) & (neg.log.padj >= 1.3)) %>% 
    mutate(Group = 'Uppers')
  Mids <- data %>% 
    filter(((avg_log2FC < 0.58) & (avg_log2FC >-0.58)) | (neg.log.padj < 1.3)) %>%
    mutate(Group = 'Mids')
  Lowers <- data %>% 
    filter((avg_log2FC <= -0.58) & (neg.log.padj >= 1.3)) %>%
    mutate(Group = 'Lowers')
  
  # Topstat for labeling volcano plot
  inUpperstopstat = character()
  inLowerstopstat = character()
  inMidstopstat = character ()
  
  
  #Uppers
  for(x in 1:nrow(Uppers)) {
    if(round(Uppers$topstat[x]) %in% round(sort(Uppers$topstat, decreasing = TRUE)[1:6])) {
      inUpperstopstat <- append(inUpperstopstat, TRUE) }
    else{
      inUpperstopstat <- append(inUpperstopstat, FALSE)
    }
  }
  
  #Lowers
  for(x in 1:nrow(Mids)){
    inMidstopstat <- append(inMidstopstat, FALSE)
  }
  for(x in 1:nrow(Lowers)) {
    if(round(Lowers$topstat[x]) %in% round(sort(Lowers$topstat, decreasing = TRUE)[1:6])) {
      inLowerstopstat <- append(inLowerstopstat, TRUE) }
    else{
      inLowerstopstat <- append(inLowerstopstat, FALSE)
    }
  }
  
  midupp <- append(inUpperstopstat, inMidstopstat)
  tops <- append(midupp, inLowerstopstat)
  Volcano_groups <- rbind(Uppers,Mids)
  Volcano_groups2 <- rbind(Volcano_groups,Lowers)
  Volcano_groups3 <- cbind(Volcano_groups2, tops)
  setDT(Volcano_groups3, keep.rownames = TRUE)[]
  Volcano_groups3
  #Volcano_groups3
  
  #Making DEGs list
  
  degs <- rbind(Uppers, Lowers)
  degs <- degs %>% mutate(statistic = (avg_log2FC/abs(avg_log2FC))*neg.log.padj)
  
  degs_ordered <- degs[order(degs$statistic, decreasing = TRUE),]
  setDT(degs_ordered, keep.rownames = TRUE)[]
  degs_ordered %>% dplyr::select(c('rn', 'statistic')) -> degsforgsea
  degs_ordered %>% dplyr::select(c('rn', 'avg_log2FC', "p_val", "p_val_adj", "neg.log.padj", "Group")) -> degslist
  colnames(degslist)[1] <- "Genes"
  
  # making gsea lists for all genes ranked either by LFC or "statistic
  allgenes <- rbind(Uppers, Mids, Lowers)
  allgenes <- allgenes %>% mutate(statistic = (avg_log2FC/abs(avg_log2FC))*neg.log.padj)
  setDT(allgenes, keep.rownames = TRUE)[]
  allgenes_ranked_by_statistic <- allgenes[order(allgenes$statistic, decreasing = TRUE),]
  allgenes_ranked_by_LFC <- allgenes[order(allgenes$avg_log2FC, decreasing = TRUE),]
  allgenes_statistic_gsea <- allgenes_ranked_by_statistic %>% dplyr::select(c('rn', 'statistic'))
  allgenes_LFC_gsea <- allgenes_ranked_by_LFC %>% dplyr::select(c('rn', 'statistic'))
  
  
  # all genes
  write.table(Volcano_groups3, file = paste0(directory, batch,"_", title_1, "vs", title_2, "_all_genes.csv"), row.names = TRUE, col.names = TRUE, sep = ',')
  
  # just degs, gsea formatted
  write.table(degsforgsea, file = paste0(directory, batch,"_", title_1, "vs", title_2, "_gseadegs.csv"), row.names = TRUE, col.names = TRUE, sep = ',')
  
  # just degs
  write.table(degslist, file = paste0(directory, batch, "_", title_1, "vs", title_2, "_DEGS.csv"), row.names = TRUE, col.names = TRUE, sep = ',')
  
  # all genes, gsea formatted, ranked by "statistic" (Aka the neg.log.fold change and direction)
  write.table(allgenes_statistic_gsea, file = paste0(directory, batch, "_", title_1, "vs", title_2, "_gsea_allgenes_rankedstatistic.csv"), row.names = TRUE, col.names = TRUE, sep = ',')
  
  # all genes, gsea formatted, ranked by LFC
  write.table(allgenes_LFC_gsea, file = paste0(directory, batch, "_", title_1, "vs", title_2, "_gsea_allgenes_rankedLFC.csv"), row.names = TRUE, col.names = TRUE, sep = ',')
  
  volcano_plot  <- ggplot(Volcano_groups3, mapping = aes(x = avg_log2FC, y=neg.log.padj, label = rn)) +
    geom_point(mapping = aes(color = Group), size = 1) +
    xlab(expression(paste(log[2], " fold change"))) +
    ylab(expression(paste(-log[10],"(q)"))) +  
    geom_text_repel(data = . %>% 
                      mutate(label = ifelse(Group %in% c("Uppers", "Lowers") & tops == TRUE,
                                            rn, "")),
                    aes(label = label), 
                    show.legend = FALSE,
                    max.overlaps = Inf,
                    force = 2,
                    box.padding = 0.2) +
    ggtitle(paste0(title_1," vs.", title_2)) +
    scale_color_discrete(name = "Groups", labels = c(paste0("Increased in ", title_2, " ", nrow(Lowers)), paste0("Unchanged: ", nrow(Mids)), paste0("Increased in ", title_1, " ", nrow(Uppers)))) +
    theme_bw()+
    theme(text = element_text(family = "Arial", color = "#2C2C2E"),
          plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
          legend.title = element_blank(),
          axis.text=element_text(size=10), 
          axis.title=element_text(size=10,face="bold"),
          legend.text = element_text(size=10),
          legend.position=c(1.02,0.5),
          legend.justification=c(0, 1), 
          legend.key.width=unit(1, "lines"), 
          legend.key.height=unit(1, "lines"), 
          plot.margin = unit(c(1, 13, 0.5, 0.5), "lines"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    geom_hline(yintercept = 1.3) +
    geom_vline(xintercept = c(-0.58,0.58))
  ggsave(file = paste0(directory,batch,"_",title_1,"vs",title_2, ".pdf"), device = cairo_pdf)
}



# DEPENDENCIES FOR SpatialDimChoose() ####
SpatialPlot <- function(
  object,
  group.by = NULL,
  features = NULL,
  images = NULL,
  cols = NULL,
  image.alpha = 1,
  crop = TRUE,
  slot = 'data',
  min.cutoff = NA,
  max.cutoff = NA,
  cells.highlight = NULL,
  cols.highlight = c('#DE2D26', 'grey50'),
  facet.highlight = FALSE,
  label = FALSE,
  label.size = 5,
  label.color = 'white',
  label.box = TRUE,
  repel = FALSE,
  ncol = NULL,
  combine = TRUE,
  pt.size.factor = 1.6,
  alpha = c(1, 1),
  stroke = 0.25,
  interactive = FALSE,
  do.identify = FALSE,
  identify.ident = NULL,
  do.hover = FALSE,
  information = NULL
) {
  if (isTRUE(x = do.hover) || isTRUE(x = do.identify)) {
    warning(
      "'do.hover' and 'do.identify' are deprecated as we are removing plotly-based interactive graphics, use 'interactive' instead for Shiny-based interactivity",
      call. = FALSE,
      immediate. = TRUE
    )
    interactive <- TRUE
  }
  if (!is.null(x = group.by) & !is.null(x = features)) {
    stop("Please specific either group.by or features, not both.")
  }
  images <- images %||% Images(object = object, assay = DefaultAssay(object = object))
  if (length(x = images) == 0) {
    images <- Images(object = object)
  }
  if (length(x = images) < 1) {
    stop("Could not find any spatial image information")
  }
  if (is.null(x = features)) {
    if (interactive) {
      return(ISpatialDimPlot(
        object = object,
        image = images[1],
        group.by = group.by,
        alpha = alpha
      ))
    }
    group.by <- group.by %||% 'ident'
    object[['ident']] <- Idents(object = object)
    data <- object[[group.by]]
    for (group in group.by) {
      if (!is.factor(x = data[, group])) {
        data[, group] <- factor(x = data[, group])
      }
    }
  } else {
    if (interactive) {
      return(ISpatialFeaturePlot(
        object = object,
        feature = features[1],
        image = images[1],
        slot = slot,
        alpha = alpha
      ))
    }
    data <- FetchData(
      object = object,
      vars = features,
      slot = slot
    )
    features <- colnames(x = data)
    # Determine cutoffs
    min.cutoff <- mapply(
      FUN = function(cutoff, feature) {
        return(ifelse(
          test = is.na(x = cutoff),
          yes = min(data[, feature]),
          no = cutoff
        ))
      },
      cutoff = min.cutoff,
      feature = features
    )
    max.cutoff <- mapply(
      FUN = function(cutoff, feature) {
        return(ifelse(
          test = is.na(x = cutoff),
          yes = max(data[, feature]),
          no = cutoff
        ))
      },
      cutoff = max.cutoff,
      feature = features
    )
    check.lengths <- unique(x = vapply(
      X = list(features, min.cutoff, max.cutoff),
      FUN = length,
      FUN.VALUE = numeric(length = 1)
    ))
    if (length(x = check.lengths) != 1) {
      stop("There must be the same number of minimum and maximum cuttoffs as there are features")
    }
    # Apply cutoffs
    data <- sapply(
      X = 1:ncol(x = data),
      FUN = function(index) {
        data.feature <- as.vector(x = data[, index])
        min.use <- SetQuantile(cutoff = min.cutoff[index], data.feature)
        max.use <- SetQuantile(cutoff = max.cutoff[index], data.feature)
        data.feature[data.feature < min.use] <- min.use
        data.feature[data.feature > max.use] <- max.use
        return(data.feature)
      }
    )
    colnames(x = data) <- features
    rownames(x = data) <- Cells(x = object)
  }
  features <- colnames(x = data)
  colnames(x = data) <- features
  rownames(x = data) <- colnames(x = object)
  facet.highlight <- facet.highlight && (!is.null(x = cells.highlight) && is.list(x = cells.highlight))
  if (do.hover) {
    if (length(x = images) > 1) {
      images <- images[1]
      warning(
        "'do.hover' requires only one image, using image ",
        images,
        call. = FALSE,
        immediate. = TRUE
      )
    }
    if (length(x = features) > 1) {
      features <- features[1]
      type <- ifelse(test = is.null(x = group.by), yes = 'feature', no = 'grouping')
      warning(
        "'do.hover' requires only one ",
        type,
        ", using ",
        features,
        call. = FALSE,
        immediate. = TRUE
      )
    }
    if (facet.highlight) {
      warning(
        "'do.hover' requires no faceting highlighted cells",
        call. = FALSE,
        immediate. = TRUE
      )
      facet.highlight <- FALSE
    }
  }
  if (facet.highlight) {
    if (length(x = images) > 1) {
      images <- images[1]
      warning(
        "Faceting the highlight only works with a single image, using image ",
        images,
        call. = FALSE,
        immediate. = TRUE
      )
    }
    ncols <- length(x = cells.highlight)
  } else {
    ncols <- length(x = images)
  }
  plots <- vector(
    mode = "list",
    length = length(x = features) * ncols
  )
  for (i in 1:ncols) {
    plot.idx <- i
    image.idx <- ifelse(test = facet.highlight, yes = 1, no = i)
    image.use <- object[[images[[image.idx]]]]
    coordinates <- GetTissueCoordinates(object = image.use)
    highlight.use <- if (facet.highlight) {
      cells.highlight[i]
    } else {
      cells.highlight
    }
    for (j in 1:length(x = features)) {
      cols.unset <- is.factor(x = data[, features[j]]) && is.null(x = cols)
      if (cols.unset) {
        cols <- hue_pal()(n = length(x = levels(x = data[, features[j]])))
        names(x = cols) <- levels(x = data[, features[j]])
      }
      plot <- SingleSpatialPlot(
        data = cbind(
          coordinates,
          data[rownames(x = coordinates), features[j], drop = FALSE]
        ),
        image = image.use,
        image.alpha = image.alpha,
        col.by = features[j],
        cols = cols,
        alpha.by = if (is.null(x = group.by)) {
          features[j]
        } else {
          NULL
        },
        pt.alpha = if (!is.null(x = group.by)) {
          alpha[j]
        } else {
          NULL
        },
        geom = if (inherits(x = image.use, what = "STARmap")) {
          'poly'
        } else {
          'spatial'
        },
        cells.highlight = highlight.use,
        cols.highlight = cols.highlight,
        pt.size.factor = pt.size.factor,
        stroke = stroke,
        crop = crop
      )
      if (is.null(x = group.by)) {
        plot <- plot +
          scale_fill_gradientn(
            name = features[j],
            colours = SpatialColors(n = 100)
          ) +
          theme(legend.position = 'top') +
          scale_alpha(range = alpha) +
          guides(alpha = FALSE)
      } else if (label) {
        plot <- LabelClusters(
          plot = plot,
          id = ifelse(
            test = is.null(x = cells.highlight),
            yes = features[j],
            no = 'highlight'
          ),
          geom = if (inherits(x = image.use, what = "STARmap")) {
            'GeomPolygon'
          } else {
            'GeomSpatial'
          },
          repel = repel,
          size = label.size,
          color = label.color,
          box = label.box,
          position = "nearest"
        )
      }
      if (j == 1 && length(x = images) > 1 && !facet.highlight) {
        plot <- plot +
          ggtitle(label = images[[image.idx]]) +
          theme(plot.title = element_text(hjust = 0.5))
      }
      if (facet.highlight) {
        plot <- plot +
          ggtitle(label = names(x = cells.highlight)[i]) +
          theme(plot.title = element_text(hjust = 0.5)) +
          NoLegend()
      }
      plots[[plot.idx]] <- plot
      plot.idx <- plot.idx + ncols
      if (cols.unset) {
        cols <- NULL
      }
    }
  }
  # if (do.identify) {
  #   return(CellSelector(
  #     plot = plot,
  #     object = identify.ident %iff% object,
  #     ident = identify.ident
  #   ))
  # } else if (do.hover) {
  #   return(HoverLocator(
  #     plot = plots[[1]],
  #     information = information %||% data[, features, drop = FALSE],
  #     axes = FALSE,
  #     # cols = c('size' = 'point.size.factor', 'colour' = 'fill'),
  #     images = GetImage(object = object, mode = 'plotly', image = images)
  #   ))
  # }
  if (length(x = images) > 1 && combine) {
    plots <- wrap_plots(plots = plots, ncol = length(x = images))
  } else if (length(x = images == 1) && combine) {
    plots <- wrap_plots(plots = plots, ncol = ncol)
  }
  return(plots)
}
DefaultImage <- function(object) {
  object <- UpdateSlots(object = object)
  images <- Images(object = object, assay = DefaultAssay(object = object))
  if (length(x = images) < 1) {
    images <- Images(object = object)
  }
  return(images[[1]])
}
UpdateSlots <- function(object) {
  object.list <- sapply(
    X = slotNames(x = object),
    FUN = function(x) {
      return(tryCatch(
        expr = slot(object = object, name = x),
        error = function(...) {
          return(NULL)
        }
      ))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  object.list <- Filter(f = Negate(f = is.null), x = object.list)
  object.list <- c('Class' = class(x = object)[1], object.list)
  object <- do.call(what = 'new', args = object.list)
  for (x in setdiff(x = slotNames(x = object), y = names(x = object.list))) {
    xobj <- slot(object = object, name = x)
    if (is.vector(x = xobj) && !is.list(x = xobj) && length(x = xobj) == 0) {
      slot(object = object, name = x) <- vector(mode = class(x = xobj), length = 1L)
    }
  }
  return(object)
}
SingleSpatialPlot <- function(
  data,
  image,
  cols = NULL,
  image.alpha = 1,
  pt.alpha = NULL,
  crop = TRUE,
  pt.size.factor = NULL,
  stroke = 0.25,
  col.by = NULL,
  alpha.by = NULL,
  cells.highlight = NULL,
  cols.highlight = c('#DE2D26', 'grey50'),
  geom = c('spatial', 'interactive', 'poly'),
  na.value = 'grey50'
) {
  geom <- match.arg(arg = geom)
  if (!is.null(x = col.by) && !col.by %in% colnames(x = data)) {
    warning("Cannot find '", col.by, "' in data, not coloring", call. = FALSE, immediate. = TRUE)
    col.by <- NULL
  }
  col.by <- col.by %iff% paste0("`", col.by, "`")
  alpha.by <- alpha.by %iff% paste0("`", alpha.by, "`")
  if (!is.null(x = cells.highlight)) {
    highlight.info <- SetHighlight(
      cells.highlight = cells.highlight,
      cells.all = rownames(x = data),
      sizes.highlight = pt.size.factor,
      cols.highlight = cols.highlight[1],
      col.base = cols.highlight[2]
    )
    order <- highlight.info$plot.order
    data$highlight <- highlight.info$highlight
    col.by <- 'highlight'
    levels(x = data$ident) <- c(order, setdiff(x = levels(x = data$ident), y = order))
    data <- data[order(data$ident), ]
  }
  plot <- ggplot(data = data, aes_string(
    x = colnames(x = data)[2],
    y = colnames(x = data)[1],
    fill = col.by,
    alpha = alpha.by
  ))
  plot <- switch(
    EXPR = geom,
    'spatial' = {
      if (is.null(x = pt.alpha)) {
        plot <- plot + geom_spatial(
          point.size.factor = pt.size.factor,
          data = data,
          image = image,
          image.alpha = image.alpha,
          crop = crop,
          stroke = stroke,
        )
      } else {
        plot <- plot + geom_spatial(
          point.size.factor = pt.size.factor,
          data = data,
          image = image,
          image.alpha = image.alpha,
          crop = crop,
          stroke = stroke,
          alpha = pt.alpha
        )
      }
      plot + coord_fixed() + theme(aspect.ratio = 1)
    },
    'interactive' = {
      plot + geom_spatial_interactive(
        data = tibble(grob = list(GetImage(object = image, mode = 'grob'))),
        mapping = aes_string(grob = 'grob'),
        x = 0.5,
        y = 0.5
      ) +
        geom_point(mapping = aes_string(color = col.by)) +
        xlim(0, ncol(x = image)) +
        ylim(nrow(x = image), 0) +
        coord_cartesian(expand = FALSE)
    },
    'poly' = {
      data$cell <- rownames(x = data)
      data[, c('x', 'y')] <- NULL
      data <- merge(
        x = data,
        y = GetTissueCoordinates(object = image, qhulls = TRUE),
        by = "cell"
      )
      plot + geom_polygon(
        data = data,
        mapping = aes_string(fill = col.by, group = 'cell')
      ) + coord_fixed() + theme_cowplot()
      
    },
    stop("Unknown geom, choose from 'spatial' or 'interactive'", call. = FALSE)
  )
  if (!is.null(x = cells.highlight)) {
    plot <- plot + scale_fill_manual(values = cols.highlight)
  }
  if (!is.null(x = cols) && is.null(x = cells.highlight)) {
    if (length(x = cols) == 1 && (is.numeric(x = cols) || cols %in% rownames(x = brewer.pal.info))) {
      scale <- scale_fill_brewer(palette = cols, na.value = na.value)
    } else if (length(x = cols) == 1 && (cols %in% c('alphabet', 'alphabet2', 'glasbey', 'polychrome', 'stepped'))) {
      colors <- DiscretePalette(length(unique(data[[col.by]])), palette = cols)
      scale <- scale_fill_manual(values = colors, na.value = na.value)
    } else {
      scale <- scale_fill_manual(values = cols, na.value = na.value)
    }
    plot <- plot + scale
  }
  plot <- plot + NoAxes() + theme(panel.background = element_blank())
  return(plot)
}
geom_spatial <-  function(
  mapping = NULL,
  data = NULL,
  image = image,
  image.alpha = image.alpha,
  crop = crop,
  stat = "identity",
  position = "identity",
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE,
  ...
) {
  layer(
    geom = GeomSpatial,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, image = image, image.alpha = image.alpha, crop = crop, ...)
  )
}
GeomSpatial <- ggproto(
  "GeomSpatial",
  Geom,
  required_aes = c("x", "y"),
  extra_params = c("na.rm", "image", "image.alpha", "crop"),
  default_aes = aes(
    shape = 21,
    colour = "black",
    point.size.factor = 1.0,
    fill = NA,
    alpha = NA,
    stroke = 0.25
  ),
  setup_data = function(self, data, params) {
    data <- ggproto_parent(Geom, self)$setup_data(data, params)
    # We need to flip the image as the Y coordinates are reversed
    data$y = max(data$y) - data$y + min(data$y)
    data
  },
  draw_key = draw_key_point,
  draw_panel = function(data, panel_scales, coord, image, image.alpha, crop) {
    # This should be in native units, where
    # Locations and sizes are relative to the x- and yscales for the current viewport.
    if (!crop) {
      y.transform <- c(0, nrow(x = image)) - panel_scales$y.range
      data$y <- data$y + sum(y.transform)
      panel_scales$x$continuous_range <- c(0, ncol(x = image))
      panel_scales$y$continuous_range <- c(0, nrow(x = image))
      panel_scales$y.range <- c(0, nrow(x = image))
      panel_scales$x.range <- c(0, ncol(x = image))
    }
    z <- coord$transform(
      data.frame(x = c(0, ncol(x = image)), y = c(0, nrow(x = image))),
      panel_scales
    )
    # Flip Y axis for image
    z$y <- -rev(z$y) + 1
    wdth <- z$x[2] - z$x[1]
    hgth <- z$y[2] - z$y[1]
    vp <- viewport(
      x = unit(x = z$x[1], units = "npc"),
      y = unit(x = z$y[1], units = "npc"),
      width = unit(x = wdth, units = "npc"),
      height = unit(x = hgth, units = "npc"),
      just = c("left", "bottom")
    )
    img.grob <- GetImage(object = image)
    
    img <- editGrob(grob = img.grob, vp = vp)
    # spot.size <- slot(object = image, name = "spot.radius")
    spot.size <- Radius(object = image)
    coords <- coord$transform(data, panel_scales)
    pts <- pointsGrob(
      x = coords$x,
      y = coords$y,
      pch = data$shape,
      size = unit(spot.size, "npc") * data$point.size.factor,
      gp = gpar(
        col = alpha(colour = coords$colour, alpha = coords$alpha),
        fill = alpha(colour = coords$fill, alpha = coords$alpha),
        lwd = coords$stroke)
    )
    vp <- viewport()
    gt <- gTree(vp = vp)
    if (image.alpha > 0) {
      if (image.alpha != 1) {
        img$raster = as.raster(
          x = matrix(
            data = alpha(colour = img$raster, alpha = image.alpha),
            nrow = nrow(x = img$raster),
            ncol = ncol(x = img$raster),
            byrow = TRUE)
        )
      }
      gt <- addGrob(gTree = gt, child = img)
    }
    gt <- addGrob(gTree = gt, child = pts)
    # Replacement for ggname
    gt$name <- grobName(grob = gt, prefix = 'geom_spatial')
    return(gt)
    # ggplot2:::ggname("geom_spatial", gt)
  }
)
InvertCoordinate <- function(x, MARGIN = 2) {
  if (!is.null(x = x)) {
    switch(
      EXPR = MARGIN,
      '1' = {
        rmin <- 'left'
        rmax <- 'right'
        cmin <- 'xmin'
        cmax <- 'xmax'
      },
      '2' = {
        rmin <- 'bottom'
        rmax <- 'top'
        cmin <- 'ymin'
        cmax <- 'ymax'
      },
      stop("'MARGIN' must be either 1 or 2", call. = FALSE)
    )
    # Fix the range so that rmin becomes rmax and vice versa
    # Needed for both points and brushes
    range <- x$range
    x$range[[rmin]] <- range[[rmax]]
    x$range[[rmax]] <- range[[rmin]]
    # Fix the cmin and cmax values, if provided
    # These are used for brush boundaries
    coords <- c(x[[cmin]], x[[cmax]])
    if (all(!is.null(x = coords))) {
      names(x = coords) <- c(cmin, cmax)
      x[[cmin]] <- quantile(
        x = x$range[[rmin]]:x$range[[rmax]],
        probs = 1 - (coords[cmax] / x$range[[rmax]]),
        names = FALSE
      )
      x[[cmax]] <- quantile(
        x = x$range[[rmin]]:x$range[[rmax]],
        probs = 1 - (coords[cmin] / x$range[[rmax]]),
        names = FALSE
      )
    }
  }
  return(x)
}
DefaultDimReduc <- function(object, assay = NULL) {
  object <- UpdateSlots(object = object)
  assay <- assay %||% DefaultAssay(object = object)
  drs.use <- c('umap', 'tsne', 'pca')
  dim.reducs <- FilterObjects(object = object, classes.keep = 'DimReduc')
  drs.assay <- Filter(
    f = function(x) {
      return(DefaultAssay(object = object[[x]]) == assay)
    },
    x = dim.reducs
  )
  if (length(x = drs.assay) > 0) {
    index <- lapply(
      X = drs.use,
      FUN = grep,
      x = drs.assay,
      ignore.case = TRUE
    )
    index <- Filter(f = length, x = index)
    if (length(x = index) > 0) {
      return(drs.assay[min(index[[1]])])
    }
  }
  index <- lapply(
    X = drs.use,
    FUN = grep,
    x = dim.reducs,
    ignore.case = TRUE
  )
  index <- Filter(f = length, x = index)
  if (length(x = index) < 1) {
    stop(
      "Unable to find a DimReduc matching one of '",
      paste(drs.use[1:(length(x = drs.use) - 1)], collapse = "', '"),
      "', or '",
      drs.use[length(x = drs.use)],
      "', please specify a dimensional reduction to use",
      call. = FALSE
    )
  }
  return(dim.reducs[min(index[[1]])])
}
FilterObjects <- function(object, classes.keep = c('Assay', 'DimReduc')) {
  object <- UpdateSlots(object = object)
  slots <- na.omit(object = Filter(
    f = function(x) {
      sobj <- slot(object = object, name = x)
      return(is.list(x = sobj) && !is.data.frame(x = sobj) && !is.package_version(x = sobj))
    },
    x = slotNames(x = object)
  ))
  slots <- grep(pattern = 'tools', x = slots, value = TRUE, invert = TRUE)
  slots <- grep(pattern = 'misc', x = slots, value = TRUE, invert = TRUE)
  slots.objects <- unlist(
    x = lapply(
      X = slots,
      FUN = function(x) {
        return(names(x = slot(object = object, name = x)))
      }
    ),
    use.names = FALSE
  )
  object.classes <- sapply(
    X = slots.objects,
    FUN = function(i) {
      return(inherits(x = object[[i]], what = classes.keep))
    }
  )
  object.classes <- which(x = object.classes, useNames = TRUE)
  return(names(x = object.classes))
}
SingleDimPlot <- function(
  data,
  dims,
  col.by = NULL,
  cols = NULL,
  pt.size = NULL,
  shape.by = NULL,
  alpha.by = NULL,
  order = NULL,
  label = FALSE,
  repel = FALSE,
  label.size = 4,
  cells.highlight = NULL,
  cols.highlight = '#DE2D26',
  sizes.highlight = 1,
  na.value = 'grey50',
  raster = NULL
) {
  pt.size <- pt.size %||% AutoPointSize(data = data, raster = raster)
  if ((nrow(x = data) > 1e5) & !isFALSE(raster)){
    message("Rasterizing points since number of points exceeds 100,000.",
            "\nTo disable this behavior set `raster=FALSE`")
  }
  raster <- raster %||% (nrow(x = data) > 1e5)
  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }
  if (!is.data.frame(x = data)) {
    data <- as.data.frame(x = data)
  }
  if (is.character(x = dims) && !all(dims %in% colnames(x = data))) {
    stop("Cannot find dimensions to plot in data")
  } else if (is.numeric(x = dims)) {
    dims <- colnames(x = data)[dims]
  }
  if (!is.null(x = cells.highlight)) {
    highlight.info <- SetHighlight(
      cells.highlight = cells.highlight,
      cells.all = rownames(x = data),
      sizes.highlight = sizes.highlight %||% pt.size,
      cols.highlight = cols.highlight,
      col.base = cols[1] %||% '#C3C3C3',
      pt.size = pt.size
    )
    order <- highlight.info$plot.order
    data$highlight <- highlight.info$highlight
    col.by <- 'highlight'
    pt.size <- highlight.info$size
    cols <- highlight.info$color
  }
  if (!is.null(x = order) && !is.null(x = col.by)) {
    if (typeof(x = order) == "logical") {
      if (order) {
        data <- data[order(!is.na(x = data[, col.by]), data[, col.by]), ]
      }
    } else {
      order <- rev(x = c(
        order,
        setdiff(x = unique(x = data[, col.by]), y = order)
      ))
      data[, col.by] <- factor(x = data[, col.by], levels = order)
      new.order <- order(x = data[, col.by])
      data <- data[new.order, ]
      if (length(x = pt.size) == length(x = new.order)) {
        pt.size <- pt.size[new.order]
      }
    }
  }
  if (!is.null(x = col.by) && !col.by %in% colnames(x = data)) {
    warning("Cannot find ", col.by, " in plotting data, not coloring plot")
    col.by <- NULL
  } else {
    # col.index <- grep(pattern = col.by, x = colnames(x = data), fixed = TRUE)
    col.index <- match(x = col.by, table = colnames(x = data))
    if (grepl(pattern = '^\\d', x = col.by)) {
      # Do something for numbers
      col.by <- paste0('x', col.by)
    } else if (grepl(pattern = '-', x = col.by)) {
      # Do something for dashes
      col.by <- gsub(pattern = '-', replacement = '.', x = col.by)
    }
    colnames(x = data)[col.index] <- col.by
  }
  if (!is.null(x = shape.by) && !shape.by %in% colnames(x = data)) {
    warning("Cannot find ", shape.by, " in plotting data, not shaping plot")
  }
  if (!is.null(x = alpha.by) && !alpha.by %in% colnames(x = data)) {
    warning(
      "Cannot find alpha variable ",
      alpha.by,
      " in data, setting to NULL",
      call. = FALSE,
      immediate. = TRUE
    )
    alpha.by <- NULL
  }
  
  plot <- ggplot(data = data)
  plot <- if (isTRUE(x = raster)) {
    plot + geom_scattermore(
      mapping = aes_string(
        x = dims[1],
        y = dims[2],
        color = paste0("`", col.by, "`"),
        shape = shape.by,
        alpha = alpha.by
      ),
      pointsize = pt.size
    )
  } else {
    plot + geom_point(
      mapping = aes_string(
        x = dims[1],
        y = dims[2],
        color = paste0("`", col.by, "`"),
        shape = shape.by,
        alpha = alpha.by
      ),
      size = pt.size
    )
  }
  plot <- plot +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    labs(color = NULL, title = col.by) +
    CenterTitle()
  if (label && !is.null(x = col.by)) {
    plot <- LabelClusters(
      plot = plot,
      id = col.by,
      repel = repel,
      size = label.size
    )
  }
  if (!is.null(x = cols)) {
    if (length(x = cols) == 1 && (is.numeric(x = cols) || cols %in% rownames(x = brewer.pal.info))) {
      scale <- scale_color_brewer(palette = cols, na.value = na.value)
    } else if (length(x = cols) == 1 && (cols %in% c('alphabet', 'alphabet2', 'glasbey', 'polychrome', 'stepped'))) {
      colors <- DiscretePalette(length(unique(data[[col.by]])), palette = cols)
      scale <- scale_color_manual(values = colors, na.value = na.value)
    } else {
      scale <- scale_color_manual(values = cols, na.value = na.value)
    }
    plot <- plot + scale
  }
  plot <- plot + theme_cowplot()
  return(plot)
}

AutoPointSize <- function(data, raster = NULL) {
  return(ifelse(
    test = isTRUE(x = raster),
    yes = 1,
    no = min(1583 / nrow(x = data), 1)
  ))
}

# FUNCTION: SpatialDimChoose() ####
# object is data, 
# image is the slide (i.e. "slice1.1), 
# group.by is the ident that is displayed on the dimplot
# alpha is transparency of spots (I think, tbh don't remember right now), just leave that as default
# new_identity is the new identity that the subset will go into. can be an identity that already exists (i.e. cell_class)
# group_name is the name of the subset for that round of clicking spots (i.e. acinar_cells)

SpatialDimChoose <- function(
  object,
  image = NULL,
  group.by = NULL,
  alpha = c(0.3, 1),
  new_identity = NULL,
  group_name = NULL
) {
  # Setup gadget UI
  ui <- miniPage(
    miniButtonBlock(miniTitleBarButton(
      inputId = 'done',
      label = 'Done',
      primary = TRUE
    )),
    miniContentPanel(
      fillRow(
        plotOutput(
          outputId = 'plot',
          height = '100%',
          click = clickOpts(id = 'click', clip = TRUE),
          hover = hoverOpts(id = 'hover', delay = 10, nullOutside = TRUE)
        ),
        height = '97%'
      ),
      verbatimTextOutput(outputId = 'info')
    )
  )
  # Get plotting data
  # Prepare plotting data
  image <- image %||% DefaultImage(object = object)
  cells.use <- Cells(x = object[[image]])
  group.by <- group.by %||% 'ident'
  group.data <- FetchData(
    object = object,
    vars = group.by,
    cells = cells.use
  )
  coords <- GetTissueCoordinates(object = object[[image]])
  plot.data <- cbind(coords, group.data)
  plot.data$selected_ <- FALSE
  Idents(object = object) <- group.by
  
  # Set up the server
  server <- function(input, output) {
    click <- reactiveValues(pt = NULL)
    plot.env <- reactiveValues(data = plot.data, alpha.by = NULL)
    
    
    clicklist <- reactiveVal(list())
    observeEvent(input$click, {
      click <- rownames(x = nearPoints(
        df = plot.data,
        coordinfo = InvertCoordinate(x = input$click),
        threshold = 10,
        maxpoints = 1))
      temp <- clicklist()
      temp[[length(temp) + 1]] <- click
      clicklist(temp)
    })
    observeEvent(
      eventExpr = input$done,
      handlerExpr = stopApp(returnValue = plot.env$plot)
    )
    # Set plot
    output$plot <- renderPlot(
      expr = {
        plot.env$plot <- SingleSpatialPlot(
          data = plot.env$data,
          image = object[[image]],
          col.by = group.by,
          crop = TRUE,
          alpha.by = plot.env$alpha.by,
          pt.size.factor = 1.6
        ) + scale_alpha_ordinal(range = alpha) + NoLegend()
        plot.env$plot
      }
    )
    output$info <- renderPrint(
      expr = {
        cell.hover <- rownames(x = nearPoints(
          df = plot.data,
          coordinfo = InvertCoordinate(x = input$click),
          threshold = 10,
          maxpoints = 1
        ))
        if(length(x = cell.hover) == 1) {
          my_list <<- c(my_list, cell.hover)
          paste(cell.hover)#, paste('Group:', plot.data[cell.hover, group.by, drop = TRUE]), collapse = '<br />')
        } else {
          NULL
        }
      }
    )
  }
  # SET LIST TO EMPTY
  my_list <- list()
  # RUN THE THANG
  runGadget(app = ui, server = server)
  # add group to identity class
  object <- SetIdent(object, value = new_identity)
  my_list <- unique(my_list)
  my_list <- unlist(my_list)
  index1 <- which(rownames(object[[]]) %in% my_list)
  Idents(object, cells = index1) <- group_name
  object@meta.data[new_identity] <- Idents(object)
  print(my_list)
  return(object)
}