library(ComplexHeatmap)

customTheme = function(sizeStripFont, xAngle, hjust, vjust, xSize,
                       ySize, xAxisSize, yAxisSize) {
  theme(strip.background = element_rect(colour = "black", fill = "white",
                                        size = 1), strip.text.x = element_text(size = sizeStripFont),
        strip.text.y = element_text(size = sizeStripFont), axis.text.x = element_text(angle = xAngle,
                                                                                      hjust = hjust, vjust = vjust, size = xSize, color = "black"),
        axis.text.y = element_text(size = ySize, color = "black"),
        axis.title.x = element_text(size = xAxisSize, color = "black"),
        axis.title.y = element_text(size = yAxisSize, color = "black"),
        panel.background = element_rect(fill = "white", color = "black"))
}

# Differential expression analysis
de_analysis = function(segments, pdata, eset, gp, lvls, contrts, tissue){
  top <- lapply(segments, function(i){
    phenodata_fullroi <- subset(pdata, segment == i)
    q3norm_fullroi <- eset[, rownames(phenodata_fullroi)]
    top <- lapply(tissue, function(skin_layer){
      cat(paste(i, skin_layer, sep=" "), fill = TRUE)
      full_roi_data <- q3norm_fullroi %>%
        t() %>%
        as.data.frame() %>%
        mutate(day = phenodata_fullroi$day,
               type = phenodata_fullroi$type,
               bacteria = phenodata_fullroi$bacteria,
               type_day = phenodata_fullroi$type_day,
               segment = phenodata_fullroi$segment,
               section = phenodata_fullroi$section,
               tissue = phenodata_fullroi$tissue) %>%
        filter(tissue == skin_layer) %>%
        group_by(segment, day, type, type_day, bacteria, section, tissue) %>%
        summarise(across(everything(), mean))
      group <- as.character(as.data.frame(full_roi_data)[, gp])
      lvls_new <- intersect(lvls, group)
      group <- factor(group, levels = lvls_new)
      design <- model.matrix(~group+0)
      colnames(design) <- gsub("group", "", colnames(design))
      lmfit <- lmFit(t(full_roi_data[, rownames(q3norm_fullroi)]), design)
      if(i == "Full ROI"){
        cts_new <- cts
      } else {
        cts_new <- cts[-grep("Placebo", cts)]
      }
      args <- as.list(cts_new)
      args$levels = design
      cont <- do.call(makeContrasts, args)
      lmfit.cont <- contrasts.fit(lmfit, cont)
      lmfit.cont.ebayes <- eBayes(lmfit.cont)
      top <- lapply(colnames(cont), function(contrast){
        topTable(lmfit.cont.ebayes, coef = contrast,
                 adjust.method = "BH", n= nrow(lmfit.cont.ebayes)) %>%
          as.data.frame() %>%
          mutate(contrast = contrast) %>%
          mutate(feature = rownames(.)) %>%
          mutate(n = 1:n())
      }) %>%
        do.call(rbind, .)
      top$tissue <- skin_layer
      top
    }) %>%
      do.call(rbind, .)
    top$contrast <- factor(top$contrast)
    top$segment <- i
    top
  }) %>%
    do.call(rbind, .)
  top
}

## plot DE results
plot_de = function(top, cts, tissue_lvls){
  my_colors <- c("#999999","#E69F00", "#56B4E9")
  names(my_colors) <- levels(top$contrast)
  contrast_lvls <- unique(top$contrast)
  contrast_keep <- cts
  p <- top %>%
    filter(contrast %in% contrast_keep) %>%
    mutate(tissue = factor(tissue, levels = tissue_lvls),
           contrast = factor(contrast, levels = contrast_keep)) %>%
    ggplot(aes(x = n, y = adj.P.Val, color = contrast)) +
    geom_line() +
    scale_x_log10() +
    ylim(c(0,1)) +
    facet_grid(tissue~segment) +
    geom_hline(yintercept = 0.1, col = 'grey', linetype = "dashed") +
    theme(legend.position = "bottom") +
    xlab("Number of signficant genes") +
    ylab("BH-FDR") +
    theme_classic() +
    theme(legend.position = "bottom") +
    scale_color_manual(name = "Contrasts", values = my_colors)
  p
}

# Geneset enrichment analysis
## GSEA function
gsea2 = function(top, segments, tissues, dbs){
  top_db <- lapply(segments, function(seg){
    top <- lapply(tissues, function(tis){
      sig_top <- subset(top, segment == seg & tissue == tis)
      ranks <- lapply(unique(sig_top$contrast), function(cts){
        x <- subset(sig_top, contrast == cts)
        ranks <- x$t
        names(ranks) <- x$feature
        ranks
      })
      names(ranks) <- unique(sig_top$contrast)
      result <- lapply(names(ranks), function(cts){
        fgseaRes <- lapply(names(dbs), function(db){
          pathways <- lapply(dbs[[db]]$data, function(i){as.character(unlist(i))})
          names(pathways) <- dbs[[db]]$pathway
          stats = ranks[[cts]]
          # pathways <- lapply(pathways, function(i){
          #   intersect(i, names(stats))
          # })
          fgsea_result <- fgsea(pathways = pathways, stats = ranks[[cts]],
                                minSize  = 5, maxSize  = 500,
                                nPermSimple=10000)
          fgsea_result$DB <- db
          fgsea_result
        }) %>%
          do.call(rbind, .)
        fgseaRes$contrast <- cts
        fgseaRes
      }) %>%
        do.call(rbind, .) %>%
        as.data.frame()
      result$tissue <- tis
      result
    }) %>%
      do.call(rbind, .)
    top$segment <- seg
    top
  }) %>%
    do.call(rbind, .)
  top_db$padj[is.na(top_db$padj)] <- 1
  top_db
}

gsea = function(top, segments, tissues, dbs){
  top_db <- lapply(segments, function(seg){
    top <- lapply(tissues, function(tis){
      sig_top <- subset(top, segment == seg & tissue == tis)
      ranks <- lapply(unique(sig_top$contrast), function(cts){
        x <- subset(sig_top, contrast == cts)
        ranks <- x$t
        names(ranks) <- x$feature
        ranks
      })
      names(ranks) <- unique(sig_top$contrast)
      result <- lapply(names(ranks), function(cts){
        fgseaRes <- lapply(names(dbs), function(db){
          fgsea_result <- fgsea(pathways = dbs[[db]], stats = ranks[[cts]],
                                minSize  = 5, maxSize  = 500)
          fgsea_result$DB <- db
          fgsea_result
        }) %>%
          do.call(rbind, .)
        fgseaRes$contrast <- cts
        fgseaRes
      }) %>%
        do.call(rbind, .) %>%
        as.data.frame()
      result$tissue <- tis
      result
    }) %>%
      do.call(rbind, .)
    top$segment <- seg
    top
  }) %>%
    do.call(rbind, .)
  top_db$padj[is.na(top_db$padj)] <- 1
  top_db
}

## plot GSEA results
plot_gsea = function(top, cts, tissue_lvls){
  my_colors <- c("#999999","#E69F00", "#56B4E9")
  names(my_colors) <- levels(top$contrast)
  top %>%
    group_by(DB, contrast, tissue, segment) %>%
    arrange(padj) %>%
    mutate(n = 1:n()) %>%
    mutate(segment = factor(segment, segments)) %>%
    mutate(tissue = factor(tissue, tissue_lvls)) %>%
    ggplot(aes(x = n, y = padj, color = contrast, linetype = DB)) +
    geom_line() +
    facet_grid(tissue~segment, scales = "free") +
    scale_x_continuous(trans="log2") +
    geom_hline(yintercept = 0.1, col = 'grey', linetype = "dashed") +
    theme(legend.position = "bottom") +
    xlab("Number of signficant genesets") +
    ylab("BH-FDR") +
    theme_classic() +
    scale_color_manual(name = "Contrasts", values = my_colors)
}

## plot single geneset
plot_geneset = function(top, seg, tissues, db, geneset){
  res <- lapply(tissues, function(tis){
    sig_top <- subset(top, segment == seg & tissue == tis)
    # compare ranks
    ranks <- lapply(as.character(unique(sig_top$contrast)), function(cts){
      x <- subset(sig_top, contrast == cts)
      ranks <- x$t
      names(ranks) <- x$feature
      ranks
    })
    names(ranks) <- unique(sig_top$contrast)
    # plot enrichment results
    lapply(unique(sig_top$contrast), function(cts){
      plotEnrichment(db[[geneset]], ranks[[cts]]) +
        ggtitle(paste(geneset, cts, tis, sep=" - "))
    })
  })
  cowplot::plot_grid(plotlist = unlist(res, recursive = FALSE))
}
plot_geneset2 = function(top, seg, tissues, db, geneset){
  res <- lapply(tissues, function(tis){
    sig_top <- subset(top, segment == seg & tissue == tis)
    # compare ranks
    ranks <- lapply(unique(sig_top$contrast), function(cts){
      x <- subset(sig_top, contrast == cts)
      ranks <- x$t
      names(ranks) <- x$feature
      ranks
    })
    names(ranks) <- unique(sig_top$contrast)
    # plot enrichment results
    lapply(unique(sig_top$contrast), function(cts){
      plotEnrichment(db[[geneset]], ranks[[cts]]) +
        ggtitle(paste(geneset, cts, tis, sep=" - "))
    })
  })
  cowplot::plot_grid(plotlist = unlist(res, recursive = FALSE))
}

## plot gsea heatmap
gsea_heatmap = function(top, db, seg, cutoff, cts_tissue_ord, tissues, cts, ht_cols, top_k=NULL){
  top_sig <- top %>%
    filter(DB == db, segment == seg, padj < cutoff)

  if(!is.null(top_k)){
    top_sig <- top_sig %>%
      group_by(contrast, tissue, segment) %>%
      dplyr::slice(1:10)
  }

  nes0 <- top %>%
    filter(DB == db, segment == seg, pathway %in% unique(top_sig$pathway)) %>%
    mutate(comp_tissue = paste(contrast, tissue, sep=",")) %>%
    arrange(tissue) %>%
    dplyr::select(pathway, NES, comp_tissue) %>%
    spread(comp_tissue, NES, fill = 0)
  nes = nes0[, cts_tissue_ord]
  rownames(nes) <- gsub("Skin Human", "", nes0$pathway)

  padj0 <- top %>%
    filter(DB == db, segment == seg, pathway %in% unique(top_sig$pathway)) %>%
    mutate(comp_tissue = paste(contrast, tissue, sep=",")) %>%
    arrange(tissue) %>%
    dplyr::select(pathway, padj, comp_tissue) %>%
    spread(comp_tissue, padj, fill = 1)
  padj = padj0[, cts_tissue_ord]
  rownames(padj) <- padj0$pathway
  all(colnames(nes) == colnames(padj))
  cont <- sapply(strsplit(colnames(nes), ","), function(i){i[[1]]})
  tis <- sapply(strsplit(colnames(nes), ","), function(i){i[[2]]})
  col_ha <- columnAnnotation(Tissue = factor(tis, levels = tissues),
                             Comparison = factor(cont, levels = cts),
                             col = ht_cols)
  colnames(nes) <- NULL

  htmp <- Heatmap(nes, name = "NES",  top_annotation = col_ha,
                  cluster_columns = FALSE,
                  width = ncol(nes)*unit(8, "mm"),
                  heatmap_legend_param = list(
                    legend_direction = "horizontal",
                    legend_width = unit(3, "cm")),
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    if(padj[i, j] < cutoff)
                      grid.text("*", x, y, gp = gpar(fontsize = 10))
                  })
  draw(htmp, heatmap_legend_side="bottom",
       annotation_legend_side="left",
       legend_grouping = "original")
}

gsea_heatmap2 = function(top, db, dbs, seg, cutoff, cts_tissue_ord, tissues, cts, ht_cols, top_k=NULL){
  top_sig <- top %>%
    filter(DB == db, segment == seg, padj < cutoff)

  ann <- dbs[[db]] %>%
    dplyr::select(-data) %>%
    as.data.frame()
  rownames(ann) <- ann$pathway

  if(!is.null(top_k)){
    top_sig <- top_sig %>%
      group_by(contrast, tissue, segment) %>%
      dplyr::slice(1:10)
  }

  nes0 <- top %>%
    filter(DB == db, segment == seg, pathway %in% unique(top_sig$pathway)) %>%
    mutate(comp_tissue = paste(contrast, tissue, sep=",")) %>%
    arrange(tissue) %>%
    dplyr::select(pathway, NES, comp_tissue) %>%
    spread(comp_tissue, NES, fill = 0)
  nes = nes0[, cts_tissue_ord]
  rownames(nes) <- nes0$pathway
  padj0 <- top %>%
    filter(DB == db, segment == seg, pathway %in% unique(top_sig$pathway)) %>%
    mutate(comp_tissue = paste(contrast, tissue, sep=",")) %>%
    arrange(tissue) %>%
    dplyr::select(pathway, padj, comp_tissue) %>%
    spread(comp_tissue, padj, fill = 1)
  padj = padj0[, cts_tissue_ord]
  rownames(padj) <- padj0$pathway
  all(colnames(nes) == colnames(padj))
  cont <- sapply(strsplit(colnames(nes), ","), function(i){i[[1]]})
  tis <- sapply(strsplit(colnames(nes), ","), function(i){i[[2]]})
  col_ha <- columnAnnotation(Tissue = factor(tis, levels = tissues),
                             Comparison = factor(cont, levels = cts),
                             col = ht_cols)
  colnames(nes) <- NULL
  if(db == "panglaodb"){
    row_ha <- rowAnnotation(organ = ann[rownames(nes),"organ"])
  } else if(db == "reactome"){
    row_ha <- rowAnnotation(level1 = ann[rownames(nes), "level1_name"])
  } else if(db == "kegg"){
    row_ha <- rowAnnotation(level1 = ann[rownames(nes), "level1"],
                            level2 = ann[rownames(nes), "level2"])
  }

  htmp <- Heatmap(nes, name = "NES",
                  top_annotation = col_ha,
                  right_annotation = row_ha,
                  cluster_columns = FALSE,
                  width = ncol(nes)*unit(8, "mm"),
                  heatmap_legend_param = list(
                    legend_direction = "horizontal",
                    legend_width = unit(3, "cm")),
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    if(padj[i, j] < cutoff)
                      grid.text("*", x, y, gp = gpar(fontsize = 10))
                  })
  draw(htmp, heatmap_legend_side="bottom",
       annotation_legend_side="left",
       legend_grouping = "original")
}
