library(dplyr)

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

# PCA analysis helper function
## join points to centroids
StatCentSeg <- ggplot2::ggproto("StatCentSeg", Stat,
                                compute_group = function(data, scales, params,
                                                         cfun=median) {
                                  data$xend <- cfun(data$x)
                                  data$yend <- cfun(data$y)
                                  return(data)
                                },
                                required_aes = c("x", "y")
)
stat_centseg <- function(mapping = NULL, data = NULL, geom = "segment",
                         position = "identity", na.rm = FALSE, show.legend = NA,
                         inherit.aes = TRUE, cfun=median, ...) {
  layer(
    stat = StatCentSeg, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, cfun = cfun, ...)
  )
}

## perform pca
compute_pcs = function(pdata, eset, segments){
  pcs <- lapply(segments, function(i){
    phenodata_fullroi <- subset(pdata, segment == i)
    q3norm_fullroi <- eset[, rownames(phenodata_fullroi)]
    pcs <- prcomp(t(q3norm_fullroi), scale. = TRUE, center = TRUE, rank. = 3)
    prop_var <- round(100*summary(pcs)$importance["Proportion of Variance", c("PC1", "PC2")], 1)
    pcs$x[, 1:2] %>%
      as.data.frame() %>%
      mutate(type_day = factor(phenodata_fullroi$type_day,
                               levels = c("Placebo", "BCG_Day1", "BCG_Day7", "BCG_Day14")),
             tissue = factor(phenodata_fullroi$tissue, levels = c("Epidermis", "Dermis+Hypodermis")),
             segment = phenodata_fullroi$segment,
             bacteria = phenodata_fullroi$bacteria,
             day = phenodata_fullroi$day,
             section = phenodata_fullroi$section)
  }) %>%
    do.call(rbind, .)
  pcs
}

## plot pcs
plot_pcs = function(pcs, segments, col_by){
  pcs %>%
    mutate(segment = factor(segment, segments)) %>%
    filter(segment %in% segments) %>%
    ggplot(aes(x = PC1, y = PC2, color = {{col_by}})) +
    geom_point() +
    stat_ellipse(level=0.95) +
    facet_grid(tissue~segment, scales = "free") +
    stat_centseg(alpha=0.3) +
    theme_classic() +
    # theme(legend.position = "bottom") +
    xlab(paste0("PC1")) +
    ylab(paste0("PC2"))
}
