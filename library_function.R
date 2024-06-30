library(readxl)

library(magrittr)
library("ggvenn")
library(gtools)
library(stringr)
library(networkD3)
library(dplyr)
library(UniProt.ws)
library(data.table)
library(tibble)
# library(NormalyzerDE)
library(Hmisc)
library(plyr)
library(colourpicker)
library(RColorBrewer)
library(grDevices)
library(reshape)
library(VennDiagram)
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
library(pheatmap)
library(ComplexHeatmap)
library(factoextra)
library(ggpubr)
library(ggrepel)
library(corrplot)
library(dendsort)
library(mclust)
library(preprocessCore)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)
library(gridExtra)
library(WGCNA)
library(NMF)
library(DEP)
library(tidyr)
library(purrr)
library(SummarizedExperiment)
library(GGally)
library(rstatix)
library(readr)
library(ropls)
library(ggplot2)
library(ggsci)
library(Cairo)
library(tidyverse)
library(extrafont)
library(RVAideMemoire)
library(mixOmics)
library(vegan)
library(ggpmisc)
library(Cairo)
library(gimme)
library(gtools)
library(lme4)
library(lmerTest)
library(MuMIn)
library(nlme)
library(GGally)
library(naturalsort)
library(xlsx)
library(ggExtra)
library(gghalves)
library(openxlsx)
loadfonts()

mytheme = theme_bw() +
  theme(axis.text = element_text(color = 'black', size = 13),
        axis.ticks.length = unit(0.2, "lines"),
        axis.ticks = element_line(size = 0.5, color = 'black'),
        title = element_text(color = 'black', size = 13),
        panel.border = element_rect(size = 1),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.position = 'right'
  )
mytheme2 = theme(axis.text = element_text(color = 'black', size = 13),
        axis.ticks.length = unit(0.2, "lines"),
        axis.ticks = element_line(size = 0.5, color = 'black'),
        title = element_text(color = 'black', size = 13),
        panel.background = element_rect(fill = '#F4F4F4'),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.position = 'right'
  )
mytheme3 = theme_classic() + theme(axis.text = element_text(color = 'black', size = 13),
                 axis.ticks.length = unit(0.2, "lines"),
                 axis.ticks = element_line(size = 0.5, color = 'black'),
                 title = element_text(color = 'black', size = 13),
                 legend.text = element_text(size = 13),
                 legend.title = element_text(size = 13),
                 legend.position = 'right'
)

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

twowise_boxplot = function(df, feature_wise = F, sample_wise = T, ylabel = 'Intensity'){
  acc = df[, 1]
  df = data.frame(colnames(df)[-1], t(df[-1]))
  colnames(df) = c('sample_id', acc)
  df = melt(df, id.vars = 'sample_id')
  if(feature_wise == T){
    p1 = ggplot(df, aes(x = variable, y = value, color = variable)) +
      geom_boxplot(outlier.size = 0.3, outlier.shape = 1, outlier.fill = NULL, size = 0.2, color = 'black', fill = '#90EE90', outlier.alpha = 0.6) +
      labs(y = ylabel, x = "Metabolites") +
      theme_bw() +
      theme(legend.position = "none", axis.text = element_text(colour = 'black', size = 5), 
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.ticks.length = unit(0.1, "lines"), panel.grid = element_line(size = 0.1),
            axis.ticks = element_line(size = 0.2, color = 'black'),
            axis.title = element_text(size = 6), plot.title = element_text(size = 8, hjust = 0.5))
    print(p1)
  }
  if(sample_wise == T){
    p2 = ggplot(df, aes(x = sample_id, y = value, fill = str_split_fixed(sample_id, '_', 2)[, 1])) +
      geom_boxplot(outlier.size = 0.3, outlier.shape = 1, outlier.fill = NULL, size = 0.2, color = 'black', outlier.alpha = 0.6, alpha = 0.3) +
      labs(y = ylabel, x = "Samples") +
      scale_fill_manual(values = c('#374E55', '#DF8F44')) +
      theme_bw() +
      theme(legend.position = "none", axis.text = element_text(colour = 'black', size = 5), 
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.ticks.length = unit(0.1, "lines"), panel.grid = element_line(size = 0.1),
            axis.ticks = element_line(size = 0.2, color = 'black'),
            axis.title = element_text(size = 6), plot.title = element_text(size = 8, hjust = 0.5))
    print(p2)
  }
}


myPCA = function(pcaMat = pcaMat, pcaGroup, plotFlag = F, filePath = './', pal = 'jama', shape = 21){
  res.pca = stats::prcomp(pcaMat, scale = T, center = T)
  imp.pca = summary(res.pca)$importance
  pcDF = data.frame(pcs = colnames(imp.pca), exp = imp.pca[2, ])
  pcDF$PC_index = as.factor(1:nrow(pcDF))
  xlabel = paste0("PC", 1, " (", round(100 * imp.pca[2, 1], 1), "%)")
  ylabel = paste0("PC", 2, " (", round(100 * imp.pca[2, 2], 1), "%)")
  if(plotFlag){pdf(filePath, height = 3.5, width = 4.42)}
  ppca = fviz_pca_ind(res.pca, geom.ind = "point", pointshape = 21, pointsize = 2,
                      fill.ind = pcaGroup, alpha.ind = 0.7,
                      palette = pal,
                      addEllipses = T, ellipse.level = 0.95, ellipse.border.remove = T,
                      legend.title = "Group") +
    labs(x = xlabel, y = ylabel, title = NULL) + 
    mytheme + theme(panel.grid = element_line(linetype = 3, size = 0.2, colour = 'gray'))
  print(ppca)
  pscree = ggplot(pcDF[1:5, ], aes(x = PC_index, y = exp, group = 1)) +
    geom_line(size = 0.5) +
    geom_col(fill = "#90EE90", color = 'black', width = 0.7) + 
    geom_point(size = 2, shape = 16) +
    geom_text(aes(label = paste0(round(exp * 100, 1), '%')), hjust = 0, vjust = -0.5, angle = 0) +
    scale_y_continuous(limits = c(0, 0.30)) +
    labs(x = 'PC index', y = 'Percentage of explained variances') + mytheme + 
    theme(panel.border = element_blank(), panel.background = element_rect(fill = '#F2F2F2'), panel.grid = element_line(color = 'white'))
  print(pscree)
  # Each gray line connects two samples from the same individual taken four years apart.
  pcx = res.pca[['x']] %>% as.data.frame()
  pcx$Group = pcaGroup
  pcx$label = str_split_fixed(row.names(pcx), '_', 2)[, 2]
  pcx$label %<>% as.factor
  pcxp1 = ggscatter(data = pcx, x = "PC1", y = "PC2", fill = 'Group', color = 'black', alpha = 0.7, shape = 21, 
                    size = 2.5, alpha = 0.7, font.label = 12, repel = F, palette = pal,
                    ellipse = T, ellipse.type = "norm", ellipse.border.remove = T, mean.point = T, 
                    ellipse.alpha = 0.1, ellipse.level = 0.95) + 
    geom_line(aes(group = label), color = 'darkgray', size = 0.1) + 
    labs(x = xlabel, y = ylabel) +
    geom_vline(xintercept = 0, size = 0.25, linetype = 2, color = 'black') +
    geom_hline(yintercept = 0, size = 0.25, linetype = 2, color = 'black') +
    mytheme + theme(panel.grid = element_blank())
  print(pcxp1)
  
  # PERMANOVA
  dist = vegdist(pcaMat, method = 'euclidean')
  site = data.frame(sample = rownames(pcaMat),
                    group = pcaGroup)
  adonis_result_dis = adonis(dist~group, site, permutations = 1000)
  adonis_result_dis 
  pcxp1_with_P = pcxp1 + annotate(x = -10, y = -10, geom = 'text', label = paste0('PERMANOVA, p=', 
                                                                                  p_format(adonis_result_dis$aov.tab$`Pr(>F)`[1])))
  print(pcxp1_with_P)
  # pcx$Stage = pcaGroup2
  # pcx = pcx[!is.na(pcx$Stage), ]
  # pcxp2 = ggscatter(data = pcx, x = "PC1", y = "PC2", fill = 'Group', color = 'Group', alpha = 0.7, shape = 'Stage', 
  #                  size = 2.5, alpha = 0.7, font.label = 15, repel = F, palette = pal,
  #                  ellipse = T, ellipse.type = "norm", ellipse.border.remove = T, mean.point = F, 
  #                  ellipse.alpha = 0.1, ellipse.level = 0.95) + 
  #   geom_line(aes(group = label), color = 'darkgray', size = 0.1) + 
  #   labs(x = xlabel, y = ylabel) +
  #   geom_vline(xintercept = 0, size = 0.25, linetype = 2, color = 'black') +
  #   geom_hline(yintercept = 0, size = 0.25, linetype = 2, color = 'black') +
  #   mytheme + theme(panel.grid = element_blank())
  # print(pcxp2)
  # 
  # pcx$Sample = row.names(pcx)
  # pcxp3 = ggscatter(data = pcx, x = "PC1", y = "PC2", fill = 'Group', color = 'Group', alpha = 0.7, shape = 'Stage', 
  #                   size = 2.5, alpha = 0.7, font.label = 15, repel = F, palette = pal,
  #                   ellipse = T, ellipse.type = "norm", ellipse.border.remove = T, mean.point = F, 
  #                   ellipse.alpha = 0.1, ellipse.level = 0.95) + 
  #   geom_line(aes(group = label), color = 'darkgray', size = 0.1) + 
  #   labs(x = xlabel, y = ylabel) +
  #   geom_vline(xintercept = 0, size = 0.25, linetype = 2, color = 'black') +
  #   geom_hline(yintercept = 0, size = 0.25, linetype = 2, color = 'black') +
  #   geom_text_repel(aes(label = Sample), size = 3.5) + 
  #   mytheme + theme(panel.grid = element_blank())
  # print(pcxp3)

  if(plotFlag){dev.off()}
}



myPLSDA = function(pcaMat = pcaMat, pcaGroup, plotFlag = F, filePath = './', pal = 'jama'){
  plsda.datatm = plsda(pcaMat, pcaGroup, ncomp = 2)
  plsda.vip = PLSDA.VIP(plsda.datatm, graph = F)
  plsda.vip = plsda.vip$tab
  a = unclass(plsda.datatm)
  plotdata = as.data.frame(a$variates$X)
  plotdata$Group = pcaGroup
  plotdata$label = str_split_fixed(row.names(plotdata), '_', 2)[,2]
  eig = a$explained_variance$X %>% round(., 3) * 100
  if(plotFlag){pdf(filePath, height = 5, width = 6)}
  pcxp = ggscatter(data = plotdata, x = "comp1", y = "comp2", fill = 'Group', color = 'black', alpha = 0.7, shape = 21, 
                   size = 2.5, alpha = 0.7, font.label = 12, repel = F, palette = pal,
                   ellipse = T, ellipse.type = "norm", ellipse.border.remove = T, mean.point = TRUE,
                   ellipse.alpha = 0.1, ellipse.level = 0.95) + 
    geom_line(aes(group = label), color = 'darkgray', size = 0.1) + 
    labs(x = paste0('P1 (', eig[1], '%)'), y = paste0('P2 (', eig[2], '%)')) +
    geom_vline(xintercept = 0, size = 0.25, linetype = 2, color = 'black') +
    geom_hline(yintercept = 0, size = 0.25, linetype = 2, color = 'black') +
    mytheme + theme(panel.grid = element_blank())
  print(pcxp)
  # segment
  plsda.vip$feature = row.names(plsda.vip)
  vipDF = plsda.vip[1:20, ]
  vipDF$feature %<>% factor(., levels = vipDF$feature)
  ff = ggplot(vipDF, aes(feature, VIP)) +
    geom_segment(aes(x = feature, xend = feature, y = 0, yend = VIP)) +
    geom_point(shape = 21, size = 3, color = '#008000' ,fill = '#008000') +
    geom_hline(yintercept = 1, linetype = 'dashed') +
    scale_y_continuous(expand = c(0,0), limits = c(0, max(vipDF$VIP) + max(vipDF$VIP) * 0.15)) +
    labs(x = '', y = 'VIP') +
    mytheme3 + theme(axis.text.x = element_text(angle = 60, hjust = 1), plot.margin = unit(c(5,1,1,25),'mm'))
  print(ff)
  if(plotFlag){dev.off()}
}




# require('CancerSubtypes')
myNMF = function(
  main_title,
  seed = 6,
  nmf_data,
  cluster_num = 2,
  nrun = 50,
  cluster_colors = brewer.pal(n = 8, name = 'Set2')
){

  # *** Execute NMF *** #
  cat("                                                     \n")
  cat("*****************************************************\n")
  cat(paste("Execute NMF: ", "Cluster=", cluster_num, "  "))
  set.seed(seed)
  system.time({
    nmf_result = ExecuteCNMF(nmf_data, clusterNum = cluster_num, nrun = nrun)
  })
  nmf_group = nmf_result$group
  nmf_distanceMatrix = nmf_result$distanceMatrix

  # *** Layout *** #
  cat("                                                     \n")
  cat("*****************************************************\n")
  cat(paste("Layout", "  ", sep = ''))
  if(cluster_num > 1) {
    cat("                                                     \n")
    cat("*****************************************************\n")
    cat(paste(main_title, "Cluster=", cluster_num, "  "))
  }else {
    cat("There is only one cluster in the group")
  }
  if(!is.null(nmf_distanceMatrix[1, 1])) {
    layout(
      matrix(c(1, 1,1,1,1,1, 2, 2,2,2), 2, 5, byrow = FALSE), 
      widths = c(2.2, 2), 
      heights = c(2, 2)
    )
  }
  
  # *** Heat map *** #
  cat("                                                     \n")
  cat("*****************************************************\n")
  cat(paste("Heat map", "  ", sep = ''))
  if(!is.null(nmf_distanceMatrix[1, 1])){
    nmf_cluster = nmf_group
    if (class(nmf_distanceMatrix) == "Similarity") {
      si = silhouette_SimilarityMatrix(nmf_cluster, nmf_distanceMatrix)
    }else {
      si = silhouette(nmf_cluster, nmf_distanceMatrix)
    }
    attr(nmf_distanceMatrix, "class") = NULL
    ind = order(nmf_cluster, -si[, "sil_width"])
    num = length(unique(nmf_cluster))
    annotation = data.frame(Cluster = as.factor(nmf_cluster))
    Var1 = cluster_colors
    names(Var1) = sort(unique(nmf_cluster))
    ann_colors = list(Cluster = Var1)
    consensusmap(
      nmf_distanceMatrix, Rowv = ind, Colv = ind, 
      main = "Clustering display", annCol = annotation, 
      annColors = ann_colors, 
      # labRow = "Sample", labCol = "Sample", 
      scale = "none"
    )
  }else{
    print(paste("nmf_distanceMatrix is NULL", "  ", sep = ''))
  }
  
  # *** silhouette *** #
  cat("                                                     \n")
  cat("*****************************************************\n")
  cat(paste("Silhouette", "  ", sep = ''))
  
  si_col = rep(
    cluster_colors[1:cluster_num],
    as.vector(table(nmf_cluster))
  )
  plot(si, col = si_col)
  par(mfrow = c(1, 1))
  
  # *** Return NMF result *** #
  cat("                                                     \n")
  cat("*****************************************************\n")
  return(nmf_cluster)
}





volcano_plot = function(res_df, x = 'pi_logFC', y = 'FDR',FC_cutoff = 2, Pval_cutoff = 0.05, keyGenes = NULL){
  colnames(res_df)
  res_df$logFC = res_df[, x]
  densityplot(res_df$logFC)
  res_df$FDR = res_df[, y]
  res_df = res_df[!is.na(res_df$FDR), ]
  res_df = res_df[!is.na(res_df$logFC), ]
  quantile(res_df$FDR)
  
  res_df$keyGenes = '0'
  res_df$keyGenes[res_df$Hugo_Symbol %in% keyGenes] = '1'
  table(res_df$keyGenes)
  
  res_df$label = NA
  res_df$label[res_df$logFC >= log2(FC_cutoff) & res_df$FDR < Pval_cutoff] = 'Sig.Up'
  res_df$label[res_df$logFC <= log2(1/FC_cutoff) & res_df$FDR < Pval_cutoff] = 'Sig.Down'
  res_df$label[is.na(res_df$label)] = paste0('Unsig [', sum(is.na(res_df$label)), ']')
  res_df$label[which(res_df$label == 'Sig.Up')] = paste0('Sig.Up [', sum(res_df$label == 'Sig.Up'), ']')
  res_df$label[which(res_df$label == 'Sig.Down')] = paste0('Sig.Down [', sum(res_df$label == 'Sig.Down'), ']')
  uniLabel = unique(res_df$label)
  res_df$label %<>% factor(., levels = c(uniLabel[grep('Sig.Down', uniLabel)], uniLabel[grep('Sig.Up', uniLabel)], uniLabel[grep('Unsig', uniLabel)]))
  table(res_df$label)
  
  mycol = c('#01BFC4', '#F8766D', 'gray')
  library(ggnewscale)
  ptt = ggplot(res_df, aes(x = logFC, y = -log10(FDR), fill = label)) +
    geom_point(aes(color = label), shape = 21, alpha = 0.4, size = 2) + 
    geom_vline(xintercept = c(-log2(FC_cutoff), log2(FC_cutoff)), linetype = 2, size = 0.3) + 
    geom_hline(yintercept = -log10(Pval_cutoff), linetype = 2, size = 0.3) + 
    scale_fill_manual(values = mycol) + 
    scale_color_manual(values = mycol) + 
    new_scale_color() +
    geom_point(data = res_df[res_df$keyGenes == '1', ], 
               aes(x = logFC, y = -log10(FDR)),
               color = 'black', shape = 21, size = 2, show.legend = F) +
    # scale_color_manual(values = c(NA, 'black')) + 
    labs(x = 'log2(Fold Change)', y = '-log10(adjusted P value)') +
    geom_text_repel(aes(label = ifelse(keyGenes == '1', Hugo_Symbol, NA)), size = 4) +
    mytheme + theme(panel.grid = element_line(size = 0.4, linetype = 'dotted'), 
                    legend.title = element_blank(), legend.margin = margin(unit(c(0, 0, -5, 0), 'mm')),
                    legend.position = 'top')
  print(ptt)
}

# 6
get_topx_df_1 <- function(df, topx){
  mat_value = df[,-1]
  mat_value_col = ncol(mat_value)
  topx_seq = seq(1, topx)
  for(i in 1:mat_value_col){
    x = as.vector(unlist(mat_value[,i]))
    x_order = order(x, decreasing = T)
    x_order_index = x_order[topx_seq]
    y = x
    y[-x_order_index] = 0
    mat_value[,i] = y
  }
  kept_index = which(rowSums(mat_value)>0)
  topx_df = data.frame(
    Symbol = as.vector(unlist(df[,1])),
    mat_value
  )[kept_index,]
  return(topx_df)
}

# 8
fot_df_correction <- function(df){
  Symbol = as.vector(unlist(df[,1]))
  Value = df[,-1]
  Value_row = nrow(Value)
  Value_col = ncol(Value)
  correction_value = matrix(0, Value_row, Value_col)
  colnames(correction_value) = colnames(Value)
  for(i in 1:Value_col){
    i_v = as.vector(unlist(Value[,i]))
    i_v_fot = i_v/sum(i_v)*1e5
    correction_value[,i] = i_v_fot
  }
  result_df = data.frame(Symbol, correction_value)
  return(result_df)
}
if(F){
global_ifot_data_1_mat_row_mean = apply(global_ifot_data_1_mat, 1, mean)
global_ifot_data_1_mat_row_mean_sort_desc_cumsum = cumsum(
  sort(global_ifot_data_1_mat_row_mean, decreasing = T)
)
global_ifot_data_1_mat_row_mean_sort_desc_cumsum_probs = global_ifot_data_1_mat_row_mean_sort_desc_cumsum/max(
  global_ifot_data_1_mat_row_mean_sort_desc_cumsum
)
plot(
  seq(1, length(global_ifot_data_1_mat_row_mean_sort_desc_cumsum_probs)), 
  global_ifot_data_1_mat_row_mean_sort_desc_cumsum_probs,
  xlab = 'Gene rank',
  ylab = 'Cumulative probability'
)
cumsum_cutoff = 0.80
abline(h = cumsum_cutoff, col = 2, lty = 3)
high_intensity_genes = names(global_ifot_data_1_mat_row_mean_sort_desc_cumsum_probs)[
  global_ifot_data_1_mat_row_mean_sort_desc_cumsum_probs < cumsum_cutoff
]
print(sort(high_intensity_genes))
print(paste('high_intensity_genes', length(high_intensity_genes)))
high_intensity_genes_index = match(high_intensity_genes, global_ifot_data_1_genes)
global_ifot_data_2 = global_ifot_data_1[-high_intensity_genes_index,]
}






PlotPathSummary_wx = function (mSetObj = NA, show.grid, imgName, format = "png", 
          dpi = 72, width = NA, x, y) 
{
  mSetObj <- .get.mSet(mSetObj)
  if (mSetObj$analSet$type == "pathora") {
    x <- mSetObj$analSet$ora.mat[, 8]
    y <- mSetObj$analSet$ora.mat[, 4]
    names(x) <- names(y) <- rownames(mSetObj$analSet$ora.mat)
    if (!.on.public.web) {
      path.nms <- rownames(mSetObj$analSet$ora.mat)
    }
  }
  else if (mSetObj$analSet$type == "pathqea") {
    x <- mSetObj$analSet$qea.mat[, 7]
    y <- mSetObj$analSet$qea.mat[, 3]
    names(x) <- names(y) <- rownames(mSetObj$analSet$qea.mat)
    if (!.on.public.web) {
      path.nms <- rownames(mSetObj$analSet$qea.mat)
    }
  }
  else if (mSetObj$analSet$type == "pathinteg") {
    x <- mSetObj$dataSet$path.mat[, 8]
    y <- mSetObj$dataSet$path.mat[, 4]
    names(x) <- names(y) <- rownames(mSetObj$dataSet$path.mat)
    if (!.on.public.web) {
      path.nms <- rownames(mSetObj$analSet$jointPAMatches)
    }
  }
  else {
    print(paste("Unknown analysis type: ", mSetObj$analSet$type))
    return(0)
  }
  y = -log10(y)
  inx <- order(y, decreasing = T)
  x <- x[inx]
  y <- y[inx]
  sqx <- sqrt(x)
  min.x <- min(sqx, na.rm = TRUE)
  max.x <- max(sqx, na.rm = TRUE)
  if (min.x == max.x) {
    max.x = 1.5 * max.x
    min.x = 0.5 * min.x
  }
  maxR <- (max.x - min.x)/40
  minR <- (max.x - min.x)/160
  radi.vec <- minR + (maxR - minR) * (sqx - min.x)/(max.x - 
                                                      min.x)
  bg.vec <- heat.colors(length(y))
  if (.on.public.web) {
    if (mSetObj$analSet$type == "pathinteg") {
      path.nms <- names(current.kegglib$path.ids)[match(names(x), 
                                                        current.kegglib$path.ids)]
    }
    else {
      path.nms <- names(current.kegglib$path.ids)[match(names(x), 
                                                        current.kegglib$path.ids)]
    }
  }
  bg = "white"
  imgName = paste(imgName, "dpi", dpi, ".", format, sep = "")
  if (is.na(width)) {
    w <- 7
  }
  else if (width == 0) {
    w <- 7
  }
  else {
    w <- width
  }
  h <- w
  mSetObj$imgSet$path.overview <- imgName
  Cairo::Cairo(file = imgName, unit = "in", dpi = dpi, width = w, 
               height = h, type = format, bg = bg)
  op <- par(mar = c(6, 5, 2, 3))
  plot(x, y, type = "n", axes = F, xlab = "Pathway Impact", 
       ylab = "-log10(p)")
  axis(1)
  axis(2)
  if (show.grid) {
    grid(col = "blue")
  }
  symbols(x, y, add = TRUE, inches = F, circles = radi.vec, 
          bg = bg.vec, xpd = T)
  if (dpi == 72) {
    width.px <- height.px <- w * dpi
    mSetObj$imgSet$circleInfo <- CalculateCircleInfo(x, 
                                                     y, radi.vec, width.px, height.px, path.nms)
  }
  par(op)
  dev.off()
  return(.set.mSet(mSetObj))
}


darken <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}


lighten <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col*factor
  col <- rgb(t(col), maxColorValue=255)
  col
}

multiplot_t = function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols), byrow = T)
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

multi.which <- function(A){
  if ( is.vector(A) ) return(which(A))
  d <- dim(A)
  T <- which(A) - 1
  nd <- length(d)
  t( sapply(T, function(t){
    I <- integer(nd)
    I[1] <- t %% d[1]
    sapply(2:nd, function(j){
      I[j] <<- (t %/% prod(d[1:(j-1)])) %% d[j]
    })
    I
  }) + 1 )
}


my_fn1 <- function(data, mapping, ...){
  ggplot(data = data, mapping = mapping) + 
    geom_point(size = 1.5, shape = 21, color = '#6BAED6', fill = '#C6DBEF') + 
    geom_smooth(method=lm, size = 0.3, se = F, color = "#FF7F00", ...) + theme_bw()
}
my_fn2_spearman = function(data, mapping, ...){
  #main correlation
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  vals = cor.test(x,y, method="spearman")[c("estimate","p.value")]
  names(vals) = c("R = ","P = ")
  mainCor =  paste0(names(vals),signif(unlist(vals),2),collapse="\n")
  corr = vals["R = "]
  
  print(mainCor)
  colFn <- colorRampPalette(c("blue", "white", "red"), interpolate ='spline')
  fill <- colFn(100)[findInterval(corr, seq(-1, 1, length=100))]
  ggplot(data = data, mapping = mapping) +
    annotate(x = 0.5, y = 0.5, label = mainCor, geom = "text", size = 4) +
    theme_bw() +
    theme(panel.background = element_rect(fill=fill))
}
my_fn2_pearson = function(data, mapping, ...){
  #main correlation
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  vals = cor.test(x,y, method="pearson")[c("estimate","p.value")]
  names(vals) = c("R = ","P = ")
  mainCor =  paste0(names(vals),signif(unlist(vals),2),collapse="\n")
  corr = vals["R = "]
  
  print(mainCor)
  colFn <- colorRampPalette(c("blue", "white", "red"), interpolate ='spline')
  fill <- colFn(100)[findInterval(corr, seq(-1, 1, length=100))]
  ggplot(data = data, mapping = mapping) +
    annotate(x = 0.5, y = 0.5, label = mainCor, geom = "text", size = 4) +
    theme_bw() +
    theme(panel.background = element_rect(fill=fill))
}
source("http://peterhaschke.com/Code/multiplot.R") 



## recursion function
traverse <- function(df,a,i,innerl){
  if(i < (ncol(df))){
    alevelinner <- as.character(unique(df[which(as.character(df[,i])==a),i+1]))
    desc <- NULL
    if(length(alevelinner) == 1) (newickout <- traverse(alevelinner,i+1,innerl))
    else {
      for(b in alevelinner) desc <- c(desc,traverse(b,i+1,innerl))
      il <- NULL; if(innerl==TRUE) il <- a
      (newickout <- paste("(",paste(desc,collapse=","),")",il,sep=""))
    }
  }
  else { (newickout <- a) }
}
traverse <- function(a,i,innerl){
  if(i < (ncol(df))){
    alevelinner <- as.character(unique(df[which(as.character(df[,i])==a),i+1]))
    desc <- NULL
    ##if(length(alevelinner) == 1) (newickout <- traverse(alevelinner,i+1,innerl))
    ##else {
    for(b in alevelinner) desc <- c(desc,traverse(b,i+1,innerl))
    il <- NULL; if(innerl==TRUE) il <- a
    (newickout <- paste("(",paste(desc,collapse=","),")",il,sep=""))
    ##}
  }
  else { (newickout <- a) }
}
## data.frame to newick function
df2newick <- function(df, innerlabel=FALSE){
  alevel <- as.character(unique(df[,1]))
  newick <- NULL
  for(x in alevel) newick <- c(newick,traverse(x,1,innerlabel))
  (newick <- paste("(",paste(newick,collapse=","),");",sep=""))
}



jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}


mytheme4 = mytheme2 + theme(legend.position = "right", 
      panel.background = element_rect(fill = '#FAF6EE'), #F5EFE2
      panel.grid.major = element_line(size = 0.3)) 
