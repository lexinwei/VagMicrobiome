source('./script/library_function.R')
load('./data/RData/refDF_mags.RData')
load('./data/RData/metaData.RData')
load('./data/RData/swabData.RData')
load('./data/RData/magAbuDF.RData')
load('./data/RData/bin2rank.RData')
load('./data/RData/bin2name.RData')
load('./data/RData/colorList.RData')
plotFlag = T
#######################################################################################
#                              associated with some metadata                          #
#######################################################################################

if(T){
  refMetData = merge(refDF[!grepl('KIT', refDF$pt_ID), ], metaData, by = 'pt_ID', all.x = T)
  dim(refMetData)
  colnames(refMetData)
  refMetData$coverage %<>% as.numeric()
  refMetData$breadth %<>% as.numeric()
  refMetData$d_prime_mean  %<>% as.numeric()
  refMetData$SNV_per_kbp = refMetData$SNV_count*1000/refMetData$length
  refMetData$Employment = refMetData$Employed
  lineageAr = read.table('./MAGs/gtdbtk.ar53.summary.tsv', sep = '\t', header = T)
  lineageBac = read.table('./MAGs/gtdbtk.bac120.summary.tsv', sep = '\t', header = T)
  lineage = rbind(lineageBac, lineageAr)
  bin2rank = data.frame('ID' = lineage$user_genome, str_split_fixed(lineage$classification, ';', 7))
  colnames(bin2rank) = c('ID', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
  refMetData$genome2 = bin2name$name[match(refMetData$genome, bin2name$genome)]
  refMetData$species = bin2rank$species[match(refMetData$genome, paste0(bin2rank$ID, '.fa'))]
  refMetData$genus = bin2rank$genus[match(refMetData$genome, paste0(bin2rank$ID, '.fa'))]
  refMetData$family = bin2rank$family[match(refMetData$genome, paste0(bin2rank$ID, '.fa'))]
  table(refMetData$family)
  table(refMetData$species)
  refMetData$last_tax = refMetData$species
  refMetData$last_tax[refMetData$last_tax == 's__'] = refMetData$genus[refMetData$last_tax == 's__']
  refMetData$last_tax[refMetData$last_tax == 'g__'] = refMetData$family[refMetData$last_tax == 'g__']
  table(refMetData$last_tax == 'f__')
  table(refMetData$last_tax)
  refMetData$last_tax[refMetData$last_tax == 'f__'] = bin2rank$class[bin2rank$ID == 'Ming_nova_VS_2_46_bwa_bin.6']
  refMetData$genome_species = paste0(refMetData$last_tax, ' (', refMetData$genome2, ')')
  table(refMetData$genome_species) %>% sort(., decreasing = T)
  refMetData$abu = magAbuDF$relative_abundance[match(paste0(refMetData$genome, refMetData$pt_ID.u),
                                                     paste0(magAbuDF$Genome, '.fa', magAbuDF$pt_ID.u))]
  refMetData$if_Lactobacillaceae = 'Non-Lactobacillaceae'
  refMetData$if_Lactobacillaceae[refMetData$family == 'f__Lactobacillaceae'] = 'Lactobacillaceae'
  refMetData$if_Lactobacillaceae %<>% factor(., levels = c('Lactobacillaceae', 'Non-Lactobacillaceae'))
  table(refMetData$if_Lactobacillaceae)
  refMetData$D_rev = 1-refMetData$d_prime_mean
}
refMetData2 = refMetData[refMetData$genome_species %in% names(table(refMetData$genome_species)[table(refMetData$genome_species) > 20]) &
                         refMetData$species != 's__',]
table(refMetData2$species) %>% sort(., decreasing = T)
refMetData2$divergent_site_count %<>% as.numeric()
sumDivergentSite = refMetData2[, c('sample', 'divergent_site_count')] %>% 
  group_by(sample) %>% summarise(sum = sum(divergent_site_count))
range(sumDivergentSite$sum)
mean(sumDivergentSite$sum)
sumDivergentSite2 = refMetData2[, c('species', 'divergent_site_count')] %>% 
  group_by(species) %>% summarise(sum = sum(divergent_site_count))
range(sumDivergentSite2$sum)
mean(sumDivergentSite2$sum)
range(refMetData2$divergent_site_count)
#### relation between coverage and nucleotide diversity
relPlot = function(df, xval, yval, zval = NULL, gval = NULL, mycolor = NULL, xLab = NULL, yLab = NULL, showSE = T, xAngle = NULL, aspect.Ratio = 1, legendPos = 'right', facet.Nrow = 1, Title = NULL){
  df$x = df[, xval]
  df$y = df[, yval]
  if(!is.null(zval)){df$z = df[, zval]}else{df$z = 'One panel';df$z %<>% as.factor()}
  if(!is.null(gval)){df$g = df[, gval]}else{df$g = 'One group';df$g %<>% as.factor()}
  for(zs in levels(df$z)){
    cat(zs, '\n\n')
    for(gs in levels(df$g)){
      lmData = df[df$g == gs & df$z == zs, ]
      if(nrow(lmData) >= 3){
        lmFit = try({lm(y ~ x, lmData)})
        slope = signif(coefficients(lmFit)[-1],3)
        if(slope < 0){symb = '-'}else{symb = '+'}
        lmFormu = paste0("lm, y = ", signif(coefficients(lmFit)[1],2), ' ', symb, ' ', abs(slope),
                         names(coefficients(lmFit)[-1]))
        cat(gs, ': ', lmFormu, '\n', sep = '')
        lmSum = summary(lmFit)
        cat('p-value:', lmSum$coefficients[2, 4] %>% signif(3), 
            'Adj. R2:', lmSum$adj.r.squared %>% signif(3), 
            'R2:', lmSum$r.squared  %>% signif(3), '\n')
        corT_Pearson = cor.test(lmData$x, lmData$y, method = 'pearson')
        cat('Pearson, p-value =', signif(corT_Pearson$p.value, 3), 'r =',signif(corT_Pearson$estimate, 3), '\n')
        corT_Spearm = cor.test(lmData$x, lmData$y, method = 'spearm')
        cat('Spearman, p-value =', signif(corT_Spearm$p.value, 3), 'r =',signif(corT_Spearm$estimate, 3), '\n')
      }
    }
  }
  wx3 = ggplot(df, aes(x = x, y = y)) +
    geom_point(aes(fill = g),shape = 21, size = 2, alpha = 0.7, color = 'black', stroke = 0.1) +
    labs(x = xLab, y = yLab, title = Title) +
    geom_smooth(method = 'lm', aes(color = g), se = showSE) +
    scale_color_manual(values = mycolor, name = 'Group') +
    scale_fill_manual(values = mycolor, name = 'Group') +
    mytheme3 + theme(legend.spacing.y = unit(0.1, 'cm'), aspect.ratio = aspect.Ratio)
  if(!is.null(zval)){
    wx3 = wx3 + facet_wrap(.~z, scales = 'free_x', nrow = facet.Nrow) + theme(strip.text = element_text(size = 13))
  }
  if(!is.null(xAngle)){
    wx3 = wx3 + theme(axis.text.x = element_text(angle = xAngle, hjust = 1), legend.position = legendPos)
  }
  print(wx3)
}
lmFitFun = function(df, xval, yval, zval = NULL, gval = NULL, xLab = NULL, yLab = NULL){
  df$x = df[, xval]
  df$y = df[, yval]
  if(!is.null(zval)){df$z = df[, zval]}else{df$z = 'One panel';df$z %<>% as.factor()}
  if(!is.null(gval)){df$g = df[, gval]}else{df$g = 'One group';df$g %<>% as.factor()}
  Zlist = c()
  Glist = c()
  SlopeS = c()
  FormuList = c()
  pS = c()
  adjR2s = c()
  R2s = c()
  Xs =  c()
  Ys = c()
  for(zs in levels(df$z)){
    cat(zs, '\n')
    for(gs in levels(df$g)){
      lmData = df[df$g == gs & df$z == zs, ]
      if(nrow(lmData) >= 3){
        Zlist = c(Zlist, zs)
        Glist = c(Glist, gs)
        Xs = c(Xs, xLab)
        Ys = c(Ys, yLab)
        slope = signif(coefficients(lmFit)[-1],3)
        SlopeS = c(SlopeS, slope)
        if(slope < 0){symb = '-'}else{symb = '+'}
        lmFormu = paste0("lm, y = ", signif(coefficients(lmFit)[1], 2), ' ', symb, ' ', abs(slope),
                         names(coefficients(lmFit)[-1]))
        FormuList = c(FormuList, lmFormu)
        lmSum = summary(lmFit)
        pS = c(pS, lmSum$coefficients[2, 4] %>% signif(3))
        adjR2s = c(adjR2s, lmSum$adj.r.squared %>% signif(3)) 
        R2s = c(R2s, lmSum$r.squared  %>% signif(3))
      }
    }
  }
  ref_df = data.frame(
    Panel = Zlist,
    Group = Glist,
    X = Xs,
    Y = Ys,
    comp = paste0(Xs, ' vs ', Ys),
    Slope = SlopeS,
    Formula = FormuList,
    Pvalue = pS,
    Adj.R2 = adjR2s,
    R2 = R2s
  )
  return(ref_df)
}
corFitFun = function(df, xval, yval, zval = NULL, gval = NULL, xLab = NULL, yLab = NULL){
  df$x = df[, xval]
  df$y = df[, yval]
  if(!is.null(zval)){df$z = df[, zval]}else{df$z = 'One panel';df$z %<>% as.factor()}
  if(!is.null(gval)){df$g = df[, gval]}else{df$g = 'One group';df$g %<>% as.factor()}
  Zlist = c()
  Glist = c()
  Xs = Ys =c()
  PearsonP = c()
  PearsonR = c()
  SpearmP = c()
  SpearmR = c()
  for(zs in levels(df$z)){
    cat(zs, '\n')
    for(gs in levels(df$g)){
      lmData = df[df$g == gs & df$z == zs, ]
      if(nrow(lmData) >= 3){
        Zlist = c(Zlist, zs)
        Glist = c(Glist, gs)
        Xs = c(Xs, xLab)
        Ys = c(Ys, yLab)
        corT_Pearson = cor.test(lmData$x, lmData$y, method = 'pearson')
        PearsonP = c(PearsonP, corT_Pearson$p.value) 
        PearsonR = c(PearsonR, corT_Pearson$estimate) 
        corT_Spearm = cor.test(lmData$x, lmData$y, method = 'spearm')
        SpearmP = c(SpearmP, corT_Spearm$p.value) 
        SpearmR = c(SpearmR, corT_Spearm$estimate) 
      }
    }
  }
  ref_df = data.frame(
    Panel = Zlist,
    Group = Glist,
    X = Xs,
    Y = Ys,
    comp = paste0(Xs, ' vs ', Ys),
    PearsonP = PearsonP,
    PearsonR = PearsonR,
    SpearmP = SpearmP,
    SpearmR = SpearmR
  )
  return(ref_df)
}
relAbuPlot = function(data0, xval, yval, zval = NULL, yLab = 'Value',  titleLab = 'title', jitterSize = 1.8,
                      y.position = NULL,step.increase = 0.12, expandMult = NULL,leName = NULL, addJitter = T,aspectRatio = NULL,
                      x.angle = 0, CLD = F, hide.ns = F, legend.position = 'right', if_log10 = F, show_Padj = F){
  data0$x = data0[, xval]
  if(if_log10){
    data0$y = log10(data0[, yval])
  }else{
    data0$y = data0[, yval]
  }
  
  if(endsWith(xval, '2')){xLab = str_remove(xval, '2$')}else{xLab = xval}
  if(!is.null(zval)){data0$z = data0[, zval]}
  if(x.angle == 0){x.hjust = 0.5}else{x.hjust = 1}
  if(is.null(leName)){
    leName = xLab
  }
  # base ggplot
  if(addJitter){
    pbO = ggboxplot(data0, x = 'x', y = 'y', xlab = xLab,
                    ylab = yLab, color = 'x', title = titleLab,
                    add = "jitter", 
                    add.params = list(shape = 21, size = jitterSize, color = 'x', fill = 'x', alpha = 0.7),
                    bxp.errorbar = T, notch = F) +
      scale_color_manual(values = get(paste0('color', xval)), name = leName) +
      scale_fill_manual(values = get(paste0('color', xval)), name = leName) +
      mytheme3 + theme(strip.text = element_text(size = 12), axis.text.x = element_text(angle = x.angle, hjust = x.hjust),
                       plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 'cm'), legend.position = legend.position)
  }else{
    pbO = ggboxplot(data0, x = 'x', y = 'y', xlab = xLab,
                    ylab = yLab, color = 'black', fill = 'x', title = titleLab,
                    bxp.errorbar = T, notch = F) +
      scale_fill_manual(values = get(paste0('color', xval)), name = leName) +
      mytheme3 + theme(strip.text = element_text(size = 12), axis.text.x = element_text(angle = x.angle, hjust = x.hjust),
                       plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 'cm'), legend.position = legend.position)
  }
  
  
  if(!CLD){
    if(hide.ns){
      if(!is.null(zval)){
        stat.test <- data0 %>%
          group_by(z) %>%
          wilcox_test(y ~ x) %>%
          adjust_pvalue(method = "BH") %>%
          add_significance("p.adj") %>%
          add_xy_position(x = xval,dodge = 0, step.increase = step.increase)
      }else{
        stat.test <- data0 %>%
          wilcox_test(y ~ x) %>%
          adjust_pvalue(method = "BH") %>%
          add_significance("p.adj") %>%
          add_xy_position(x = xval,dodge = 0, step.increase = step.increase)
      }
      
      if(!show_Padj){
        stat.test = stat.test[stat.test$p<=0.05,]
        stat.test = BBmisc::dropNamed(stat.test, c('p.adj', 'p.adj.signif'))
        Plab = 'p'
      }else{
        stat.test = stat.test[stat.test$p.adj<=0.05,]
        Plab = 'p.adj'
      }
      
      if(!is.null(y.position)){
        stat.test$y.position = y.position
      }
      pbO = pbO + 
        # scale_y_continuous(expand = expansion(mult = c(0.05, 0.5))) +
        stat_pvalue_manual(stat.test, label = Plab, tip.length = 0.02, 
                           remove.bracket = F, hide.ns = T)
    }else{
      CN = combn(na.omit(droplevels(unique(data0$x))) %>% levels, 2, simplify = FALSE)
      pbO = pbO + stat_compare_means(comparisons = CN, method = 'wilcox.test', label = '..p..')
    }
  }else{
    if(!is.null(zval)){
      Lmodel <- lm(y ~ x * z, data = data0)
      model_means <- emmeans(object = Lmodel, specs = ~ x | z) 
    }else{
      Lmodel <- lm(y ~ x, data = data0)
      model_means <- emmeans(object = Lmodel, specs = ~ x) 
    }
    model_means_cld <- cld(object = model_means,
                           adjust = "Tukey",
                           Letters = letters,
                           alpha = 0.05) %>% as.data.frame()
    pbO = pbO +
      geom_text(aes(x, y = max(data0$y) + (max(data0$y) - min(data0$y)) *0.025, label= str_trim(.group)),
                data = model_means_cld, vjust = -0.5)
  }
  if(!is.null(expandMult)){
    pbO = pbO + scale_y_continuous(expand = expandMult)
  }
  if((!is.null(zval))){
    pbO = pbO + facet_grid(.~z)
  }
  if(!is.null(aspectRatio)){
    pbO = pbO + theme(aspect.ratio  = aspectRatio)
  }
  print(pbO)
}
if(F){
  ####   overview: x = SNVs/kbp , y = D', colored by pN/pS, size by coverage, shape by species group ####
  if(T){
    a = refMetData2[refMetData2$SNV_per_kbp < 10 & refMetData2$D_rev > 0.2,]
    refMetData2$SampleID = swabData$SampleID.u[match(refMetData2$sample, swabData$fqID)]
    paper_var = c('SampleID', 'genome_species', 'family', 
                  'length', 'coverage', 'breadth', 'breadth_minCov', 
                  'SNV_per_kbp','divergent_site_count',  'nucl_diversity', 'pnps', 'D_rev')
    paper_col = c('SampleID', 'Genome', 'Family',
                  'Length', 'Coverage', 'Breadth', 'Breadth_minCov',
                  'SNVs/kbp', 'SNPs count', 'π', 'pN/pS', "1-D'")
    refMetData2_paper = refMetData2[, paper_var]
    colnames(refMetData2_paper) = paper_col
    refMetData2_paper$Family %<>% str_remove(., 'f__')
    refMetData2_paper$Genome %<>% str_remove(., 's__')
    write.xlsx(refMetData2_paper, 
               file = './table/evoluionary_metrices_this_study_MAGs.xlsx')
    cairo_pdf('./plot_mags/SNVs_vs_D_vs_coverage_vs_pNpS_2.pdf', width = 8, height = 3, onefile = T)
    wx1 = ggplot(refMetData2, aes(x = SNV_per_kbp, y = 1-d_prime_mean, size = coverage, shape = if_Lactobacillaceae,
                                  fill = pnps)) +
      geom_point(alpha = 0.7, color = 'black', stroke = 0.1) +
      labs(x = 'SNVs/kbp', y = "1 - D'") +
      gradient_fill(c("#00AFBB", "#E7B800", "#FC4E07")) +
      # gradient_fill(c("white", "#FC4E07")) +
      guides(fill = guide_colourbar(title = 'pN/pS')) +
      scale_shape_manual(values = c(21, 22), name = 'Group', guide = guide_legend(override.aes = list(size=5))) +
      scale_size_continuous(range = c(1,6), name = 'Coverage') +
      mytheme3 + theme(legend.spacing.y = unit(0.1, 'cm'), 
                        aspect.ratio = 1)
    print(wx1)
    wx2= ggplot(refMetData2, aes(x = nucl_diversity, y = 1-d_prime_mean, size = coverage, shape = if_Lactobacillaceae,
                                 fill = pnps)) +
      geom_point(alpha = 0.7, color = 'black', stroke = 0.1) +
      labs(x = expression(pi), y = "1 - D'") +
      gradient_fill(c("#00AFBB", "#E7B800", "#FC4E07")) +
      # gradient_fill(c("white", "#FC4E07")) +
      guides(fill = guide_colourbar(title = 'pN/pS')) +
      scale_shape_manual(values = c(21, 22), name = 'Group', guide = guide_legend(override.aes = list(size=5))) +
      scale_size_continuous(range = c(1,6), name = 'Coverage') +
      mytheme3 + theme(legend.spacing.y = unit(0.1, 'cm'), 
                       aspect.ratio = 1)
    print(wx2)
    lmLac = lm(d_prime_mean ~ SNV_per_kbp, refMetData2[refMetData2$if_Lactobacillaceae == 'Lactobacillaceae', ])
    summary(lmLac)
    lmNonLac = lm(d_prime_mean ~ SNV_per_kbp, refMetData2[refMetData2$if_Lactobacillaceae == 'Non-Lactobacillaceae', ])
    summary(lmNonLac)
    wx3 = ggplot(refMetData2, aes(x = SNV_per_kbp, y = d_prime_mean)) +
      geom_point(refMetData2,mapping = aes(size = coverage, shape = if_Lactobacillaceae,
                                           fill = pnps), alpha = 0.7, color = 'black', stroke = 0.1) +
      labs(x = 'SNVs/kbp', y = "D'") +
      gradient_fill(c("#00AFBB", "#E7B800", "#FC4E07")) +
      # gradient_fill(c("white", "#FC4E07")) +
      guides(fill = guide_colourbar(title = 'pN/pS')) +
      scale_shape_manual(values = c(21, 22), name = 'Group', guide = guide_legend(override.aes = list(size=5))) +
      scale_size_continuous(range = c(1,6), name = 'Coverage') +
      geom_smooth(method = 'lm', aes(color = if_Lactobacillaceae)) +
      scale_color_manual(values = colorLactobacillaceae, name = 'Group') +
      mytheme3 + theme(legend.spacing.y = unit(0.1, 'cm'), 
                       aspect.ratio = 1)
    print(wx3)
    dev.off()
  }
  if(T){
    cairo_pdf('./plot_mags/pi_vs_pNpS_vs_coverage_vs_D.pdf', width = 8, height = 3, onefile = T)
    wx4 = ggplot(refMetData2, aes(x = nucl_diversity, y = pnps, size = coverage, shape = if_Lactobacillaceae,
                                  fill = 1-d_prime_mean)) +
      geom_point(alpha = 0.7, color = 'black', stroke = 0.1) +
      labs(x = expression(pi), y = "pN/pS") +
      gradient_fill(c("#00AFBB", "#E7B800", "#FC4E07")) +
      # gradient_fill(c("white", "#FC4E07")) +
      guides(fill = guide_colourbar(title = "1-D'")) +
      scale_shape_manual(values = c(21, 22), name = 'Group', guide = guide_legend(override.aes = list(size=5))) +
      scale_size_continuous(range = c(1,6), name = 'Coverage') +
      mytheme3 + theme(legend.spacing.y = unit(0.1, 'cm'), legend.position = 'right',
                       legend.box.spacing = unit(0.05, 'cm'), aspect.ratio = 1)
    print(wx4)
    dev.off()
  }
  if(T){ # overview: x = SNVs/kbp , y = D', colored by pN/pS, size by coverage, shape by species group
    cairo_pdf('./plot_mags/SNVs_vs_D_vs_coverage_vs_pNpS_facetBy_Ethnicity.pdf', width = 10.28, height = 4.48)
    wx2 = ggplot(refMetData2, aes(x = SNV_per_kbp, y = d_prime_mean, size = coverage, shape = if_Lactobacillaceae,
                                  fill = pnps)) +
      geom_point(alpha = 0.7, color = 'black', stroke = 0.1) +
      labs(x = 'SNVs/kbp', y = "D'") +
      gradient_fill(c("#00AFBB", "#E7B800", "#FC4E07")) +
      # gradient_fill(c("white", "#FC4E07")) +
      guides(fill = guide_colourbar(title = 'pN/pS')) +
      scale_shape_manual(values = c(21, 22), name = 'Group', guide = guide_legend(override.aes = list(size=5))) +
      scale_size_continuous(range = c(1,6), name = 'Coverage') +
      mytheme3 + theme(legend.spacing.y = unit(0.1, 'cm'), legend.box.spacing = unit(0.2, 'cm'), aspect.ratio = 1, 
                       legend.position = 'right', strip.text = element_text(size = 13)) +
      facet_grid(.~Ethnicity2)
    print(wx2)
    Eth = 'Other'
    lmLac = lm(d_prime_mean ~ SNV_per_kbp, refMetData2[refMetData2$if_Lactobacillaceae == 'Lactobacillaceae' & refMetData2$Ethnicity2 == Eth, ])
    summary(lmLac)
    lmNonLac = lm(d_prime_mean ~ SNV_per_kbp, refMetData2[refMetData2$if_Lactobacillaceae == 'Non-Lactobacillaceae' & refMetData2$Ethnicity2 == Eth, ])
    summary(lmNonLac)
    wx4 = ggplot(refMetData2, aes(x = SNV_per_kbp, y = d_prime_mean)) +
      geom_point(refMetData2,mapping = aes(size = coverage, shape = if_Lactobacillaceae,
                                           fill = pnps), alpha = 0.7, color = 'black', stroke = 0.1) +
      labs(x = 'SNVs/kbp', y = "D'") +
      gradient_fill(c("#00AFBB", "#E7B800", "#FC4E07")) +
      # gradient_fill(c("white", "#FC4E07")) +
      guides(fill = guide_colourbar(title = 'pN/pS')) +
      scale_shape_manual(values = c(21, 22), name = 'Group', guide = guide_legend(override.aes = list(size=5))) +
      scale_size_continuous(range = c(1,6), name = 'Coverage') +
      geom_smooth(method = 'lm', aes(color = if_Lactobacillaceae)) +
      scale_color_manual(values = colorLactobacillaceae, name = 'Group') +
      mytheme3 + theme(legend.spacing.y = unit(0.1, 'cm'), legend.box.spacing = unit(0.2, 'cm'), aspect.ratio = 1, 
                       legend.position = 'right', strip.text = element_text(size = 13)) +
      facet_grid(.~Ethnicity2, scales = 'free_x')
    print(wx4)
    dev.off()
  }
  table(refMetData2$genome_species) %>% length()
  pnpsDF = refMetData2 %>% group_by(genome) %>% 
    summarise(mean_pNpS = mean(pnps),
              mean_D = mean(D_rev, na.rm = T),
              mean_pi = mean(nucl_diversity),
              mean_coverage = mean(coverage))
  pnpsDF$Species =  refMetData2$genome_species[match(pnpsDF$genome, refMetData2$genome)]
  pnpsDF$Species = paste0(str_sub(pnpsDF$Species, 4, 4),
                          str_replace_all(pnpsDF$Species %>% str_remove(., 's__'), '^[A-Z][a-z]+', '.'))
  pnpsDF$Family = bin2rank$family[match(pnpsDF$genome %>% str_remove(., '\\.fa'), bin2rank$ID)] %>% str_remove('f__')
  # pnpsDF2 = pnpsDF[pnpsDF$Family %in% names(colorFamily), ]
  # unclassified Coriobacteriales
  
  table(pnpsDF$Family)
  # pnpsDF$Family[pnpsDF$Family == 'Atopobiaceae'] = 'Atopobiaceae | Actinomycetaceae'
  # pnpsDF$Family[pnpsDF$Family == 'Dialisteraceae'] = 'Dialisteraceae | Veillonellaceae'
  # pnpsDF$Family[pnpsDF$Family == 'Megasphaeraceae'] = 'Megasphaeraceae | Veillonellaceae'
  # pnpsDF$Family[pnpsDF$Family == 'Fastidiosipilaceae'] = 'Fastidiosipilaceae | Oscillospiraceae'
  # colorFamily_Mag = colorFamily
  # names(colorFamily_Mag)[names(colorFamily_Mag) == 'Actinomycetaceae'] = 'Atopobiaceae | Actinomycetaceae'
  # names(colorFamily_Mag)[names(colorFamily_Mag) == 'Veillonellaceae'] = 'Dialisteraceae | Veillonellaceae'
  # colorFamily_Mag = c(colorFamily_Mag, 'Megasphaeraceae | Veillonellaceae' = '#F994DA',
  #                     'Eggerthellaceae' = '#79D2D6',
  #                     'Fastidiosipilaceae | Oscillospiraceae' = '#75C4FD')\
  colorFamily_Mag = colorFamily
  names(colorFamily_Mag)[names(colorFamily_Mag) == 'Actinomycetaceae'] = 'Atopobiaceae'
  names(colorFamily_Mag)[names(colorFamily_Mag) == 'Veillonellaceae'] = 'Dialisteraceae'
  colorFamily_Mag = c(colorFamily_Mag, 'Megasphaeraceae' = '#F994DA',
                      'Eggerthellaceae' = '#E20023',
                      'Fastidiosipilaceae' = '#250000')
  
  wx5 = ggplot(pnpsDF, aes(x = mean_pi, y = mean_pNpS, size = mean_coverage, fill = Family)) +
    geom_vhlines(xintercept = median(pnpsDF$mean_pi), yintercept = median(pnpsDF$mean_pNpS), 
                 linetype = 'dashed', color = 'gray') +
    geom_point(alpha = 0.7, color = 'black', stroke = 0.3, shape = 21) +
    labs(x = expression(~~Avg. ~~ pi), y = "Avg. pN/pS") +
    scale_fill_manual(values = colorFamily_Mag, name = 'Family') +
    scale_size_continuous(range = c(2,6), name = 'Normalized abundance') +
    geom_text_repel(aes(label = Species), size = 2) +
    mytheme3 + theme(legend.spacing.y = unit(0.1, 'cm'), 
                     legend.box.spacing = unit(0.05, 'cm'), aspect.ratio = 1)
  cairo_pdf('./plot/MAGs_evo_type3.pdf', width = 6.84, height = 3.72, onefile = T)
  print(wx5)
  dev.off()
  #### one vs one ####
  if(T){ # overview: x = Pi , y = D', colored by pN/pS, size by coverage, shape by species group
    cairo_pdf('./plot_mags/evolutionary_metrics_one_vs_one.pdf', width = 8, height = 3, onefile = T)
    relPlot(refMetData2, xval = 'SNV_per_kbp', yval = 'd_prime_mean', xLab = 'SNVs/kbp', y = "D'", 
            gval = 'if_Lactobacillaceae', mycolor = colorLactobacillaceae)
    relPlot(refMetData2, xval = 'SNV_per_kbp', yval = 'D_rev', xLab = 'SNVs/kbp', y = "1 - D'", 
            gval = 'if_Lactobacillaceae', mycolor = colorLactobacillaceae)
    relPlot(refMetData2, xval = 'nucl_diversity', yval = 'd_prime_mean', xLab = expression(italic(pi)), yLab = "D'",
            gval = 'if_Lactobacillaceae', mycolor = colorLactobacillaceae)
    relPlot(refMetData2, xval = 'nucl_diversity', yval = 'D_rev', xLab = expression(italic(pi)), yLab = "1 - D'",
            gval = 'if_Lactobacillaceae', mycolor = colorLactobacillaceae)
    relPlot(refMetData2, xval = 'coverage', yval = 'SNV_per_kbp', xLab = 'Coverage', yLab = 'SNVs/kbp', 
            gval = 'if_Lactobacillaceae', mycolor = colorLactobacillaceae)
    relPlot(refMetData2, xval = 'coverage', yval = 'nucl_diversity', xLab = 'Coverage', yLab = expression(italic(pi)), 
            gval = 'if_Lactobacillaceae', mycolor = colorLactobacillaceae)
    relPlot(refMetData2[refMetData2$coverage < 500, ], xval = 'coverage', yval = 'nucl_diversity', xLab = 'Coverage', yLab = expression(italic(pi)), 
            gval = 'if_Lactobacillaceae', mycolor = colorLactobacillaceae)
    relPlot(refMetData2, xval = 'coverage', yval = 'd_prime_mean', xLab = 'Coverage', yLab = "D'",
            gval = 'if_Lactobacillaceae', mycolor = colorLactobacillaceae)
    relPlot(refMetData2, xval = 'coverage', yval = 'D_rev', xLab = 'Coverage', yLab = "1 - D'",
            gval = 'if_Lactobacillaceae', mycolor = colorLactobacillaceae)
    relPlot(refMetData2, xval = 'SNV_per_kbp', yval = 'pnps', xLab = 'SNVs/kbp', yLab = "pN/pS",
            gval = 'if_Lactobacillaceae', mycolor = colorLactobacillaceae)
    relPlot(refMetData2, xval = 'nucl_diversity', yval = 'pnps', xLab = expression(italic(pi)),  yLab = "pN/pS",
            gval = 'if_Lactobacillaceae', mycolor = colorLactobacillaceae)
    relPlot(refMetData2, xval = 'coverage', yval = 'pnps', xLab = 'Coverage', yLab = "pN/pS",
            gval = 'if_Lactobacillaceae', mycolor = colorLactobacillaceae)
    dev.off()
  }
  cairo_pdf('./plot_mags/species_pi_vs_pnps.pdf', width = 20, height = 10)
  relPlot(refMetData2, xval = 'nucl_diversity', zval = 'genome_species',xAngle = 30,
          yval = 'pnps', xLab = expression(italic(pi)),  yLab = "pN/pS",facet.Nrow = 3,aspect.Ratio = 1/1.5,
          gval = 'if_Lactobacillaceae', mycolor = colorLactobacillaceae)
  dev.off()
  cairo_pdf('./plot_mags/species_pi_vs_pnps_iners.pdf', width = 8, height = 3, onefile = T)
  for(sp in unique(refMetData2$genome_species)){
    print(sp)
    relPlot(refMetData2[refMetData2$genome_species == sp,], 
            xval = 'nucl_diversity', xAngle = 30,Title = sp,
            yval = 'pnps', xLab = expression(italic(pi)),  yLab = "pN/pS",facet.Nrow = 3,aspect.Ratio = 1,
            gval = 'if_Lactobacillaceae', mycolor = colorLactobacillaceae, showSE = F)
    relPlot(refMetData2[refMetData2$genome_species == sp,], 
            xval = 'nucl_diversity', xAngle = 30,Title = sp,
            yval = 'pnps', xLab = expression(italic(pi)),  yLab = "pN/pS",facet.Nrow = 3,aspect.Ratio = 1,
            gval = 'Ethnicity2', mycolor = colorEthnicity2, showSE = F)
  }
  dev.off()
  #### one vs one facet by ethnicity ####
  if(T){
    cairo_pdf('./plot_mags/evolutionary_metrics_one_vs_one_facetBy_Ethnicity.pdf', width = 11.48, height = 2.52, onefile = T)
    relPlot(refMetData2, xval = 'SNV_per_kbp', yval = 'd_prime_mean', zval = 'Ethnicity2', xLab = 'SNVs/kbp', yLab = "D'", 
            gval = 'if_Lactobacillaceae', mycolor = colorLactobacillaceae)
    relPlot(refMetData2, xval = 'SNV_per_kbp', yval = 'D_rev', zval = 'Ethnicity2', xLab = 'SNVs/kbp', yLab = "1 - D'", 
            gval = 'if_Lactobacillaceae', mycolor = colorLactobacillaceae)
    relPlot(refMetData2, xval = 'nucl_diversity', yval = 'd_prime_mean', zval = 'Ethnicity2',xLab = expression(italic(pi)), yLab = "D'", 
            gval = 'if_Lactobacillaceae', mycolor = colorLactobacillaceae, xAngle = 45)
    relPlot(refMetData2, xval = 'nucl_diversity', yval = 'D_rev', zval = 'Ethnicity2',xLab = expression(italic(pi)), yLab = "1 - D'", 
            gval = 'if_Lactobacillaceae', mycolor = colorLactobacillaceae, xAngle = 45)
    relPlot(refMetData2, xval = 'coverage', yval = 'SNV_per_kbp', zval = 'Ethnicity2', xLab = 'Coverage', yLab = 'SNVs/kbp', 
            gval = 'if_Lactobacillaceae', mycolor = colorLactobacillaceae, xAngle = 45)
    relPlot(refMetData2, xval = 'coverage', yval = 'nucl_diversity', zval = 'Ethnicity2', xLab = 'Coverage', yLab =  expression(italic(pi)), 
            gval = 'if_Lactobacillaceae', mycolor = colorLactobacillaceae, xAngle = 45)
    relPlot(refMetData2, xval = 'coverage', yval = 'd_prime_mean', zval = 'Ethnicity2', xLab = 'Coverage', yLab = "D'", 
            gval = 'if_Lactobacillaceae', mycolor = colorLactobacillaceae, xAngle = 45)
    relPlot(refMetData2, xval = 'coverage', yval = 'D_rev', zval = 'Ethnicity2', xLab = 'Coverage', yLab = "1 - D'", 
            gval = 'if_Lactobacillaceae', mycolor = colorLactobacillaceae, xAngle = 45)
    relPlot(refMetData2, xval = 'SNV_per_kbp', yval = 'pnps', zval = 'Ethnicity2', 'SNVs/kbp', yLab = "pN/pS",
            gval = 'if_Lactobacillaceae', mycolor = colorLactobacillaceae)
    relPlot(refMetData2, xval = 'nucl_diversity', yval = 'pnps', zval = 'Ethnicity2',xLab = expression(italic(pi)), yLab = "pN/pS",
            gval = 'if_Lactobacillaceae', mycolor = colorLactobacillaceae, xAngle = 45)
    relPlot(refMetData2, xval = 'coverage', yval = 'pnps', zval = 'Ethnicity2',xLab = 'Coverage', yLab = "pN/pS",
            gval = 'if_Lactobacillaceae', mycolor =colorLactobacillaceae, xAngle = 45)
    dev.off()
  }
}
cairo_pdf(file = './plot_mags/evolutionary_metrics_one_vs_one_facetBy_Ethnicity_pi_vs_pnps.pdf', width = 11.51, height = 4)
relPlot(refMetData2, xval = 'nucl_diversity', yval = 'pnps', zval = 'Ethnicity2',xLab = expression(italic(pi)), yLab = "pN/pS", xAngle = 30, aspect.Ratio = 1/1.5,
        gval = 'if_Lactobacillaceae', mycolor = colorLactobacillaceae)
dev.off()

###### correlation among evolutionary metrics by each species ######
if(T){
  refMetData2$genome_species %<>% as.factor()
  # SNVs.D.DF = corFitFun(df = refMetData2, xval = 'SNV_per_kbp', yval = 'D_rev', gval = 'genome_species',
  #                      xLab = 'SNVs/kbp', yLab = "1 - D'")
  # SNVs.pNpS.DF = corFitFun(refMetData2, xval = 'SNV_per_kbp', yval = 'pnps', gval = 'genome_species',
  #                         xLab = 'SNVs/kbp', yLab = "pN/pS")
  # Coverage.SNVs.DF = corFitFun(refMetData2, xval = 'coverage', yval = 'SNV_per_kbp', gval = 'genome_species',
  #                             xLab = 'Coverage', yLab = "SNVs/kbp")
  # Coverage.D.DF = corFitFun(refMetData2, xval = 'coverage', yval = 'D_rev', gval = 'genome_species',
  #                          xLab = 'Coverage', yLab = "1 - D'")
  SNVs.D.DF = corFitFun(df = refMetData2, xval = 'nucl_diversity', yval = 'D_rev', gval = 'genome_species',
                        xLab = "π", yLab = "1 - D'")
  SNVs.pNpS.DF = corFitFun(refMetData2, xval = 'nucl_diversity', yval = 'pnps', gval = 'genome_species',
                           xLab = 'π', yLab = "pN/pS")
  Coverage.SNVs.DF = corFitFun(refMetData2, xval = 'coverage', yval = 'nucl_diversity', gval = 'genome_species',
                               xLab = 'Coverage', yLab = 'π')
  Coverage.D.DF = corFitFun(refMetData2, xval = 'coverage', yval = 'D_rev', gval = 'genome_species',
                            xLab = 'Coverage', yLab = "1 - D'")
  compDF = rbind(SNVs.D.DF, SNVs.pNpS.DF, Coverage.SNVs.DF, Coverage.D.DF)
  compDF$comp %<>% factor(., levels = c("π vs 1 - D'", "π vs pN/pS", "Coverage vs π", "Coverage vs 1 - D'"))
  compDF$sig = 'Non-Significant'
  compDF$sig[compDF$PearsonP< 0.05] = 'Significant'
  compDF$sig %<>% factor(., levels = c('Significant', 'Non-Significant'))
  compDF$change[compDF$PearsonR > 0] = 'Positive'
  compDF$change[compDF$PearsonR < 0] = 'Negative'
  compDF$Group %<>% factor(., levels = unique(SNVs.D.DF$Group[order(SNVs.D.DF$PearsonR, decreasing = F)]))
  compDF$if_Lactobacillaceae = refMetData2$if_Lactobacillaceae[match(compDF$Group, refMetData2$genome_species)]
  myP = ggplot(compDF, aes(x=Group, y=PearsonR)) +
    geom_segment(aes(y=Group, yend=Group, x=0, xend=PearsonR), color="gray", size = 1.8) +
    geom_point(aes(x = PearsonR, y = Group, shape = sig, fill = sig, color = change), stroke = 1, size = 3.5) +
    facet_grid(if_Lactobacillaceae~comp, scales = 'free', space = 'free_y') +
    # scale_size_continuous(name = expression('Adj.' ~ italic(R)^2), range = c(0.5,5)) +
    scale_shape_manual(name = NULL, values = c(19, 21)) +
    scale_color_manual(name = 'Correlation', values = c(Positive = '#D99F16',Negative = '#2597A2')) +
    scale_fill_manual(name = '', values = c(NA, 'white')) +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'gray') +
    scale_x_continuous(expand = expansion(mult = c(0.08,0.08))) +
    # labs(x = "Beta coefficient", y = 'Species') +
    labs(x = "Pearson's r", y = 'Species') +
    guides(fill = 'none',
           shape = guide_legend(order = 3),
           color = guide_legend(order = 1),
           size = guide_legend(order = 2)) +
    geom_vline(xintercept = 0, linetype = 'dashed', color = 'gray', size = 0.3) +
    mytheme3 + theme(strip.text = element_text(size = 13), strip.background = element_rect(fill = '#E2E5E9', color = NA),
                     # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                     legend.position = 'right')
  cairo_pdf('./plot_mags/relationship among evolutionary metrics by each species_2.pdf', width = 13.51, height = 5.17)
  print(myP)
  dev.off()
}
######## ggpair ######
if(F){
  cor.test(refMetData$coverage %>% log10, refMetData$nucl_diversity)
  cor.test(refMetData$coverage %>% log10, refMetData$dn_ds_common)
  ggplot(refMetData, aes(x = coverage, y = nucl_diversity)) +
    geom_point(color = '#2271B5', fill = '#9ECAE1', shape = 21, alpha = 0.9, size = 2) +
    geom_smooth(method = 'lm') +
    labs(x = 'Coverage', y = 'Nucleotide diversiry') +
    mytheme3
  ggplot(refMetData, aes(x = coverage, y = dn_ds_common)) +
    geom_point(color = '#2271B5', fill = '#9ECAE1', shape = 21, alpha = 0.9, size = 2) +
    geom_smooth(method = 'lm') +
    labs(x = 'Coverage', y = 'dN/dS') +
    mytheme3
  
  refMetData_sub = refMetData[, c('nucl_diversity', 'pnps', 'dnds','dn_ds_common',  'coverage', 'breadth',
                                  'breadth_minCov')]
  colnames(refMetData_sub) = c('pi', 'pN/pS', 'dN/dS', 'dN/dS (common)', 'Coverage', 'Breath', 'minCov breath')
  refMetData_sub %<>% add_column(`log10(coverage)` = log10(refMetData_sub$Coverage), .after = 'Coverage')
  cairo_pdf('./plot_mags/corr_of_diversity_vs_coverage_mags.pdf', height = 12, width = 12)
  ggpairs(refMetData_sub,
          diag =  list(continuous = wrap("barDiag", fill="#6BAED6")),
          lower = list(continuous = my_fn1),
          upper = list(continuous = my_fn2_pearson)) +
    theme(panel.grid = element_blank(), strip.text = element_text(colour = 'black', size = 12),
          axis.text = element_text(size = 10),
          axis.text.x = element_text(size = 10, colour = 'black', angle = 45, hjust = 1))
  dev.off()
  
}
###### comparsion of Lactobacillus vs Non-Lactobacillus #####
if(T){
  ys = c('genome_species', 'if_Lactobacillaceae', 'Ethnicity2', 
         'SNV_per_kbp', 'D_rev', 'nucl_diversity', 'pnps')
  useData = refMetData2[, ys]
  useData2 = melt(useData, id.vars = c('genome_species', 'if_Lactobacillaceae', 'Ethnicity2'))
  colnames(useData2) = c('Species', 'Group', 'Ethnicity', 'Metrics', 'Value')
  useData2$Metrics %<>% as.character()
  useData2$Metrics[useData2$Metrics == 'nucl_diversity'] = 'Pi'
  useData2$Metrics[useData2$Metrics == 'SNV_per_kbp'] = 'SNVs/kbp'
  useData2$Metrics[useData2$Metrics == 'dnds'] = 'dN/dS'
  useData2$Metrics[useData2$Metrics == 'pnps'] = 'pN/pS'
  useData2$Metrics[useData2$Metrics == 'D_rev'] = "1 - D'"
  useData2$Metrics %<>% factor(., levels = c('SNVs/kbp', 'Pi', "1 - D'", 'pN/pS'))
}
if(T){
  cairo_pdf('./plot_mags/evolution_metrics_nonLac_vs_Lac.pdf', width = 10.47, height = 1.78)
  stat.test <- useData2 %>% group_by(Metrics) %>% 
    wilcox_test(Value ~ Group) %>% 
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj") %>%
    add_xy_position(x = 'Group', dodge = 0, fun = 'max')
  stat.test$y.position = c(52, 0.025, 0.42, 1.9)
  stat.test$p.adj %<>% signif(., digits = 3)
  p = ggboxplot(useData2, x = 'Group', y = 'Value', fill = 'Group')
  pp = p %>% facet(nrow = 1, facet.by = 'Metrics', scales = 'free_x') +
    stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.01,coord.flip = T) +
    mytheme + theme(strip.text = element_text(size = 13)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.13))) + 
    scale_x_discrete(limits=rev)+
    # labs(title = paste0(str_replace(sp, '.fna', ''), ' | ', x), 
    #      color = x, y = ysNorm[[y]], x = x) + 
    scale_fill_manual(values = colorLactobacillaceae) +
    coord_flip()
  print(pp)
  dev.off()
}
if(T){
  cairo_pdf('./plot_mags/evolution_metrics_nonLac_vs_Lac_facetBy_Ethnicity.pdf', width = 10.85, height = 4.69)
  stat.test2 <- useData2 %>% group_by(Metrics, Ethnicity) %>% 
    wilcox_test(Value ~ Group) %>% 
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj") %>%
    add_xy_position(x = 'Group', dodge = 0)
  stat.test2$y.position = c(rep(52, 5),
                            rep(0.025, 5), 
                            rep(0.42, 5), 
                            rep(1.85, 5))
  stat.test2$p.adj %<>% signif(., digits = 3)
  p = ggboxplot(useData2, x = 'Group', y = 'Value', 
                # add = "jitter", add.params = list(shape = 1, size = 0.5),
                fill = 'Group')
  pp = p %>% facet(nrow = 1, facet.by = c('Ethnicity', 'Metrics'), scales = 'free') +
    stat_pvalue_manual(stat.test2, label = "p.adj", tip.length = 0.01,coord.flip = T) +
    mytheme + theme(strip.text = element_text(size = 13)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.13))) + 
    scale_x_discrete(limits=rev)+
    # labs(title = paste0(str_replace(sp, '.fna', ''), ' | ', x), 
    #      color = x, y = ysNorm[[y]], x = x) + 
    scale_fill_manual(values = colorLactobacillaceae) +
    coord_flip()
  print(pp)
  dev.off()
}

if(T){
  xs = c(#'Insulin', 'B6',
    'abu',
    'Ethnicity2','Employment', 'Housing',
    'Term_char', # 'Abortions',
    'BMI', 'Marriage','Age',  'trimester',
    'FOB', 'Depression')  # interest priority
  ys = c('nucl_diversity', 'SNV_per_kbp', 'D_rev', 'pnps')
  lmeRes = data.frame()
  for(y in ys){
    for(gp in c('Lactobacillaceae', 'Non-Lactobacillaceae')){
      formula_char = paste0(y, ' ~ (Intercept) + ', str_c(capitalize(xs), collapse = ' + '), ' + (1|Subject)')
      dat = refMetData2[refMetData2$if_Lactobacillaceae == gp, c(xs, 'pt_ID', y)] %>% na.omit
      lmeTry = try({lme(as.formula(paste0(y, ' ~ ', str_c(xs, collapse = '+'))), random = ~1|pt_ID,
                        data = dat, method = 'ML')})
      if(class(lmeTry) != 'try-error'){
        lme.Multi = lmeTry
        lme.Multi.R2 = r.squaredGLMM(lme.Multi)[1]
        lme.Multi.aov = anova(lme.Multi)
        vv_aov_df = data.frame(Group = gp,
                               y = y,
                               x = row.names(lme.Multi.aov), 
                               lme.multi.Model = formula_char,
                               lme.multi.F_value = lme.Multi.aov$`F-value`,
                               lme.multi.P_value = lme.Multi.aov$`p-value`,
                               lme.multi.R2 = lme.Multi.R2)
        lmeRes = rbind(lmeRes, vv_aov_df)
      }
    }
  }
  
  
  
  cairo_pdf('./plot_mags/evolution_metrics_Ethnicity_By_nonLac_vs_Lac_2.pdf', width = 5.09, height = 3.46, onefile = T)
  relAbuPlot(data0 = refMetData2[refMetData2$if_Lactobacillaceae == 'Lactobacillaceae',], 
             xval = 'Ethnicity2', yLab = expression(pi),titleLab = 'Lactobacillaceae',addJitter = F,
             yval = 'nucl_diversity', if_log10 = F, hide.ns = T, jitterSize = 1.2,
             y.position = seq(0.013, 0.013+0.003*7, 0.003),
             expandMult = expansion(mult = c(0.03, 0.1)))
  relAbuPlot(data0 = refMetData2[refMetData2$if_Lactobacillaceae == 'Non-Lactobacillaceae',], 
             xval = 'Ethnicity2', yLab = expression(pi),titleLab = 'Non-Lactobacillaceae',
             yval = 'nucl_diversity', if_log10 = F, hide.ns = T, addJitter = F,
             y.position = seq(0.026, 0.026 + 0.007*9, 0.007),
             expandMult = expansion(mult = c(0.03, 0.1)))
  relAbuPlot(data0 = refMetData2[refMetData2$if_Lactobacillaceae == 'Lactobacillaceae',], 
             xval = 'Ethnicity2', yLab = 'SNVs/kbp',titleLab = 'Lactobacillaceae',addJitter = F,
             yval = 'SNV_per_kbp', if_log10 = F, hide.ns = T, jitterSize = 1.2,
             y.position = seq(22, 22+3*2, 2),
             expandMult = expansion(mult = c(0.03, 0.1)))
  relAbuPlot(data0 = refMetData2[refMetData2$if_Lactobacillaceae == 'Non-Lactobacillaceae',], 
             xval = 'Ethnicity2', yLab = 'SNVs/kbp',titleLab = 'Non-Lactobacillaceae',addJitter = F,
             yval = 'SNV_per_kbp',  if_log10 = F, hide.ns = T, jitterSize = 1.2,
             y.position = seq(55, 55+12*7, 12),
             expandMult = expansion(mult = c(0.03, 0.1)))
  relAbuPlot(data0 = refMetData2[refMetData2$if_Lactobacillaceae == 'Lactobacillaceae',], 
             xval = 'Ethnicity2', yLab = "1 - D'",titleLab = 'Lactobacillaceae',addJitter = F,
             yval = 'D_rev', if_log10 = F, hide.ns = T, jitterSize = 1.2,
             y.position = seq(0.4, 0.4 + 5*0.08, 0.08)+0.02,
             expandMult = expansion(mult = c(0.03, 0.1)))
  relAbuPlot(data0 = refMetData2[refMetData2$if_Lactobacillaceae == 'Non-Lactobacillaceae',], 
             xval = 'Ethnicity2', yLab = "1 - D'",titleLab = 'Non-Lactobacillaceae',addJitter = F,
             yval = 'D_rev',if_log10 = F, hide.ns = T, jitterSize = 1.2,
             y.position = seq(0.22, 0.22+0.04*5, 0.04),
             expandMult = expansion(mult = c(0.03, 0.1)))
  relAbuPlot(data0 = refMetData2[refMetData2$if_Lactobacillaceae == 'Lactobacillaceae',], 
             xval = 'Ethnicity2', yLab = "pN/pS",titleLab = 'Lactobacillaceae',addJitter = F,
             yval = 'pnps',if_log10 = F, hide.ns = T, jitterSize = 1.2,
             y.position = seq(2, 2+0.25*3, 0.25),
             expandMult = expansion(mult = c(0.03, 0.1)))
  relAbuPlot(data0 = refMetData2[refMetData2$if_Lactobacillaceae == 'Non-Lactobacillaceae',], 
             xval = 'Ethnicity2', yLab = "pN/pS",titleLab = 'Non-Lactobacillaceae',addJitter = F,
             yval = 'pnps', if_log10 = F, hide.ns = T, jitterSize = 1.2,
             y.position = seq(1.5, 1.5+0.3*6, 0.3),
             expandMult = expansion(mult = c(0.03, 0.1)))
  dev.off()
}
###### comparsion of anaerobic or aerobic #####
if(T){
  ys = c('genome_species', 'if_lac', 
         'SNV_per_kbp', 'nucl_diversity', 'pnps', 'dnds')
  ysNorm = list(SNV_per_kbp = 'SNVs per Kbp', nucl_diversity = 'Nucleotide diversity',
                dnds = 'dN/dS', pnps = 'pN/pS')
  useData = refMetData[, ys]
  useData$genome_species %>% table() %>% sort(., decreasing = T)
  table(useData$genome_species) %>% sort(., decreasing = T)
  openxlsx::write.xlsx(data.frame(
    species = names(table(useData$genome_species) %>% sort(., decreasing = T)),
    sampleNum = table(useData$genome_species) %>% sort(., decreasing = T) %>% as.numeric(),
    anaerobic_or_aerobic = NA
  ), file = './data/anaerobic or aerobic.xlsx',rowNames = F)
  useData1 = useData[useData$genome_species %in% 
                       names(table(useData$genome_species)[table(useData$genome_species) > 10]), ]
  table(useData1$genome_species)
  useData1$genome_species
  useData2 = melt(useData1, id.vars = c('genome_species', 'if_lac'))
  colnames(useData2) = c('Species', 'Group', 'Metrics', 'Value')
  useData2$Metrics %<>% as.character()
  useData2$Metrics[useData2$Metrics == 'nucl_diversity'] = 'Pi'
  useData2$Metrics[useData2$Metrics == 'SNV_per_kbp'] = 'SNVs per Kbp'
  useData2$Metrics[useData2$Metrics == 'dnds'] = 'dN/dS'
  useData2$Metrics[useData2$Metrics == 'pnps'] = 'pN/pS'
  useData2$Metrics %<>% factor(., levels = c('SNVs per Kbp', 'Pi', 'pN/pS', 'dN/dS'))
  stat.test <- useData2 %>% group_by(Metrics) %>% 
    wilcox_test(Value ~ Group) %>% 
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj") %>%
    add_xy_position(x = 'Group', dodge = 0, fun = 'max')
  stat.test$y.position = c(52, 0.025, 2.3, 1.25)
  stat.test$p.adj %<>% signif(., digits = 3)
  p = ggboxplot(useData2, x = 'Group', y = 'Value', fill = 'Group')
  pp = p %>% facet(nrow = 1, facet.by = 'Metrics', scales = 'free_x') +
    stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.01,coord.flip = T) +
    mytheme + theme(strip.text = element_text(size = 13)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) + 
    # labs(title = paste0(str_replace(sp, '.fna', ''), ' | ', x), 
    #      color = x, y = ysNorm[[y]], x = x) + 
    scale_fill_manual(values = c('#DFA800', '#22A2AE')) +
    coord_flip()
  cairo_pdf('./plot_mags/evolution_metrics_nonLac_vs_Lac.pdf', width = 9.78, height = 1.78)
  print(pp)
  dev.off()
}

plotFlag = T
#### single variable comparison lm & lme #####
if(T){
  # ys = c('nucl_diversity', 'dnds', 'pnps', 'dn_ds_common')
  ys = c('nucl_diversity', 'dnds', 'pnps', 'SNV_per_kbp', 'D_rev')
  # ysNorm = list(nucl_diversity = 'Nucleotide diversity', dnds = 'dN/dS', pnps = 'pN/pS', dn_ds_common = 'dN/dS (common)')
  ysNorm = list(nucl_diversity = 'Nucleotide diversity', dnds = 'dN/dS', pnps = 'pN/pS', SNV_per_kbp = 'SNVs/kbp', D_rev = "1 - D'")
  xs = c('abu','trimester', 'Age_char', 'BMI_char', 'Ethnicity', 'Ethnicity2', 'Language', 'Marriage', 'FOB', 'Employment', 'Housing', 
         'Term_char', 'Term_char2', 'Gravida_char', 'Abortions_char', 'Living_children_char', 
         'Abnormal_pap', 'Depression', 'PTSD', 'PPD_pos', 'DM2', 'Asthma', 'PCN_allergy', 'Sulfa_allergy', 'NKDA',
         'PNV', 'Fe', 'Insulin', 'B6', 'Progesterone',
         'IPV_hx', 'TED_hx', 'Tobacco_hx', 'EtOH_hx', 
         'GBS', 'PCN', 'Antibiotic',
         'Induction', 'Baby_gender', 'Baby_weight_char', 
         # repeated numeric 
         'Age', 'BMI', 'Term', 'Gravida', 'Abortions', 'Living_children', 'Baby_weight')
  comp_DF = data.frame()
  uniRes = data.frame()
  # sp = unique(refMetData$genome_species)[1]
  # x = xs[1]
  for (sp in unique(refMetData2$genome_species)) {
    cat('sp =', sp, '\n')
    tmpDF = refMetData2[refMetData2$genome_species == sp, ]
    for(y in ys){
      cat('  y =', y, '\n')
      for (x in xs) {
        cat('    x =', x, '\n')
        tmpDF2 = data.frame(
          'y' = tmpDF[, y],
          'x' = tmpDF[, x],
          trimester = tmpDF$trimester,
          pt_ID = tmpDF$pt_ID
        )
        tmpDF2 = na.omit(tmpDF2)
        tmpDF2 = tmpDF2[!is.infinite(tmpDF2$y), ]
        if(length(na.omit(unique(tmpDF2$x))) >= 2){
          # lm
          lmUni = lm(y ~ x, tmpDF2)
          lmUni.summary = summary(lmUni)
          lmUni.aov = anova(lmUni)
          # lmer
          formula_char = paste0(y, ' ~ (Intercept) + ', x, ' + (1|Subject)')
          lmeTry = try({lmer(as.formula(paste0('y ~ x + (1|pt_ID)')),
                             data = tmpDF2, REML = F)}, silent = F)
          if(class(lmeTry) != 'try-error'){
            lmeUni = lmeTry
            lmeUni.R2 = r.squaredGLMM(lmeUni)[1]
            lmeUni.aov = anova(lmeUni)
            sv = c(sp, y, x, (lmUni.aov['x', c('F value', 'Pr(>F)')]) %>% unlist(), lmUni.summary$r.squared %>% as.numeric(),
                   formula_char, lmeUni.aov$`F value`[1], lmeUni.aov$`Pr(>F)`[1], lmeUni.R2)
          }else{
            sv = c(sp, y, x, (lmUni.aov['x', c('F value', 'Pr(>F)')]) %>% unlist(), lmUni.summary$r.squared %>% as.numeric(),
                   formula_char, NA, NA, NA)
          }
          uniRes = rbind(uniRes, sv)
        }
        if(class(tmpDF2$x) != 'numeric' & 
           length(na.omit(unique(tmpDF2$x))) >= 2 & 
           nrow(tmpDF2) != length(unique(tmpDF2$x))){
          mm = aggregate(y ~ x, tmpDF2, mean)
          mmed = aggregate(y ~ x, tmpDF2, median)
          size = table(tmpDF2$x) %>% as.data.frame()
          # do pairwise comparison
          tmp_comp_DF = compare_means(y ~ x, tmpDF2, p.adjust.method = 'BH') %>% as.data.frame()
          tmp_comp_DF$y = y
          tmp_comp_DF$x = x
          tmp_comp_DF$species = sp
          tmp_comp_DF$size1 = size$Freq[match(tmp_comp_DF$group1, size$Var1)]
          tmp_comp_DF$size2 = size$Freq[match(tmp_comp_DF$group2, size$Var1)]
          
          tmp_comp_DF$mean1 = mm$y[match(tmp_comp_DF$group1, mm$x)]
          tmp_comp_DF$mean2 = mm$y[match(tmp_comp_DF$group2, mm$x)]
          
          tmp_comp_DF$log2_FC = log2(tmp_comp_DF$mean1/tmp_comp_DF$mean2)
          
          tmp_comp_DF$med1 = mmed$y[match(tmp_comp_DF$group1, mmed$x)]
          tmp_comp_DF$med2 = mmed$y[match(tmp_comp_DF$group2, mmed$x)]
          
          tmp_comp_DF = tmp_comp_DF[, c('species', 'y', 'x', 'group1', 'group2', 'size1', 'size2', 'mean1', 'mean2', 'log2_FC', 'med1', 'med2',
                                        'p', 'p.adj', 'p.format', 'p.signif')]
          comp_DF = rbind(comp_DF, tmp_comp_DF)
          
        }
      }
    }
  }
  colnames(uniRes) = c('species', 'y', 'x', 'lm.uni.F_value', 'lm.uni.P_value', 'lm.uni.R2',
                       'lme.uni.Model', 'lme.uni.F_value', 'lme.uni.P_value', 'lme.uni.R2')
  uniRes$lm.uni.F_value %<>% as.numeric()
  uniRes$lm.uni.P_value %<>% as.numeric()
  uniRes$lm.uni.R2 %<>% as.numeric()
  uniRes$lme.uni.F_value %<>% as.numeric()
  uniRes$lme.uni.P_value %<>% as.numeric()
  uniRes$lme.uni.R2 %<>% as.numeric()
  if(plotFlag){
    openxlsx::write.xlsx(comp_DF, file = './data_mags/microdiversity_pairwise_group_comp_2.xlsx')  
    openxlsx::write.xlsx(uniRes, file = './data_mags/microdiversity_uni_regression_v2.xlsx')
  }
}

##### lme regression fixed model #####
# nucl div associations to metadata variables
colnames(refMetData2)
table(refMetData2$species)
xs = c(#'Insulin', 'B6',
  'abu',
  'Ethnicity2','Employment', 'Housing',
  'Term_char', # 'Abortions',
  'BMI', 'Marriage','Age',  'trimester',
  'FOB', 'Depression')  # interest priority
# xs = c(#'Insulin', 'B6', 
#   'abu', 
#   'Employment', 'Housing',
#   'Depression',
#   'FOB', 'Marriage',
#   'Term_char', # 'Abortions', 
#   'Ethnicity2', 'trimester',
#   'BMI', # 'BMI_char', # 'trimester'
#   'Age')  # interest priority
# ys = c('nucl_diversity', 'dnds', 'pnps', 'dn_ds_common')
ys = c('nucl_diversity', 'dnds', 'pnps', 'SNV_per_kbp', 'D_rev')
x_random = c('pt_ID')

# sp = "Lactobacillus_vaginalis (NCBI)" 
# y = 'nucl_diversity'
lmeRes = data.frame()
# refMetData = refMetData[refMetData$coverage > 10, ]
for (sp in unique(refMetData2$genome_species)) {
  cat('  sp =', sp, '\n')
  for (y in ys) {
    cat('  y =', y, '\n')
    dat = refMetData2[refMetData2$genome_species == sp, c(xs, x_random, y)]
    dat = na.omit(dat)
    sind = which(sapply(dat, function(col) length(unique(col))) > 1)
    if(nrow(dat) > 3 & length(sind) > 1){
      for (c in 1:ncol(dat)){
        if(is.factor(dat[, c])) {
          dat[,c] = droplevels(dat[,c])
        }else if(is.character(dat[,c])){
          dat[,c] = factor(dat[,c])
        }
      }
      valid_x = colnames(dat)[!colnames(dat) %in% c(ys, x_random)]
      formula_char = paste0(y, ' ~ (Intercept) + ', str_c(capitalize(valid_x), collapse = ' + '), ' + (1|Subject)')
      
      lmeTry = try({lme(as.formula(paste0(y, ' ~ ', str_c(valid_x, collapse = '+'))), random = ~1|pt_ID,
                        data = dat, method = 'ML')}, silent = F)
      if(class(lmeTry) != 'try-error'){
        lme.Multi = lmeTry
        lme.Multi.R2 = r.squaredGLMM(lme.Multi)[1]
        lme.Multi.aov = anova(lme.Multi)
        vv_aov_df = data.frame(species = sp,
                               y = y,
                               x = row.names(lme.Multi.aov), 
                               lme.multi.Model = formula_char,
                               lme.multi.F_value = lme.Multi.aov$`F-value`,
                               lme.multi.P_value = lme.Multi.aov$`p-value`,
                               lme.multi.R2 = lme.Multi.R2)
        lmeRes = rbind(lmeRes, vv_aov_df)
      }
    }
  }
}

lmeRes = lmeRes[lmeRes$x != '(Intercept)', ]
table(lmeRes$x)
table(lmeRes$y)
table(lmeRes$species)
# merge
uniRes = readxl::read_xlsx('./data_mags/microdiversity_uni_regression_v2.xlsx')
res = merge(uniRes, lmeRes, by = c('species', 'y', 'x'), all = T)
# openxlsx::write.xlsx(res, './data/microdiveristy_lme_regression_fixed_v1.xlsx')
res_mulisig = res[which(res$lme.multi.P_value < 0.05), ]
nrow(res_mulisig)
# res_mulisig = res[which(res$lme.multi.P_value < 0.05 & res$lme.uni.P_value < 0.05),]
# nrow(res_mulisig)
res_mulisig = arrange(res_mulisig, species, y, x, lme.multi.P_value)
openxlsx::write.xlsx(res_mulisig, './data_mags/microdiveristy_lme_regression_fixed_sig_v2.xlsx')

##### plot boxplot for candidate species and variable all references and mags #####
if(plotFlag){
  pdf('./plot_mags/microdiveristy_boxplot_all_cov10.pdf', width = 10, height = 7)
}
x = 'Ethnicity2'
y = 'dn_ds_common' # 'nucl_diversity'
for (sp in unique(refMetData$genome_species)) {
  tmpDF = refMetData[refMetData$genome_species == sp, ]
  tmpDF2 = data.frame(
    'fq_ID' = tmpDF$sample,
    'pt_ID' = tmpDF$pt_ID,
    'pt_ID.u' = tmpDF$pt_ID.u,
    'genome_species' = tmpDF$genome_species,
    y = tmpDF[, y],
    x = tmpDF[, x]
  )
  #pairs I want to compare in list format for stat_compare_means
  if(plotFlag){
    if(class(tmpDF2$x) == 'factor' & length(droplevels(unique(tmpDF2$x))) > 1){
      #all pairs I want to compare
      if(is.factor(tmpDF2$x)){
        CN <- combn(na.omit(droplevels(unique(tmpDF2$x))) %>% levels, 2, simplify = FALSE)
      }else{
        CN <- combn(na.omit(unique(tmpDF2$x)), 2, simplify = FALSE)
      }
      p = ggboxplot(na.omit(tmpDF2), x = 'x', y = 'y', color = 'x', add = 'jitter',
                    add.paramsshape = 'pt_ID') +
        mytheme3 + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        stat_compare_means(comparisons=CN, method = "wilcox.test", p.adjust.method = "BH", aes(label=..p.adj..)) +
        labs(title = paste0(str_replace(sp, '.fna', ''), ' | ', x), 
             color = x,
             y = ysNorm[[y]], x = x) + 
        scale_color_brewer(palette = 'Dark2')
      print(p)
    }else{
      if(nrow(na.omit(tmpDF2)) > 0){
        p2 = ggplot(na.omit(tmpDF2), aes(x = x, y = y)) +
          geom_smooth(method = 'lm', color = '#969696', size = 0.5, fill = "#D5D5D5") +
          geom_point(aes(color = 'Blue', fill = 'Blue'), shape = 21, size = 2, alpha = 0.5) +
          labs(x = x, y = ysNorm[[y]], title = paste0(str_replace(sp, '_', ' '), ' | ', x)) +
          stat_poly_eq(formula = y ~ x,
                       aes(group=1, label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), 
                       parse = TRUE) +
          mytheme + theme(strip.text = element_text(size = 12), legend.position = 'none')
        print(p2)
      }
    }
  }
}
dev.off()

##### plot boxplot for candidate aka significant species and variable #####
cand_sp_var = res_mulisig[, c('species', 'y', 'x')]
nrow(cand_sp_var)
plotFlag = T
ysNorm = list(nucl_diversity = 'Nucleotide diversity', dnds = 'dN/dS', pnps = 'pN/pS', dn_ds_common = 'dN/dS (common)')

if(plotFlag){
  pdf('./plot/microdiveristy_super_sig_2_v2_cov10.pdf', width = 10, height = 7)
}
for (pr in 1:nrow(cand_sp_var)) {
  sp = cand_sp_var[pr, ][1] %>% unlist
  y = cand_sp_var[pr, ][2] %>% unlist
  x = cand_sp_var[pr, ][3] %>% unlist
  tmpDF = refMetData[refMetData$genome_species == sp, ]
  tmpDF2 = data.frame(
    'fq_ID' = tmpDF$sample,
    'pt_ID' = tmpDF$pt_ID,
    'pt_ID.u' = tmpDF$pt_ID.u,
    'genome_species' = tmpDF$genome_species,
    y = tmpDF[, y],
    x = tmpDF[, x]
  )
  #pairs I want to compare in list format for stat_compare_means
  if(plotFlag){
    if(class(tmpDF2$x) == 'factor'){
      #all pairs I want to compare
      if(is.factor(tmpDF2$x)){
        CN <- combn(na.omit(droplevels(unique(tmpDF2$x))) %>% levels, 2, simplify = FALSE)
      }else{
        CN <- combn(na.omit(unique(tmpDF2$x)), 2, simplify = FALSE)
      }
      p = ggboxplot(na.omit(tmpDF2), x = 'x', y = 'y', color = 'x', add = 'jitter',
                    add.paramsshape = 'pt_ID') +
        mytheme3 + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        stat_compare_means(comparisons=CN, method = "wilcox.test", p.adjust.method = "BH", aes(label=..p.adj..)) +
        labs(title = paste0(str_replace(sp, '.fna', ''), ' | ', x), 
             color = x,
             y = ysNorm[[y]], x = x) + 
        scale_color_brewer(palette = 'Dark2')
      print(p)
    }else{
      p2 = ggplot(na.omit(tmpDF2), aes(x = x, y = y)) +
        geom_smooth(method = 'lm', color = '#969696', size = 0.5, fill = "#D5D5D5") +
        geom_point(aes(color = 'Blue', fill = 'Blue'), shape = 21, size = 2, alpha = 0.5) +
        labs(x = x, y = ysNorm[[y]], title = paste0(str_replace(sp, '_', ' '), ' | ', x)) +
        stat_poly_eq(formula = y ~ x,
                     aes(group=1, label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), 
                     parse = TRUE) +
        mytheme + theme(strip.text = element_text(size = 12), legend.position = 'none')
      print(p2)
    }
  }
}
dev.off()

###### careful plot selected boxplot for paper ####

cairo_pdf('./plot/microdiversity_super_v3.pdf', width = 4.34,height = 3.55, onefile = T)
if(T){
  sp = 'Lactobacillus gasseri'
  if(T){
    tmpDF = refMetData[refMetData$genome_species == sp, ]
    tmpDF$Term = tmpDF$Term_char
    relAbuPlot(data0 = tmpDF, xval = 'Term', yLab = expression(pi),titleLab = sp,
               yval = 'nucl_diversity', x.angle = 45, if_log10 = F, hide.ns = F, 
               # y.position = c(0.14, 0.148),
               expandMult = expansion(mult = c(0.03, 0.1)))
  }
  sp = 'Lactobacillus iners' 
  if(T){
    tmpDF = refMetData[refMetData$genome_species == sp, ]
    relAbuPlot(data0 = tmpDF, xval = 'Ethnicity2', yLab = expression(pi),titleLab = sp,
               yval = 'nucl_diversity', x.angle = 45, if_log10 = F, hide.ns = T, 
               y.position = c(0.012, 0.014,0.016,0.018,0.020),
               expandMult = expansion(mult = c(0.03, 0.1)))
    relAbuPlot(data0 = tmpDF, xval = 'Ethnicity2', yLab = 'pN/pS',titleLab = sp,
               yval = 'pnps', x.angle = 45, if_log10 = F, hide.ns = T, 
               y.position = c(1.5,1.7,1.9,2.1,2.3,2.5),
               expandMult = expansion(mult = c(0.03, 0.1)))
    relAbuPlot(data0 = tmpDF, xval = 'FOB', yLab = 'pN/pS',titleLab = sp,
               yval = 'pnps', x.angle = 45, if_log10 = F, hide.ns = T, 
               # y.position = c(1.5,1.7,1.9,2.1,2.3,2.5),
               expandMult = expansion(mult = c(0.03, 0.1)))
    relAbuPlot(data0 = tmpDF, xval = 'Ethnicity2', yLab = 'dN/dS',titleLab = sp,
               yval = 'dnds', x.angle = 45, if_log10 = F, hide.ns = T, 
               y.position = c(1.5,1.7,1.9,2.1,2.3,2.5)/10,
               expandMult = expansion(mult = c(0.03, 0.1)))
    relAbuPlot(data0 = tmpDF, xval = 'Housing', yLab = 'dN/dS',titleLab = sp,
               yval = 'dnds', x.angle = 45, if_log10 = F, hide.ns = T, 
               # y.position = c(1.5,1.7,1.9,2.1,2.3,2.5)/10,
               expandMult = expansion(mult = c(0.03, 0.1)))
    relAbuPlot(data0 = tmpDF, xval = 'Ethnicity2', yLab = 'dN/dS (common)',titleLab = sp,
               yval = 'dn_ds_common', x.angle = 45, if_log10 = F, hide.ns = T, 
               y.position = c(1.5,1.65,1.8,1.95)/4 - 0.02,
               expandMult = expansion(mult = c(0.03, 0.1)))
  }
  sp = 'Fannyhessea vaginae'
  if(T){
    tmpDF = refMetData[refMetData$genome_species == sp, ]
    relAbuPlot(data0 = tmpDF, xval = 'Ethnicity2', yLab = expression(pi),titleLab = sp,
               yval = 'nucl_diversity', x.angle = 45, if_log10 = F, hide.ns = T, 
               y.position = c(0.02,0.025,0.03,0.035,0.04,0.045,0.05),
               expandMult = expansion(mult = c(0.03, 0.1)))
    relAbuPlot(data0 = tmpDF, xval = 'Ethnicity2', yLab = 'pN/pS',titleLab = sp,
               yval = 'pnps', x.angle = 45, if_log10 = F, hide.ns = T, 
               y.position = c(1.5,1.7,1.9,2.1,2.3,2.5)-0.7,
               expandMult = expansion(mult = c(0.03, 0.1)))
    relAbuPlot(data0 = tmpDF, xval = 'Ethnicity2', yLab = 'dN/dS',titleLab = sp,
               yval = 'dnds', x.angle = 45, if_log10 = F, hide.ns = T, 
               expandMult = expansion(mult = c(0.03, 0.1)))
    relAbuPlot(data0 = tmpDF, xval = 'Ethnicity2', yLab = 'dN/dS (common)',titleLab = sp,
               yval = 'dn_ds_common', x.angle = 45, if_log10 = F, hide.ns = T, 
               y.position = c(1.5,1.65,1.8)/10,
               expandMult = expansion(mult = c(0.03, 0.1)))
    tmpDF$Trimester = tmpDF$trimester
    relAbuPlot(data0 = tmpDF, xval = 'Trimester', yLab = 'pN/pS',titleLab = sp,
               yval = 'pnps', if_log10 = F, hide.ns = T, 
               y.position = c(0.75, 0.82),
               expandMult = expansion(mult = c(0.03, 0.1)))
  }
  sp = 'Gardnerella vaginalis'
  if(T){
    tmpDF = refMetData[refMetData$genome_species == sp, ]
    relAbuPlot(data0 = tmpDF, xval = 'Ethnicity2', yLab = expression(pi),titleLab = sp,
               yval = 'nucl_diversity', x.angle = 45, if_log10 = F, hide.ns = T, 
               y.position = c(0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.055)+0.003,
               expandMult = expansion(mult = c(0.03, 0.1)))
    relAbuPlot(data0 = tmpDF, xval = 'Ethnicity2', yLab = 'pN/pS',titleLab = sp,
               yval = 'pnps', x.angle = 45, if_log10 = F, hide.ns = T, 
               y.position = c(1.5,1.7,1.9,2.1,2.3,2.5)-0.3,
               expandMult = expansion(mult = c(0.03, 0.1)))
    relAbuPlot(data0 = tmpDF, xval = 'Ethnicity2', yLab = 'dN/dS',titleLab = sp,
               yval = 'dnds', x.angle = 45, if_log10 = F, hide.ns = T, 
               y.position = c(0.14,0.15,0.16,0.17)-0.01,
               expandMult = expansion(mult = c(0.03, 0.1)))
    relAbuPlot(data0 = tmpDF, xval = 'Ethnicity2', yLab = 'dN/dS (common)',titleLab = sp,
               yval = 'dn_ds_common', x.angle = 45, if_log10 = F, hide.ns = T, 
               y.position = c(0.24, 0.27),
               expandMult = expansion(mult = c(0.03, 0.1)))
  }
  sp = 'Bifidobacteriaceae bacterium NR047'
  if(T){
    tmpDF = refMetData[refMetData$genome_species == sp, ]
    relAbuPlot(data0 = tmpDF, xval = 'Ethnicity2', yLab = expression(pi),titleLab = sp,
               yval = 'nucl_diversity', x.angle = 45, if_log10 = F, hide.ns = T, 
               y.position = c(0.02,0.025,0.03,0.035,0.04,0.045)-0.003,
               expandMult = expansion(mult = c(0.03, 0.1)))
    relAbuPlot(data0 = tmpDF, xval = 'Ethnicity2', yLab = 'pN/pS',titleLab = sp,
               yval = 'pnps', x.angle = 45, if_log10 = F, hide.ns = T, 
               y.position = c(1.5,1.7,1.9,2.1),
               expandMult = expansion(mult = c(0.03, 0.1)))
    relAbuPlot(data0 = tmpDF, xval = 'Ethnicity2', yLab = 'dN/dS',titleLab = sp,
               yval = 'dnds', x.angle = 45, if_log10 = F, hide.ns = T, 
               y.position = c(0.14,0.15,0.16,0.17)-0.02,
               expandMult = expansion(mult = c(0.03, 0.1)))
    relAbuPlot(data0 = tmpDF, xval = 'Ethnicity2', yLab = 'dN/dS (common)',titleLab = sp,
               yval = 'dn_ds_common', x.angle = 45, if_log10 = F, hide.ns = T, 
               y.position = c(0.116),
               expandMult = expansion(mult = c(0.03, 0.1)))
  }
  sp = 'Lactobacillus crispatus'
  if(T){
    tmpDF = refMetData[refMetData$genome_species == sp, ]
    relAbuPlot(data0 = tmpDF, xval = 'Ethnicity2', yLab = expression(pi),titleLab = sp,
               yval = 'nucl_diversity', x.angle = 45, if_log10 = F, hide.ns = T, 
               y.position = c(0.02,0.025,0.03,0.035,0.04,0.045,0.05)-0.01,
               expandMult = expansion(mult = c(0.03, 0.1)))
    relAbuPlot(data0 = tmpDF, xval = 'Ethnicity2', yLab = 'pN/pS',titleLab = sp,
               yval = 'pnps', x.angle = 45, if_log10 = F, hide.ns = T, 
               y.position = c(1.5,1.7,1.9)-0.17,
               expandMult = expansion(mult = c(0.03, 0.1)))
    relAbuPlot(data0 = tmpDF, xval = 'Ethnicity2', yLab = 'dN/dS',titleLab = sp,
               yval = 'dnds', x.angle = 45, if_log10 = F, hide.ns = T, 
               y.position = c(0.13),
               expandMult = expansion(mult = c(0.03, 0.1)))
    relAbuPlot(data0 = tmpDF, xval = 'Ethnicity2', yLab = 'dN/dS (common)',titleLab = sp,
               yval = 'dn_ds_common', x.angle = 45, if_log10 = F, hide.ns = T, 
               y.position = c(0.4,0.44,0.48,0.52,0.56),
               expandMult = expansion(mult = c(0.03, 0.1)))
    relAbuPlot(data0 = tmpDF, xval = 'Employment', yLab = 'dN/dS (common)',titleLab = sp,
               yval = 'dn_ds_common', x.angle = 45, if_log10 = F, hide.ns = T, 
               y.position = c(0.4,0.44,0.48),
               expandMult = expansion(mult = c(0.03, 0.1)))
  }
  sp = 'Megasphaera genomosp type 1'
  if(T){
    tmpDF = refMetData[refMetData$genome_species == sp, ]
    relAbuPlot(data0 = tmpDF, xval = 'Employment', yLab = 'dN/dS (common)',titleLab = sp,
               yval = 'dn_ds_common', x.angle = 45, if_log10 = F, hide.ns = T, 
               y.position = c(0.44,0.48),
               expandMult = expansion(mult = c(0.03, 0.1)))
    relAbuPlot(data0 = tmpDF, xval = 'Employment', yLab = 'dN/dS',titleLab = sp,
               yval = 'dnds', x.angle = 45, if_log10 = F, hide.ns = T, 
               y.position = c(0.14, 0.148),
               expandMult = expansion(mult = c(0.03, 0.1)))
    relAbuPlot(data0 = tmpDF, xval = 'Employment', yLab = 'pN/pS',titleLab = sp,
               yval = 'pnps', x.angle = 45, if_log10 = F, hide.ns = T, 
               y.position = c(0.8, 0.9, 1)+0.05,
               expandMult = expansion(mult = c(0.03, 0.1)))
  }
}
dev.off()

##### carefully!! ####
cairo_pdf('./plot/microdiversity_super_good_results.pdf', width = 6.05,height = 3.55, onefile = T)
if(T){
  sp = 'Fannyhessea vaginae'
  tmpDF = refMetData[refMetData$genome_species == sp, ]
  tmpDF$Trimester = tmpDF$trimester
  relAbuPlot(data0 = tmpDF, xval = 'Trimester', yLab = 'pN/pS',titleLab = sp,
             yval = 'pnps', if_log10 = F, hide.ns = T, 
             y.position = c(0.75, 0.82),aspectRatio = 1.2/1,
             expandMult = expansion(mult = c(0.03, 0.1)))
  relAbuPlot(data0 = tmpDF, xval = 'Trimester', yLab = "1 - D'",titleLab = sp,
             yval = 'D_rev', if_log10 = F, hide.ns = T, 
             y.position = 0.25, 
             aspectRatio = 1.2/1,
             expandMult = expansion(mult = c(0.03, 0.1)))
  sp = 'Gardnerella vaginalis'
  tmpDF = refMetData[refMetData$genome_species == sp, ]
  tmpDF$Trimester = tmpDF$trimester
  relAbuPlot(data0 = tmpDF, xval = 'Trimester', yLab = "1 - D'",titleLab = sp,
             yval = 'D_rev', if_log10 = F, hide.ns = T, 
             y.position =  seq(0.16, 0.16+0.02*3, 0.02),  aspectRatio = 1.2/1,
             expandMult = expansion(mult = c(0.03, 0.1)))
  
  sp = 'Lactobacillus iners'
  tmpDF = refMetData[refMetData$genome_species == sp, ]
  relAbuPlot(data0 = tmpDF, xval = 'FOB', yLab = 'pN/pS',titleLab = sp,
             yval = 'pnps',if_log10 = F, hide.ns = T, aspectRatio = 1.2/1,
             y.position = 1.5,
             expandMult = expansion(mult = c(0.03, 0.1)))
  
  sp = 'Lactobacillus crispatus'
  tmpDF = refMetData[refMetData$genome_species == sp, ]
  relAbuPlot(data0 = tmpDF, xval = 'Ethnicity2', yLab = "Ethnicity",titleLab = sp,
             yval = 'D_rev', if_log10 = F, hide.ns = T, 
             y.position = seq(0.29, 0.29 + 0.04*5, 0.04), 
             aspectRatio = 1/1.2,
             expandMult = expansion(mult = c(0.03, 0.1)))
}
dev.off()

##### count associations res1 #####
# number of significant univariate association
dim(res)
table(res$y)
res_filter = res[res$y %in% c('nucl_diversity', 'pnps', 'D_rev') &
                   # res$species %in% res$species[which(res$lme.multi.P_value < 0.05)] &
                   res$x %in% xs[!(xs %in% c('abu'))], ] # dn_ds_common, pnps, dnds
nrow(res_filter)
# uniVarSigDF = res_filter[which(res_filter$lme.uni.P_value <= 0.05), ]
uniVarSigDF = res_filter[which(res_filter$lme.uni.P_value <= 0.05 | res_filter$lme.multi.P_value <= 0.05), ]
nrow(uniVarSigDF)
uniVarSigDF$x = str_replace_all(uniVarSigDF$x, '_char', '')
uniVarSigDF$x %>% table()
uniVarSigDF$x[uniVarSigDF$x == 'Ethnicity2'] = 'Ethnicity'
uniVarSigDF$x[uniVarSigDF$x == 'Term2'] = 'Term'
uniVarSigDF$x[uniVarSigDF$x == 'trimester'] = 'Trimester'
uniVarSigDF$lm.uni.P_value_log = -log10(uniVarSigDF$lm.uni.P_value)
uniVarSigDF$lme.uni.P_value_log = -log10(uniVarSigDF$lme.uni.P_value)
uniVarSigDF$lme.multi.P_value_log = -log10(uniVarSigDF$lme.multi.P_value)
uniqueVar = uniVarSigDF[, c('species', 'y', 'x', 'lm.uni.P_value_log', 'lme.uni.P_value_log',
                            'lme.multi.P_value_log', 'lme.multi.P_value', 'lme.uni.P_value')]
uniqueVar$y[uniqueVar$y == 'pnps'] = 'pN/pS'
uniqueVar$y[uniqueVar$y == 'dnds'] = 'dN/dS'
uniqueVar$y[uniqueVar$y == 'dn_ds_common'] = 'common dN/dS'
uniqueVar$y[uniqueVar$y == 'D_rev'] = "1 - D'"
uniqueVar$y[uniqueVar$y == 'SNV_per_kbp'] = "SNVs/kbp"
uniqueVar$y[uniqueVar$y == 'nucl_diversity'] = 'π'
table(uniqueVar$y)
table(uniqueVar$species) %>% length()
uniqueVar$sp_y = paste0(uniqueVar$species %>% str_sub(1,1), '.',
                        uniqueVar$species %>% str_replace('^[A-Z]\\w+', ''), 
                        ' (', uniqueVar$y, ')')
YsName = sort(table(uniqueVar$sp_y)) %>% names()
ord = c(sort(table(uniqueVar$x)) %>% names(),
        YsName[grep('(π)', YsName)],
        YsName[grep("(1 - D')", YsName)],
        # YsName[grep('(dN/dS)', YsName)],
        YsName[grep('(pN/pS)', YsName)])
grid.col = c(rep('#E5C494', length(YsName[grep('(π)', YsName)])),
             rep('#80B1D3', length(YsName[grep("(1 - D')", YsName)])),
             rep('#9C97F6', length(YsName[grep('(pN/pS)', YsName)])),
             rep('#CDC5CF', length(table(uniqueVar$x))))
names(grid.col) = c(YsName[grep('(π)', YsName)], 
                    YsName[grep("(1 - D')", YsName)],
                    YsName[grep('(pN/pS)', YsName)],
                    table(uniqueVar$x) %>% names())
uniqueVar = uniqueVar[order(uniqueVar$lme.multi.P_value, decreasing = T), ]
col_fun = data.frame(
  uniqueVar$sp_y,
  uniqueVar$x
)
col_fun$col = '#ECF0D9'
col_fun$col[uniqueVar$lme.multi.P_value < 0.05] = '#FB8072' # '#DECBE4' 
col_fun$col[uniqueVar$lme.multi.P_value < 0.05 & uniqueVar$lme.uni.P_value < 0.05] ='#BC80BD'  # '#FFA18F' 
table(col_fun$col)
library(circlize)
cairo_pdf('./plot_mags/circos_microdiveristy_v3.pdf', width = 6.25, height = 8.52)
{chordDiagram(uniqueVar[, c('x', 'sp_y')], order = ord, annotationTrack = "grid", link.sort = T,
              grid.col = grid.col, col = col_fun$col, transparency = 0.4,
              annotationTrackHeight = mm_h(2.3),
              preAllocateTracks = list(track.height = 0.3))
  circos.track(track.index = 1,
               panel.fun = function(x, y) {
                 circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                             facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.6)
               }, bg.border = NA) # here set bg.border to NA is important
}
dev.off()
###### longitudinal lme , Susan #####
# function
lmePlot = function(lmeData, yLab, if_logged = T, P.pos = c(10,2.5), expandMulti = NULL,
                   yPos = 2,groupV = NULL, addSE = T, linewith = 0.8, titleLab = NULL){
  if(is.null(groupV)){
    lmeFit <- lme(y ~ Sample_GA, random = ~ 1 + Sample_GA | pt_ID, data = lmeData, control = lmeControl(maxIter = 500, msMaxIter = 500))
    new.dat <- data.frame(Sample_GA = 10:40)
    new.dat$pred <- predict(lmeFit, newdata=new.dat,level=0)
    #create design matrix
    Designmat <- model.matrix(eval(eval(lmeFit$call$fixed)[-2]), new.dat[-ncol(new.dat)])
    #compute standard error for predictions
    predvar <- diag(Designmat %*% lmeFit$varFix %*% t(Designmat))
    new.dat$SE <- sqrt(predvar)
    new.dat$hi2SE <- new.dat$pred + 2*new.dat$SE
    new.dat$lo2SE <- new.dat$pred - 2*new.dat$SE
    # new.dat$hi2SE <- exp(new.dat$pred + 2*new.dat$SE)
    # new.dat$lo2SE <- exp(new.dat$pred - 2*new.dat$SE)
    # new.dat$pred <- exp(new.dat$pred)
    # Reverse log-transform for the vaginal community
    pval = summary(lmeFit)$tTable["Sample_GA", "p-value"] %>% signif(., 2)
    adoL_p = paste0("~~italic(P) == ", pval)
    if(if_logged){
      # lmeData$y %<>% exp()
    }
    lP = ggplot(new.dat, aes(x = Sample_GA, y = pred)) + 
      geom_line(color = "red", size = linewith) +
      geom_point(data = lmeData, aes(x = Sample_GA, y = y), size=1, shape = 21, color = '#8FB7D3', fill = '#8FB7D3', alpha = 0.5) +
      geom_ribbon(aes(ymin = lo2SE,ymax = hi2SE), alpha=0.2, fill="red") +
      labs(x = "Gestational Weeks", y = yLab, title = titleLab) +
      annotate("text", x = P.pos[1], y = P.pos[2], size=4, label = adoL_p, parse = T, hjust = 0) +
      mytheme + theme(aspect.ratio = 1)
  }else{
    pval = groupR = c()
    new.dat = data.frame()
    for(i in levels(lmeData[, groupV])){
      trRes = try({
        lmeFitTmp <- lme(y ~ Sample_GA, random = ~ 1 + Sample_GA | pt_ID, data = lmeData[lmeData[, groupV]==i, ], 
                         control = lmeControl(maxIter = 500, msMaxIter = 500))
      })
      if (class(trRes) == "try-error") {
        cat('! can not perform lme for', i)
      }else{
        new.dat_tmp <- data.frame(Sample_GA = 10:40)
        new.dat_tmp$pred <- predict(lmeFitTmp, newdata=new.dat_tmp,level=0)
        #create design matrix
        Designmat <- model.matrix(eval(eval(lmeFitTmp$call$fixed)[-2]), new.dat_tmp[-ncol(new.dat_tmp)])
        #compute standard error for predictions
        predvar <- diag(Designmat %*% lmeFitTmp$varFix %*% t(Designmat))
        new.dat_tmp$SE <- sqrt(predvar)
        new.dat_tmp$hi2SE <- new.dat_tmp$pred + 2*new.dat_tmp$SE
        new.dat_tmp$lo2SE <- new.dat_tmp$pred - 2*new.dat_tmp$SE
        # new.dat_tmp$hi2SE <- exp(new.dat_tmp$pred + 2*new.dat_tmp$SE)
        # new.dat_tmp$lo2SE <- exp(new.dat_tmp$pred - 2*new.dat_tmp$SE)
        # new.dat_tmp$pred <- exp(new.dat_tmp$pred)
        new.dat_tmp$Group = i
        # Reverse log-transform for the vaginal community
        pval = c(pval, summary(lmeFitTmp)$tTable["Sample_GA", "p-value"] %>% signif(., 2))
        new.dat = rbind(new.dat, new.dat_tmp)
        groupR = c(groupR, i)
        
      }
      
    }
    if(if_logged){
      # lmeData$y %<>% exp()
    }
    if(endsWith(groupV, '2')){lname = str_remove(groupV, '2$')}else{lname = groupV}
    adoL_p = paste0('"', groupR, '"', "~~italic(P) == ", pval)
    lmeData$groupV = lmeData[, groupV]
    lP = ggplot(new.dat, aes(x = Sample_GA, y = pred, fill = Group)) +
      geom_point(data = lmeData, aes(x = Sample_GA, y = y, color = groupV, fill = groupV), size=1, shape = 21, alpha = 0.5) +
      scale_color_manual(values = get(paste0('color', groupV)), name = lname) +
      scale_fill_manual(values = get(paste0('color', groupV)), name = lname)
    if(addSE){
      lP = lP + geom_ribbon(aes(ymin = lo2SE,ymax = hi2SE), alpha=0.2)
    }
    if(yPos == 'auto'){
      yPos = 0.1*length(adoL_p)
    }
    lP = lP + geom_line(aes(color = Group), size = linewith) +
      new_scale_fill() + scale_color_manual(values = get(paste0('color', groupV)), name = lname) +
      new_scale_color() + scale_fill_manual(values = get(paste0('color', groupV)), name = lname) +
      labs(x = "Gestational Weeks", y = yLab, title = titleLab) +
      annotate("text", x = P.pos[1], y = yPos, size=4, label = adoL_p, parse = T, hjust = 0) +
      mytheme + theme(aspect.ratio = 1)
  }
  if(!is.null(expandMulti)){
    lP = lP + scale_y_continuous(expand = expandMulti)
  }
  print(lP)
}
cairo_pdf('./plot/lme_microdiversity.pdf', width = 4.08, height = 3.01, onefile = T)
if(T){
  refMetData$Term2 = refMetData$Term_char2
  refMetData$Sample_GA = swabData$Sample_GA[match(refMetData$pt_ID.u, swabData$pt_ID.u)]
  ys = c('nucl_diversity', 'dnds', 'pnps', 'dn_ds_common')
  for(sp in unique(refMetData$genome_species)){
    for(yv in ys){
      lmeDF = refMetData[refMetData$genome_species == sp & refMetData$trimester != 'P', ]
      lmeDF$y = lmeDF[, yv]
      titleLab = sp
      if(yv == 'dnds'){
        yLab = 'dN/dS'
      }else if(yv == 'pnps'){
        yLab = 'pN/pS'
      }else if(yv == 'nucl_diversity'){
        yLab = expression(pi)
      }else if(yv == 'dn_ds_common'){
        yLab = 'dN/dS (common)'
      }
      nrow(lmeDF)
      try({
        lmePlot(lmeData = lmeDF, yLab = yLab, titleLab = titleLab, if_logged = F, P.pos = c(10, max(  lmeDF$y)))
      })
      try({
        lmePlot(lmeData = lmeDF, yLab = yLab, titleLab = titleLab, if_logged = F, P.pos = c(10, 1),
                groupV = 'Term2', yPos = c(max(  lmeDF$y)-0.08,max(  lmeDF$y)))
      })
      try({
        lmePlot(lmeData = lmeDF, yLab = yLab, titleLab = titleLab, if_logged = F, P.pos = c(10, 1),
                groupV = 'Ethnicity2',yPos ='auto')
      })
    }
  }
}
dev.off()
###### longitudinal - GAMM #####
library(voxel)
library(gamm4)

gammPlot = function(gammDF, groupCovs, x.Pos = 10, y.Pos = 0, titleLab = NULL, yLab = NULL){
  if(endsWith(groupCovs,'2')){lname = groupCovs %>% str_remove(., '2$')}else{lname = groupCovs}
  gammDF$g = gammDF[, groupCovs]
  gamm2 <- gamm4::gamm4(y ~ Term2 + Ethnicity2 + BMI +  s(Sample_GA, by = g), 
                        data=gammDF, random =  ~ (1 | pt_ID))
  pval = summary(gamm2$gam)$s.table[, 4] %>% signif(., 2)
  names(pval) %<>% str_remove(., 's\\(.*\\)\\:') %>% str_remove('^g')
  adoL_p = paste0('"', names(pval), '"',  "~~italic(P) == ", pval)
  gammP2 <- plotGAMM_2(gammFit <- gamm2, smooth.cov <- "Sample_GA", groupCovs = groupCovs,
                       plotCI <- T, rawOrFitted = "raw", grouping = "pt_ID", plotCI = F) +
    scale_color_manual(values = get(paste0('color', groupCovs)), name = lname) +
    scale_fill_manual(values = get(paste0('color', groupCovs)), name = lname) +
    labs(x = 'Gestational Weeks', y = yLab, title = titleLab) +
    annotate(x =x.Pos, y = y.Pos, label = adoL_p, geom = "text", size = 4.5, parse = TRUE, hjust = 0) +
    mytheme
  print(gammP2)
}

cairo_pdf('./plot/gamm_microdiveristy.pdf', width = 4.66, height = 3.26, onefile = T)
if(T){
  refMetData$Term2 = refMetData$Term_char2
  refMetData$Sample_GA = swabData$Sample_GA[match(refMetData$pt_ID.u, swabData$pt_ID.u)]
  ys = c('nucl_diversity', 'dnds', 'pnps', 'dn_ds_common')
  for(sp in unique(refMetData$genome_species)){
    for(yv in ys){
      gammDF = refMetData[refMetData$genome_species == sp & refMetData$trimester != 'P', ]
      gammDF$y = gammDF[, yv]
      titleLab = sp
      if(yv == 'dnds'){
        yLab = 'dN/dS'
      }else if(yv == 'pnps'){
        yLab = 'pN/pS'
      }else if(yv == 'nucl_diversity'){
        yLab = expression(pi)
      }else if(yv == 'dn_ds_common'){
        yLab = 'dN/dS (common)'
      }
      nrow(gammDF)
      try({
        gammPlot(gammDF = gammDF,groupCovs = 'Term2', x.Pos = 10, y.Pos = max(gammDF$y) - seq(1,2,1)*max(gammDF$y)/10, yLab = yLab,titleLab = titleLab)
        gammPlot(gammDF = gammDF,groupCovs = 'Ethnicity2', x.Pos = 10, y.Pos = max(gammDF$y) - seq(1,5,1)*max(gammDF$y)/10, yLab = yLab,titleLab = titleLab)
      })
    }
  }
}
dev.off()  

##### lme regression selection model #####
# nucl div associations to metadata variables
colnames(allMetData)

xs = c('BMI_char', 'Ethnicity', 'Language', 'Marriage', 'FOB', 'Employed', 'Housing', 
       'Term_char', 
       'Abnormal_pap', 'Depression', 'PTSD', 'PPD_pos', 'DM2', 'Asthma', 'PCN_allergy', 'Sulfa_allergy', 'NKDA',
       'PNV', 'Fe', 'Insulin', 'B6', 'Progesterone',
       'IPV_hx', 'TED_hx', 'Tobacco_hx', 'EtOH_hx', 
       'GBS', 'PCN', 'Antibiotic', 
       'Induction', 'Baby_gender',
       # repeated numeric 
       'Age', 'Gravida', 'Abortions', 'Living_children', 'Baby_weight')
ys = c('nucl_diversity', 'dnds', 'pnps', 'dn_ds_common')
x_random = c('pt_ID', 'trimester')

sp = "Lactobacillus_vaginalis (NCBI)" 
y = 'nucl_diversity'
lmeBestRes = data.frame()
for (sp in unique(allMetData$genome_species)) {
  cat('  sp =', sp, '\n')
  for (y in ys) {
    cat('  y =', y, '\n')
    dat = allMetData[allMetData$genome_species == sp, c(xs, x_random, y)]
    dat = dat[!is.na(dat[, y]),]
    sind = which(sapply(dat, function(col) length(unique(col) %>% na.omit)) > 1)
    if(nrow(dat) > 3 & length(sind) > 1 & sum(c('trimester', 'pt_ID', y) %in% names(sind)) == 3){
      dat = dat[, sind]
      for (c in 1:ncol(dat)){
        if(is.factor(dat[, c])) {
          dat[,c] = droplevels(dat[,c])
        }else if(is.character(dat[,c])){
          dat[,c] = factor(dat[,c])
        }
      }
      addFx = colnames(dat)[!colnames(dat) %in% c(ys, x_random)]
      lmeNull=try({lme(as.formula(paste0(y, ' ~ trimester')), 
                       random = ~1|pt_ID,
                       data = dat, method = "ML")})
      if(class(lmeNull) != 'try-error'){
        scopeList = paste0('~.+', paste(addFx, collapse = '+'))
        lmeBest = stepAIC(lmeNull, dir = "forward", scope = scopeList)
        lmeBest.R2 = r.squaredGLMM(lmeBest)[1]
        lmeBest.aov = anova(lmeBest)
        lmeBest.res = data.frame(species = sp,
                                 y = y, 
                                 x = row.names(lmeBest.aov), 
                                 lme.best.Model = paste0(y , ' ~ $(Intercept) + ', str_c(row.names(lmeBest.aov), collapse = ' + '),' + (1|pt_ID)'),
                                 lme.best.F_value = lmeBest.aov$`F-value`,
                                 lme.best.P_value = lmeBest.aov$`p-value`,
                                 lme.best.R2 = lmeBest.R2)
        lmeBestRes = rbind(lmeBestRes, lmeBest.res)
      }
    }
  }
}

lmeBestRes = lmeBestRes[(lmeBestRes$x != '(Intercept)') & (!is.na(lmeBestRes$lme.best.P_value)), ]
table(lmeBestRes$x)
table(lmeBestRes$y)
table(lmeBestRes$species)
# merge
res2 = merge(uniRes, lmeBestRes, by = c('species', 'y', 'x'), all = T)
res2_mulisig = res2[which(res2$lme.best.P_value < 0.05 & res2$lm.uni.P_value < 0.05),]
res2_mulisig = arrange(res2_mulisig, species, y, x, lme.best.P_value)
xlsx::write.xlsx(res2_mulisig, './Routput/nucl_div_dnds_lme_regression_best.xlsx')


##### count associations res2 #####
# number of significant univariate association
dim(res2)
res2_filter = res2[res2$y %in% c('nucl_diversity', 'dn_ds_common') &
                     res2$x %in% xs, ] # dn_ds_common, pnps, dnds
uniVarSigDF = res2_filter[which(res2_filter$lme.uni.P_value <= 0.05), ]
uniVarSigDF$x = str_replace_all(uniVarSigDF$x, '_char', '')
uniVarSigDF$lm.uni.P_value_log = -log10(uniVarSigDF$lm.uni.P_value)
uniVarSigDF$lme.uni.P_value_log = -log10(uniVarSigDF$lme.uni.P_value)
uniVarSigDF$lme.best.P_value_log = -log10(uniVarSigDF$lme.best.P_value)
uniqueVar = uniVarSigDF[, c('species', 'y', 'x', 'lm.uni.P_value_log', 'lme.uni.P_value_log',
                            'lme.best.P_value_log', 'lme.best.P_value')]
uniqueVar$sp_y = paste0(uniqueVar$species, ': ', uniqueVar$y)
ord = c(sort(table(uniqueVar$x)) %>% names(),
        sort(table(uniqueVar$sp_y)) %>% names())
grid.col = c(rep('#80B1D3', length(table(uniqueVar$sp_y))), rep('#E5C494', length(table(uniqueVar$x))))
names(grid.col) = c(table(uniqueVar$sp_y) %>% names(), table(uniqueVar$x) %>% names())
col_fun = data.frame(
  uniqueVar$sp_y,
  uniqueVar$x
)
col_fun$col = '#DBDBDB'
col_fun$col[uniqueVar$lme.best.P_value < 0.05] = 'red' 
table(col_fun$col)
library(circlize)
pdf('./plot/circos_association_genitic_diversity_best.pdf', width = 17, height = 15)
chordDiagram(uniqueVar[, c('x', 'sp_y')], order = ord, annotationTrack = "grid", 
             grid.col = grid.col, col = col_fun$col, transparency = 0.5,
             annotationTrackHeight = mm_h(1.5),
             preAllocateTracks = list(track.height = 0.3))
circos.track(track.index = 1,
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                           facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.6)
             }, bg.border = NA) # here set bg.border to NA is important
dev.off()

##### look at the intersect #####

df1 = res_filter[which(res_filter$lme.uni.P_value <= 0.05), ]
df1$id = paste(df1$species, df1$y, df1$x, sep = '|')
dim(df1)
df2 = res2_filter[which(res2_filter$lme.uni.P_value <= 0.05), ]
df2$id = paste(df2$species, df2$y, df2$x, sep = '|')
dim(df2)

intersect(df1$id, df2$id) %>% length()

df1 = df1[which(df1$lme.multi.P_value < 0.05),]
df2 = df2[which(df2$lme.best.P_value < 0.05),]


intersect(df1$id, df2$id) %>% length()