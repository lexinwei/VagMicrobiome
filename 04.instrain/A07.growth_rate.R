source('./script/library_function.R')
# load('./data/RData/refDF0.RData')
load('./data/RData/refDF_minCovBreath10.RData')
load('./data/RData/metaData.RData')
load('./data/RData/swabData.RData')
load('./data/RData/metaphlanDF.RData')
load('./data/RData/colorList.RData')
load('./data/RData/abuList.RData')
s_df = abuList[['s']]
sp2lineage = s_df[, 1:2]
refDF0 = refDF
refDF0$coverage %<>% as.numeric()
refDF = data.frame()
for(samp in unique(refDF0$sample)){
  subf = refDF0[refDF0$sample == samp, ]
  subf$abu = subf$coverage/sum(subf$coverage)
  refDF = rbind(refDF, subf)
}

baseDir = './irep/bPTR_tsv'
statList = list.files(baseDir)
ptrDF = data.frame()
for(s in statList){
  if(s == "Ming_nova_VS_2_16_bwa.bPTR.tsv"){
    next
  }
  statDF = read.table(paste0(baseDir, '/', s), header = F)
  s2 = str_extract(s, '[0-9]_[0-9][0-9]')
  colnames(statDF) = c('genome', 'ori', 'ter', 'ptr')
  statDF$genome %<>% str_remove(., '/share/home/jianglab/weixin/workspace/vag/irep/NCBI_reference/s__') %>% str_remove(., '\\.fasta') %>% str_replace_all(., '_', ' ')
  statDF$ptr[statDF$ptr == 'n/a'] = NA
  statDF$sample = s2
  ptrDF = rbind(ptrDF, statDF)
}
table(ptrDF$genome) %>% names()
ptrDF2 = ptrDF[!is.na(ptrDF$ptr), ]
ptrDF2$ori %<>% str_replace_all(., ',', '')
ptrDF2$ter %<>% str_replace_all(., ',', '')
ptrDF2$ptr %<>% as.numeric()
table(ptrDF2$genome) %>% length()
ptrDF2$SampleID.u = swabData$SampleID.u[match(ptrDF2$sample, swabData$seq_ID)]

#######################################################################################
#                              associated with some metadata                          #
#######################################################################################

if(T){
  refMetData = merge(refDF[!grepl('KIT', refDF$pt_ID), ], metaData, by = 'pt_ID', all.x = T)
  refMetData$Marriage %<>% as.character
  refMetData$Marriage[is.na(refMetData$Marriage)] = 'Unknown'
  refMetData$Marriage %<>% factor(., levels = c('Married', 'Unmarried', 'Unknown'))
  refMetData$Employment = refMetData$Employed
  dim(refMetData)
  colnames(refMetData)
  refMetData$coverage %<>% as.numeric()
  refMetData$breadth %<>% as.numeric()
  refMetData$genus = str_split_fixed(refMetData$genome, '_', 2)[, 1]
  table(refMetData$genus)
  
  refMetData$genome_species = refMetData$species %>% str_remove_all('^s__|\\.fasta$') %>% str_replace_all('_', ' ')
  table(refMetData$genome_species) %>% sort(., decreasing = T)
  setdiff(refMetData$genome_species, metaphlanDF$Species)
  refMetData$abu = metaphlanDF$relative_abundance[match(
    paste0(refMetData$pt_ID.u, '|', refMetData$genome_species),
    paste0(metaphlanDF$sample_id, '|', metaphlanDF$Species)
  )]
  refMetData$SNS_per_kbp = as.numeric(refMetData$SNS_count)*1000/refMetData$length
  refMetData$SNV_per_kbp = refMetData$SNV_count*1000/refMetData$length
  refMetData$d_prime_mean %<>% as.numeric()
  refMetData$r2_mean %<>% as.numeric()
  refMetData$lineage = sp2lineage$clade_name[match(refMetData$genome_species, sp2lineage$taxa %>% str_replace_all(., '_', ' '))]
  lineageDF = str_split_fixed(refMetData$lineage, '\\|', 7) %>% as.data.frame()
  table(lineageDF$V5) %>% length()
  lineageDF2 = lineageDF[!duplicated(lineageDF), ]
  refMetData$Phylum = lineageDF$V2 %>% str_remove('p__')
  refMetData$Family = lineageDF$V5 %>% str_remove('f__')
  table(  refMetData$Family )
  refMetData$if_Lactobacillaceae[lineageDF$V5 == 'f__Lactobacillaceae'] = 'Lactobacillaceae'
  refMetData$if_Lactobacillaceae[lineageDF$V5 != 'f__Lactobacillaceae'] = 'Non-Lactobacillaceae'
  refMetData$if_Lactobacillaceae %<>% as.factor()
  table(refMetData$if_Lactobacillaceae)
  aerobicDF = openxlsx::read.xlsx('./data/anaerobic or aerobic.xlsx')
  table(aerobicDF$anaerobic_or_aerobic2)
  refMetData$if_aerobic = aerobicDF$anaerobic_or_aerobic2[match(refMetData$genome_species, aerobicDF$species)]
  table(refMetData$if_aerobic)
  refMetData$D_rev = 1-refMetData$d_prime_mean
  
  refMetData2$ptr = ptrDF2$ptr[match(refMetData2$SampleID, ptrDF2$SampleID.u)]
}
refMetData2 = refMetData[refMetData$genome_species %in% names(table(refMetData$genome_species)[table(refMetData$genome_species) > 20]),]
table(refMetData2$genome_species) %>% length()
table(refMetData2$Family)
refMetData2$Family %<>% as.factor()
table(refMetData2$if_aerobic)
refMetData2$divergent_site_count %<>% as.numeric()
range(refMetData2$divergent_site_count)
mean(refMetData2$divergent_site_count)
sum(refMetData2$divergent_site_count)

range(refMetData2$SNPs_per_bp*1000)
mean(refMetData2$SNPs_per_bp*1000)

sumDivergentSite = refMetData2[, c('sample', 'divergent_site_count')] %>% 
  group_by(sample) %>% summarise(sum = sum(divergent_site_count))
range(sumDivergentSite$sum)
mean(sumDivergentSite$sum)
sd(sumDivergentSite$sum)
sumDivergentSite2 = refMetData2[, c('species', 'divergent_site_count')] %>% 
  group_by(species) %>% summarise(sum = sum(divergent_site_count))
range(sumDivergentSite2$sum)
mean(sumDivergentSite2$sum)
######### relationship among evolutionary metrics #####
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
    wx3 = wx3 + theme(axis.text.x = element_text(angle = xAngle, hjust = 1))
  }
  if(!is.null(legendPos)){
    wx3 = wx3 + theme(legend.position = legendPos)
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
                    color = 'x', title = titleLab,
                    add = "jitter", 
                    add.params = list(shape = 21, size = jitterSize, color = 'x', fill = 'x', alpha = 0.7),
                    bxp.errorbar = T, notch = F) +
      labs(y = yLab) +
      scale_color_manual(values = get(paste0('color', xval)), name = leName) +
      scale_fill_manual(values = get(paste0('color', xval)), name = leName) +
      mytheme3 + theme(strip.text = element_text(size = 12), axis.text.x = element_text(angle = x.angle, hjust = x.hjust),
                       plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 'cm'), legend.position = legend.position)
  }else{
    pbO = ggboxplot(data0, x = 'x', y = 'y', xlab = xLab,
                    color = 'black', fill = 'x', title = titleLab,
                    bxp.errorbar = T, notch = F) +
      labs(y = yLab) +
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
  #### one vs one ####
  if(T){ # overview: x = Pi , y = D', colored by pN/pS, size by coverage, shape by species group
    cairo_pdf('./plot/evolutionary_metrics_one_vs_one_ptr.pdf', width = 8, height = 3.63, onefile = T)
    relPlot(refMetData2, xval = 'nucl_diversity', yval = 'ptr', xLab = expression(italic(pi)), yLab = "PTR",
            gval = 'if_Lactobacillaceae', mycolor = colorLactobacillaceae)
    relPlot(refMetData2, xval = 'coverage', yval = 'ptr', xLab = 'Coverage',yLab = "PTR",
            gval = 'if_Lactobacillaceae', mycolor = colorLactobacillaceae)
    relPlot(refMetData2, xval = 'pnps', yval = 'ptr', xLab = 'pN/pS', yLab = "PTR",
            gval = 'if_Lactobacillaceae', mycolor = colorLactobacillaceae)
    relPlot(refMetData2, xval = 'D_rev', yval = 'ptr', xLab = "1 - D'",yLab = "PTR",
            gval = 'if_Lactobacillaceae', mycolor = colorLactobacillaceae)
    dev.off()
  }
  if(T){ # overview: x = Pi , y = D', colored by pN/pS, size by coverage, shape by species group
    cairo_pdf('./plot/evolutionary_metrics_one_vs_one_ptr2.pdf', width = 8, height = 3.63, onefile = T)
    relPlot(refMetData2, xval = 'nucl_diversity', yval = 'ptr', xLab = expression(italic(pi)), yLab = "PTR",
            mycolor = '#85B1D2')
    relPlot(refMetData2, xval = 'coverage', yval = 'ptr', xLab = 'Coverage', yLab = "PTR",
            mycolor = '#85B1D2')
    relPlot(refMetData2, xval = 'pnps', yval = 'ptr', xLab = 'pN/pS',  yLab = "PTR",
            mycolor = '#85B1D2')
    relPlot(refMetData2, xval = 'D_rev', yval = 'ptr', xLab = "1 - D'", yLab = "PTR",
            mycolor = '#85B1D2')
    dev.off()
  }
}

###### comparsion of Lactobacillus vs Non-Lactobacillus #####
if(T){
  ys = c('genome_species', 'if_Lactobacillaceae', 'Ethnicity2', 
         'SNV_per_kbp', 'D_rev', 'nucl_diversity', 'pnps', 'ptr')
  useData = refMetData2[, ys]
  useData2 = melt(useData, id.vars = c('genome_species', 'if_Lactobacillaceae', 'Ethnicity2'))
  colnames(useData2) = c('Species', 'Group', 'Ethnicity', 'Metrics', 'Value')
  useData2$Metrics %<>% as.character()
  useData2$Metrics[useData2$Metrics == 'nucl_diversity'] = 'Pi'
  useData2$Metrics[useData2$Metrics == 'SNV_per_kbp'] = 'SNVs/kbp'
  useData2$Metrics[useData2$Metrics == 'dnds'] = 'dN/dS'
  useData2$Metrics[useData2$Metrics == 'pnps'] = 'pN/pS'
  useData2$Metrics[useData2$Metrics == 'D_rev'] = "1 - D'"
  useData2$Metrics[useData2$Metrics == 'ptr'] = "PTR"
  useData2$Metrics %<>% factor(., levels = c('SNVs/kbp', 'Pi', "1 - D'", 'pN/pS', "PTR"))
}
if(T){
  cairo_pdf('./plot/evolution_metrics_nonLac_vs_Lac_ptr.pdf', width = 12.86, height = 1.78, onefile = T)
  stat.test <- useData2 %>% group_by(Metrics) %>% 
    wilcox_test(Value ~ Group) %>% 
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj") %>%
    add_xy_position(x = 'Group', dodge = 0, fun = 'max')
  stat.test$y.position = c(52, 0.025, 0.42, 1.9, 4)
  stat.test$p.adj %<>% signif(., digits = 3)
  p = ggboxplot(useData2, x = 'Group', y = 'Value', color = 'Group',
                add = 'jitter', add.params = list(size = 0.1))
  pp = p %>% facet(nrow = 1, facet.by = 'Metrics', scales = 'free_x') +
    stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.01,coord.flip = T) +
    mytheme + theme(strip.text = element_text(size = 13)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.13))) + 
    scale_x_discrete(limits=rev)+
    # labs(title = paste0(str_replace(sp, '.fna', ''), ' | ', x), 
    #      color = x, y = ysNorm[[y]], x = x) + 
    scale_color_manual(values = colorLactobacillaceae) +
    coord_flip()
  print(pp)
  dev.off()
}
##### comparison among species #####
if(T){
  cairo_pdf('./plot/evolution_metrics_species_ptr.pdf', width = 12.36, height = 3.85, onefile = T)
  piMean = useData2[useData2$Metrics == 'Pi', c('Species', 'Value')] %>% 
    group_by(Species) %>% summarise(mean = median(Value))
  useData2$Species = factor(useData2$Species, levels = piMean$Species[order(piMean$mean, decreasing = T)])
  p = ggboxplot(useData2, x = 'Species', y = 'Value', color = 'Group',
                add = 'jitter', add.params = list(size = 0.1))
  pp = p %>% facet(nrow = 1, facet.by = 'Metrics', scales = 'free_x') +
    # stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.01,
    #                    coord.flip = T,hide.ns = T) +
    mytheme + theme(strip.text = element_text(size = 13)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.13))) + 
    scale_x_discrete(limits=rev)+
    # labs(title = paste0(str_replace(sp, '.fna', ''), ' | ', x), 
    #      color = x, y = ysNorm[[y]], x = x) + 
    scale_color_manual(values = colorLactobacillaceae) +
    # geom_text(aes(Species, y = emmean,
    #               label= str_trim(.group)),
    #           data = model_means_cld, vjust = -0.5) +
    coord_flip()
  print(pp)
  
  for(m in c(unique(useData2$Metrics))){
    useData3 = useData2[useData2$Metrics == m,]
    Lmodel <- lm(Value ~ Species, data = useData3)
    model_means <- emmeans(object = Lmodel, specs = ~ Species, ) 
    model_means_cld <- cld(object = model_means,
                           adjust = "Tukey",
                           Letters = letters,
                           alpha = 0.05) %>% as.data.frame()
    p = ggboxplot(useData3, x = 'Species', y = 'Value', color = 'Group',
                  add = 'jitter', add.params = list(size = 0.1))
    pp = p %>% facet(nrow = 1, facet.by = 'Metrics', scales = 'free_x') +
      # stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.01,
      #                    coord.flip = T,hide.ns = T) +
      mytheme + theme(strip.text = element_text(size = 13)) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.13))) + 
      scale_x_discrete(limits=rev)+
      # labs(title = paste0(str_replace(sp, '.fna', ''), ' | ', x), 
      #      color = x, y = ysNorm[[y]], x = x) + 
      scale_color_manual(values = colorLactobacillaceae) +
      geom_text(aes(Species, y = emmean,
                    label= str_trim(.group)),
                data = model_means_cld, vjust = 0.5) +
      coord_flip()
    print(pp)
  }
}
dev.off()

if(T){
  cairo_pdf('./plot/evolution_metrics_nonLac_vs_Lac_facetBy_Ethnicity2_ptr.pdf', width = 12.85, height = 4.69, onefile = T)
  stat.test2 <- useData2 %>% group_by(Metrics, Ethnicity) %>% 
    wilcox_test(Value ~ Group) %>% 
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj") %>%
    add_xy_position(x = 'Group', dodge = 0)
  stat.test2$y.position = c(rep(52, 5),
                            rep(0.025, 5), 
                            rep(0.42, 5), 
                            rep(1.85, 5),
                            rep(4, 5))
  stat.test2$p.adj %<>% signif(., digits = 3)
  p = ggboxplot(useData2, x = 'Group', y = 'Value', 
                add = "jitter", add.params = list(shape = 1, size = 0.1),
                color = 'Group')
  pp = p %>% facet(nrow = 1, facet.by = c('Ethnicity', 'Metrics'), scales = 'free') +
    stat_pvalue_manual(stat.test2, label = "p.adj", tip.length = 0.01,coord.flip = T) +
    mytheme + theme(strip.text = element_text(size = 13)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.13))) + 
    scale_x_discrete(limits=rev)+
    # labs(title = paste0(str_replace(sp, '.fna', ''), ' | ', x), 
    #      color = x, y = ysNorm[[y]], x = x) + 
    scale_color_manual(values = colorLactobacillaceae) +
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
  ys = c('nucl_diversity', 'SNV_per_kbp', 'D_rev', 'pnps', 'ptr')
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
  lmeRes_ptr = lmeRes[lmeRes$y == 'ptr', ]
  cairo_pdf('./plot/evolution_metrics_Ethnicity_By_nonLac_vs_Lac_ptr.pdf', width = 4.45, height = 4.18, onefile = T)
  relAbuPlot(data0 = refMetData2[refMetData2$if_Lactobacillaceae == 'Lactobacillaceae',], 
             xval = 'Ethnicity2', yLab = 'PTR',titleLab = 'Lactobacillaceae',addJitter = T,
             yval = 'ptr', if_log10 = F, hide.ns = T, jitterSize = 0.1, 
             y.position = seq(3.5, 3.5+0.4*4, 0.4),
             aspectRatio = 1.2/1,x.angle = 30,
             expandMult = expansion(mult = c(0.03, 0.1)))
  relAbuPlot(data0 = refMetData2[refMetData2$if_Lactobacillaceae == 'Non-Lactobacillaceae',], 
             xval = 'Ethnicity2', yLab = 'PTR',titleLab = 'Non-Lactobacillacea',addJitter = T,
             yval = 'ptr', if_log10 = F, hide.ns = T, jitterSize = 0.1, 
             y.position = seq(3.5, 3.5+0.4*6, 0.4)-0.4,
             aspectRatio = 1.2/1,x.angle = 30,
             expandMult = expansion(mult = c(0.03, 0.1)))
  relAbuPlot(data0 = refMetData2[refMetData2$if_Lactobacillaceae == 'Lactobacillaceae',], 
             xval = 'Housing', yLab = 'PTR',titleLab = 'Lactobacillaceae',addJitter = T,
             yval = 'ptr', if_log10 = F, hide.ns = T, jitterSize = 0.1, 
             y.position = seq(3.5, 3.5+0.4*2, 0.4),
             aspectRatio = 1.2/1,x.angle = 30,
             expandMult = expansion(mult = c(0.03, 0.1)))
  relAbuPlot(data0 = refMetData2[refMetData2$if_Lactobacillaceae == 'Non-Lactobacillaceae',], 
             xval = 'Housing', yLab = 'PTR',titleLab = 'Non-Lactobacillacea',addJitter = T,
             yval = 'ptr', if_log10 = F, hide.ns = T, jitterSize = 0.1, 
             y.position = seq(3.5, 3.5+0.4*1, 0.4)-0.4,
             aspectRatio = 1.2/1,x.angle = 30,
             expandMult = expansion(mult = c(0.03, 0.1)))
  relAbuPlot(data0 = refMetData2[refMetData2$if_Lactobacillaceae == 'Lactobacillaceae',], 
             xval = 'Employment', yLab = 'PTR',titleLab = 'Lactobacillaceae',addJitter = T,
             yval = 'ptr', if_log10 = F, hide.ns = T, jitterSize = 0.1, 
             y.position = seq(3.5, 3.5+0.4*3, 0.4),
             aspectRatio = 1.2/1,x.angle = 30,
             expandMult = expansion(mult = c(0.03, 0.1)))
  relAbuPlot(data0 = refMetData2[refMetData2$if_Lactobacillaceae == 'Non-Lactobacillaceae',], 
             xval = 'Employment', yLab = 'PTR',titleLab = 'Non-Lactobacillacea',addJitter = T,
             yval = 'ptr', if_log10 = F, hide.ns = T, jitterSize = 0.1, 
             y.position = seq(3.5, 3.5+0.4*1, 0.4)-0.4,
             aspectRatio = 1.2/1,x.angle = 30,
             expandMult = expansion(mult = c(0.03, 0.1)))
  dev.off()
}
