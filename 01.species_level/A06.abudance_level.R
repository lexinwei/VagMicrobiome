source('./script/library_function.R')
load('./data/RData/swabData.RData')
load('./data/RData/metaData.RData')
load('./data/RData/abuList.RData')
load('./data/RData/clade2taxidDF.RData')
load('./data/RData/colorList.RData')
load('./data/RData/sampleVagitypes.RData')
swabData$Term = metaData$Term_char[match(swabData$pt_ID, metaData$pt_ID)]
swabData$Term2 = metaData$Term_char2[match(swabData$pt_ID, metaData$pt_ID)]
swabData$Ethnicity = metaData$Ethnicity[match(swabData$pt_ID, metaData$pt_ID)]
swabData$Ethnicity2 = metaData$Ethnicity2[match(swabData$pt_ID, metaData$pt_ID)]
if(T){ # prepare abundance matrix
  abuS = abuList[['s']]
  if(F){
    abuS2 = abuS[, c(1:3, match(swabData$seq_ID[!startsWith(swabData$pt_ID.u, 'KIT')], 
                                colnames(abuS)))]
    newColName = swabData$SampleID.u[match(colnames(abuS2)[4:322], swabData$seq_ID)]
    colnames(abuS2) = c(colnames(abuS2)[1:3], newColName)
    abuS2 = abuS2[, c(1:3,   naturalorder(colnames(abuS2)[4:322]) + 3)]
  
    write.xlsx(abuS2, file = './table/abuS2.xlsx')
  }
  
  abuS_mat = as.matrix(abuS[abuS$taxa != 'unclassified', 4:ncol(abuS)])
  row.names(abuS_mat) = abuS$taxa[abuS$taxa != 'unclassified']
  range(abuS_mat)
  # abuS_mat = apply(abuS_mat, 2, function(x){x/sum(x)})
  abuS_mat = abuS_mat/100
  range(abuS_mat)
  range(rowSums(abuS_mat))
  colSums(abuS_mat)
  abuS_mat_T = t(abuS_mat)
  interSample = intersect(rownames(abuS_mat_T), swabData$seq_ID[!grepl('^KIT', swabData$sample_ID)])
  sampleOrder = sampleVagitypes$seq_ID[order(sampleVagitypes$vagitypes_Susan,
                                             sampleVagitypes$vagitypes_Fettweis2, 
                                             sampleVagitypes$Term_char, sampleVagitypes$Ethnicity)]
  sampleOrder = sampleOrder[sampleOrder %in% interSample]
  abuMat = abuS_mat[, match(sampleOrder, colnames(abuS_mat))]
  dim(abuMat)
  abuMat_T = t(abuMat)
  table(colSums(abuMat_T)==0)
  # we retained taxa that either (1) 5% of the profiles exhibited an abundance of at least 1%, or (2) at least 15% of profiles exhibited an abundance of at least 0.1%. 
  criteria1 = which(rowSums(abuMat >= 0.01) >= ncol(abuMat)*0.05)
  criteria2 = which(rowSums(abuMat >= 0.001) >= ncol(abuMat)*0.15)
  setdiff(criteria2, criteria1)
  abuMat2 = abuMat[union(criteria1, criteria2), ]
  dim(abuMat2)
  abuMatOther = abuMat[-union(criteria1, criteria2), ]
  dim(abuMatOther)
  unkeepSp = row.names(abuMatOther)
  colSums(abuMatOther)
  
  abuMat3 = rbind(abuMat2, colSums(abuMatOther))
  abuMat3 = rbind(abuMat3, 1-colSums(abuMat3)) %>% as.matrix()
  colSums(abuMat3)
  row.names(abuMat3)[(nrow(abuMat2)+1):nrow(abuMat3)] = c('Other', 'Unclassified')
  range(abuMat3)
  abuMat3[abuMat3 < 0] = 0
  # shonnon diversity
  shannonDF0 = read.table('./metaphlan4_results/diversity/alpha_T_shannon.tsv', sep = '\t')
  shannonDF0$sample = str_remove(row.names(shannonDF0), 'Ming_nova_VS_') %>% str_remove(., '_bwa')
}
if(T){
  refDF0 = t(abuMat2*100) %>% as.data.frame()
  colnames(refDF0) %<>% str_replace_all(., '_', ' ')
  refDF0$seq_ID = row.names(refDF0)
  refDF0$Trimester = swabData$trimester[match(refDF0$seq_ID, swabData$seq_ID)]
  refDF0$Term = swabData$Term[match(refDF0$seq_ID, swabData$seq_ID)]
  refDF0$Term2 = swabData$Term2[match(refDF0$seq_ID, swabData$seq_ID)]
  refDF0$Ethnicity = swabData$Ethnicity[match(refDF0$seq_ID, swabData$seq_ID)]
  refDF0$Ethnicity2 = swabData$Ethnicity2[match(refDF0$seq_ID, swabData$seq_ID)]
  refDF0$pt_ID = swabData$pt_ID[match(refDF0$seq_ID, swabData$seq_ID)]
  refDF0$SampleID = swabData$SampleID[match(refDF0$seq_ID, swabData$seq_ID)]
  refDF0$Sample_GA = swabData$Sample_GA[match(refDF0$seq_ID, swabData$seq_ID)]
  refDF0$Housing = metaData$Housing[match(refDF0$pt_ID, metaData$pt_ID)]
  refDF0$Marriage = metaData$Marriage[match(refDF0$pt_ID, metaData$pt_ID)]
  refDF0$Employment = metaData$Employed[match(refDF0$pt_ID, metaData$pt_ID)]
  refDF0$Depression = metaData$Depression[match(refDF0$pt_ID, metaData$pt_ID)]
  refDF0$FOB = metaData$FOB[match(refDF0$pt_ID, metaData$pt_ID)]
  refDF0$BMI =  metaData$BMI[match(refDF0$pt_ID, metaData$pt_ID)]
  refDF0$Age =  metaData$Age[match(refDF0$pt_ID, metaData$pt_ID)]
  refDF0$Vagitype = sampleVagitypes$Vagitypes[match(refDF0$seq_ID, sampleVagitypes$seq_ID)]
  refDF0$CST = sampleVagitypes$CST[match(refDF0$seq_ID, sampleVagitypes$seq_ID)]
}
if(T){
  dim(abuMat2)
  colSums(abuMat2) %>% range()
  rowSums(abuMat2) %>% sort(.,decreasing = T)
  metaphlanDF = melt(abuMat2)
  colnames(metaphlanDF) = c('Species', 'seq_ID', 'relative_abundance')
  range(metaphlanDF$relative_abundance)
  metaphlanDF$patient_id = swabData$pt_ID[match(metaphlanDF$seq_ID, swabData$seq_ID)]
  metaphlanDF$sample_id = swabData$pt_ID.u[match(metaphlanDF$seq_ID, swabData$seq_ID)]
  metaphlanDF$Species  %<>% as.character() %>% str_replace_all(., '_', ' ')
  table(metaphlanDF$Species)
  metaphlanDF$Species %<>% factor(., levels = names(colorSpecies)[1:15] %>% rev())
  metaphlanDF$pt_ID = str_split_fixed(metaphlanDF$sample_id, '_', 2)[, 1]
  
  metaphlanDF$Pregnancy.period = metaData$Pregnancy.period[match(metaphlanDF$patient_id, metaData$pt_ID)]
  metaphlanDF$trimester = swabData$trimester[match(metaphlanDF$sample_id, swabData$pt_ID.u)]
  metaphlanDF$week = swabData$Sample_GA[match(metaphlanDF$sample_id, swabData$pt_ID.u)]
  metaphlanDF$week_PP = metaphlanDF$week
  metaphlanDF$week_PP[metaphlanDF$trimester == 'P'] = 'P'
  metaphlanDF$Term_char = metaData$Term_char[match(metaphlanDF$patient_id, metaData$pt_ID)]
  metaphlanDF$Term_char2 = metaData$Term_char2[match(metaphlanDF$patient_id, metaData$pt_ID)]
  metaphlanDF$Ethnicity = metaData$Ethnicity[match(metaphlanDF$patient_id, metaData$pt_ID)] %>% as.character()
  metaphlanDF$Ethnicity[metaphlanDF$Ethnicity == 'Pacific Islander'] = 'Pacific'
  metaphlanDF$Ethnicity %<>% factor(., levels = c('White', 'Asian', 'Latina', 'Black', 'Pacific', 'Other'))
  metaphlanDF$Ethnicity2 = metaData$Ethnicity2[match(metaphlanDF$patient_id, metaData$pt_ID)]
}

###### function ######
relAbuPlot = function(data0, xval, yval, zval = NULL, y.position = NULL,step.increase = 0.12, expandMult = NULL, hide.x.label = F,
                        x.angle = 0, CLD = F, hide.ns = F, legend.position = 'right', if_log10 = F, show_Padj = F,
                      legend.direction = 'vertical'){
  data0$x = data0[, xval]
  if(if_log10){
    data0$y = log10(data0[, yval]+0.01)
    yLab = paste0('log10(Relative abundance)')
  }else{
    data0$y = data0[, yval]
    yLab = 'Relative abundance'
  }
  titleLab = yval
  if(endsWith(xval, '2')){xLab = str_remove(xval, '2$')}else{xLab = xval}
  if(!is.null(zval)){data0$z = data0[, zval]}
  if(x.angle == 0){x.hjust = 0.5}else{x.hjust = 1}
  # base ggplot
  pbO = ggboxplot(data0, x = 'x', y = 'y', xlab = xLab,
            ylab = yLab, color = 'x', title = titleLab,
            add = "jitter", 
            add.params = list(shape = 21, size = 1.8, color = 'x', fill = 'x', alpha = 0.7),
                              bxp.errorbar = T, notch = F) +
    scale_color_manual(values = get(paste0('color', xval)), name = xLab) +
    scale_fill_manual(values = get(paste0('color', xval)), name = xLab) +
    mytheme3 + theme(strip.text = element_text(size = 12), axis.text.x = element_text(angle = x.angle, hjust = x.hjust),
                     plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 'cm'), legend.position = legend.position)

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
    pbO = pbO + facet_grid(.~z) +
      theme(legend.position = legend.position, legend.direction = legend.direction, 
              legend.spacing.x = unit(0.1, 'cm'), panel.spacing.x = unit(0.1, 'cm'))
  }
  if(hide.x.label){
    pbO = pbO + theme(axis.text.x = element_blank())
  }
  print(pbO)
}
##### lme vs Ethnicity#####
if(T){
  metaphlanDF$Employment = metaData$Employed[match(metaphlanDF$pt_ID, metaData$pt_ID)]
  metaphlanDF$Housing = metaData$Housing[match(metaphlanDF$pt_ID, metaData$pt_ID)]
  metaphlanDF$Depression = metaData$Depression[match(metaphlanDF$pt_ID, metaData$pt_ID)]
  metaphlanDF$BMI = metaData$BMI[match(metaphlanDF$pt_ID, metaData$pt_ID)]
  metaphlanDF$Age = metaData$Age[match(metaphlanDF$pt_ID, metaData$pt_ID)]
  metaphlanDF$Sample_GA = swabData$Sample_GA[match(metaphlanDF$pt_ID, swabData$pt_ID)]
  metaphlanDF$FOB = metaData$FOB[match(metaphlanDF$pt_ID, metaData$pt_ID)]
  metaphlanDF$Marriage = metaData$Marriage[match(metaphlanDF$pt_ID, metaData$pt_ID)]
  metaphlanDF$Marriage %<>% as.character()
  metaphlanDF$Marriage[is.na(metaphlanDF$Marriage)] = 'Unknown'
  metaphlanDF$Marriage %<>% factor(., levels = c('Married', 'Unmarried', 'Unknown'))
  # xs = c(#'Insulin', 'B6', 
  #   'Ethnicity2','Employment', 'Housing',
  #   'Term_char', # 'Abortions', 
  #   'BMI', 'Marriage','Age',  'trimester',
  #   'FOB', 'Depression')  # interest priorit
  xs = c(#'Insulin', 'B6', 
    'Employment', 'Housing',
    'Depression',
    'FOB', 'Marriage',
    'Term_char', # 'Abortions', 
    'Ethnicity2', 'trimester',
    'BMI', # 'BMI_char', # 'trimester'
    'Age')  # interest priority

  y = 'relative_abundance'
  lme.Multi.df = data.frame()
  for(s in unique(metaphlanDF$Species)){
    dat =  metaphlanDF[metaphlanDF$Species == s, ]
    lmeTry = lme(as.formula(paste0(y, ' ~ ', str_c(xs, collapse = '+'))), random = ~ 1 | pt_ID,
                                   data = dat, method = 'ML')
    library(multcomp)
    
    lme.Multi = lmeTry
    # summary(glht(lme.Multi, linfct = mcp(Ethnicity2 = 'Tukey')), test = adjusted("BH"))
    lme.Multi.R2 = r.squaredGLMM(lme.Multi)[1]
    lme.Multi.aov = anova(lme.Multi)
    formula_char = paste0(y, ' ~ (Intercept) + ', str_c(capitalize(xs), collapse = ' + '), ' + (1|Subject)')
    lme.Multi.df_tmp = data.frame(species = s,
                           y = 'relative_abundance',
                           x = row.names(lme.Multi.aov), 
                           lme.multi.Model = formula_char,
                           lme.multi.F_value = lme.Multi.aov$`F-value`,
                           lme.multi.P_value = lme.Multi.aov$`p-value`,
                           lme.multi.R2 = lme.Multi.R2)
    lme.Multi.df = rbind(lme.Multi.df, lme.Multi.df_tmp)
  }
  lme.Multi.df = lme.Multi.df[lme.Multi.df$x != '(Intercept)', ]
  openxlsx::write.xlsx(lme.Multi.df, file = './data/abundance_lme.Multi.df_v2.xlsx', rowNames =F)
  lme.Multi.df_sig = lme.Multi.df[lme.Multi.df$lme.multi.P_value < 0.05, ]
  # facet by trimester
  lme.Multi.df2 = data.frame()
  for(s in unique(metaphlanDF$Species)){
    for(t in c('T1', 'T2', 'T3', 'P')){
      dat =  metaphlanDF[metaphlanDF$Species == s & metaphlanDF$trimester == t, ]
      lmeTry = try({lme(as.formula(paste0(y, ' ~ ', str_c(setdiff(xs, 'trimester'), collapse = '+'))), random = ~ 1 | pt_ID,
                   data = dat, method = 'ML')})
      if(class(lmeTry) != 'try-error'){
        lme.Multi = lmeTry
        lme.Multi.R2 = r.squaredGLMM(lme.Multi)[1]
        lme.Multi.aov = anova(lme.Multi)
        formula_char = paste0(y, ' ~ (Intercept) + ', str_c(capitalize(xs), collapse = ' + '), ' + (1|Subject)')
        lme.Multi.df2_tmp = data.frame(species = s,
                                       trimester = t,
                                       y = 'relative_abundance',
                                       x = row.names(lme.Multi.aov), 
                                       lme.multi.Model = formula_char,
                                       lme.multi.F_value = lme.Multi.aov$`F-value`,
                                       lme.multi.P_value = lme.Multi.aov$`p-value`,
                                       lme.multi.R2 = lme.Multi.R2)
        lme.Multi.df2 = rbind(lme.Multi.df2, lme.Multi.df2_tmp)
      }
    }
  }
  lme.Multi.df2 = lme.Multi.df2[lme.Multi.df2$x != '(Intercept)', ]
  openxlsx::write.xlsx(lme.Multi.df2, file = './data/abundance_lme.Multi.df2_facet_by_trimester_v2.xlsx', rowNames =F)
  lme.Multi.df2_sig = lme.Multi.df2[lme.Multi.df2$lme.multi.P_value < 0.05, ]
}
# overview plot of ethnicity
if(T){
  lme.Multi.df = readxl::read_excel(path = './data/abundance_lme.Multi.df_v2.xlsx')
  lme.Multi.df_eth = lme.Multi.df[lme.Multi.df$x == 'Ethnicity2', ]
  lme.Multi.df_eth$species = factor(lme.Multi.df_eth$species, 
                                    levels = lme.Multi.df_eth$species[order(lme.Multi.df_eth$lme.multi.P_value, decreasing = T)])
  lme.Multi.df_eth$label = 'Significant'
  lme.Multi.df_eth$label[lme.Multi.df_eth$lme.multi.P_value > 0.05] = 'Non-significant'
  bbp = ggplot(lme.Multi.df_eth, aes(x=-log10(lme.multi.P_value), y=species)) +
    geom_col(fill="gray", width = 0.25, alpha = 0.8) +
    geom_point(aes(fill = label), color = 'black', stroke = 0.2, shape = 21, size = 4) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
    scale_fill_manual(values = c('Significant' = '#F8766D', 'Non-significant' = '#01BFC4'), 
                      name = 'Difference in Ethnicity',
                      guide = guide_legend(title.position = 'left', title.hjust = 1)) +
    geom_vline(xintercept = -log10(0.05), linetype = 'dashed', color = 'gray') +
    mytheme3 + theme(strip.background = element_rect(fill = '#E2E5E9', color = NA), 
                     # axis.title = element_text(size = 12), axis.text = element_text(size = 12), 
                     # legend.title = element_text(size = 12), 
                     # legend.text =  element_text(size = 12), 
                     legend.position = 'bottom', legend.box.margin = margin(-10,100,0,-120),
                     legend.margin = margin(0,0,0,0), legend.spacing.x = unit(0.02, 'cm')
                     ) +
    xlab(expression('-log10(' ~ italic(P) ~ 'value)')) + ylab("")
  cairo_pdf('./plot/ethnicity_lme_pvalue_v2.pdf',height = 3.83, width = 4.72)
  print(bbp)
  dev.off()
}

##### plot relative abundance boxplot ####
# one for all, by term2
cairo_pdf('./plot/abundance_boxplot_by_term2.pdf', width = 5.14, height = 5.53, onefile = T)
if(T){

  stat.test <- metaphlanDF %>%
    group_by(Species) %>%
    wilcox_test(relative_abundance ~ Term_char2) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj")

  pgW = ggplot(metaphlanDF, aes(x = log10(relative_abundance*100 + 0.01), y = Species)) + 
    labs(y = 'Species', x = 'log10(Relative abundance)', title = NULL) +
    geom_boxplot(aes(color = Term_char2), outlier.shape = NA, width = 0.8) +
    geom_point(aes(color = Term_char2), position=position_jitterdodge(jitter.width = 0.3,dodge.width = 0.8), size = 0.3) + 
    scale_color_manual(name = 'Term', values = colorTerm2) +
    mytheme + theme(text = element_text(size = 15, colour = 'black'), legend.position = 'top',
                    legend.spacing.x = unit(0.1,'cm'),
                    legend.title = element_blank(), legend.margin = margin(0,0,-7,0))
  print(pgW + geom_text(data = stat.test, aes(x = 3.5, y = Species, label = signif(p, 2)), hjust = 1, size = 4)+
          theme(legend.margin = margin(0,30,-5,-7)))
  # print(pgW + geom_text(data = stat.test, aes(x = 3, y = Species, label = signif(p.adj, 2)), hjust = 1, size = 4.5))
}
dev.off()
# one by one bxplot, roughly
cairo_pdf('./plot/abundance_boxplot.pdf', width = 5.89, height = 4.73, onefile = T)
for(yDiv in names(colorSpecies)[1:15]){
  relAbuPlot(data0 = refDF0, xval = 'Term', yval = yDiv, if_log10 = T)
  relAbuPlot(data0 = refDF0, xval = 'Term2', yval = yDiv, if_log10 = T)
  relAbuPlot(data0 = refDF0, xval = 'Ethnicity2', yval = yDiv, if_log10 = T, hide.ns = T, step.increase = 20.)
  relAbuPlot(data0 = refDF0, xval = 'Ethnicity2', yval = yDiv, if_log10 = T, CLD = T)
}
dev.off()
cairo_pdf('./plot/abundance_boxplot_by_trimester.pdf', width = 6.38, height = 4.73, onefile = T)
for(yDiv in names(colorSpecies)[1:15]){
  try({
    relAbuPlot(data0 = refDF0, xval = 'Ethnicity2', yval = yDiv, z = 'Trimester', x.angle = 45, if_log10 = T, hide.ns = T, 
               step.increase = 0.15, 
               expandMult = expansion(mult = c(0.03, 0.1)))
  })
}
dev.off()

##### careful re plot good results #####
if(T){
  lme.Multi.df = read.xlsx('./data/abundance_lme.Multi.df_v2.xlsx')
  lme.Multi.df.sig = lme.Multi.df[lme.Multi.df$lme.multi.P_value < 0.05, ]
  lme.Multi.df.sig_sel = lme.Multi.df.sig[lme.Multi.df.sig$x %in% c('Marriage', 'Housing', 'Employment', 'FOB', 'Depression'), ]
  cairo_pdf('./plot/abundance_boxplot_by_sig.pdf', width = 4.95, height = 3.59, onefile = T)
  for(r in 1:nrow(lme.Multi.df.sig_sel)){
    xv = lme.Multi.df.sig_sel$x[r]
    yv = lme.Multi.df.sig_sel$species[r]
   
    relAbuPlot(data0 = refDF0, xval = xv, yval = yv, if_log10 = T, 
               expandMult = expansion(mult = c(0.03, 0.1)), legend.position = 'none')
  }
  dev.off()
  cairo_pdf('./plot/abundance_boxplot_by_Marriage.pdf', width = 3.57, height = 3.59, onefile = T)
  refDF0$Marriage %<>% as.character()
  refDF0$Marriage[is.na(refDF0$Marriage)] = 'Unknown'
  refDF0$Marriage %<>% factor(., levels = c('Married', 'Unmarried', 'Unknown'))
  relAbuPlot(data0 = refDF0, xval = 'Marriage', yval = 'Fannyhessea vaginae',  if_log10 = T, 
             hide.ns = T, CLD = F,
             expandMult = expansion(mult = c(0.03, 0.1)), legend.position = 'none')
  dev.off()
  
  cairo_pdf('./plot/abundance_boxplot_by_Depression.pdf', width = 2.84, height = 3.59, onefile = T)
  relAbuPlot(data0 = refDF0, xval = 'Depression', yval = 'Gardnerella vaginalis',  if_log10 = T, 
             hide.ns = T, CLD = F,
             expandMult = expansion(mult = c(0.03, 0.1)), legend.position = 'none')
  dev.off()
}
if(T){
  cairo_pdf('./plot/abundance_boxplot_by_term.pdf', width = 2.97, height = 3.59, onefile = T)
  relAbuPlot(data0 = refDF0, xval = 'Term', yval = 'Lactobacillus crispatus', if_log10 = T, 
             expandMult = expansion(mult = c(0.03, 0.1)), legend.position = 'none')
  relAbuPlot(data0 = refDF0, xval = 'Term', yval = 'Gardnerella vaginalis', if_log10 = T, 
             expandMult = expansion(mult = c(0.03, 0.1)), legend.position = 'none')
  relAbuPlot(data0 = refDF0, xval = 'Term', yval = 'Lactobacillus iners', if_log10 = T, 
             expandMult = expansion(mult = c(0.03, 0.1)), legend.position = 'none')
  dev.off()
}
cairo_pdf('./plot/fanny.pdf', width = 1.82, height = 3.33)
relAbuPlot(data0 = refDF0, xval = 'Term2', yval = 'Fannyhessea vaginae', if_log10 = T, 
           expandMult = expansion(mult = c(0.03, 0.1)), legend.position = 'none')
dev.off()
if(T){
  cairo_pdf('./plot/abundance_boxplot_by_Ethnicity2_v2.pdf', width = 3.82, height = 3.80, onefile = T)
  relAbuPlot(data0 = refDF0, xval = 'Ethnicity2', yval = "Megasphaera genomosp type 1", if_log10 = T, hide.ns = T, 
             step.increase = 20, legend.position = 'none', 
             expandMult = expansion(mult = c(0.03, 0.1)))
  relAbuPlot(data0 = refDF0, xval = 'Ethnicity2', yval = 'Lactobacillus crispatus', if_log10 = T, hide.ns = T, 
             y.position = c(3,4,5,6,7,8,9), step.increase = 20, legend.position = 'none', 
             expandMult = expansion(mult = c(0.03, 0.1)))
  relAbuPlot(data0 = refDF0, xval = 'Ethnicity2', yval = 'Gardnerella vaginalis', if_log10 = T, hide.ns = T, 
             y.position = c(3,4,5,6,7)/1.5+0.3, step.increase = 20, legend.position = 'none', 
             expandMult = expansion(mult = c(0.03, 0.1)))
  relAbuPlot(data0 = refDF0, xval = 'Ethnicity2', yval = 'Lactobacillus iners', if_log10 = T, hide.ns = T, 
             y.position = c(3,4,5,6)/2 + 0.7, step.increase = 20, legend.position = 'none', 
             expandMult = expansion(mult = c(0.03, 0.1)))
  relAbuPlot(data0 = refDF0, xval = 'Ethnicity2', yval = 'Fannyhessea vaginae', if_log10 = T, hide.ns = T, 
             y.position = c(3,4,5,6)/2.2 + 0.8, step.increase = 20, legend.position = 'none', 
             expandMult = expansion(mult = c(0.03, 0.1)))
  relAbuPlot(data0 = refDF0, xval = 'Ethnicity2', yval = 'Lactobacillus jensenii', if_log10 = T, hide.ns = T, 
             y.position = c(3,4,5)/2+0.7, step.increase = 20, legend.position = 'none',
             expandMult = expansion(mult = c(0.03, 0.1)))
  relAbuPlot(data0 = refDF0, xval = 'Ethnicity2', yval = 'Lactobacillus gasseri', if_log10 = T, hide.ns = T, 
             y.position = c(3,4,5,6)/2 + 0.7, step.increase = 20, legend.position = 'none', 
             expandMult = expansion(mult = c(0.03, 0.1)))
  relAbuPlot(data0 = refDF0, xval = 'Ethnicity2', yval = 'Actinomycetaceae SGB989', if_log10 = T, hide.ns = T, 
             y.position = c(3,4,5,6,7,8)/1.5-0.2, step.increase = 20, legend.position = 'none', 
             expandMult = expansion(mult = c(0.03, 0.1)))
  relAbuPlot(data0 = refDF0, xval = 'Ethnicity2', yval = 'Veillonellaceae bacterium DNF00626', if_log10 = T, hide.ns = T, 
             y.position = c(3,4,5,6,7,8)/2-0.2, step.increase = 20, legend.position = 'none', 
             expandMult = expansion(mult = c(0.03, 0.1)))
  # relAbuPlot(data0 = refDF0, xval = 'Ethnicity2', yval = 'Veillonellaceae bacterium DNF00626', if_log10 = T, hide.ns = T, 
  #            y.position = c(3,4,5,6), step.increase = 20, legend.position = 'none', show_Padj = T,
  #            expandMult = expansion(mult = c(0.03, 0.1)))
  relAbuPlot(data0 = refDF0, xval = 'Ethnicity2', yval = 'Coriobacteriales bacterium DNF00809', if_log10 = T, hide.ns = T, 
             y.position = c(3,4,5,6,7,8)/2.5-0.5, step.increase = 20, legend.position = 'none', 
             expandMult = expansion(mult = c(0.03, 0.1)))
  relAbuPlot(data0 = refDF0, xval = 'Ethnicity2', yval = 'Aerococcus christensenii', if_log10 = T, hide.ns = T, 
             y.position = c(3,4,5,6,7,8)/2.1-0.6, step.increase = 10, legend.position = 'none', 
             expandMult = expansion(mult = c(0.03, 0.1)))
  relAbuPlot(data0 = refDF0, xval = 'Ethnicity2', yval = 'Candida albicans', if_log10 = T, hide.ns = T, 
             y.position = c(3,4,5,6,7)/2, step.increase = 20, legend.position = 'none', 
             expandMult = expansion(mult = c(0.03, 0.1)))
  relAbuPlot(data0 = refDF0, xval = 'Ethnicity2', yval = 'Hungateiclostridiaceae SGB4003', if_log10 = T, hide.ns = T, 
             y.position = c(3,4,5,6)/2-0.5, step.increase = 20, legend.position = 'none', 
             expandMult = expansion(mult = c(0.03, 0.1)))
  dev.off()
}
relAbuPlot(data0 = refDF0, xval = 'Ethnicity2', yval = "Megasphaera genomosp type 1", if_log10 = T, hide.ns = T, 
           step.increase = 20, legend.position = 'top', 
           expandMult = expansion(mult = c(0.03, 0.1))) + theme(legend.spacing.x = unit(0.08,'cm'))
if(T){
  cairo_pdf('./plot/abundance_boxplot_by_Ethnicity2_facetBy_trimester_v3.pdf', width = 7.73, height = 3.33, onefile = T)
  relAbuPlot(data0 = refDF0, xval = 'Ethnicity2', yval = 'Lactobacillus crispatus', z = 'Trimester', x.angle = 45, if_log10 = T, hide.ns = T, 
             y.position = c(2.5:5.5, 2.5:7.5)/1.2+0.3, step.increase = 20,  
             expandMult = expansion(mult = c(0.03, 0.1)), hide.x.label = T)
  relAbuPlot(data0 = refDF0, xval = 'Ethnicity2', yval = 'Gardnerella vaginalis', z = 'Trimester', x.angle = 45, if_log10 = T, hide.ns = T, 
             y.position = c(2.5:5.5, 2.5)/1.7 + 0.7, step.increase = 20, 
             expandMult = expansion(mult = c(0.03, 0.1)), hide.x.label = T)
  relAbuPlot(data0 = refDF0, xval = 'Ethnicity2', yval = 'Lactobacillus iners', z = 'Trimester', x.angle = 45, if_log10 = T, hide.ns = T, 
             y.position = c(2.5:4.5, 2.5:6.5)/1.5 + 0.5, step.increase = 20, 
             expandMult = expansion(mult = c(0.03, 0.1)), hide.x.label = T)
  relAbuPlot(data0 = refDF0, xval = 'Ethnicity2', yval = 'Lactobacillus jensenii', z = 'Trimester', x.angle = 45, if_log10 = T, hide.ns = T, 
             y.position = c(2.5, 2.5:6.5)/1.5+0.5, step.increase = 20, 
             expandMult = expansion(mult = c(0.03, 0.1)), hide.x.label = T)
  relAbuPlot(data0 = refDF0, xval = 'Ethnicity2', yval = 'Lactobacillus gasseri', z = 'Trimester', x.angle = 45, if_log10 = T, hide.ns = T, 
             y.position = c(2.5:5.5)/1.5 + 0.5, step.increase = 20, 
             expandMult = expansion(mult = c(0.03, 0.1)), hide.x.label = T)
  relAbuPlot(data0 = refDF0, xval = 'Ethnicity2', yval = 'Fannyhessea vaginae', z = 'Trimester', x.angle = 45, if_log10 = T, hide.ns = T, 
             y.position = c(2.5,3.5,2.5:5.5)/1.5 + 0.5, step.increase = 20, 
             expandMult = expansion(mult = c(0.03, 0.1)), hide.x.label = T)
  relAbuPlot(data0 = refDF0, xval = 'Ethnicity2', yval = 'Lactobacillus vaginalis', z = 'Trimester', x.angle = 45, if_log10 = T, hide.ns = T, 
             y.position = c(2.5, 2.5:7.5)/1.3+0.3, step.increase = 20, 
             expandMult = expansion(mult = c(0.03, 0.1)), hide.x.label = T)
  relAbuPlot(data0 = refDF0, xval = 'Ethnicity2', yval = 'Megasphaera genomosp type 1', z = 'Trimester', x.angle = 45, if_log10 = T, hide.ns = T, 
             y.position = c(1.9), step.increase = 20, 
             expandMult = expansion(mult = c(0.03, 0.1)), hide.x.label = T)
  relAbuPlot(data0 = refDF0, xval = 'Ethnicity2', yval = "Candida albicans", z = 'Trimester', x.angle = 45, if_log10 = T, hide.ns = F, 
             y.position = c(1.9), step.increase = 20, 
             expandMult = expansion(mult = c(0.03, 0.1)), hide.x.label = T)
  relAbuPlot(data0 = refDF0, xval = 'Ethnicity2', yval = "Candida albicans", z = 'Trimester', x.angle = 45, if_log10 = T, hide.ns = F, 
             y.position = c(1.9), step.increase = 20, CLD = T,
             expandMult = expansion(mult = c(0.03, 0.4)), hide.x.label = T)
  relAbuPlot(data0 = refDF0, xval = 'Ethnicity2', yval = "Actinomycetaceae SGB989", z = 'Trimester', x.angle = 45, if_log10 = T, hide.ns = F, 
             expandMult = expansion(mult = c(0.03, 0.08)), 
             hide.x.label = T)
  relAbuPlot(data0 = refDF0, xval = 'Ethnicity2', yval = "Actinomycetaceae SGB989", z = 'Trimester', x.angle = 45, if_log10 = T, hide.ns = F, 
             expandMult = expansion(mult = c(0.03, 0.4)), CLD = T,
             hide.x.label = T)
  relAbuPlot(data0 = refDF0, xval = 'Ethnicity2', yval = "Aerococcus christensenii", z = 'Trimester', x.angle = 45, if_log10 = T, hide.ns = T, 
             expandMult = expansion(mult = c(0.03, 0.1)),  y.position = c(2.5:7.5,2.5:6.5)/2-0.3, step.increase = 20, 
             hide.x.label = T)
  dev.off()
}

if(T){
  cairo_pdf('./plot/abundance_boxplot_by_Term2_facetBy_trimester.pdf', width = 3.74, height = 3.25, onefile = T)
  relAbuPlot(data0 = refDF0, xval = 'Term2', yval = 'Lactobacillus crispatus', z = 'Trimester', x.angle = 45, if_log10 = T, 
             hide.ns = F, expandMult = expansion(mult = c(0.03, 0.1)), hide.x.label = T,
             legend.position = 'none', legend.direction = 'horizontal')
  relAbuPlot(data0 = refDF0, xval = 'Term2', yval = 'Gardnerella vaginalis', z = 'Trimester', x.angle = 45, if_log10 = T, 
             hide.ns = F, expandMult = expansion(mult = c(0.03, 0.1)), hide.x.label = T,
             legend.position = 'none', legend.direction = 'horizontal')
  relAbuPlot(data0 = refDF0, xval = 'Term2', yval = 'Lactobacillus iners', z = 'Trimester', x.angle = 45, if_log10 = T, 
             hide.ns = F, expandMult = expansion(mult = c(0.03, 0.1)), hide.x.label = T,
             legend.position = 'none', legend.direction = 'horizontal')
  relAbuPlot(data0 = refDF0, xval = 'Term2', yval = 'Lactobacillus jensenii', z = 'Trimester', x.angle = 45, if_log10 = T, 
             hide.ns = F, expandMult = expansion(mult = c(0.03, 0.1)), hide.x.label = T,
             legend.position = 'none', legend.direction = 'horizontal')
  relAbuPlot(data0 = refDF0, xval = 'Term2', yval = 'Lactobacillus gasseri', z = 'Trimester', x.angle = 45, if_log10 = T, 
             hide.ns = F, expandMult = expansion(mult = c(0.03, 0.1)), hide.x.label = T,
             legend.position = 'none', legend.direction = 'horizontal')
  relAbuPlot(data0 = refDF0, xval = 'Term2', yval = 'Fannyhessea vaginae', z = 'Trimester', x.angle = 45, if_log10 = T, 
             hide.ns = F, expandMult = expansion(mult = c(0.03, 0.1)), hide.x.label = T,
             legend.position = 'none', legend.direction = 'horizontal')
  relAbuPlot(data0 = refDF0, xval = 'Term2', yval = 'Lactobacillus vaginalis', z = 'Trimester', x.angle = 45, if_log10 = T, 
             hide.ns = F, expandMult = expansion(mult = c(0.03, 0.1)), hide.x.label = T,
             legend.position = 'none', legend.direction = 'horizontal')
  relAbuPlot(data0 = refDF0, xval = 'Term2', yval = "Megasphaera genomosp type 1", z = 'Trimester', x.angle = 45, if_log10 = T, 
             hide.ns = F, expandMult = expansion(mult = c(0.03, 0.1)), hide.x.label = T,
             legend.position = 'none', legend.direction = 'horizontal')
  relAbuPlot(data0 = refDF0, xval = 'Term2', yval = "Aerococcus christensenii", z = 'Trimester', x.angle = 45, if_log10 = T, 
             hide.ns = F, expandMult = expansion(mult = c(0.03, 0.1)), hide.x.label = T,
             legend.position = 'none', legend.direction = 'horizontal')
  dev.off()
}
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
cairo_pdf('./plot/lme_abundance.pdf', width = 4.08, height = 3.01, onefile = T)
if(T){
  lmeDF = refDF0[refDF0$Trimester != 'P', ]
  yDiv = 'Lactobacillus crispatus'
  lmeDF$y = log10(lmeDF[, yDiv] + 0.01)
  lmePlot(lmeData = lmeDF, yLab = 'log10(Relative abundance)', titleLab = yDiv, if_logged = T, P.pos = c(10,1))
  lmePlot(lmeData = lmeDF, yLab = 'log10(Relative abundance)', titleLab = yDiv,
          if_logged = T, P.pos = c(10,1), groupV = 'Term2',addSE = T,
          expandMulti = expansion(mult = c(0.05, 0.1)), yPos = c(1, NA))
  lmePlot(lmeData = lmeDF, yLab = 'log10(Relative abundance)', titleLab = yDiv,
          if_logged = T, P.pos = c(10,1), groupV = 'Ethnicity2',addSE = T,
          expandMulti = expansion(mult = c(0.05, 0.1)), yPos = c(0.7,1))
  
  yDiv = 'Gardnerella vaginalis'
  lmeDF$y = log10(lmeDF[, yDiv] + 0.01)
  lmePlot(lmeData = lmeDF, yLab = 'log10(Relative abundance)', titleLab = yDiv, if_logged = T, P.pos = c(10, -2))
  lmePlot(lmeData = lmeDF, yLab = 'log10(Relative abundance)', titleLab = yDiv,
          if_logged = T, P.pos = c(10,1), groupV = 'Term2',addSE = T,
          expandMulti = expansion(mult = c(0.05, 0.1)), yPos = c(-2,-1.7))
  lmePlot(lmeData = lmeDF, yLab = 'log10(Relative abundance)', titleLab = yDiv,
          if_logged = T, P.pos = c(10,1), groupV = 'Ethnicity2',addSE = F,
          expandMulti = expansion(mult = c(0.05, 0.1)), yPos = c(2.0,2.3,2.6,2.9)-4)
  yDiv = 'Lactobacillus iners'
  lmeDF$y = log10(lmeDF[, yDiv] + 0.01)
  lmePlot(lmeData = lmeDF, yLab = 'log10(Relative abundance)', titleLab = yDiv, if_logged = T, P.pos = c(10, 0.5))
  lmePlot(lmeData = lmeDF, yLab = 'log10(Relative abundance)', titleLab = yDiv,
          if_logged = T, P.pos = c(10,1), groupV = 'Term2',addSE = T,
          expandMulti = expansion(mult = c(0.05, 0.1)), yPos = c(-2,-1.7)+4)
  lmePlot(lmeData = lmeDF, yLab = 'log10(Relative abundance)', titleLab = yDiv,
          if_logged = T, P.pos = c(10,1), groupV = 'Ethnicity2',addSE = F,
          expandMulti = expansion(mult = c(0.05, 0.1)), yPos = c(2.0,2.3,2.6))
  yDiv = 'Lactobacillus jensenii'
  lmeDF$y = log10(lmeDF[, yDiv] + 0.01)
  lmePlot(lmeData = lmeDF, yLab = 'log10(Relative abundance)', titleLab = yDiv, if_logged = T, P.pos = c(10, 1))
  lmePlot(lmeData = lmeDF, yLab = 'log10(Relative abundance)', titleLab = yDiv,
          if_logged = T, P.pos = c(10,1), groupV = 'Term2',addSE = T,
          expandMulti = expansion(mult = c(0.05, 0.1)), yPos = c(1,1.3))
  lmePlot(lmeData = lmeDF, yLab = 'log10(Relative abundance)', titleLab = yDiv,
          if_logged = T, P.pos = c(10,1), groupV = 'Ethnicity2',addSE = F,
          expandMulti = expansion(mult = c(0.05, 0.1)), yPos = c(2.0,2.3,2.6))
  yDiv = 'Lactobacillus jensenii'
  lmeDF$y = log10(lmeDF[, yDiv] + 0.01)
  lmePlot(lmeData = lmeDF, yLab = 'log10(Relative abundance)', titleLab = yDiv, if_logged = T, P.pos = c(10, 1))
  lmePlot(lmeData = lmeDF, yLab = 'log10(Relative abundance)', titleLab = yDiv,
          if_logged = T, P.pos = c(10,1), groupV = 'Term2',addSE = T,
          expandMulti = expansion(mult = c(0.05, 0.1)), yPos = c(1,1.3))
  lmePlot(lmeData = lmeDF, yLab = 'log10(Relative abundance)', titleLab = yDiv,
          if_logged = T, P.pos = c(10,1), groupV = 'Ethnicity2',addSE = F,
          expandMulti = expansion(mult = c(0.05, 0.1)), yPos = c(2.0,2.3,2.6))
  yDiv = 'Fannyhessea vaginae'
  lmeDF$y = log10(lmeDF[, yDiv] + 0.01)
  lmePlot(lmeData = lmeDF, yLab = 'log10(Relative abundance)', titleLab = yDiv, if_logged = T, P.pos = c(10, 2))
  lmePlot(lmeData = lmeDF, yLab = 'log10(Relative abundance)', titleLab = yDiv,
          if_logged = T, P.pos = c(10,1), groupV = 'Term2',addSE = T,
          expandMulti = expansion(mult = c(0.05, 0.1)), yPos = c(1.7,2))
  lmePlot(lmeData = lmeDF, yLab = 'log10(Relative abundance)', titleLab = yDiv,
          if_logged = T, P.pos = c(10,1), groupV = 'Ethnicity2',addSE = T,
          expandMulti = expansion(mult = c(0.05, 0.1)), yPos = c(2.0,2.3))
  yDiv = 'Lactobacillus vaginalis'
  lmeDF$y = log10(lmeDF[, yDiv] + 0.01)
  lmePlot(lmeData = lmeDF, yLab = 'log10(Relative abundance)', titleLab = yDiv, if_logged = T, P.pos = c(10, 2))
  lmePlot(lmeData = lmeDF, yLab = 'log10(Relative abundance)', titleLab = yDiv,
          if_logged = T, P.pos = c(10,1), groupV = 'Term2',addSE = T,
          expandMulti = expansion(mult = c(0.05, 0.1)), yPos = c(1.7,2))
  lmePlot(lmeData = lmeDF, yLab = 'log10(Relative abundance)', titleLab = yDiv,
          if_logged = T, P.pos = c(10,1), groupV = 'Ethnicity2',addSE = T,
          expandMulti = expansion(mult = c(0.05, 0.1)), yPos = c(2.0))
}
dev.off()
###### longitudinal - GAMM #####
library(voxel)
library(gamm4)

gammPlot = function(gammDF, groupCovs, x.Pos = 10, y.Pos = 0, titleLab = NULL){
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
    labs(x = 'Gestational Weeks', y = 'log10(Relative abundance)', title = titleLab) +
    annotate(x =x.Pos, y = y.Pos, label = adoL_p, geom = "text", size = 4.5, parse = TRUE, hjust = 0) +
    mytheme
  print(gammP2)
}

cairo_pdf('./plot/gamm_abundance.pdf', width = 4.66, height = 3.26, onefile = T)
for(yDiv in names(colorSpecies)[1:15]){
  gammDF = refDF0[refDF0$Trimester != 'P', ]
  gammDF$y = log10(gammDF[, yDiv] + 0.01)
  gammPlot(gammDF = gammDF,groupCovs = 'Term2', x.Pos = 10, y.Pos = c(2.2, 2.6), titleLab = yDiv)
  gammPlot(gammDF = gammDF,groupCovs = 'Ethnicity2', x.Pos = 10, y.Pos = seq(2, 3.6,0.4)+0.2, titleLab = yDiv)
}
dev.off()

###### longitudinal - GBMT #####
refDF1 = refDF0[refDF0$Trimester != 'P', ]
refDF1$Trimester = factor(refDF1$Trimester, levels = c('T1', 'T2', 'T3'))
colnames(refDF1) %<>% str_replace_all(., ' ', '_')
varNames = c(colnames(refDF1[, 1:15]))
if(T){
  ng2 <- gbmt(x.names = varNames, unit="SampleID", time="Sample_GA", 
               d=2, ng=2, data=refDF1, scaling=0)
  ng3 <- gbmt(x.names = varNames, unit="SampleID", time="Sample_GA", 
               d=2, ng=3, data=refDF1, scaling=0)
  ng4 <- gbmt(x.names = varNames, unit="SampleID", time="Sample_GA", 
               d=2, ng=4, data=refDF1, scaling=0)
  ng5 <- gbmt(x.names = varNames, unit="SampleID", time="Sample_GA", 
              d=2, ng=5, data=refDF1, scaling=0)
  ng6 <- gbmt(x.names = varNames, unit="SampleID", time="Sample_GA", 
              d=2, ng=6, data=refDF1, scaling=0)
  ng7 <- gbmt(x.names = varNames, unit="SampleID", time="Sample_GA", 
              d=2, ng=7, data=refDF1, scaling=0)
  ng8 <- gbmt(x.names = varNames, unit="SampleID", time="Sample_GA", 
              d=2, ng=8, data=refDF1, scaling=0)
  ngIC = rbind(ng2$ic, ng3$ic, ng4$ic, ng5$ic, ng6$ic, ng7$ic, ng8$ic) %>% as.data.frame()
  row.names(ngIC) = paste0('ng', 2:8)
  plot(2:8, ngIC$bic, xlab = 'Number of groups to create', ylab = 'BIC', type = 'ol')
}
ng4$assign.list
gbmtDF = data.frame()
for(i in names(ng4$assign.list)){
  tmp = data.frame(
    Group = i,
    Person = ng4$assign.list[[i]]
  )
  gbmtDF = rbind(gbmtDF, tmp)
}

plot(ng4, 
     n.ahead=4,
     # transparency=80,
     # unit = "Italy",
     # group=1,
     # equal.scale=TRUE, trim=0.05, # same scale to ease comparisons 
     mar = c(3.1,2.55,3.1,1.2), # margin- bottom, left, top, right
) ## overlapped groups

refDF1$gbmtGroup = gbmtDF$Group[match(refDF1$SampleID, gbmtDF$Person)]
table(refDF1$CST, refDF1$gbmtGroup) 
table(refDF1$Vagitype, refDF1$gbmtGroup) 
table(refDF1$Ethnicity, refDF1$gbmtGroup) 
table(refDF1$Term, refDF1$gbmtGroup) 

pMat = refDF1[, 1:15]
annoRow = data.frame(
  Group = refDF1$gbmtGroup %>% as.character()
)
row.names(annoRow) = row.names(pMat)
pheatmap(pMat, scale = 'none', show_rownames = F, annotation_row = annoRow)


