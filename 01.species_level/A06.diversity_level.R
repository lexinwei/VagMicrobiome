source('./script/library_function.R')
load('./data/RData/swabData.RData')
load('./data/RData/metaData.RData')
load('./data/RData/colorList.RData')
load('./data/RData/sampleVagitypes.RData')
##### diversity dir #####
divDir = './metaphlan4_results/diversity'
divFiles = list.files(divDir)
######## alpha diversity ######
if(T){
  alphaFileList = divFiles[startsWith(divFiles, 'alpha_')]
  alphaList = c('Dominance gini index', 'Richness', 'Shannon index', 'Simpson index')
  alphaDF0 = data.frame()
  for(a in alphaFileList){
    tmp0 = read.table(paste0('./metaphlan4_results/diversity/', a), sep = '\t')
    tmp0$seq_ID = str_remove(row.names(tmp0), 'Ming_nova_VS_') %>% str_remove(., '_bwa')
    if(nrow(alphaDF0) == 0){alphaDF0 = tmp0}else{alphaDF0 = merge(alphaDF0, tmp0)}
  }
  colnames(alphaDF0) = c('seq_ID', alphaList)
  alphaDF0$Trimester = swabData$trimester[match(alphaDF0$seq_ID, swabData$seq_ID)]
  alphaDF0$Term = swabData$Term[match(alphaDF0$seq_ID, swabData$seq_ID)]
  alphaDF0$Term2 = swabData$Term2[match(alphaDF0$seq_ID, swabData$seq_ID)]
  alphaDF0$Ethnicity = swabData$Ethnicity[match(alphaDF0$seq_ID, swabData$seq_ID)]
  alphaDF0$Ethnicity2 = swabData$Ethnicity2[match(alphaDF0$seq_ID, swabData$seq_ID)]
  alphaDF0$pt_ID = swabData$pt_ID[match(alphaDF0$seq_ID, swabData$seq_ID)]
  alphaDF0$Sample_GA = swabData$Sample_GA[match(alphaDF0$seq_ID, swabData$seq_ID)]
  alphaDF0$Housing = metaData$Housing[match(alphaDF0$pt_ID, metaData$pt_ID)]
  alphaDF0$Employed = metaData$Employed[match(alphaDF0$pt_ID, metaData$pt_ID)]
  alphaDF0$Depression = metaData$Depression[match(alphaDF0$pt_ID, metaData$pt_ID)]
  alphaDF0$FOB = metaData$FOB[match(alphaDF0$pt_ID, metaData$pt_ID)]
  alphaDF0$Marriage = metaData$Marriage[match(alphaDF0$pt_ID, metaData$pt_ID)]
  alphaDF0$BMI =  metaData$BMI[match(alphaDF0$pt_ID, metaData$pt_ID)]
  alphaDF0$Age =  metaData$Age[match(alphaDF0$pt_ID, metaData$pt_ID)]
  alphaDF0$Vagitype = sampleVagitypes$Vagitypes[match(alphaDF0$seq_ID, sampleVagitypes$seq_ID)]
  alphaDF0$CST = sampleVagitypes$CST[match(alphaDF0$seq_ID, sampleVagitypes$seq_ID)]
  swabData_sub = swabData %>% filter(Trimester != 'NC') 
  alphaDF1 = alphaDF0[alphaDF0$seq_ID %in% swabData_sub$seq_ID, ]
  alphaDF1$SampleID = swabData$SampleID.u[match(alphaDF1$seq_ID, swabData$seq_ID)]
  alphaDF1$Vagitype_full = sampleVagitypes$vagitypes_Fettweis[match(alphaDF1$seq_ID, sampleVagitypes$seq_ID)]
  table_alphaDF1 = alphaDF1[, c('SampleID', "Shannon index", 'CST','Vagitype', 'Vagitype_full')]
  # write.xlsx(table_alphaDF1, './table/alpha_CST_Vagitype.xlsx')
}

###### function ######
alphaDivPlot = function(data0, x, y, y.position = NULL, z = NULL, CN = NULL, cn = T, xLab, cldSize = 5, if_coord_flip = F,yPos = NULL,
                        x.angle = 0, CLD = F, hide.ns = F, legend.position = 'right', expandMulti = NULL, marG = NULL){
  data0$x = data0[, x]
  data0$y = data0[, y]
  data0$z = data0[, z]
  if(x.angle == 0){x.hjust = 0.5}else{x.hjust = 1}
  if(endsWith(x, '2')){lname = str_remove(x, '2$')}else{lname = x}
  if(is.null(CN)){CN = combn(na.omit(droplevels(unique(data0$x))) %>% levels, 2, simplify = FALSE)}
  if(!cn){
    stat.test <- data0 %>%
      group_by(z) %>%
      wilcox_test(y ~ x) %>%
      adjust_pvalue(method = "BH") %>%
      add_significance("p.adj")
    stat.test <- stat.test %>%
      add_xy_position(x = "x", dodge = 0)
    if(hide.ns){
      stat.test = stat.test[stat.test$p < 0.05, ]
      if(!is.null(y.position)){
        stat.test$y.position = y.position 
      }
    }
    stat.test = BBmisc::dropNamed(stat.test, c('p.adj', 'p.adj.signif'))
    pbO = ggboxplot(data0, x = x, y = y, xlab = lname,
                    ylab = y, color = x, 
                    add = "jitter", 
                    add.params = list(shape = 21, size = 1.8, color = x, fill = x, alpha = 0.7),
                    bxp.errorbar = T, notch = F) +
      # geom_smooth(method = "lm", se=T, aes(group=1), size = 1, color = '#F8B469') + 
      stat_pvalue_manual(
        stat.test,  label = "p", tip.length = 0.02, remove.bracket = F, hide.ns = hide.ns
      ) +
      scale_color_manual(values = get(paste0('color', x)), name = lname) +
      scale_fill_manual(values = get(paste0('color', x)), name = lname) +
      mytheme3 + theme(strip.text = element_text(size = 12),axis.text.x = element_text(angle = x.angle, hjust = x.hjust),
                       legend.position = legend.position)
  }else{
    pbO = ggboxplot(data0, x = x, y = y, xlab = lname,
                    ylab = y, color = x, 
                    add = "jitter", 
                    add.params = list(shape = 21, size = 1.8, color = x,
                                      fill = x, alpha = 0.7),
                    bxp.errorbar = T, notch = F) +
      # geom_smooth(method = "lm", se=T, aes(group=1), size = 1, color = '#F8B469') + 
      stat_compare_means(comparisons = CN, method = 'wilcox.test', label = '..p..') +
      scale_color_manual(values = get(paste0('color', x)), name = lname) +
      scale_fill_manual(values = get(paste0('color', x)), name = lname) +
      mytheme3 + theme(strip.text = element_text(size = 12),axis.text.x = element_text(angle = x.angle, hjust = x.hjust),
                      legend.position = legend.position)
  }
  if(!is.null(expandMulti)){
    pbO = pbO + scale_y_continuous(expand = expandMulti)
  }
  if(!is.null(marG)){
    pbO = pbO + theme(plot.margin = marG)
  }
  if((!is.null(z)) & (!CLD)){
    pbO = pbO + facet_grid(.~z)
    print(pbO)
  }
  if(is.null(z) & (!CLD)){
    print(pbO)
  }
  if(CLD){
    if(!is.null(z)){
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
    if(is.null(yPos)){
      yPos = max(data0$y) + (max(data0$y) - min(data0$y)) *0.025
    }

    pbO2 = ggboxplot(data0, x = x, y = y, xlab = lname,
                    ylab = y, color = x, 
                    add = "jitter", 
                    add.params = list(shape = 21, size = 1.8, color = x,
                                      fill = x, alpha = 0.7),
                    bxp.errorbar = T, notch = F) +
      # geom_smooth(method = "lm", se=T, aes(group=1), size = 1, color = '#F8B469') + 
      # stat_compare_means(comparisons = CN) +
      scale_color_manual(values = get(paste0('color', x)), name = lname) +
      scale_fill_manual(values = get(paste0('color', x)), name = lname) +
      # facet_grid(. ~ z) +
      geom_text(aes(x, y = yPos, 
                    label= str_trim(.group)),vjust = 0.5,
                data = model_means_cld, size = cldSize) +
      mytheme3 + theme(strip.text = element_text(size = 12), axis.text.x = element_text(angle = x.angle, hjust = x.hjust),
                       legend.position = legend.position)
    if((!is.null(z))){
      pbO2 = pbO2 + facet_grid(.~z)
    }
    if(!is.null(marG)){
      pbO2 = pbO2 + theme(plot.margin = marG)
    }
    if(!is.null(expandMulti)){
      pbO2 = pbO2 + scale_y_continuous(expand = expandMulti)
    }
    if(if_coord_flip){
      pbO2 = pbO2 + coord_flip()
    }
    print(pbO2)
  }
}
###### rough plot for investigate their relationshop ######
if(F){ 
  cairo_pdf('./plot/alpha_diversity_boxplot.pdf', width = 9.27, height = 5.75, onefile = T)
  alphaList = c('Richness', 'Shannon index', 'Simpson index')
  for(yDiv in alphaList){
    alphaDivPlot(data0 = alphaDF1, x = 'Term', y = yDiv)
    alphaDivPlot(data0 = alphaDF1, x = 'Term2', y = yDiv)
    alphaDivPlot(data0 = alphaDF1, x = 'Ethnicity', y = yDiv)
    alphaDivPlot(data0 = alphaDF1, x = 'Ethnicity2', y = yDiv)
    alphaDivPlot(data0 = alphaDF1, x = 'Trimester', y = yDiv)
    alphaDivPlot(data0 = alphaDF1, x = 'Housing', y = yDiv)
    alphaDivPlot(data0 = alphaDF1, x = 'Employed', y = yDiv)
  
    alphaDivPlot(data0 = alphaDF1, x = 'Term', y = yDiv, z = 'Trimester')
    alphaDivPlot(data0 = alphaDF1, x = 'Term2', y = yDiv, z = 'Trimester')
    alphaDivPlot(data0 = alphaDF1, x = 'Ethnicity', y = yDiv, z = 'Trimester', x.angle = 45)
    alphaDivPlot(data0 = alphaDF1, x = 'Ethnicity2', y = yDiv, z = 'Trimester', x.angle = 45)
    alphaDivPlot(data0 = alphaDF1, x = 'Housing', y = yDiv, z = 'Trimester', x.angle = 45)
    alphaDivPlot(data0 = alphaDF1, x = 'Employed', y = yDiv, x.angle = 45, z = 'Trimester')
    
    alphaDivPlot(data0 = alphaDF1, x = 'Trimester', y = yDiv, z = 'Term')
    alphaDivPlot(data0 = alphaDF1, x = 'Trimester', y = yDiv, z = 'Term2')
    alphaDivPlot(data0 = alphaDF1, x = 'Trimester', y = yDiv, z = 'Ethnicity')
    alphaDivPlot(data0 = alphaDF1, x = 'Trimester', y = yDiv, z = 'Ethnicity2')
    alphaDivPlot(data0 = alphaDF1, x = 'Trimester', y = yDiv, z = 'Housing')
    alphaDivPlot(data0 = alphaDF1, x = 'Trimester', y = yDiv, z = 'Employed')
    
    alphaDivPlot(data0 = alphaDF1, x = 'Ethnicity', y = yDiv, z = 'Term', x.angle = 45)
    alphaDivPlot(data0 = alphaDF1, x = 'Ethnicity', y = yDiv, z = 'Term2', x.angle = 45)
    alphaDivPlot(data0 = alphaDF1, x = 'Ethnicity2', y = yDiv, z = 'Term', x.angle = 45)
    alphaDivPlot(data0 = alphaDF1, x = 'Ethnicity2', y = yDiv, z = 'Term2', x.angle = 45)
  }
  dev.off()
}
##### lme #####
if(T){
  alphaDF1$Employment = metaData$Employed[match(alphaDF1$pt_ID, metaData$pt_ID)]
  alphaDF1$Housing = metaData$Housing[match(alphaDF1$pt_ID, metaData$pt_ID)]
  alphaDF1$Depression = metaData$Depression[match(alphaDF1$pt_ID, metaData$pt_ID)]
  alphaDF1$BMI = metaData$BMI[match(alphaDF1$pt_ID, metaData$pt_ID)]
  alphaDF1$Age = metaData$Age[match(alphaDF1$pt_ID, metaData$pt_ID)]
  alphaDF1$Sample_GA = swabData$Sample_GA[match(alphaDF1$pt_ID, swabData$pt_ID)]
  alphaDF1$trimester = swabData$trimester[match(alphaDF1$seq_ID, swabData$seq_ID)]
  alphaDF1$Marriage = metaData$Marriage[match(alphaDF1$pt_ID, metaData$pt_ID)]
  alphaDF1$Marriage %<>% as.character()
  alphaDF1$Marriage[is.na(alphaDF1$Marriage)] = 'Unknown'
  alphaDF1$Marriage %<>% factor(.,levels = c('Married', 'Unmarried', 'Unknown'))
  alphaDF1$FOB = metaData$FOB[match(alphaDF1$pt_ID, metaData$pt_ID)]
  xs = c(#'Insulin', 'B6', 
    'Employment', 'Housing',
    'Depression',
    'Term', # 'Abortions', 
    'Ethnicity2', 'trimester','FOB','Marriage',
    'BMI', # 'BMI_char', # 'trimester'
    'Age')  # interest priority
  lme.Multi.df = data.frame()
  for(y in c('Simpson index', 'Shannon index',"Richness")){
    dat = alphaDF1
    dat$y = alphaDF1[, y]
    lmeTry = lme(as.formula(paste0('y ~ ', str_c(xs, collapse = '+'))), random = ~ 1 | pt_ID,
                 data = dat, method = 'ML')
    lme.Multi = lmeTry
    lme.Multi.R2 = r.squaredGLMM(lme.Multi)[1]
    lme.Multi.aov = anova(lme.Multi)
    formula_char = paste0(y, ' ~ (Intercept) + ', str_c(capitalize(xs), collapse = ' + '), ' + (1|Subject)')
    lme.Multi.df_tmp = data.frame(y = y,
                                  x = row.names(lme.Multi.aov), 
                                  lme.multi.Model = formula_char,
                                  lme.multi.F_value = lme.Multi.aov$`F-value`,
                                  lme.multi.P_value = lme.Multi.aov$`p-value`,
                                  lme.multi.R2 = lme.Multi.R2)
    lme.Multi.df = rbind(lme.Multi.df, lme.Multi.df_tmp)
  }
  lme.Multi.df = lme.Multi.df[lme.Multi.df$x != '(Intercept)', ]
  openxlsx::write.xlsx(lme.Multi.df, './data/alpha_diversity_lme.Multi.df_v2.xlsx', rowNames =F)
  lme.Multi.df_sig = lme.Multi.df[lme.Multi.df$lme.multi.P_value < 0.05, ]
  # facet by trimester
  lme.Multi.df2 = data.frame()
  for(y in c('Simpson index', 'Shannon index',"Richness")){
    dat = alphaDF1
    dat$y = alphaDF1[, y]
    for(t in c('T1', 'T2', 'T3', 'P')){
      dat2 =  dat[dat$trimester == t, ]
      lmeTry = try({lme(as.formula(paste0('y ~ ', str_c(setdiff(xs, 'trimester'), collapse = '+'))), random = ~ 1 | pt_ID,
                    data = dat2, method = 'ML')})
      if(class(lmeTry) != 'try-error'){
        lme.Multi = lmeTry
        lme.Multi.R2 = r.squaredGLMM(lme.Multi)[1]
        lme.Multi.aov = anova(lme.Multi)
        formula_char = paste0(y, ' ~ (Intercept) + ', str_c(capitalize(xs), collapse = ' + '), ' + (1|Subject)')
        lme.Multi.df2_tmp = data.frame(trimester = t,
                                       y = y,
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
  openxlsx::write.xlsx(lme.Multi.df2, './data/alpha_diverisity_lme.Multi.df2_facet_by_trimester_v2.xlsx', rowNames =F)
  lme.Multi.df2_sig = lme.Multi.df2[lme.Multi.df2$lme.multi.P_value < 0.05, ]
  # facet by ethnicity
  lme.Multi.df3 = data.frame()
  for(y in c('Simpson index', 'Shannon index',"Richness")){
    dat = alphaDF1
    dat$y = alphaDF1[, y]
    for(t in c('White', 'Asian', 'Latina', 'Black', 'Other')){
      dat2 =  dat[dat$Ethnicity2 == t, ]
      lmeTry = try({lme(as.formula(paste0('y ~ ', str_c(setdiff(xs, c('Term', 'Housing', "FOB", "Marriage" ,'Ethnicity2')), collapse = '+'))), random = ~ 1 | pt_ID,
                        data = dat2, method = 'ML')})
      if(class(lmeTry) != 'try-error'){
        lme.Multi = lmeTry
        lme.Multi.R2 = r.squaredGLMM(lme.Multi)[1]
        lme.Multi.aov = anova(lme.Multi)
        formula_char = paste0(y, ' ~ (Intercept) + ', str_c(capitalize(xs), collapse = ' + '), ' + (1|Subject)')
        lme.Multi.df3_tmp = data.frame(trimester = t,
                                       y = y,
                                       x = row.names(lme.Multi.aov), 
                                       lme.multi.Model = formula_char,
                                       lme.multi.F_value = lme.Multi.aov$`F-value`,
                                       lme.multi.P_value = lme.Multi.aov$`p-value`,
                                       lme.multi.R2 = lme.Multi.R2)
        lme.Multi.df3 = rbind(lme.Multi.df3, lme.Multi.df3_tmp)
      }
    }
  }
  lme.Multi.df3 = lme.Multi.df3[lme.Multi.df3$x != '(Intercept)', ]
  openxlsx::write.xlsx(lme.Multi.df3, './data/alpha_diverisity_lme.Multi.df3_facet_by_ethnicity_v2.xlsx', rowNames =F)
  lme.Multi.df3_sig = lme.Multi.df3[lme.Multi.df3$lme.multi.P_value < 0.05, ]
  # facet by term2
  lme.Multi.df4 = data.frame()
  for(y in c('Simpson index', 'Shannon index',"Richness")){
    dat = alphaDF1
    dat$y = alphaDF1[, y]
    for(t in c('Preterm', 'Full-term')){
      dat2 =  dat[dat$Term2 == t, ]
      lmeTry = try({lme(as.formula(paste0('y ~ ', str_c(setdiff(xs, c('FOB','Marriage', 'Depression', 'Term2', 'Term')), collapse = '+'))), random = ~ 1 | pt_ID,
                        data = dat2, method = 'ML')})
      if(class(lmeTry) != 'try-error'){
        lme.Multi = lmeTry
        lme.Multi.R2 = r.squaredGLMM(lme.Multi)[1]
        lme.Multi.aov = anova(lme.Multi)
        formula_char = paste0(y, ' ~ (Intercept) + ', str_c(capitalize(xs), collapse = ' + '), ' + (1|Subject)')
        lme.Multi.df4_tmp = data.frame(trimester = t,
                                       y = y,
                                       x = row.names(lme.Multi.aov), 
                                       lme.multi.Model = formula_char,
                                       lme.multi.F_value = lme.Multi.aov$`F-value`,
                                       lme.multi.P_value = lme.Multi.aov$`p-value`,
                                       lme.multi.R2 = lme.Multi.R2)
        lme.Multi.df4 = rbind(lme.Multi.df4, lme.Multi.df4_tmp)
      }
    }
  }
  lme.Multi.df4 = lme.Multi.df4[lme.Multi.df4$x != '(Intercept)', ]
  openxlsx::write.xlsx(lme.Multi.df4, './data/alpha_diverisity_lme.Multi.df4_facet_by_term2_v2.xlsx', rowNames =F)
  lme.Multi.df4_sig = lme.Multi.df4[lme.Multi.df4$lme.multi.P_value < 0.05, ]
}
###### re plot for good results ######
if(T){
  lme.Multi.df = read.xlsx('./data/alpha_diversity_lme.Multi.df_v2.xlsx')
  lme.Multi.df_sig = lme.Multi.df[lme.Multi.df$lme.multi.P_value < 0.05, ]
  cairo_pdf('./plot/alpha_diversity_boxplot_Employment.pdf', width = 4.59, height = 3.49, onefile = T)
  alphaDivPlot(data0 = alphaDF1, x = 'Employment', y = 'Shannon index', 
               legend.position = 'none',cn=F, hide.ns = T, 
               y.position = seq(3,0.3*4+3, 0.3))
  dev.off()
  cairo_pdf('./plot/alpha_diversity_boxplot_FOB.pdf', width = 2.66, height = 3.59, onefile = T)
  alphaDivPlot(data0 = alphaDF1, x = 'FOB', y = 'Shannon index', legend.position = 'none')
  # alphaDivPlot(data0 = alphaDF1, x = 'Ethnicity2', y = 'Shannon index', legend.position = 'none')
  dev.off()
}
if(T){
  cairo_pdf('./plot/alpha_diversity_boxplot_by_term.pdf', 
            width = 2.66, height = 3.59, onefile = T)
  alphaDivPlot(data0 = alphaDF1, x = 'Term', y = 'Shannon index', legend.position = 'none')
  alphaDivPlot(data0 = alphaDF1, x = 'Term', y = 'Simpson index', legend.position = 'none')
  alphaDivPlot(data0 = alphaDF1, x = 'Term2', y = 'Shannon index', legend.position = 'none')
  alphaDivPlot(data0 = alphaDF1, x = 'Term2', y = 'Simpson index', legend.position = 'none')
  dev.off()
}
if(T){
  cairo_pdf('./plot/alpha_diversity_boxplot_by_Ethnicity2.pdf', width = 3.47, height = 3.59, onefile = T)
  alphaDivPlot(data0 = alphaDF1, x = 'Ethnicity2', y = 'Shannon index', legend.position = 'none', cn=F, hide.ns = T, y.position = c(3,3.3,3.6,3.9))
  alphaDivPlot(data0 = alphaDF1, x = 'Ethnicity2', y = 'Simpson index', legend.position = 'none', cn=F, hide.ns = T, y.position = c(1,1.1,1.2,1.3))
  alphaDivPlot(data0 = alphaDF1, x = 'Ethnicity2', y = 'Shannon index', legend.position = 'none', CLD = T)
  alphaDivPlot(data0 = alphaDF1, x = 'Ethnicity2', y = 'Simpson index', legend.position = 'none', CLD = T)
  dev.off()
}

if(T){
  cairo_pdf('./plot/alpha_diversity_boxplot_by_term_facetBy_trimester.pdf', width = 5.03, height = 3.54, onefile = T)
  alphaDivPlot(data0 = alphaDF1, x = 'Term', y = 'Simpson index', z = 'Trimester', x.angle = 45, cn = F,
               expandMulti = expansion(mult = c(0.03, 0.1)))
  alphaDivPlot(data0 = alphaDF1, x = 'Term', y = 'Shannon index', z = 'Trimester', x.angle = 45, cn = F,
               expandMulti = expansion(mult = c(0.03, 0.1)))
  alphaDivPlot(data0 = alphaDF1, x = 'Term2', y = 'Shannon index', z = 'Trimester', x.angle = 45, cn = F, 
               expandMulti = expansion(mult = c(0.03, 0.1)))
  alphaDivPlot(data0 = alphaDF1, x = 'Term2', y = 'Simpson index', z = 'Trimester', x.angle = 45, cn = F,
               expandMulti = expansion(mult = c(0.03, 0.1)))
  dev.off()
}
if(T){
  cairo_pdf('./plot/alpha_diversity_boxplot_by_ethnicity_facetBy_trimester.pdf', width = 6.43, height = 3.54, onefile = T)
  alphaDivPlot(data0 = alphaDF1, x = 'Ethnicity2', y = 'Shannon index', z = 'Trimester', 
               x.angle = 45, cn = F, hide.ns = T, expandMulti = expansion(mult = c(0.03, 0.1)),
               y.position = c(3,3.4,3,3.4,3.8,4.2))
  alphaDivPlot(data0 = alphaDF1, x = 'Ethnicity2', y = 'Simpson index', z = 'Trimester', 
               x.angle = 45, cn = F, hide.ns = T, expandMulti = expansion(mult = c(0.03, 0.1)),
               y.position = c(1,1,1.15,1.3,1.45))
  alphaDivPlot(data0 = alphaDF1, x = 'Ethnicity2', y = 'Richness', z = 'Trimester', 
               x.angle = 45, cn = F, hide.ns = T,expandMulti = expansion(mult = c(0.03, 0.1)),
               y.position = c(90,110,130, 110, 130, 150, 170, 140, 160))
  dev.off()
}
if(T){
  cairo_pdf('./plot/alpha_diversity_boxplot_by_ethnicity_facetBy_term.pdf', width = 6.43, height = 3.54, onefile = T)
  alphaDivPlot(data0 = alphaDF1, x = 'Ethnicity2', y = 'Shannon index', z = 'Term', x.angle = 45, cn = F, hide.ns = T,
               y.position = c(2.8, 3.1, 3.4), expandMulti = expansion(mult = c(0.03,0.08)))
  alphaDivPlot(data0 = alphaDF1, x = 'Ethnicity2', y = 'Simpson index', z = 'Term', x.angle = 45, cn = F, hide.ns = T,
               y.position = c(1, 1.1), expandMulti = expansion(mult = c(0.03,0.08)))
  alphaDivPlot(data0 = alphaDF1, x = 'Ethnicity2', y = 'Shannon index', z = 'Term2', x.angle = 45, cn = F, hide.ns = T,
               y.position = c(2.8, 3.1, 3.4), expandMulti = expansion(mult = c(0.03,0.08)))
  alphaDivPlot(data0 = alphaDF1, x = 'Ethnicity2', y = 'Simpson index', z = 'Term2', x.angle = 45, cn = F, hide.ns = T,
               y.position = c(1, 1.1), expandMulti = expansion(mult = c(0.03,0.08)))
  dev.off()
}
if(T){
  cairo_pdf('./plot/alpha_diversity_boxplot_by_term_facetBy_ethnicity.pdf', width = 6.43, height = 3.54, onefile = T)
  alphaDivPlot(data0 = alphaDF1, x = 'Term2', y = 'Shannon index', z = 'Ethnicity2', x.angle = 45, 
               y.position = c(2.8, 3.1, 3.4), expandMulti = expansion(mult = c(0.03,0.08)))
  alphaDivPlot(data0 = alphaDF1, x = 'Term2', y = 'Simpson index', z = 'Ethnicity2', x.angle = 45, 
               y.position = c(1, 1.1), expandMulti = expansion(mult = c(0.03,0.08)))
  dev.off()
}

if(T){
  cairo_pdf('./plot/alpha_diversity_boxplot_by_trimester_facetBy_term.pdf', width = 6.43, height = 3.54, onefile = T)
  alphaDivPlot(data0 = alphaDF1, x = 'Trimester', y = 'Shannon index', z = 'Term2', 
               y.position = c(2.8, 3.1, 3.4), expandMulti = expansion(mult = c(0.03,0.08)))
  alphaDivPlot(data0 = alphaDF1, x = 'Trimester', y = 'Simpson index', z = 'Term2', 
               y.position = c(1, 1.1), expandMulti = expansion(mult = c(0.03,0.08)))
  dev.off()
}
if(T){
  cairo_pdf('./plot/alpha_diversity_boxplot_by_trimester_facetBy_ethnicity.pdf', width = 6.43, height = 3.54, onefile = T)
  alphaDivPlot(data0 = alphaDF1, x = 'Trimester', y = 'Shannon index', z = 'Ethnicity2', 
               y.position = c(2.8, 3.1, 3.4), expandMulti = expansion(mult = c(0.03,0.08)))
  alphaDivPlot(data0 = alphaDF1, x = 'Trimester', y = 'Simpson index', z = 'Ethnicity2', 
               y.position = c(1, 1.1), expandMulti = expansion(mult = c(0.03,0.08)))
  dev.off()
}
##### alpha_diversity_boxplot_by_vagitype ######
if(T){
  cairo_pdf('./plot/alpha_diversity_boxplot_by_vagitype.pdf', width = 5.19, height = 4.20, onefile = T)
  alphaDivPlot(data0 = alphaDF1, x = 'Vagitype', y = 'Shannon index', x.angle = 45,CLD = T,z = NULL, legend.position = 'none',
               marG = margin(0.5,0.5,0.5,2, 'cm'), cldSize = 4.5, expandMulti = expansion(mult = c(0.05, 0.2)))
  alphaDivPlot(data0 = alphaDF1, x = 'Vagitype', y = 'Simpson index', x.angle = 45,CLD = T,z = NULL, legend.position = 'none',
               marG = margin(0.5,0.5,0.5,2, 'cm'), cldSize = 4.5, expandMulti = expansion(mult = c(0.05, 0.2)))
  alphaDivPlot(data0 = alphaDF1, x = 'Vagitype', y = 'Richness', x.angle = 45,CLD = T,z = NULL, legend.position = 'none',
               marG = margin(0.5,0.5,0.5,2, 'cm'), cldSize = 4.5, expandMulti = expansion(mult = c(0.05, 0.2)))
  dev.off()
}
##### alpha_diversity_boxplot_by_vagitype t ######
if(T){
  cairo_pdf('./plot/alpha_diversity_boxplot_by_vagitype.pdf', width = 5.19, height = 4.20, onefile = T)
  cairo_pdf('./plot/alpha_diversity_boxplot_by_vagitype_T.pdf', width = 5.70, height = 2.64, onefile = T)
  alphaDF1$Vagitype = factor(alphaDF1$Vagitype, levels = rev(levels(alphaDF1$Vagitype)))
  alphaDivPlot(data0 = alphaDF1, x = 'Vagitype', y = 'Shannon index',CLD = T,z = NULL, legend.position = 'none',if_coord_flip = T,yPos = 3.1,
               marG = margin(0.5,0.5,0.5,2, 'cm'), cldSize = 4.5, expandMulti = expansion(mult = c(0.05, 0.1)))
  dev.off()
}

if(T){
  cairo_pdf('./plot/alpha_diversity_boxplot_by_cst.pdf', width = 2.73, height = 2.81, onefile = T)
  alphaDivPlot(data0 = alphaDF1, x = 'CST', y = 'Shannon index',CLD = T,z = NULL, legend.position = 'none',
               cldSize = 4.5, expandMulti = expansion(mult = c(0.05, 0.2)))
  alphaDivPlot(data0 = alphaDF1, x = 'CST', y = 'Simpson index',CLD = T,z = NULL, legend.position = 'none',
               cldSize = 4.5, expandMulti = expansion(mult = c(0.05, 0.2)))
  alphaDivPlot(data0 = alphaDF1, x = 'CST', y = 'Richness',CLD = T,z = NULL, legend.position = 'none',
               cldSize = 4.5, expandMulti = expansion(mult = c(0.05, 0.2)))
  dev.off()
}
###### alpha_diversity one sample each trimester #####
if(F){
  cairo_pdf('./plot/alpha_diversity_boxplot_sub.pdf', width = 9.27, height = 5.75, onefile = T)
  alphaList = c('Richness', 'Shannon index', 'Simpson index')
  swabData_sub = swabData %>% filter(Trimester != 'NC') %>% 
    group_by(pt_ID, trimester) %>% arrange(Sample_GA) %>% filter(row_number()==1)
  alphaDF1 = alphaDF0[alphaDF0$seq_ID %in% swabData_sub$seq_ID, ]
  for(yDiv in alphaList){
    alphaDivPlot(data0 = alphaDF1, x = 'Term', y = yDiv)
    alphaDivPlot(data0 = alphaDF1, x = 'Term2', y = yDiv)
    alphaDivPlot(data0 = alphaDF1, x = 'Ethnicity', y = yDiv)
    alphaDivPlot(data0 = alphaDF1, x = 'Ethnicity2', y = yDiv)
    alphaDivPlot(data0 = alphaDF1, x = 'Trimester', y = yDiv)
    alphaDivPlot(data0 = alphaDF1, x = 'Housing', y = yDiv)
    alphaDivPlot(data0 = alphaDF1, x = 'Employed', y = yDiv)
    
    alphaDivPlot(data0 = alphaDF1, x = 'Term', y = yDiv, z = 'Trimester')
    alphaDivPlot(data0 = alphaDF1, x = 'Term2', y = yDiv, z = 'Trimester')
    alphaDivPlot(data0 = alphaDF1, x = 'Ethnicity', y = yDiv, z = 'Trimester', x.angle = 45)
    alphaDivPlot(data0 = alphaDF1, x = 'Ethnicity2', y = yDiv, z = 'Trimester', x.angle = 45)
    alphaDivPlot(data0 = alphaDF1, x = 'Housing', y = yDiv, z = 'Trimester', x.angle = 45)
    alphaDivPlot(data0 = alphaDF1, x = 'Employed', y = yDiv, x.angle = 45, z = 'Trimester')
  }
  dev.off()
}
table(metaData$Ethnicity, metaData$Term_char)
table(metaData$Ethnicity, metaData$Term_char2)
table(metaData$Ethnicity2, metaData$Term_char)
table(metaData$Ethnicity, metaData$Term_char2)

###### longitudinal - lm #####
if(T){
  cairo_pdf('./plot/alpha_diversity_boxplot_lm.pdf', width = 4.29, height = 2.69, onefile = T)
  gammDF = alphaDF1[alphaDF1$Trimester != 'P', ]
  gammDF$diversity = gammDF$`Shannon index`
  gammDF$diversity_log = log(gammDF$diversity)
  
  ap1 = ggplot(alphaDF1[alphaDF1$Trimester != 'P', ], aes(x = Sample_GA, y = `Shannon index`)) +
    geom_point(shape = 21, aes(color = Ethnicity2, fill = Ethnicity2), alpha = 0.5) +
    geom_smooth(aes(color = Ethnicity2), method = 'lm', se = F, name = 'Ethnicity') +
    scale_color_manual(values = colorEthnicity2, name = 'Ethnicity') +
    scale_fill_manual(values = colorEthnicity2, name = 'Ethnicity') +
    labs(x = 'Gestational Weeks') +
    mytheme + theme(aspect.ratio = 1)
  ap2 = ggplot(alphaDF1[alphaDF1$Trimester != 'P', ], aes(x = Sample_GA, y = `Shannon index`)) +
    geom_point(shape = 21, aes(color = Term2, fill = Term2), alpha = 0.5) +
    geom_smooth(aes(color = Term2, fill = Term2), method = 'lm', se = T, name = 'Term', alpha = 0.2) +
    scale_color_manual(values = colorTerm2, name = 'Term') +
    scale_fill_manual(values = colorTerm2, name = 'Term') +
    labs(x = 'Gestational Weeks') +
    mytheme + theme(aspect.ratio = 1)
  ap3 = ggplot(alphaDF1[alphaDF1$Trimester != 'P', ], aes(x = Sample_GA, y = `Simpson index`)) +
    geom_point(shape = 21, aes(color = Ethnicity2, fill = Ethnicity2), alpha = 0.5) +
    geom_smooth(aes(color = Ethnicity2), method = 'lm', se = F, name = 'Ethnicity') +
    scale_color_manual(values = colorEthnicity2, name = 'Ethnicity') +
    scale_fill_manual(values = colorEthnicity2, name = 'Ethnicity') +
    labs(x = 'Gestational Weeks') +
    mytheme + theme(aspect.ratio = 1)
  ap4 = ggplot(alphaDF1[alphaDF1$Trimester != 'P', ], aes(x = Sample_GA, y = `Simpson index`)) +
    geom_point(shape = 21, aes(color = Term2, fill = Term2), alpha = 0.5) +
    geom_smooth(aes(color = Term2, fill = Term2), method = 'lm', se = T, name = 'Term', alpha = 0.2) +
    scale_color_manual(values = colorTerm2, name = 'Term') +
    scale_fill_manual(values = colorTerm2, name = 'Term') +
    labs(x = 'Gestational Weeks') +
    mytheme + theme(aspect.ratio = 1)
  print(ap1)
  print(ap2)
  print(ap3)
  print(ap4)
  dev.off()
}
###### longitudinal - GAMM #####
library(voxel)
library(gamm4)

if(T){
  alphaDF1$Ethnicity2 %<>% as.character()
  alphaDF1$Ethnicity2[  alphaDF1$Ethnicity2 == 'Latina'] = 'Latino'
  alphaDF1$Ethnicity2 %<>% factor(., levels = names(colorEthnicity2))
  colorTerm2 = c('Early-term' = "#DF8F44", 'Full-term' = "#374E55")
  alphaDF1$Term2 %<>% as.character()
  alphaDF1$Term2[  alphaDF1$Term2 == 'Preterm'] = 'Early-term' 
  alphaDF1$Term2 %<>% factor(., levels = names(colorTerm2))
  table(  alphaDF1$Term2 )
  gammDF = alphaDF1[alphaDF1$Trimester != 'P', ]
  gammDF$diversity = gammDF$`Shannon index`
  gammDF$diversity = gammDF$`Simpson index`
  gammDF$diversity_log = log(gammDF$diversity)
  gammDF$y = gammDF$diversity_log
  
  gamm1 <- gamm4::gamm4(y ~ Ethnicity2 + Employed + Housing + Term2 + BMI + Marriage + Age +  FOB + Depression + s(Sample_GA, by = Term2), 
                        data=gammDF, random =  ~ (1 | pt_ID))
  pval = summary(gamm1$gam)$s.table[, 4] %>% signif(., 2)
  names(pval) %<>% str_remove(., 's\\(.*\\)\\:') %>% str_remove('Term2')
  adoL_p = paste0('"', names(pval), '"',  "~~italic(P) == ", pval)
  gammP1 <- plotGAMM_2(gammFit <- gamm1, smooth.cov <- "Sample_GA", groupCovs = "Term2",
                     plotCI <- T, rawOrFitted = "raw", grouping = "pt_ID") +
    scale_color_manual(values = colorTerm2, name = 'Term') + 
    scale_fill_manual(values = colorTerm2, name = 'Term') +
    scale_y_continuous(breaks = seq(-8,4,2)) +
    labs(x = 'Gestational Weeks', y = 'log(Shannon index)', title = NULL) +
    annotate(x =c(40, 40), y = c(-7, -8), label = adoL_p, geom = "text", size = 4.5, parse = TRUE, hjust = 1) +
    mytheme + theme(legend.position = 'top', legend.margin = margin(0,0,-7,0))
  cairo_pdf('./plot/gamm_alpha_diversity_term_v3.pdf', width = 3.39, height =3.58, onefile = T)
  print(gammP1)
  dev.off()
  gamm2 <- gamm4::gamm4(y ~Ethnicity2 + Employed + Housing + Term + BMI + Marriage + Age + FOB + Depression + s(Sample_GA, by = Ethnicity2), 
                        data=gammDF, random =  ~ (1 | pt_ID))
  pval = summary(gamm2$gam)$s.table[, 4] %>% signif(., 2)
  names(pval) %<>% str_remove(., 's\\(.*\\)\\:') %>% str_remove('Ethnicity2')
  adoL_p = paste0('"', names(pval), '"',  "~~italic(P) == ", pval)
  gammP2 <- plotGAMM_2(gammFit <- gamm2, smooth.cov <- "Sample_GA", groupCovs = "Ethnicity2",
                      plotCI <- T, rawOrFitted = "raw", grouping = "pt_ID", plotCI = F) +
    scale_color_manual(values = colorEthnicity2, name = 'Ethnicity') +
    scale_fill_manual(values = colorEthnicity2, name = 'Ethnicity') +
    labs(x = 'Gestational Weeks', y = 'log(Shannon index)', title = NULL) +
    annotate(x =rep(40, 5), y = seq(-8, -8+4*0.7, 0.7) %>% rev(), label = adoL_p, geom = "text", size = 4.5, parse = TRUE, hjust = 1) +
    mytheme + theme(legend.position = 'right', legend.margin = margin(0,0,-7,0), legend.spacing.x = unit(0.05, 'cm'))
  cairo_pdf('./plot/gamm_alpha_diversity_ethnicity_v3.pdf', width = 4.32, height =3.41, onefile = T)
  print(gammP2)
  dev.off()
 
}
dev.off()
###### longitudinal lme , Susan #####
# function
lmePlot = function(lmeData, yLab, if_logged = T, P.pos = c(10,2.5), expandMulti = NULL,yPos = 2,groupV = NULL, addSE = T, linewith = 0.8){
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
    new.dat$hi2SE <- exp(new.dat$pred + 2*new.dat$SE)
    new.dat$lo2SE <- exp(new.dat$pred - 2*new.dat$SE)
    new.dat$pred <- exp(new.dat$pred)
    # Reverse log-transform for the vaginal community
    pval = summary(lmeFit)$tTable["Sample_GA", "p-value"] %>% signif(., 2)
    adoL_p = paste0("~~italic(P) == ", pval)
    if(if_logged){
      lmeData$y %<>% exp()
    }
    lP = ggplot(new.dat, aes(x = Sample_GA, y = pred)) + 
      geom_line(color = "red", size = linewith) +
      geom_point(data = lmeData, aes(x = Sample_GA, y = y), size=1, shape = 21, color = '#8FB7D3', fill = '#8FB7D3', alpha = 0.5) +
      geom_ribbon(aes(ymin = lo2SE,ymax = hi2SE), alpha=0.2, fill="red") +
      labs(x = "Gestational Weeks", y = yLab) +
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
        new.dat_tmp$hi2SE <- exp(new.dat_tmp$pred + 2*new.dat_tmp$SE)
        new.dat_tmp$lo2SE <- exp(new.dat_tmp$pred - 2*new.dat_tmp$SE)
        new.dat_tmp$pred <- exp(new.dat_tmp$pred)
        new.dat_tmp$Group = i
        # Reverse log-transform for the vaginal community
        pval = c(pval, summary(lmeFitTmp)$tTable["Sample_GA", "p-value"] %>% signif(., 2))
        new.dat = rbind(new.dat, new.dat_tmp)
        groupR = c(groupR, i)

      }
 
    }
    if(if_logged){
      lmeData$y %<>% exp()
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
      labs(x = "Gestational Weeks", y = yLab) +
      annotate("text", x = P.pos[1], y = yPos, size=4, label = adoL_p, parse = T, hjust = 0) +
      mytheme + theme(aspect.ratio = 1)
  }
  if(!is.null(expandMulti)){
    lP = lP + scale_y_continuous(expand = expandMulti)
  }
  print(lP)
}

cairo_pdf('./plot/lme_alpha_diversity.pdf', width = 4.29, height = 2.69, onefile = T)
if(T){
  lmeDF = alphaDF1[alphaDF1$Trimester != 'P', ]
  lmeDF$y = log(lmeDF$`Shannon index`)
  lmePlot(lmeData = lmeDF, yLab = 'Shannon index', if_logged = T)
  lmePlot(lmeData = lmeDF, yLab = 'Shannon index', if_logged = T, P.pos = c(10,1), groupV = 'Term2',addSE = T,
          expandMulti = expansion(mult = c(0.05, 0.1)), yPos = c(2.7, 3))
  lmePlot(lmeData = lmeDF, yLab = 'Shannon index', if_logged = T,addSE = F, P.pos = c(10,1), groupV = 'Ethnicity2',
          expandMulti = expansion(mult = c(0.05, 0.1)), yPos = c(2.0,2.3,2.6,2.9))
  
  lmeDF$y = log(lmeDF$`Simpson index`)
  lmePlot(lmeData = lmeDF, yLab = 'Simpson index', if_logged = T, P.pos = c(10,1), expandMulti = expansion(mult = c(0.05, 0.1)))
  lmePlot(lmeData = lmeDF, yLab = 'Simpson index', if_logged = T, P.pos = c(10,1), groupV = 'Term2',
          expandMulti = expansion(mult = c(0.05, 0.1)), yPos = c(1, 1.1))
  lmePlot(lmeData = lmeDF, yLab = 'Simpson index', if_logged = T,addSE = F, P.pos = c(10,1), groupV = 'Ethnicity2',
          expandMulti = expansion(mult = c(0.05, 0.1)), yPos = c(0.76, 0.84, 0.92, 1))
  
}
dev.off()





