source('~/workspace/library_function.R')
load('./data/metaData.RData')
load('./data/refDF.RData')
load('./data/magDF.RData')
load('./data/spName.RData')


#######################################################################################
#                              associated with some metadata                          #
#######################################################################################
#### single variable comparison lm & lme #####
refMetData = merge(refDF[!grepl('KIT', refDF$pt_ID), ], metaData, by = 'pt_ID', all.x = T)
magMetData = merge(magDF[!grepl('KIT', magDF$pt_ID), ], metaData, by = 'pt_ID', all.x = T)
allMetData = rbind(refMetData, magMetData)
allMetData$genome_species[allMetData$type == 'NCBI'] = paste0(
  allMetData$species, ' (NCBI)'
)[allMetData$type == 'NCBI']
allMetData$genome_species[allMetData$type == 'MAG'] = paste0(
  allMetData$species, ' (', allMetData$genome, ')'
)[allMetData$type == 'MAG']
plotFlag = F
if(T){
  ys = c('nucl_diversity', 'dnds', 'pnps', 'dn_ds_common')
  ysNorm = list(nucl_diversity = 'Nucleotide diversity', dnds = 'dN/dS', pnps = 'pN/pS', dn_ds_common = 'dN/dS (common)')
  xs = c('trimester', 'Age_char', 'BMI_char', 'Ethnicity', 'Language', 'Marriage', 'FOB', 'Employed', 'Housing', 
           'Term_char', 'Gravida_char', 'Abortions_char', 'Living_children_char', 
           'Abnormal_pap', 'Depression', 'PTSD', 'PPD_pos', 'DM2', 'Asthma', 'PCN_allergy', 'Sulfa_allergy', 'NKDA',
           'PNV', 'Fe', 'Insulin', 'B6', 'Progesterone',
           'IPV_hx', 'TED_hx', 'Tobacco_hx', 'EtOH_hx', 
           'GBS', 'PCN', 'Antibiotic',
           'Induction', 'Baby_gender', 'Baby_weight_char', 
           # repeated numeric 
           'Age', 'BMI', 'Term', 'Gravida', 'Abortions', 'Living_children', 'Baby_weight')
  comp_DF = data.frame()
  uniRes = data.frame()
  sp = unique(allMetData$genome_species)[1]
  x = xs[1]
  for (sp in unique(allMetData$genome_species)) {
    cat('sp =', sp, '\n')
    tmpDF = allMetData[allMetData$genome_species == sp, ]
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
        if(length(na.omit(unique(tmpDF2$x))) >= 2){
          # lm
          lmUni = lm(y ~ x, tmpDF2)
          lmUni.summary = summary(lmUni)
          lmUni.aov = anova(lmUni)
          # lmer
          formula_char = paste0(y, ' ~ (Intercept) + ', x, ' + (1|Subject)')
          lmeTry = try({lmer(as.formula(paste0('y ~ x + trimester + (1|pt_ID)')),
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
    xlsx::write.xlsx(comp_DF, file = './Routput/genome_wise_nucl_div_dnds_metadata_pairwise_comp_2.xlsx')  
    xlsx::write.xlsx(uniRes, file = './Routput/genome_wise_nucl_div_dnds_metadata_uni_model_2.xlsx')
  }
}

##### lme regression fixed model #####
# nucl div associations to metadata variables
colnames(allMetData)

xs = c(#'Insulin', 'B6', 'Employed',
       'Depression',
       'Term_char', 'Abortions', 
       'Ethnicity', 'BMI_char',  'Age', 'trimester')  # interest priority
ys = c('nucl_diversity', 'dnds', 'pnps', 'dn_ds_common')
x_random = c('pt_ID')

sp = "Lactobacillus_vaginalis (NCBI)" 
y = 'nucl_diversity'
lmeRes = data.frame()
for (sp in unique(allMetData$genome_species)) {
  cat('  sp =', sp, '\n')
  for (y in ys) {
    cat('  y =', y, '\n')
    dat = allMetData[allMetData$genome_species == sp, c(xs, x_random, y)]
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
res = merge(uniRes, lmeRes, by = c('species', 'y', 'x'), all = T)
res_mulisig = res[which(res$lme.multi.P_value < 0.05 & res$lme.uni.P_value < 0.05),]
res_mulisig = arrange(res_mulisig, species, y, x, lme.multi.P_value)
xlsx::write.xlsx(res_mulisig, './Routput/nucl_div_dnds_lme_regression_fixed.xlsx')


##### plot boxplot for candidate species and variable #####
cand_sp_var = res_mulisig[, c('species', 'y', 'x')]
if(plotFlag){
  pdf('./plot/genome_wise_nucl_div_dnds_metadata_super_sig.pdf', width = 10, height = 7)
}
for (pr in 1:nrow(cand_sp_var)) {
  sp = cand_sp_var[pr, ][1] %>% unlist
  y = cand_sp_var[pr, ][2] %>% unlist
  x = cand_sp_var[pr, ][3] %>% unlist
  tmpDF = allMetData[allMetData$genome_species == sp, ]
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

##### count associations res1 #####
# number of significant univariate association
dim(res)
res_filter = res[res$y %in% c('nucl_diversity', 'dn_ds_common') &
                   res$x %in% xs, ] # dn_ds_common, pnps, dnds
uniVarSigDF = res_filter[which(res_filter$lme.uni.P_value <= 0.05), ]
uniVarSigDF$x = str_replace_all(uniVarSigDF$x, '_char', '')
uniVarSigDF$lm.uni.P_value_log = -log10(uniVarSigDF$lm.uni.P_value)
uniVarSigDF$lme.uni.P_value_log = -log10(uniVarSigDF$lme.uni.P_value)
uniVarSigDF$lme.multi.P_value_log = -log10(uniVarSigDF$lme.multi.P_value)
uniqueVar = uniVarSigDF[, c('species', 'y', 'x', 'lm.uni.P_value_log', 'lme.uni.P_value_log',
                            'lme.multi.P_value_log', 'lme.multi.P_value')]
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
col_fun$col[uniqueVar$lme.multi.P_value < 0.05] = 'red' 
table(col_fun$col)
library(circlize)
pdf('./plot/circos_association_genitic_diversity_2.pdf', width = 14, height = 13)
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
