source('~/workspace/library_function.R')
load('./data/metaData.RData')
##### mags sample-wise profile ##### 
ided_magsDF = data.frame()
for(ss in swabData$fqID){
  tmpDat = read.delim2(paste0('./instrain_res/mags_profile/genome_info/', ss, '_bwa.IS_genome_info.tsv'))
  tmpDat = tmpDat[, c("genome", 'length', "coverage", "breadth", "nucl_diversity", "length",
                      "true_scaffolds", "detected_scaffolds", "coverage_median",
                      "coverage_std", "coverage_SEM", "breadth_minCov",
                      "breadth_expected", "nucl_diversity_rarefied",
                      "conANI_reference", "popANI_reference", "SNS_count",
                      "SNV_count", "reads_mean_PID", "reads_unfiltered_reads")]
  tmpDat$breadth_minCov %<>% as.numeric()
  ided_magsDF_tmp = tmpDat[tmpDat$breadth_minCov >= 0.1, ]
  if(nrow(ided_magsDF_tmp) > 0){
    ided_magsDF_tmp$sample = ss
    ided_magsDF = rbind(ided_magsDF, ided_magsDF_tmp)
  }
}
ided_magsDF$pt_ID = swabData$pt_ID[match(ided_magsDF$sample, swabData$fqID)]
ided_magsDF$pt_ID.u = swabData$pt_ID.u[match(ided_magsDF$sample, swabData$fqID)]
ided_magsDF$genome %>% unique() %>% length()
ided_magsDF$genome %>% table()
ided_magsDF$sample %>% unique() %>% length()
ided_magsDF$pt_ID %>% unique() %>% length()
ided_magsDF$pt_ID %>% table()
ided_magsDF$sample %>% table()
ided_magsDF$nucl_diversity %<>% as.numeric()
ided_magsDF$pt_ID.u %<>% factor(., levels = mixedsort(ided_magsDF$pt_ID.u %>% unique))
ided_magsDF$genome %<>% str_replace(., '.fna', '')
spName = read.delim2('./tax_of_mags_co.txt', sep = '\t')
ided_magsDF$species = spName$lca_species[match(ided_magsDF$genome, spName$Bin)]
if(F){
  pdf('./plot/mags_genome_wise_nucl_div.pdf', width = 12, height = 4)
  for (sp in unique(ided_magsDF$genome)) {
    tmpDF = ided_magsDF[ided_magsDF$genome == sp, ]
    p = ggplot(tmpDF, aes(x = pt_ID.u, y = nucl_diversity, color = pt_ID)) +
      geom_line(aes(group = pt_ID)) + geom_point() +
      mytheme2 + theme(axis.text.x = element_text(angle = 45, size = 5, hjust = 1)) +
      labs(title = paste0(sp, ' | ', ided_magsDF$species[ided_magsDF$genome == sp]), y = paste0('Nucleotide diversity'), x = '')
    print(p)
  }
  dev.off()
}
#######################################################################################
#                              associated with some metadata                          #
#######################################################################################
#### single variable comparison #####
nucMetData = merge(ided_magsDF[!grepl('KIT', ided_magsDF$pt_ID), ], metaData, by = 'pt_ID', all.x = T)
plotFlag = F
if(T){
  if(plotFlag){
    pdf('./plot/ref_genome_wise_nucl_div_metadata_sig.pdf', width = 10, height = 7)
  }
  vars = c('Age_char', 'BMI_char', 'Ethinicity', 'Language', 'Marriage', 'FOB', 'Employed', 'Housing', 
           'Term_char', 'Gravida_char', 'Abortions_char', 'Living_children_char', 
           'Abnormal_pap', 'Depression', 'PTSD', 'PPD_pos', 'DM2', 'Asthma', 'PCN_allergy', 'Sulfa_allergy', 'NKDA',
           'PNV', 'Fe', 'Insulin', 'B6', 'Progesterone',
           'IPV_hx', 'TED_hx', 'Tobacco_hx', 'EtOH_hx', 
           'GBS', 'PCN', 'Antibiotic', 
           'Induction', 'Baby_gender', 'Baby_weight_char',
           # repeated numeric 
           'Age', 'BMI', 'Term', 'Gravida', 'Abortions', 'Living_children', 'Baby_weight')
  comp_DF = data.frame()
  aovRes = data.frame()
  sp = unique(nucMetData$genome)[1]
  va = vars[1]
  for (sp in unique(nucMetData$genome)) {
    tmpDF = nucMetData[nucMetData$genome == sp, ]
    for (va in vars) {
      tmpDF2 = data.frame(
        'fq_ID' = tmpDF$sample,
        'pt_ID' = tmpDF$pt_ID,
        'pt_ID.u' = tmpDF$pt_ID.u,
        'nucl_diversity' = tmpDF$nucl_diversity,
        'var' = tmpDF[, va]
      )
      if(class(tmpDF2$var) == 'numeric' & length(na.omit(unique(tmpDF2$var))) >= 2 & length(na.omit(unique(tmpDF2$var)))<length(na.omit(tmpDF2$var))){
        aM = lm(nucl_diversity ~ var, tmpDF2)
        asm = summary(aM)
        aMano = anova(aM)
        sv = c(sp, va, (aMano['var', c('F value', 'Pr(>F)')]) %>% unlist(), asm$r.squared %>% as.numeric())
        aovRes = rbind(aovRes, sv)
      }else if(class(tmpDF2$var) == 'factor' & 
               length(na.omit(unique(tmpDF2$var))) >= 2 & 
               nrow(tmpDF2) != length(unique(tmpDF2$var))){
        mm = aggregate(nucl_diversity ~ var, tmpDF2, mean)
        mmed = aggregate(nucl_diversity ~ var, tmpDF2, median)
        size = table(tmpDF2$var) %>% as.data.frame()
        aM = lm(nucl_diversity ~ var, tmpDF2)
        asm = summary(aM)
        aMano = anova(aM)
        sv = c(sp, va, (aMano['var', c('F value', 'Pr(>F)')]) %>% unlist(), asm$r.squared %>% as.numeric())
        aovRes = rbind(aovRes, sv)
        # do pairwise comparison
        tmp_comp_DF = compare_means(nucl_diversity ~ var, tmpDF2, p.adjust.method = 'BH') %>% as.data.frame()
        tmp_comp_DF$var = va
        tmp_comp_DF$species = sp
        tmp_comp_DF$size1 = size$Freq[match(tmp_comp_DF$group1, size$Var1)]
        tmp_comp_DF$size2 = size$Freq[match(tmp_comp_DF$group2, size$Var1)]
        
        tmp_comp_DF$mean1 = mm$nucl_diversity[match(tmp_comp_DF$group1, mm$var)]
        tmp_comp_DF$mean2 = mm$nucl_diversity[match(tmp_comp_DF$group2, mm$var)]
        
        tmp_comp_DF$log2_FC = log2(tmp_comp_DF$mean1/tmp_comp_DF$mean2)
        
        tmp_comp_DF$med1 = mmed$nucl_diversity[match(tmp_comp_DF$group1, mmed$var)]
        tmp_comp_DF$med2 = mmed$nucl_diversity[match(tmp_comp_DF$group2, mmed$var)]
        
        tmp_comp_DF = tmp_comp_DF[, c('species', 'var', 'group1', 'group2', 'size1', 'size2', 'mean1', 'mean2', 'log2_FC', 'med1', 'med2',
                                      'p', 'p.adj', 'p.format', 'p.signif')]
        comp_DF = rbind(comp_DF, tmp_comp_DF)
        #all pairs I want to compare
        if(is.factor(tmpDF2$var)){
          CN <- combn(na.omit(droplevels(unique(tmpDF2$var))) %>% levels, 2, simplify = FALSE)
        }else{
          CN <- combn(na.omit(unique(tmpDF2$var)), 2, simplify = FALSE)
        }
        #pairs I want to compare in list format for stat_compare_means
        if(plotFlag){
          if(sum(tmp_comp_DF$p.signif != 'ns') > 0){
            p = ggboxplot(na.omit(tmpDF2), x = 'var', y = 'nucl_diversity', color = 'var', add = 'jitter',
                          add.paramsshape = 'pt_ID') +
              mytheme3 + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
              stat_compare_means(comparisons=CN, method = "wilcox.test", p.adjust.method = "BH", aes(label=..p.adj..)) +
              labs(title = paste0(str_replace(sp, '.fna', ''), ' | ', va), 
                   color = va,
                   y = paste0('Nucleotide diversity'), x = va) + 
              scale_color_brewer(palette = 'Dark2')
            print(p)
          }
        }
      }
    }
  }
  if(plotFlag){
    dev.off()
    write.xlsx(comp_DF, file = './Routput/nucl_div_pairwise_comp.xlsx')  
  }
  colnames(aovRes) = c('species', 'var', 'lm.uni.F_value', 'lm.uni.P_value', 'lm.uni.R2')
  aovRes$lm.uni.F_value %<>% as.numeric()
  aovRes$lm.uni.P_value %<>% as.numeric()
}

##### lme regression #####
library(lme4)
library(lmerTest)
library(MuMIn)
# nucl div associations to metadata varibles
colnames(nucMetData)
vars_cand = c('genome', vars[!grepl('_char', vars)], 'pt_ID', 'nucl_diversity') # numeric version
vars_cand = c('genome', vars[1:36], 'pt_ID', 'nucl_diversity') # char version
lmeBestRes = data.frame()
lmeSimpleRes = data.frame()
for (sp in unique(nucMetData$genome)) {
  dat = nucMetData[nucMetData$genome == sp, vars_cand]
  sind = which(sapply(dat, function(col) length(unique(col))) > 1)
  if(nrow(dat) > 3 & length(sind) > 1){
    dat = dat[, sind]
    for (c in 1:ncol(dat)){
      if(is.factor(dat[, c])) {
        dat[,c] = droplevels(dat[,c])
      }else if(is.character(dat[,c])){
        dat[,c] = factor(dat[,c]) 
      }
    }
    lmemod=lme(nucl_diversity ~ NULL, random = ~1|pt_ID,
               data = dat, method = "ML")
    addFx = colnames(dat)[!colnames(dat) %in% c('genome', 'pt_ID', 'nucl_diversity')]
    for (vv in addFx){
      if(sum(is.na(dat[vv])) == 0){
        vvFormula = paste0('nucl_diversity ~ ', vv)
        result <- try({
          lme(as.formula(vvFormula), random=~1|pt_ID,
              data=dat,method="ML")
        }, silent = TRUE)
        if(class(result) == 'try-error'){
          print(paste('Error in MEEM', sp, vv))
        }else{
          vvM = result
          r_Sq1 = r.squaredGLMM(vvM)[1]
          vvaov = anova(vvM)
          vv_aov_df = data.frame(species = sp,
                                 var = row.names(vvaov), 
                                 lme.uni.Model = paste0('nucl_diversity ~ (Intercept)+', vv, '+(1|pt_ID)'),
                                 lme.uni.F_value = vvaov$`F-value`,
                                 lme.uni.P_value = vvaov$`p-value`,
                                 lme.uni.R2 = r_Sq1)
          lmeSimpleRes = rbind(lmeSimpleRes, vv_aov_df)
        }
      }
    }
    scopeList = paste0('~.+', paste(addFx, collapse = '+'))
    themod = stepAIC(lmemod,dir="forward", scope = scopeList)
    r_Sq2 = r.squaredGLMM(themod)[1]
    themod_aov = anova(themod)
    themod_aov_df = data.frame(species = sp,
                               var = row.names(themod_aov), 
                               lme.multi.Model = paste0('nucl_diversity ~ ', paste(row.names(themod_aov), collapse = '+'),'+(1|pt_ID)'),
                               lme.multi.F_value = themod_aov$`F-value`,
                               lme.multi.P_value = themod_aov$`p-value`,
                               lme.multi.R2 = r_Sq2)
    lmeBestRes = rbind(lmeBestRes, themod_aov_df)
  }
}

lmeBestRes_num = lmeBestRes[lmeBestRes$var != '(Intercept)' & (!is.na(lmeBestRes$lme.multi.P_value)), ]
lmeSimpleRes_num = lmeSimpleRes[lmeSimpleRes$var != '(Intercept)' & (!is.na(lmeSimpleRes$lme.uni.P_value)), ]

lmeBestRes_char = lmeBestRes[lmeBestRes$var != '(Intercept)' & (!is.na(lmeBestRes$lme.multi.P_value)), ]
lmeSimpleRes_char = lmeSimpleRes[lmeSimpleRes$var != '(Intercept)'& (!is.na(lmeSimpleRes$lme.uni.P_value)), ]

# merge
lmeBestResBind = rbind(lmeBestRes_num, lmeBestRes_char)
lmeSimpleResBind = rbind(lmeSimpleRes_num, lmeSimpleRes_char)
lmeResMerge = merge(lmeBestResBind, lmeSimpleResBind, by = c('species', 'var'), all.x = T)
lmeResMerge = lmeResMerge[!duplicated(lmeResMerge), ]
x = merge(lmeResMerge, aovRes, by = c('species', 'var'), all.x = T)
x_mulisig = x[which(x$lme.multi.P_value < 0.05 & x$lme.uni.P_value < 0.05),]
x_mulisig = arrange(x_mulisig, species, var, lme.multi.P_value)
write.xlsx(x_mulisig, './Routput/nucl_div_lme_regression.xlsx')

##### plot boxplot for candidate species and variable #####
cand_sp_var = x_mulisig[!duplicated(x_mulisig[, 1:2]), c(1, 2)]
if(plotFlag){
  pdf('./plot/ref_genome_wise_nucl_div_metadata_super_sig.pdf', width = 10, height = 7)
}
for (pr in 1:nrow(cand_sp_var)) {
  sp = cand_sp_var[pr, ][1] %>% unlist
  tmpDF = nucMetData[nucMetData$genome == sp, ]
  va = cand_sp_var[pr, ][2] %>% unlist
  tmpDF2 = data.frame(
    'fq_ID' = tmpDF$sample,
    'pt_ID' = tmpDF$pt_ID,
    'pt_ID.u' = tmpDF$pt_ID.u,
    'nucl_diversity' = tmpDF$nucl_diversity,
    'var' = tmpDF[, va]
  )
  
  #pairs I want to compare in list format for stat_compare_means
  if(plotFlag){
    if(class(tmpDF2$var) == 'factor'){
      #all pairs I want to compare
      if(is.factor(tmpDF2$var)){
        CN <- combn(na.omit(droplevels(unique(tmpDF2$var))) %>% levels, 2, simplify = FALSE)
      }else{
        CN <- combn(na.omit(unique(tmpDF2$var)), 2, simplify = FALSE)
      }
      p = ggboxplot(na.omit(tmpDF2), x = 'var', y = 'nucl_diversity', color = 'var', add = 'jitter',
                    add.paramsshape = 'pt_ID') +
        mytheme3 + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        stat_compare_means(comparisons=CN, method = "wilcox.test", p.adjust.method = "BH", aes(label=..p.adj..)) +
        labs(title = paste0(str_replace(sp, '.fna', ''), ' | ', va), 
             color = va,
             y = paste0('Nucleotide diversity'), x = va) + 
        scale_color_brewer(palette = 'Dark2')
      print(p)
    }else{
      p2 = ggplot(na.omit(tmpDF2), aes(x = var, y = nucl_diversity)) +
        geom_smooth(method = 'lm', color = '#969696', size = 0.5, fill = "#D5D5D5") +
        geom_point(aes(color = 'Blue', fill = 'Blue'), shape = 21, size = 2, alpha = 0.5) +
        labs(x = va, y = 'Nucleotide diversity', title = paste0(str_replace(sp, '_', ' '), ' | ', va)) +
        stat_poly_eq(formula = y ~ x,
                     aes(group=1, label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), 
                     parse = TRUE) +
        mytheme + theme(strip.text = element_text(size = 12), legend.position = 'none')
      print(p2)
      
    }
  }
}
dev.off()

#####  #####











