source('./script/library_function.R')
load('./data/RData/abuList.RData')
load('./data/RData/swabData.RData')
load('./data/RData/metaData.RData')
load('./data/RData/sampleVagitypes.RData')
load('./data/RData/colorList.RData')

if(T){ # prepare abundance matrix
  abuS = abuList[['s']]
  # discard NC samples
  removeSamples = swabData$seq_ID[grepl('^KIT', swabData$pt_ID)]
  length(removeSamples)
  abuS1 = abuS[, !colnames(abuS) %in% removeSamples]
  range(rowSums(abuS1[, 4:ncol(abuS1)])) 
  # discard unclassified
  abuS2 = abuS1[abuS1$taxa != 'unclassified', ]
  # keep samples with swab info
  abuS_mat = abuS2[, colnames(abuS2) %in% swabData$seq_ID] 
  row.names(abuS_mat) = abuS2$taxa
  range(abuS_mat)
  # abuS_mat = apply(abuS_mat, 2, function(x){x/sum(x)})
  abuS_mat = abuS_mat/100
  range(abuS_mat)
  colSums(abuS_mat)
  abuS_mat_T = t(abuS_mat)
  dim(abuS_mat)
  abuMat = abuS_mat
  abuMat_T = t(abuMat)
}

if(T){ # calculate bray curtis distance and prepare the metadata
  bray_curtis_dist <- vegan::vegdist(abuMat_T, method = "bray")
  pcoa <- cmdscale(bray_curtis_dist, k = (nrow(abuMat) - 1), eig = TRUE)
  efit = envfit(pcoa, abuMat_T)
  efitSig = data.frame(as.data.frame(scores(efit, "vectors")) * ordiArrowMul(efit),
                       r2 = efit$vectors$r,
                       Pval = efit$vectors$pvals)
  range(efitSig$Dim1)
  range(efitSig$Dim2)
  range(efitSig$r2)
  efitSig = efitSig[efitSig$Pval < 0.05, ]
  
  points = as.data.frame(pcoa$points)
  eig = pcoa$eig
  bray_curtis_pcoa_df = data.frame(seqID = row.names(abuMat_T),
                                    pcoa1 = points[,1],
                                    pcoa2 = points[,2])
  bray_curtis_pcoa_df$pt_ID.u = swabData$pt_ID.u[match(bray_curtis_pcoa_df$seqID, swabData$seq_ID)]
  bray_curtis_pcoa_df$pt_ID = swabData$pt_ID[match(bray_curtis_pcoa_df$seqID, swabData$seq_ID)]
  bray_curtis_pcoa_df$Age = metaData$Age_char[match(bray_curtis_pcoa_df$pt_ID, metaData$pt_ID)]
  bray_curtis_pcoa_df$BMI = metaData$BMI_char[match(bray_curtis_pcoa_df$pt_ID, metaData$pt_ID)]
  bray_curtis_pcoa_df$Age_num = metaData$Age[match(bray_curtis_pcoa_df$pt_ID, metaData$pt_ID)]
  bray_curtis_pcoa_df$BMI_num = metaData$BMI[match(bray_curtis_pcoa_df$pt_ID, metaData$pt_ID)]
  bray_curtis_pcoa_df$Employment = metaData$Employed[match(bray_curtis_pcoa_df$pt_ID, metaData$pt_ID)]
  bray_curtis_pcoa_df$Housing = metaData$Housing[match(bray_curtis_pcoa_df$pt_ID, metaData$pt_ID)]
  bray_curtis_pcoa_df$Term = metaData$Term_char[match(bray_curtis_pcoa_df$pt_ID, metaData$pt_ID)]
  bray_curtis_pcoa_df$Term2 = metaData$Term_char2[match(bray_curtis_pcoa_df$pt_ID, metaData$pt_ID)]
  bray_curtis_pcoa_df$Ethnicity = metaData$Ethnicity[match(bray_curtis_pcoa_df$pt_ID, metaData$pt_ID)]
  bray_curtis_pcoa_df$Ethnicity2 = metaData$Ethnicity2[match(bray_curtis_pcoa_df$pt_ID, metaData$pt_ID)]
  bray_curtis_pcoa_df$Trimester = swabData$trimester[match(bray_curtis_pcoa_df$pt_ID.u, swabData$pt_ID.u)]
  bray_curtis_pcoa_df$Gestational_Week = swabData$Sample_GA[match(bray_curtis_pcoa_df$pt_ID.u, swabData$pt_ID.u)]
  bray_curtis_pcoa_df$Vagitype = sampleVagitypes$vagitypes_Fettweis2[match(bray_curtis_pcoa_df$pt_ID.u, sampleVagitypes$pt_ID.u)]
  bray_curtis_pcoa_df$CST = sampleVagitypes$vagitypes_Susan[match(bray_curtis_pcoa_df$pt_ID.u, sampleVagitypes$pt_ID.u)]
  
}

PCoA_plot = function(matT, dataDF, Var, colorVar, alpha = 0.7, pointSize = 1.8, plotMarginal = T, if_addP = T,
                     if_coTest = F,R2val = R2val, Pval = Pval,pPos.x = -0.5, efitDF = NULL){
  if(Var == 'Term2'){
    dataDF$Term = dataDF$Term2
    Var = 'Term'
  }else if(Var == 'Ethnicity2'){
    dataDF$Ethnicity = dataDF$Ethnicity2
    Var = 'Ethnicity'
  }else{
    Var = Var 
  }
  dataDF$Var = dataDF[, Var]
  if(if_coTest){
    adoL_r = paste0("italic(R)^2 == ", round(R2val,2))
    adoL_p = paste0( "italic(P) == ", Pval)
  }else{
    ado1 = adonis2(matT~Var, data = dataDF, method="bray", by = 'margin', na.action = 'na.omit')
    adoL_r = paste0("italic(R)^2 == ", round(ado1$R2[1],2))
    adoL_p = paste0( "italic(P) == ", ado1$`Pr(>F)`[1])
  }
  if(plotMarginal){
    p2 = ggpubr::ggscatterhist(
            dataDF, x = "pcoa1", y = "pcoa2", shape = 21, 
            color = 'black', stroke = 0.1, fill = Var, palette = colorVar,
            margin.params = list(color = Var, fill = NA, size = 0.5),
            size = pointSize, alpha = alpha, print = FALSE, 
            main.plot.size = 2, margin.plot.size = 1,
            xlab = paste0("PCo1 (", round(eig[1], 1), '%)'), 
            ylab =  paste0("PCo2 (", round(eig[2], 1), '%)'),
            margin.plot = "boxplot", legend = 'right',
            ggtheme = mytheme + theme(aspect.ratio=1)
         )
    if(if_addP){
      p2$sp = p2$sp + 
        annotate(x =pPos.x, y = -0.425, label = 'Adonis test', geom = "text", size = 4.5, hjust = 0) +
        annotate(x =pPos.x, y = -0.5, label = adoL_r, geom = "text", size = 4.5, parse = TRUE, hjust = 0) +
        annotate(x =pPos.x, y = -0.575, label = adoL_p, geom = "text", size = 4.5, parse = TRUE, hjust = 0)
    }
    print(p2)
  }else{
    bray_curtis_plot1 = ggplot(data = dataDF, aes(x=pcoa1, y=pcoa2, fill = Var, color = Var)) +
      geom_point(shape = 21, alpha = alpha, stroke = 0.1, size = pointSize, color = 'black') +
      labs(x = paste0("PCo1 (", round(eig[1], 1), '%)'), 
           y =  paste0("PCo2 (", round(eig[2], 1), '%)'),
           title = "Bray-Curtis PCoA") + 
      scale_fill_manual(name = Var, values = colorVar) +
      mytheme + theme(aspect.ratio=1)
    if(!is.null(efitDF)){
      efitDF$var = row.names(efitDF)
      efitDF = efitDF[str_replace_all(row.names(efitDF), '_', ' ') %in% names(colorSpecies), ]
      efitDF$var = paste0(str_sub(efitDF$var, 1,1), '. ', str_split_fixed(efitDF$var, '_', 2)[, 2] %>% str_replace_all(., '_', ' '))
      bray_curtis_plot1 =  bray_curtis_plot1 +
        new_scale_fill() +
        geom_point(data = efitDF, aes(x = Dim1, y = Dim2), shape = 18, size = 3, color = '#FDAE61') +
        geom_text_repel(efitDF, size = 4,seed = 12,color = '#FDAE61',
                        mapping = aes(x = Dim1, y = Dim2, label = var), max.overlaps = 25)
        # gradient_color(c("#00AFBB", "#E7B800", "#FC4E07")) +
        # guides(color = guide_colourbar(title = 'Contribution')) +
        # theme(axis.text = element_text(size = 18), axis.title= element_text(size = 18),
        #       legend.text= element_text(size = 18), legend.title= element_text(size = 18),
        #       legend.spacing.x = unit(0.2,'cm'), title = element_text(size = 18), aspect.ratio = 1)
        
    }
    if(if_addP){
      bray_curtis_plot1 = bray_curtis_plot1 + annotate(x =pPos.x, y = -0.425, label = 'Adonis test', geom = "text", size = 4.5, hjust = 0) +
        annotate(x =pPos.x, y = -0.5, label = adoL_r, geom = "text", size = 4.5, parse = TRUE, hjust = 0) +
        annotate(x =pPos.x, y = -0.575, label = adoL_p, geom = "text", size = 4.5, parse = TRUE, hjust = 0)
    }
    print(bray_curtis_plot1)
  }
}


##### multiple covariate premonova#####
if(T){
  bray_curtis_pcoa_df$FOB = metaData$FOB[match(bray_curtis_pcoa_df$pt_ID, metaData$pt_ID)]
  bray_curtis_pcoa_df$Marriage = metaData$Marriage[match(bray_curtis_pcoa_df$pt_ID, metaData$pt_ID)]
  bray_curtis_pcoa_df$Marriage %<>% as.character()
  bray_curtis_pcoa_df$Marriage[is.na(bray_curtis_pcoa_df$Marriage)] = 'Unknown'
  bray_curtis_pcoa_df$Marriage %<>% factor(., levels = c('Married', 'Unmarried', 'Unknown'))
  bray_curtis_pcoa_df$Depression = metaData$Depression[match(bray_curtis_pcoa_df$pt_ID, metaData$pt_ID)]
  set.seed(233)
  ado1 = adonis2(abuMat_T~Ethnicity2+Employment+Housing+Term+BMI_num+Marriage+Age_num+Trimester+FOB+Depression, 
                 data = bray_curtis_pcoa_df, method="bray", by = 'margin', na.action = 'na.omit')
  fit1DF = as.data.frame(ado1)
  fit1DF$Factor = row.names(fit1DF) %>% str_remove('2$|_num$')
  fit1DF$Pvalue = fit1DF$`Pr(>F)`
  fit1DF = fit1DF[1:(nrow(fit1DF) -2), ]
  fit1DF$Factor %<>% factor(., levels = fit1DF$Factor[order(fit1DF$R2)])
  # openxlsx::write.xlsx(fit1DF, file = './data/adonis_cofactor_v2.xlsx',rowNames = T)
  fit1DF = openxlsx::read.xlsx('./data/adonis_cofactor_v2.xlsx', rowNames = TRUE)
  fit1DF$Factor %<>% factor(., levels = fit1DF$Factor[order(fit1DF$R2)])
  bbp = ggplot(fit1DF, aes(x=R2, y=Factor)) +
    geom_col(fill="gray", width = 0.25, alpha = 0.8) +
    geom_point(aes(fill = Pvalue), stroke = 0.2, shape = 21, size = 5) +
    scale_fill_distiller(palette = 4, name = expression(~ italic(P) ~ 'value')) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.2)), breaks = seq(0,0.1,0.02)) +
    mytheme3 + theme(strip.text = element_text(size = 13), strip.background = element_rect(fill = '#E2E5E9', color = NA),
                     # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                     legend.position = 'right') +
    xlab(expression(italic(R)^2)) + ylab("")
  cairo_pdf('./plot/adonis_cofactor_v2.pdf', height = 3.92, width = 4.5)
  print(bbp)
  dev.off()
  # add vagitype
  set.seed(233)
  ado2 = adonis2(abuMat_T~Vagitype+Ethnicity2+Employment+Housing+Term+BMI_num+Marriage+Age_num+Trimester+FOB+Depression, 
                 data = bray_curtis_pcoa_df, method="bray", by = 'margin', na.action = 'na.omit')
  fit2DF = as.data.frame(ado2)
  fit2DF$Factor = row.names(fit2DF) %>% str_remove('2$|_num$')
  fit2DF$Pvalue = fit2DF$`Pr(>F)`
  fit2DF = fit2DF[1:(nrow(fit2DF) -2), ]
  fit2DF$Factor %<>% factor(., levels = fit2DF$Factor[order(fit2DF$R2)])
  openxlsx::write.xlsx(fit2DF, file = './data/adonis_cofactor_add_Vagitype_v2.xlsx',rowNames = T)
  
  # add CST
  set.seed(233)
  ado3 = adonis2(abuMat_T~CST+Ethnicity2+Employment+Housing+Term+BMI_num+Marriage+Age_num+Trimester+FOB+Depression, 
                 data = bray_curtis_pcoa_df, method="bray", by = 'margin', na.action = 'na.omit')
  fit3DF = as.data.frame(ado3)
  fit3DF$Factor = row.names(fit3DF) %>% str_remove('2$|_num$')
  fit3DF$Pvalue = fit3DF$`Pr(>F)`
  fit3DF = fit3DF[1:(nrow(fit3DF) -2), ]
  fit3DF$Factor %<>% factor(., levels = fit3DF$Factor[order(fit3DF$R2)])
  openxlsx::write.xlsx(fit3DF, file = './data/adonis_cofactor_add_CST_v2.xlsx',rowNames = T)
}

# colored by other variables add pvalue adjusted by cofactors
if(T){
  cairo_pdf('./plot/bray_curtis_pcoa_by_others_co_v4.pdf', width = 6, height = 3.55, onefile = T)
  fit1DF = openxlsx::read.xlsx('./data/adonis_cofactor_v2.xlsx', rowNames = T)
  for(v in c('Ethnicity2', 'Employment', 'Housing', 'Term',  'BMI','Marriage', 'Age','Trimester', 'FOB', 'Depression')){
    print(v)
    R2val = fit1DF[row.names(fit1DF) == v, 'R2']
    Pval = fit1DF[row.names(fit1DF) == v, 'Pr(>F)']
    # PCoA_plot(matT = abuMat_T, dataDF = bray_curtis_pcoa_df, R2val = R2val[1],Pval = Pval[1],if_coTest = T,
    #           Var = v, colorVar = get(paste0("color", v)))
    PCoA_plot(matT = abuMat_T, dataDF = bray_curtis_pcoa_df, R2val = R2val[1],Pval = Pval[1],if_coTest = T,
              alpha = 0.9, pointSize = 2,if_addP = F,
              Var = v, colorVar = get(paste0("color", v)))
    CN = combn(na.omit(droplevels(unique(bray_curtis_pcoa_df[, v]))) %>% levels, 2, simplify = FALSE)
    pc1 = ggboxplot(bray_curtis_pcoa_df, x = v, y = 'pcoa1', color = v) +
      stat_compare_means(comparisons = CN, method = 'wilcox.test', label = '..p..') +
      scale_color_manual(values =  get(paste0("color", v)), name = v) +
      labs(x = v, y = 'PCo1')
    print(pc1) 
    pc2 = ggboxplot(bray_curtis_pcoa_df, x = v, y = 'pcoa2', color = v) +
      stat_compare_means(comparisons = CN, method = 'wilcox.test', label = '..p..') +
      scale_color_manual(values =  get(paste0("color", v)), name = v) +
      labs(x = v, y = 'PCo2')
    print(pc2) 
  }
  dev.off()
  cairo_pdf('./plot/bray_curtis_pcoa_by_Ethnicity.pdf', width = 5, height = 5, onefile = T)
  if(T){
    v = 'Ethnicity2'
    PCoA_plot(matT = abuMat_T, dataDF = bray_curtis_pcoa_df, R2val = R2val[1],Pval = Pval[1],if_coTest = T,
              alpha = 0.9,pointSize = 2,if_addP = F,
              Var = v, colorVar = get(paste0("color", v)))
    CN = combn(na.omit(droplevels(unique(bray_curtis_pcoa_df$Ethnicity2))) %>% levels, 2, simplify = FALSE)
    pc1 = ggboxplot(bray_curtis_pcoa_df, x = 'Ethnicity2', y = 'pcoa1', color = 'Ethnicity2') +
      stat_compare_means(comparisons = CN, method = 'wilcox.test', label = '..p..') +
      scale_color_manual(values = colorEthnicity2, name = 'Ethnicity') +
      labs(x = 'Ethnicity', y = 'PCo1')
    print(pc1) 
    pc2 = ggboxplot(bray_curtis_pcoa_df, x = 'Ethnicity2', y = 'pcoa2', color = 'Ethnicity2') +
      stat_compare_means(comparisons = CN, method = 'wilcox.test', label = '..p..') +
      scale_color_manual(values = colorEthnicity2, name = 'Ethnicity') +
      labs(x = 'Ethnicity', y = 'PCo2')
    print(pc2) 
  }                              
  dev.off()           
  fit2DF = openxlsx::read.xlsx('./data/adonis_cofactor_add_Vagitype_v2.xlsx',rowNames = T)
  for(v in c('Vagitype')){
    print(v)
    R2val = fit2DF[row.names(fit2DF) == v, 'R2']
    Pval = fit2DF[row.names(fit2DF) == v, 'Pr(>F)']
    PCoA_plot(matT = abuMat_T, dataDF = bray_curtis_pcoa_df,R2val = R2val[1],Pval = Pval[1],if_coTest = T,
              Var = v, colorVar = get(paste0("color", v))
    )
  }
  fit3DF = openxlsx::read.xlsx('./data/adonis_cofactor_add_CST_v2.xlsx',rowNames = T)
  for(v in c('CST')){
    print(v)
    R2val = fit3DF[row.names(fit3DF) == v, 'R2']
    Pval = fit3DF[row.names(fit3DF) == v, 'Pr(>F)']
    PCoA_plot(matT = abuMat_T, dataDF = bray_curtis_pcoa_df,R2val = R2val[1],Pval = Pval[1],if_coTest = T,
              Var = v, colorVar = get(paste0("color", v))
    )
  }
  dev.off()
}
if(T){ # colored by Vagitype
  cairo_pdf('./plot/bray_curtis_pcoa_by_vagitypes_2.pdf', width = 8.15, height = 3.78, onefile = T)
  PCoA_plot(matT = abuMat_T, dataDF = bray_curtis_pcoa_df, pointSize = 2, alpha = 0.9,
            Var = 'Vagitype', colorVar = colorVagitype, plotMarginal = F)
  PCoA_plot(matT = abuMat_T, dataDF = bray_curtis_pcoa_df, efit = efitSig,pointSize = 2, alpha = 0.9,
            Var = 'Vagitype', colorVar = colorVagitype, plotMarginal = F)

  dev.off()
}
if(T){ # colored by CST
  cairo_pdf('./plot/bray_curtis_pcoa_by_CST_2.pdf', width = 4.53, height = 3.78)
  PCoA_plot(matT = abuMat_T, dataDF = bray_curtis_pcoa_df,pointSize = 2, alpha = 0.9,
            Var = 'CST', colorVar = colorCST, plotMarginal = F)
  dev.off()
}
# colored by other variables
cairo_pdf('./plot/bray_curtis_pcoa_by_others.pdf', width = 6.57, height = 3.99, onefile = T)
for(v in c('Term2','Term', 'Ethnicity', 'Trimester', 'Age', 'BMI', 'Housing', 'Employed')){
  print(v)
  PCoA_plot(matT = abuMat_T, dataDF = bray_curtis_pcoa_df,
            Var = v, colorVar = get(paste0("color", v)))
}
dev.off()


cairo_pdf('./plot/bray_curtis_anosim_by_others.pdf', width = 3.67, height = 4.3, onefile = T)
for(v in c('Term2','Term', 'Ethnicity2','Ethnicity', 'Trimester', 'Age', 'BMI', 'Housing', 'Employed')){
  print(v)
  dataDF = bray_curtis_pcoa_df[!is.na(bray_curtis_pcoa_df[, v]), ]
  dataMat = abuMat_T[!is.na(bray_curtis_pcoa_df[, v]), ]
  anosim_result = anosim(dataMat,dataDF[, v], permutations = 999)
  par(mar = c(8, 5, 2.1, 1.1),  mgp=c(3.5,0.5,0))
  plot(anosim_result, col = c('gray', get(paste0("color", v))), 
       xlab = "", ylab = "dis.rank", xaxt="n", yaxt="n")
  axis(1, labels = c('Between', levels(bray_curtis_pcoa_df[, v])), 
       at = 1:(nlevels(bray_curtis_pcoa_df[, v])+1), las=2, tck = -0.025)
  axis(2, las=1,tck = -0.025)
}
dev.off()
######## t-sne #######
library(Rtsne)
cairo_pdf('./plot/tSNE.pdf', width = 4.79, height = 4.81)
if(T){
  set.seed(1997)
  tSNE_fit <- Rtsne(abuMat_T[, colSums(abuMat_T) != 0], pca = T, pca_scale = T)
  library(tidyverse)
  tSNE_df <- tSNE_fit$Y %>% as.data.frame()
  tSNE_df$Var = row.names(abuMat_T)
  tSNE_df$Participant = swabData$pt_ID[match(tSNE_df$Var, swabData$seq_ID)]
  tSNE_df$Term = swabData$Term[match(tSNE_df$Var, swabData$seq_ID)]
  tSNE_df$Term2 = swabData$Term2[match(tSNE_df$Var, swabData$seq_ID)]
  tSNE_df$Ethnicity = swabData$Ethnicity[match(tSNE_df$Var, swabData$seq_ID)]
  tSNE_df$Trimester = swabData$trimester[match(tSNE_df$Var, swabData$seq_ID)]
  tSNE_df$Participant %<>% as.factor()
  metaData$ptID = order(metaData$pt_ID) 
  tSNE_df$ptID = metaData$ptID[match(tSNE_df$Participant, metaData$pt_ID)]
  tSNE_df %>%
    ggplot(aes(x = V1, y = V2, fill = Participant))+
    geom_point(shape = 21, alpha = 0.2, stroke = 0.2, size = 3) +
    geom_text(aes(label = ptID, color = Participant), size = 2) +
    scale_fill_manual(values = colorParticipant, name = 'Participant') +
    scale_color_manual(values = colorParticipant, name = 'Participant') +
    labs(x = 't-SNE 1', y = 't-SNE 2', title = 't-SNE') +
    # guides(fill = guide_legend(nrow = , title.position = 'top')) +
    mytheme + theme(legend.position="none") 
}
dev.off()
######## PCA #######
if(T){
  mad.5 = abuMat[apply(abuMat,1,mad) > 0, ]
  mad.5 = abuMat[which(rowSums(abuMat) != 0), ]
  row.names(mad.5) %<>% str_replace_all(., '_', ' ')
  mad.5 = mad.5[row.names(mad.5) %in% names(colorSpecies), ]
  featureS = paste0(row.names(mad.5) %>% str_sub(., 1,1), '.',
                    str_replace(row.names(mad.5) , '^[A-Za-z]+', ''))
  row.names(mad.5) = featureS
  res.pca <- prcomp(t(mad.5))
  res.pca$center
  library(ggbiplot)
  ggb = ggbiplot(res.pca, obs.scale = 1, var.scale = 1,
                 groups = bray_curtis_pcoa_df$Vagitype,
                 varname.abbrev = F,alpha = 0.7,
                 varname.adjust = 1.1, varname.size = 3,
                 ellipse = T, var.axes = T) +
    scale_color_manual(values = colorVagitype, name = 'Vagitype') +
    mytheme
  cairo_pdf('./plot/bray_curtis_pca.pdf', width = 6.63, height = 3.05)
  print(ggb)
  dev.off()
}
if(T){
  mad.5 = abuMat[apply(abuMat,1,mad) > 0, ]
  mad.5 = abuMat[which(rowSums(abuMat) != 0), ]
  row.names(mad.5) %<>% str_replace_all(., '_', ' ')
  mad.5 = mad.5[row.names(mad.5) %in% names(colorSpecies), ]
  featureS = paste0(row.names(mad.5) %>% str_sub(., 1,1), '.',
                    str_replace(row.names(mad.5) , '^[A-Za-z]+', ''))
  row.names(mad.5) = featureS
  res.pca <- prcomp(t(mad.5))

  PCAfeature = row.names(res.pca[["rotation"]])
  cairo_pdf('./plot/pca_species.pdf', width = 9.00, height = 5.83, onefile = T)
  Ppca = fviz_pca_ind(res.pca, geom.ind = "point", pointshape = 21, pointsize = 2.5, stroke = 0.2,
                      fill.ind = bray_curtis_pcoa_df$Vagitype, alpha.ind = 0.7, 
                      addEllipses = T, ellipse.level = 0.95, ellipse.border.remove = T,
                      legend.title = "Vagitype", mean.point = F) +
    scale_fill_manual(values = colorVagitype) + 
    mytheme + theme(panel.grid = element_line(linetype = 3, size = 0.2, colour = 'gray'))
  print(Ppca)
  coordcontrib = facto_summarize(res.pca, result = c('coord', "contrib"), 
                  element = "var")
  library(colorspace)
  library(shades)
  library(ggsci)
  ppa = Ppca + 
    new_scale_fill() +
    geom_point(coordcontrib, 
               mapping = aes(x = Dim.1, y = Dim.2, color = contrib),
               shape = 18, size = 4) +
    geom_text_repel(coordcontrib, size = 4,seed = 12,
                    mapping = aes(x = Dim.1, y = Dim.2, label = name, color = contrib), max.overlaps = 25) +
    gradient_color(c("#00AFBB", "#E7B800", "#FC4E07")) +
    guides(color = guide_colourbar(title = 'Contribution')) +
    theme(axis.text = element_text(size = 18), axis.title= element_text(size = 18),
          legend.text= element_text(size = 18), legend.title= element_text(size = 18),
          legend.spacing.x = unit(0.2,'cm'), title = element_text(size = 18), aspect.ratio = 1)
  print(ppa)
  dev.off()
  
  
  
  
  
  Peig = fviz_eig(res.pca, addlabels = TRUE, barfill = '#8FB7D3',
           barcolor = 'black', linecolor = 'orange', ncp = 10) +
    scale_y_continuous(limits = c(0, 50), expand = c(0,0)) +
    theme(axis.text = element_text(size = 13, color = 'black'), axis.title = element_text(size = 13),
          title = element_text(size = 13), plot.margin = unit(c(0.5, 2.5, 0.5, 0.5),'cm'))
  print(Peig)
  dev.off()
  cairo_pdf('./plot/pca_species_contrib.pdf', width = 5.16, height = 4.20, onefile = T)
  Pcontrib = fviz_contrib(res.pca,  choice = 'var', axes = 1:3, top = 20, fill = '#8FB7D3',
                          color = 'black') +
    theme(axis.text = element_text(size = 13, color = 'black'), axis.title = element_text(size = 13),
          title = element_text(size = 13), plot.margin = unit(c(0.5, 0.5, 0.5, 2),'cm'))
  print(Pcontrib)
  Pvar = fviz_pca_var(res.pca, col.var = "contrib", 
               geom = c("point", "text"),
               max.overlaps = 200, 
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
               repel = TRUE) +
    theme(axis.text = element_text(size = 13, color = 'black'), axis.title = element_text(size = 13),
          title = element_text(size = 13))
  print(Pvar)
  dev.off()

}

