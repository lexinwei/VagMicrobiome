source('./script/library_function.R')
load('./data/RData/colorList.RData')
baseDir = './instrain_res/PRJEB34536/gene'
statList = list.files(baseDir,pattern = '^inStrain_gene_statistics')
statDF_pi_upDF = statDF_pi_downDF = statDF_dNdS_upDF = statDF_dNdS_downDF = data.frame()
statDF_dNdSDF = statDF_pNpS_upDF = statDF3_all = data.frame()
statDF_pNpS_downDF = statDF_pNpSDF = statDF_dNdSc_upDF = statDF_dNdSc_downDF = statDF_dNdScDF = data.frame()
statDF_dNdSc_big = statDF_dNdS_big = statDF_pNpS_big = data.frame()
for(spP in statList){
  statDF = read.csv(paste0(baseDir, '/', spP))
  sp = str_extract(spP, 's__.*') %>% str_remove('\\.fasta.csv') %>% str_remove('s__') %>%
    str_replace_all('_', ' ')
  statDF$pvalCov.adj = p.adjust(statDF$pvalCov, 'BH')
  statDF$pvalPi.adj = p.adjust(statDF$pvalPi, 'BH')
  statDF$pval_dNdS.adj = p.adjust(statDF$pval_dNdS, 'BH')
  statDF$pval_pNpS.adj = p.adjust(statDF$pval_pNpS, 'BH')
  statDF$pval_dNdS_common.adj = p.adjust(statDF$pval_dNdS_common, 'BH')
  statDF$fcPi_log2 = log2(statDF$fcPi)
  statDF$fc_dNdS_log2 = log2(statDF$fc_dNdS)
  statDF$fc_pNpS_log2 = log2(statDF$fc_pNpS)
  statDF$fc_dNdS_common_log2 = log2(statDF$fc_dNdS_common)
  # *** exclude gene less than 3 samples ***
  # hist(statDF$sampleNum)
  table(statDF$sampleNum >= 3)
  statDF1 = statDF[statDF$sampleNum >= 3, ]
  # *** exclude gene without annotation or annotated as DUF ***
  statDF2 = statDF1[!grepl('^DUF|Uncharacterized', statDF1$gene_name), ]
  # *** exclude gene with significant difference of coverage ***
  table(statDF2$pvalCov > 0.05)
  # statDF3 = statDF2[statDF2$pvalCov > 0.05, ]
  statDF3 = statDF2
  statDF3$Species = sp
  statDF3_all = rbind(statDF3_all, statDF3)
  # *** cat the genes with FC > fcCutoff, p<0.05 ***
  plot(log2(statDF3$fcPi), -log10(statDF3$pvalPi))
  fcCutoff = 2
  statDF_pi_up = statDF3[statDF3$fcPi >= fcCutoff & statDF3$pvalPi.adj < 0.05, ]
  statDF_pi_up = statDF_pi_up[order(statDF_pi_up$pvalPi.adj), ]
  statDF_pi_down = statDF3[statDF3$fcPi <= 1/fcCutoff & statDF3$pvalPi.adj < 0.05, ]
  statDF_pi_down = statDF_pi_down[order(statDF_pi_down$pvalPi.adj), ]
  
  # *** cat the genes with FC > fcCutoff, p<0.05 ***
  plot(log2(statDF3$fc_dNdS), -log10(statDF3$pval_dNdS.adj))
  fcCutoff = 2
  statDF_dNdS_up = statDF3[(statDF3$fc_dNdS >= fcCutoff & !is.infinite(statDF3$fc_dNdS)& !is.na(statDF3$fc_dNdS))  & 
                             statDF3$pval_dNdS.adj < 0.05, ]
  statDF_dNdS_up = statDF_dNdS_up[order(statDF_dNdS_up$pval_dNdS.adj), ]
  statDF_dNdS_down = statDF3[(statDF3$fc_dNdS <= 1/fcCutoff & !is.infinite(statDF3$fc_dNdS)& !is.na(statDF3$fc_dNdS))  & 
                               statDF3$pval_dNdS.adj < 0.05, ]
  statDF_dNdS_down = statDF_dNdS_down[order(statDF_dNdS_down$pval_dNdS.adj), ]

  # *** cat the genes with FC > fcCutoff, p<0.05 ***
  plot(log2(statDF3$fc_pNpS), -log10(statDF3$pval_pNpS.adj))
  fcCutoff = 2
  statDF_pNpS_up = statDF3[(statDF3$fc_pNpS >= fcCutoff & !is.infinite(statDF3$fc_pNpS)& !is.na(statDF3$fc_pNpS))  & 
                             statDF3$pval_pNpS.adj < 0.05, ]
  statDF_pNpS_up = statDF_pNpS_up[order(statDF_pNpS_up$pval_pNpS.adj), ]
  statDF_pNpS_down = statDF3[(statDF3$fc_pNpS <= 1/fcCutoff & !is.infinite(statDF3$fc_pNpS)& !is.na(statDF3$fc_pNpS))  & 
                               statDF3$pval_pNpS.adj < 0.05, ]
  statDF_pNpS_down = statDF_pNpS_down[order(statDF_pNpS_down$pval_pNpS.adj), ]
  # *** cat the genes with FC > fcCutoff, p<0.05 ***
  plot(log2(statDF3$fc_dNdS_common), -log10(statDF3$pval_dNdS_common.adj))
  fcCutoff = 2
  statDF_dNdSc_up = statDF3[(statDF3$fc_dNdS_common >= fcCutoff & !is.infinite(statDF3$fc_dNdS_common)& !is.na(statDF3$fc_dNdS_common))  & 
                             statDF3$pval_dNdS_common.adj < 0.05, ]
  statDF_dNdSc_up = statDF_dNdSc_up[order(statDF_dNdSc_up$pval_dNdS_common.adj), ]
  statDF_dNdSc_down = statDF3[(statDF3$fc_dNdS_common <= 1/fcCutoff &  !is.infinite(statDF3$fc_dNdS_common)& !is.na(statDF3$fc_dNdS_common))  & 
                                statDF3$pval_dNdS_common.adj < 0.05, ]
  statDF_dNdSc_down = statDF_dNdSc_down[order(statDF_dNdSc_down$pval_dNdS_common.adj), ]
  # *** cat the genes with dN/dS>1 frequency >= 0.5
  hist(statDF3$NSd_bigger_1_freq)
  statDF_dNdS = statDF3[statDF3$NSd_bigger_1_freq >= 0.5, ]
  statDF_dNdS = statDF_dNdS[order(statDF_dNdS$NSd_bigger_1_freq), ]
  # *** cat the genes with dN/dS>1 frequency >= 0.5
  hist(statDF3$NSp_bigger_1_freq)
  statDF_pNpS = statDF3[statDF3$NSp_bigger_1_freq >= 0.5, ]
  statDF_pNpS = statDF_pNpS[order(statDF_pNpS$NSp_bigger_1_freq), ]
  # *** cat the genes with dN/dS>1 frequency >= 0.5
  hist(statDF3$NSdc_bigger_1_freq)
  statDF_dNdSc = statDF3[statDF3$NSdc_bigger_1_freq >= 0.5, ]
  statDF_dNdSc = statDF_dNdSc[order(statDF_dNdSc$NSdc_bigger_1_freq), ]
  
  statDF_pi_upDF = rbind(statDF_pi_upDF,statDF_pi_up)
  statDF_pi_downDF = rbind(statDF_pi_downDF,statDF_pi_down)
  statDF_dNdS_upDF = rbind(statDF_dNdS_upDF,statDF_dNdS_up)
  statDF_dNdS_downDF = rbind(statDF_dNdS_downDF,statDF_dNdS_down)
  statDF_dNdSDF = rbind(statDF_dNdSDF,statDF_dNdS)
  statDF_pNpS_upDF = rbind(statDF_pNpS_upDF,statDF_pNpS_up)
  statDF_pNpS_downDF = rbind(statDF_pNpS_downDF,statDF_pNpS_down)
  statDF_pNpSDF = rbind(statDF_pNpSDF,statDF_pNpS)
  statDF_dNdSc_upDF = rbind(statDF_dNdSc_upDF,statDF_dNdSc_up)
  statDF_dNdSc_downDF = rbind(statDF_dNdSc_downDF,statDF_dNdSc_down)
  statDF_dNdScDF = rbind(statDF_dNdScDF,statDF_dNdSc)
  
  statDF_dNdS_big = rbind(statDF_dNdS_big, statDF3[statDF3$meanGen_dNdS > 1 & !is.na(statDF3$meanGen_dNdS) & !is.infinite(statDF3$meanGen_dNdS), ])
  statDF_pNpS_big = rbind(statDF_pNpS_big, statDF3[statDF3$meanGen_pNpS > 1 & !is.na(statDF3$meanGen_pNpS) & !is.infinite(statDF3$meanGen_pNpS), ])
  statDF_dNdSc_big = rbind(statDF_dNdSc_big, statDF3[statDF3$meanGen_dNdS_common > 1 & !is.na(statDF3$meanGen_dNdS_common) & !is.infinite(statDF3$meanGen_dNdS_common), ])
}
l = list('pi_Up' = statDF_pi_upDF,
         'pi_Down' = statDF_pi_downDF,
         'dNdS_Up' = statDF_dNdS_upDF,
         'dNdS_Down' = statDF_dNdS_downDF,
         'dNdS>1' = statDF_dNdS_big,
         'dNdS>1 Freq >= 0.5' = statDF_dNdSDF,
         'pNpS_Up' = statDF_pNpS_upDF,
         'pNpS_Down' = statDF_pNpS_downDF,
         'pNpS>1' = statDF_pNpS_big,
         'pNpS>1 Freq >= 0.5' = statDF_pNpSDF,
         'dNdSc_Up' = statDF_dNdSc_upDF,
         'dNdSc_Down' = statDF_dNdSc_downDF,
         'common dNdS>1' = statDF_dNdSc_big,
         'common dNdS>1 Freq >= 0.5' = statDF_dNdScDF)

openxlsx::write.xlsx(l, paste0('./data/instrain_gene/PRJEB34536_gene_vs_genome.adj.xlsx'))
# volcano_plot(statDF3, x = 'pi_logFC', y = 'pi_pvalue.adj', FC_cutoff = 1.2)
save(statDF3_all, file = './data/RData/PRJEB34536_statDF3_all.RData')
save(l, file = './data/RData/PRJEB34536_ll.RData')
###### scatter plot for each species, SNVs vs SNS, dN/dS vs pNpS ######
load('./data/RData/PRJEB34536_statDF3_all.RData')
cairo_pdf(file = './plot/PRJEB34536_scatterplot_dNdS_vs_pNpS.pdf', width = 4.47, height = 4.3, onefile = T)
for(sp in unique(statDF3_all$Species)){
  tmpDF = statDF3_all[statDF3_all$Species == sp, ]
  ps = ggplot(tmpDF, aes(x =SNVperKbp, y = SNSperKbp)) +
    geom_point(shape = 21) +
    geom_abline(slope = 1, linetype = 'dashed') +
    labs(title = sp, x = 'SNV per Kbp', y = 'SNS per Kbp') + 
    scale_x_continuous(expand = expansion(mult = c(0, 0.02)), 
                       limits = c(0, max(tmpDF$SNVperKbp))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02)), 
                       limits = c(0, max(tmpDF$SNSperKbp))) +
    mytheme3 + theme(aspect.ratio = 1)
  print(ps)
  ps1 = ggplot(tmpDF, aes(x = meanGen_pNpS, y = meanGen_dNdS)) +
    geom_point(shape = 21, aes(fill = 'a'), color = 'black', stroke = 0.2, size = 2.5) +
    geom_abline(slope = 1, linetype = 'dashed') +
    scale_fill_manual(values = alpha('#4292C6', 0.3)) + 
    labs(title = sp, x = 'pNpS', y = 'dN/dS') + 
    geom_hline(yintercept = 1, color = 'gray', linetype = 'dashed') +
    geom_vline(xintercept = 1, color = 'gray', linetype = 'dashed') +
    scale_x_continuous(expand = expansion(mult = c(0, 0.02)), limits = c(0, max(tmpDF$meanGen_pNpS))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02)), limits = c(0, max(tmpDF$meanGen_dNdS))) +
    geom_text_repel(aes(label = ifelse(meanGen_pNpS > 1 | meanGen_dNdS > 1, gene_name, ""))) +
    mytheme3 + theme(aspect.ratio = 1, legend.position = 'none')
  print(ps1)
  ps2 = ggplot(tmpDF, aes(x = meanGen_pi, y = meanGen_dNdS_common)) +
    geom_point(shape = 21, aes(fill = 'a'), color = 'black', stroke = 0.2, size = 2.5) +
    scale_fill_manual(values = alpha('#4292C6', 0.3)) + 
    labs(title = sp, x = expression(pi), y = 'common dN/dS') + 
    geom_hline(yintercept = 1, color = 'gray', linetype = 'dashed') +
    scale_x_continuous(expand = expansion(mult = c(0, 0.02)), limits = c(0, max(tmpDF$meanGen_pi))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02)), limits = c(0, max(tmpDF$meanGen_dNdS_common))) +
    geom_text_repel(aes(label = ifelse(meanGen_dNdS_common > 1, gene_name, ""))) +
    mytheme3 + theme(aspect.ratio = 1, legend.position = 'none')
  print(ps2)
  ps3 = ggplot(tmpDF, aes(x = meanGeno_pi, y = meanGen_pi)) +
    geom_point(shape = 21, aes(fill = 'a'), color = 'black', stroke = 0.2, size = 2.5) +
    geom_abline(slope = 1, linetype = 'dashed') +
    scale_fill_manual(values = alpha('#4292C6', 0.3)) + 
    labs(title = sp, x = 'Genome pi', y = 'Gene pi') + 
    mytheme3 + theme(aspect.ratio = 1, legend.position = 'none')
  print(ps3)
  ps4 = ggplot(tmpDF, aes(x = meanGeno_pNpS, y = meanGen_pNpS)) +
    geom_point(shape = 21, aes(fill = 'a'), color = 'black', stroke = 0.2, size = 2.5) +
    geom_abline(slope = 1, linetype = 'dashed') +
    scale_fill_manual(values = alpha('#4292C6', 0.3)) + 
    labs(title = sp, x = 'Genome pN/pS', y = 'Gene pN/pS') + 
    geom_hline(yintercept = 1, color = 'gray', linetype = 'dashed') +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02)), limits = c(0, max(tmpDF$meanGen_pNpS, na.rm = T))) +
    geom_text_repel(aes(label = ifelse(meanGen_pNpS > 1, gene_name, ""))) +
    mytheme3 + theme(aspect.ratio = 1, legend.position = 'none')
  print(ps4)
  ps5 = ggplot(tmpDF, aes(x = meanGeno_dNdS, y = meanGen_dNdS)) +
    geom_point(shape = 21, aes(fill = 'a'), color = 'black', stroke = 0.2, size = 2.5) +
    geom_abline(slope = 1, linetype = 'dashed') +
    scale_fill_manual(values = alpha('#4292C6', 0.3)) + 
    labs(title = sp, x = 'Genome dN/dS', y = 'Gene dN/dS') + 
    geom_hline(yintercept = 1, color = 'gray', linetype = 'dashed') +
      scale_y_continuous(expand = expansion(mult = c(0, 0.02)), limits = c(0, max(tmpDF$meanGen_dNdS, na.rm = T))) +
    geom_text_repel(aes(label = ifelse(meanGen_dNdS > 1, gene_name, ""))) +
    mytheme3 + theme(aspect.ratio = 1, legend.position = 'none')
  print(ps5)
  ps6 = ggplot(tmpDF, aes(x = meanGeno_dNdS_common, y = meanGen_dNdS_common)) +
    geom_point(shape = 21, aes(fill = 'a'), color = 'black', stroke = 0.2, size = 2.5) +
    geom_abline(slope = 1, linetype = 'dashed') +
    scale_fill_manual(values = alpha('#4292C6', 0.3)) + 
    labs(title = sp, x = 'Genome dN/dS (common)', y = 'Gene dN/dS (common)') + 
    geom_hline(yintercept = 1, color = 'gray', linetype = 'dashed') +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02)), limits = c(0, max(tmpDF$meanGen_dNdS_common, na.rm = T))) +
    geom_text_repel(aes(label = ifelse(meanGen_dNdS_common > 1, gene_name, ""))) +
    mytheme3 + theme(aspect.ratio = 1, legend.position = 'none')
  print(ps6)
}
dev.off()
###### statistics for mean pnps and dN/dS distribution ####
load('./data/RData/PRJEB34536_statDF3_all.RData')
paper_var = c('Species', 'gene', 'gene_name', 'gene_description', 
              'sampleNum', 'meanCoverage', 'meanBreadth', 'meanBreadth_minCov', 
              'meanGeno_pi', 'meanGen_pi', 'fcPi', 
              'meanGen_pNpS', 'meanGen_dNdS')
paper_col = c('Species', 'Gene ID', 'Gene name', 'Gene description', 
              'Sample count', 'Mean coverage', 'Mean breadth', 'Mean Breadth_minCov',
              'Mean π of genome', 'Mean π of gene', 'Fold change of gene/genome π', 
              'Mean pN/pS of gene', 'Mean dN/dS of gene')
geneDF_to_paper = statDF3_all[, paper_var]
colnames(geneDF_to_paper) = paper_col
write.xlsx(geneDF_to_paper, file = './table/gene_evoluionary_metrices_validation.xlsx')


load('./data/RData/PRJEB34536_ll.RData')
# statDF_dNdS = l[["dNdS>1" ]]
# dNdSgene = statDF_dNdS[statDF_dNdS$sampleNum > 10 & statDF_dNdS$delFreq < 0.75 &
#                           statDF_dNdS$NSd_bigger_1_freq > 0.25, ]
# statDF_pNpS = l[["pNpS>1" ]]
# pNpSgene = statDF_pNpS[statDF_pNpS$sampleNum > 10 & statDF_pNpS$delFreq < 0.75 &
#                          statDF_pNpS$NSd_bigger_1_freq > 0.25, ]

dNdSgene = statDF3_all[statDF3_all$sampleNum > 10, ]
pNpSgene = statDF3_all[statDF3_all$sampleNum > 10, ]
wdata = rbind(
  data.frame(
    Group = 'dN/dS',
    Ratio = dNdSgene$meanGen_dNdS
  ),
  data.frame(
    Group = 'pN/pS',
    Ratio = pNpSgene$meanGen_pNpS
  )
)
if(T){
  # 1. Create the histogram plot
  phist <- gghistogram(
    wdata, x = "Ratio", ylab = 'Count',
    add = "none", rug = TRUE, bins = 30,
    fill = "Group", palette = c("#00AFBB", "#E7B800")
  ) +    mytheme +theme(legend.position = 'top', legend.direction = 'vertical') 
  # 2. Create the density plot with y-axis on the right
  # Remove x axis elements
  pdensity <- ggdensity(
    wdata, x = "Ratio", ylab = 'Density',
    color= "Group", palette = c("#00AFBB", "#E7B800"),
    alpha = 0
  ) +
    scale_y_continuous(position = "right")  +
    mytheme2 +theme(legend.position = 'top', legend.direction = 'vertical') +
    cowplot::theme_half_open(13, rel_small = 1) +
    rremove("x.axis")+
    rremove("xlab") +
    rremove("x.text") +
    rremove("x.ticks") +
    rremove("legend") 
  # 3. Align the two plots and then overlay them.
  aligned_plots <- cowplot::align_plots(phist, pdensity, align="hv", axis="tblr")
  cairo_pdf('./plot/PRJEB34536_distribution_of_dnds_pnps_genes.pdf', width = 3.94, height = 4.39)
  cowplot::ggdraw(aligned_plots[[1]]) + cowplot::draw_plot(aligned_plots[[2]])
  dev.off()
}

###### statistics for pnps>1 dN/dS>1 gene numbers ####
load('./data/RData/PRJEB34536_statDF3_all.RData')
load('./data/RData/PRJEB34536_ll.RData')
# statDF_dNdS = l[["dNdS>1" ]]
# dNdSgene = statDF_dNdS[statDF_dNdS$sampleNum > 10 & statDF_dNdS$delFreq < 0.75 &
#                           statDF_dNdS$NSd_bigger_1_freq > 0.25, ]
# statDF_pNpS = l[["pNpS>1" ]]
# pNpSgene = statDF_pNpS[statDF_pNpS$sampleNum > 10 & statDF_pNpS$delFreq < 0.75 &
#                          statDF_pNpS$NSd_bigger_1_freq > 0.25, ]
dNdSgene = statDF3_all[statDF3_all$sampleNum > 10 , ]
pNpSgene = statDF3_all[statDF3_all$sampleNum > 10, ]
# dNdSgene = statDF3_all[statDF3_all$sampleNum > 10 &
#                          statDF3_all$NSd_bigger_1_freq > 0, ]
# pNpSgene = statDF3_all[statDF3_all$sampleNum > 10 &
#                          statDF3_all$NSp_bigger_1_freq > 0, ]
wdata = rbind(
  data.frame(
    Group = 'dN/dS>1',
    Frequency = dNdSgene$NSd_bigger_1_freq
  ),
  data.frame(
    Group = 'pN/pS>1',
    Frequency = pNpSgene$NSp_bigger_1_freq
  )
)
if(T){
# 1. Create the histogram plot
  phist <- gghistogram(
    wdata, x = "Frequency", ylab = 'Count',
    add = "mean", rug = TRUE, bins = 30,
    fill = "Group", palette = c("#00AFBB", "#E7B800")
  ) +    mytheme +theme(legend.position = 'top', legend.direction = 'vertical') 
  # 2. Create the density plot with y-axis on the right
  # Remove x axis elements
  pdensity <- ggdensity(
    wdata, x = "Frequency", ylab = 'Density',
    color= "Group", palette = c("#00AFBB", "#E7B800"),
    alpha = 0
  ) +
    scale_y_continuous(position = "right")  +
    mytheme2 +theme(legend.position = 'top', legend.direction = 'vertical') +
    cowplot::theme_half_open(13, rel_small = 1) +
    rremove("x.axis")+
    rremove("xlab") +
    rremove("x.text") +
    rremove("x.ticks") +
    rremove("legend") 
  # 3. Align the two plots and then overlay them.
  aligned_plots <- cowplot::align_plots(phist, pdensity, align="hv", axis="tblr")
  cairo_pdf('./plot/PRJEB34536_stat_of_dnds_pnps_genes_include_0.pdf', width = 3.94, height = 4.39)
  cowplot::ggdraw(aligned_plots[[1]]) + cowplot::draw_plot(aligned_plots[[2]])
  dev.off()
}

###### statistics for pnps>1 dN/dS>1 fraction for each species ####
load('./data/RData/PRJEB34536_statDF3_all.RData')
load('./data/RData/PRJEB34536_ll.RData')
# statDF_dNdS = l[["dNdS>1" ]]
# dNdSgene = statDF_dNdS[statDF_dNdS$sampleNum > 10 & statDF_dNdS$delFreq < 0.75 &
#                           statDF_dNdS$NSd_bigger_1_freq > 0.25, ]
# statDF_pNpS = l[["pNpS>1" ]]
# pNpSgene = statDF_pNpS[statDF_pNpS$sampleNum > 10 & statDF_pNpS$delFreq < 0.75 &
#                          statDF_pNpS$NSd_bigger_1_freq > 0.25, ]
dNdSgene = statDF3_all[statDF3_all$sampleNum > 10, ]
pNpSgene = statDF3_all[statDF3_all$sampleNum > 10, ]

Species = dNdSfreq = pNpSfreq = c()
for(sp in unique(statDF3_all$Species)){
  statDF_sub = statDF3_all[statDF3_all$Species == sp, ]
  dNdSfreq = c(dNdSfreq, sum((!is.na(statDF_sub$meanGen_dNdS)) & statDF_sub$meanGen_dNdS > 1 & (!is.infinite(statDF_sub$meanGen_dNdS)))/nrow(statDF_sub))
  pNpSfreq = c(pNpSfreq, sum((!is.na(statDF_sub$meanGen_pNpS)) & statDF_sub$meanGen_pNpS > 1 & (!is.infinite(statDF_sub$meanGen_pNpS)))/nrow(statDF_sub))
  Species = c(Species, sp)
}
freqDF = rbind(data.frame(
    Species = Species,
    freq = c(dNdSfreq, pNpSfreq),
    Metrics = c(rep('dN/dS', length(dNdSfreq)), rep('pN/pS', length(dNdSfreq))),
    Group = 'Avg. Ratio > 1'
  ),
  data.frame(
    Species = Species,
    freq = 1 - c(dNdSfreq, pNpSfreq),
    Metrics = c(rep('dN/dS', length(dNdSfreq)), rep('pN/pS', length(dNdSfreq))),
    Group = 'Avg. Ratio < 1'
  )
)
ind = which(freqDF$Group == 'Avg. Ratio > 1' & freqDF$Metrics == 'pN/pS')
freqDF$Species %<>% factor(., levels = freqDF$Species[ind][order(freqDF$freq[ind])])
gs = ggplot(freqDF, aes(fill=Group, y=freq, x=Species)) + 
  geom_bar(position="fill", stat="identity", width = 0.7, color = 'black', size = 0.3) +
  labs(y = 'Percentage of gene (%)', x = 'Species') +
  scale_fill_manual(values = c('#A6CEE3', '#1F78B4')) +
  scale_y_continuous(labels = scales::percent_format(suffix = ''), expand = c(0,0)) +
  facet_grid(.~Metrics) +
  mytheme2 + theme(strip.text = element_text(size = 13), axis.text.x = element_text(size = 8),
                   legend.position = 'top',aspect.ratio = 1.5/1,
                   axis.ticks.length.x = unit(0.05 , 'cm'),
                   legend.margin = margin(5,5,-5,10), panel.spacing = unit(0.7, "lines")) + coord_flip()
cairo_pdf('./plot/percentage_of_pnps_bigger1.pdf', width = 10, height = 4)
print(gs)
dev.off()
###### circular packing to show up and down pi genes ######
load('./data/RData/PRJEB34536_statDF3_all.RData')
load('./data/RData/PRJEB34536_ll.RData')
statDF_pi_Up = l[["pi_Up" ]]
statDF_pi_Up$Label = 'pi Up'
statDF_pi_Down = l[["pi_Down" ]]
statDF_pi_Down$Label = 'pi Down'
statDF_pi_sig = rbind(statDF_pi_Up, statDF_pi_Down)
keykeyGene = c('Gram_pos_anchor', 'Muc_B2', 'MucBP_2',  'MucBP')

packingDF = statDF_pi_sig[statDF_pi_sig$Species == "Gardnerella vaginalis", ]
packingDF = statDF_pi_Up
packingDF = statDF_pi_sig
packingDF = statDF_pi_sig[statDF_pi_sig$sampleNum > 10, ]
packingDF = packingDF[order(-abs(packingDF$fcPi_log2)), ]
# packingDF = packingDF[packingDF$gene_name %in%keykeyGene, ]
packingDF$gene %<>% str_replace_all(., '\\.', '-')
if(T){
  edgesDF = data.frame(
    from = paste0(packingDF$Species), 
    to = paste0(packingDF$Species, '.', packingDF$gene)
  )
  verticesDF = data.frame(
    name = c(unique(edgesDF$from), unique(edgesDF$to)),
    # size = c(rep(0.5, length(unique(edgesDF$from))), rep(0.1, length(unique(edgesDF$to)))),
    Gene = c(unique(edgesDF$from), str_split_i(unique(edgesDF$to), '\\.', -1))
  )
  verticesDF$shortName = packingDF$gene_name[match(verticesDF$Gene, packingDF$gene)]
  verticesDF$shortName[is.na(verticesDF$shortName)] = verticesDF$name[is.na(verticesDF$shortName)]
  verticesDF$FoldChange_gene = packingDF$fcPi[match(verticesDF$name, paste0(packingDF$Species, '.',  packingDF$gene))]
  verticesDF$FoldChange_gene_log2 = log2(verticesDF$FoldChange_gene)
  verticesDF$size = packingDF$sampleNum[match(verticesDF$name, paste0(packingDF$Species, '.',  packingDF$gene))]
  verticesDF$size = verticesDF$size/100
  verticesDF$size[is.na(verticesDF$size)] = 1
  # verticesDF$FoldChange_gene_log2[is.na(  verticesDF$FoldChange_gene_log2)] = 0
  verticesDF$Change[verticesDF$FoldChange_gene_log2 > 0] = 'Up'
  verticesDF$Change[verticesDF$FoldChange_gene_log2 < 0] = 'Down'
  verticesDF$Change %<>% factor(., levels = c('Up', 'Down'))
  table(verticesDF$Change)
  verticesDF$Description = packingDF$gene_description[match(verticesDF$Gene, packingDF$gene)]
  verticesDF$Key = 'A:Other'
  verticesDF$Key[grepl('Ribosomal', verticesDF$Description, ignore.case = T)] = 'B:Ribosome'
  verticesDF$Key[grepl('toxin', verticesDF$Description, ignore.case = T)] = 'Toxin-antitoxin system'
  verticesDF$Key[grepl('Integrase', verticesDF$Description, ignore.case = T)] = 'Integrase'
  verticesDF$Key[grepl('Transposase', verticesDF$Description, ignore.case = T)] = 'Transposase'
  verticesDF$Key[(verticesDF$shortName %in% keykeyGene)] = 'Z:Bacterial adhesion component'
  table(verticesDF$Key)
  keyCategory= c('Z:Bacterial adhesion component', 'Integrase', 'Transposase', 'Toxin-antitoxin system','B:Ribosome', 'A:Other')
  verticesDF$Key %<>% factor(.,levels = keyCategory)
  keyColor = c('#E4211C', '#377EB8', '#33A02B', '#984EA3', '#A6761D', 'gray')
  names(keyColor) = keyCategory
  # verticesDF$shortName[duplicated(verticesDF$shortName)]
  mygraph = graph_from_data_frame(edgesDF, vertices = verticesDF )
  # set.seed(86786)
  # set.seed(6567)
  set.seed(17358)
  cP = ggraph(mygraph, layout = 'circlepack') +  # weight = size
    geom_node_circle(aes(color = Key, fill = FoldChange_gene_log2),  alpha = 0.7, linetype = 1) +
    geom_node_text(aes(label = shortName), size = 1.5) +
    # scale_color_manual(values = c('Up' = '#E4211C', 'Down' = '#377EB8'), na.value = 'gray') +
    scale_color_manual(values = keyColor) +
    colorspace::scale_fill_continuous_divergingx(name = 'log2(Fold Change)', rev = F, palette = 'BrBG', mid = 0, na.value = 'white') + 
    theme_void() + theme(legend.text = element_text(size = 13), legend.key.size = unit(0.5, 'cm'), 
                         legend.title = element_text(size = 13), plot.margin = unit(c(0,1,0,0), 'cm'))
  cairo_pdf('./plot/PRJEB34536_pi_high_low_gene_circle_packing_2.pdf', width = 11.97 , height = 8.35)
  print(cP)
  dev.off()
}

###### circular packing to show pN/pS > 1 genes ######
load('./data/RData/PRJEB34536_statDF3_all.RData')
load('./data/RData/PRJEB34536_ll.RData')
# x
statLow = statDF3_all[statDF3_all$sampleNum > 20 & statDF3_all$delFreq < 0.5 & 
                        statDF3_all$meanGen_pNpS < 1 & (!is.na(statDF3_all$meanGen_pNpS)) & (!is.infinite(statDF3_all$meanGen_pNpS)) &
                        statDF3_all$NSp_bigger_1_freq == 0 &
                        statDF3_all$pval_pNpS.adj < 0.05 & statDF3_all$fc_pNpS_log2 < -1, ]
packingDF = statLow
packingDF$freq = packingDF$NSp_bigger_1_freq

statDF_pNpS = l[["pNpS>1" ]]  # meanGen pN/dS > 1
# 0,2
packingDF = statDF_pNpS[statDF_pNpS$sampleNum > 20 & statDF_pNpS$delFreq < 0.5 & statDF_pNpS$NSp_bigger_1_freq > 0.35 &
                          statDF_pNpS$pval_pNpS.adj < 0.05 & statDF_pNpS$fc_pNpS_log2 > 1, ]
# 3
packingDF = statDF_pNpS[statDF_pNpS$sampleNum > 10 & 
                          statDF_pNpS$meanGen_pNpS > 2 &
                          statDF_pNpS$NSp_bigger_1_freq > 0.25, ]
packingDF$freq = packingDF$NSp_bigger_1_freq

table(packingDF$Species)
packingDF = packingDF[order(-packingDF$freq), ]
keykeyGene = c('Gram_pos_anchor', 'Muc_B2', 'MucBP_2',  'MucBP')
# packingDF = packingDF[packingDF$gene_name %in%keykeyGene, ]
packingDF$gene %<>% str_replace_all(., '\\.', '-')
if(T){
  edgesDF = data.frame(
    from = paste0(packingDF$Species), 
    to = paste0(packingDF$Species, '.', packingDF$gene)
  )
  verticesDF = data.frame(
    name = c(unique(edgesDF$from), unique(edgesDF$to)),
    # size = c(rep(0.5, length(unique(edgesDF$from))), rep(0.1, length(unique(edgesDF$to)))),
    Gene = c(unique(edgesDF$from), str_split_i(unique(edgesDF$to), '\\.', -1))
  )
  verticesDF$shortName = packingDF$gene_name[match(verticesDF$Gene, packingDF$gene)]
  verticesDF$shortName[is.na(verticesDF$shortName)] = verticesDF$name[is.na(verticesDF$shortName)]
  verticesDF$FoldChange_gene = packingDF$fcPi[match(verticesDF$name, paste0(packingDF$Species, '.',  packingDF$gene))]
  verticesDF$FoldChange_gene_log2 = log2(verticesDF$FoldChange_gene)
  verticesDF$size = packingDF$sampleNum[match(verticesDF$name, paste0(packingDF$Species, '.',  packingDF$gene))]
  verticesDF$size = verticesDF$size/100
  verticesDF$size[is.na(verticesDF$size)] = 1
  verticesDF$Frequency = packingDF$freq[match(verticesDF$name, paste0(packingDF$Species, '.',  packingDF$gene))]
  # verticesDF$FoldChange_gene_log2[is.na(  verticesDF$FoldChange_gene_log2)] = 0
  verticesDF$Description = packingDF$gene_description[match(verticesDF$Gene, packingDF$gene)]
  verticesDF$Key = 'A:Other'
  verticesDF$Key[grepl('Ribosomal', verticesDF$Description, ignore.case = T)] = 'B:Ribosome'
  verticesDF$Key[grepl('toxin', verticesDF$Description, ignore.case = T)] = 'Toxin-antitoxin system'
  verticesDF$Key[grepl('Integrase', verticesDF$Description, ignore.case = T)] = 'Integrase'
  verticesDF$Key[grepl('Transposase', verticesDF$Description, ignore.case = T)] = 'Transposase'
  verticesDF$Key[(verticesDF$shortName %in% keykeyGene)] = 'Z:Bacterial adhesion component'
  table(verticesDF$Key)
  keyCategory= c('Z:Bacterial adhesion component', 'Integrase', 'Transposase', 'Toxin-antitoxin system','B:Ribosome', 'A:Other')
  verticesDF$Key %<>% factor(.,levels = keyCategory)
  keyColor = c('#E4211C', '#377EB8', '#33A02B', '#984EA3', '#A6761D', 'gray')
  names(keyColor) = keyCategory
  # verticesDF$shortName[duplicated(verticesDF$shortName)]
  mygraph = graph_from_data_frame(edgesDF, vertices = verticesDF )
}
set.seed(789)
cP = ggraph(mygraph, layout = 'circlepack') +  # weight = size
  geom_node_circle(aes(color = Key, fill = Frequency),  alpha = 0.7, linetype = 1) +
  geom_node_text(aes(label = shortName), size = 1.5) +
  # scale_color_manual(values = c('Up' = '#E4211C', 'Down' = '#377EB8'), na.value = 'gray') +
  scale_color_manual(values = keyColor) +
  # colorspace::scale_fill_continuous_divergingx(name = 'pN/pS>1 frequency', palette = 'RdBu', na.value = 'white') +
  scale_fill_gradient(name = 'pN/pS>1 frequency', # Earth Spectral Roma
                      low = '#FFFEE5', high = '#42AB5D', na.value = 'white') +
  theme_void() + theme(legend.text = element_text(size = 13), legend.key.size = unit(0.5, 'cm'), 
                       legend.title = element_text(size = 13), plot.margin = unit(c(0,1,0,0), 'cm'))
cairo_pdf('./plot/PRJEB34536_pNpS_gene_circle_packing_3.pdf', width = 8.81 , height = 4.76)
print(cP)
dev.off()
###### circular packing to show dNdS > 1 genes ######
load('./data/RData/PRJEB34536_statDF3_all.RData')
load('./data/RData/PRJEB34536_ll.RData')
# 0,2(0.3),3
statDF_dNdS = l[["dNdS>1"]]  # meanGen dN/dS > 1
packingDF = statDF_dNdS[statDF_dNdS$sampleNum > 10 & 
                          statDF_dNdS$NSd_bigger_1_freq > 0.25, ]

packingDF$freq = packingDF$NSd_bigger_1_freq
keykeyGene = c('Gram_pos_anchor', 'Muc_B2', 'MucBP_2',  'MucBP')
# packingDF = packingDF[packingDF$gene_name %in%keykeyGene, ]
packingDF$gene %<>% str_replace_all(., '\\.', '-')
if(T){
  edgesDF = data.frame(
    from = paste0(packingDF$Species), 
    to = paste0(packingDF$Species, '.', packingDF$gene)
  )
  verticesDF = data.frame(
    name = c(unique(edgesDF$from), unique(edgesDF$to)),
    # size = c(rep(0.5, length(unique(edgesDF$from))), rep(0.1, length(unique(edgesDF$to)))),
    Gene = c(unique(edgesDF$from), str_split_i(unique(edgesDF$to), '\\.', -1))
  )
  verticesDF$shortName = packingDF$gene_name[match(verticesDF$Gene, packingDF$gene)]
  verticesDF$shortName[is.na(verticesDF$shortName)] = verticesDF$name[is.na(verticesDF$shortName)]
  verticesDF$FoldChange_gene = packingDF$fcPi[match(verticesDF$name, paste0(packingDF$Species, '.',  packingDF$gene))]
  verticesDF$FoldChange_gene_log2 = log2(verticesDF$FoldChange_gene)
  verticesDF$size = packingDF$sampleNum[match(verticesDF$name, paste0(packingDF$Species, '.',  packingDF$gene))]
  verticesDF$size = verticesDF$size/100
  verticesDF$size[is.na(verticesDF$size)] = 1
  verticesDF$Frequency = packingDF$freq[match(verticesDF$name, paste0(packingDF$Species, '.',  packingDF$gene))]
  # verticesDF$FoldChange_gene_log2[is.na(  verticesDF$FoldChange_gene_log2)] = 0
  verticesDF$Description = packingDF$gene_description[match(verticesDF$Gene, packingDF$gene)]
  verticesDF$Key = 'A:Other'
  verticesDF$Key[grepl('Ribosomal', verticesDF$Description, ignore.case = T)] = 'B:Ribosome'
  verticesDF$Key[grepl('toxin', verticesDF$Description, ignore.case = T)] = 'Toxin-antitoxin system'
  verticesDF$Key[grepl('Integrase', verticesDF$Description, ignore.case = T)] = 'Integrase'
  verticesDF$Key[grepl('Transposase', verticesDF$Description, ignore.case = T)] = 'Transposase'
  verticesDF$Key[(verticesDF$shortName %in% keykeyGene)] = 'Z:Bacterial adhesion component'
  table(verticesDF$Key)
  keyCategory= c('Z:Bacterial adhesion component', 'Integrase', 'Transposase', 'Toxin-antitoxin system','B:Ribosome', 'A:Other')
  verticesDF$Key %<>% factor(.,levels = keyCategory)
  keyColor = c('#E4211C', '#377EB8', '#33A02B', '#984EA3', '#A6761D', 'gray')
  names(keyColor) = keyCategory
  # verticesDF$shortName[duplicated(verticesDF$shortName)]
  mygraph = graph_from_data_frame(edgesDF, vertices = verticesDF )
}
set.seed(35353)
set.seed(345)
cP = ggraph(mygraph, layout = 'circlepack') +  # weight = size
  geom_node_circle(aes(color = Key, fill = Frequency),  alpha = 1, linetype = 1) +
  geom_node_text(aes(label = shortName), size = 2) +
  # scale_color_manual(values = c('Up' = '#E4211C', 'Down' = '#377EB8'), na.value = 'gray') +
  scale_color_manual(values = keyColor) +
  scale_fill_gradient(name = 'dN/dS>1 frequency', # Earth Spectral Roma
                       low = '#FFFDBF', high = '#66C2A5', na.value = 'white') +
  theme_void() + theme(legend.text = element_text(size = 13), legend.key.size = unit(0.5, 'cm'), 
                       legend.title = element_text(size = 13), plot.margin = unit(c(0,1,0,0), 'cm'))
cairo_pdf('./plot/PRJEB34536_dNdS_gene_circle_packing_3.pdf', width = 6.92 , height = 4.36)
print(cP)
dev.off()

######### evolutionary metrics for keyGenes, color,fill,shape,size #####
load('./data/RData/PRJEB34536_statDF3_all.RData')
load('./data/RData/PRJEB34536_ll.RData')
keyGenes = c('Muc_B2', 'MucBP', 'MucBP_2', 'Gram_pos_anchor')
myPal = c('#FB8072', '#FDB462', '#B3DE68', '#80B1D3')
names(myPal) = keyGenes
statKey = statDF3_all[statDF3_all$gene_name %in% keyGenes & statDF3_all$sampleNum > 5,]
nrow(statKey)
myShape = c(22:25)
names(myShape) = keyGenes
statKey$Species %<>% str_replace_all(., 'Lactobacillus', 'L.')
names(colorSpecies)  %<>% str_replace_all(., 'Lactobacillus', 'L.')
statKey$Species %<>% factor(., levels = intersect(names(colorSpecies), statKey$Species))
statKey$gene_name %<>% factor(., levels = keyGenes)
sP = ggplot(statKey, aes(x = meanGen_dNdS, y = meanGen_pNpS)) +
  geom_point(aes(color = Species, fill = meanGen_pi, shape = gene_name, size = NSp_bigger_1_freq), stroke = 1) +
  geom_hline(yintercept = 1, linewidth = 0.3, linetype = 'dashed', color = 'gray') +
  geom_vline(xintercept = 1, linewidth = 0.3, linetype = 'dashed', color = 'gray') +
  scale_shape_manual(values = myShape, name = 'Gene') +
  # scale_fill_manual(values = colorSpecies) +
  scale_fill_gradient(low = 'white', high = 'red', name = expression(pi)) +
  scale_color_manual(values = colorSpecies, name = 'Species') +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 3)), 
         color = guide_legend(override.aes = list(shape = 21, size = 3)),
         size = guide_legend(title = 'pN/pS>1 frequency')) +
  labs(x = 'Avg. dN/dS', y = 'Avg. pN/pS') +
  mytheme3 + theme(legend.spacing.y = unit(0.1, 'cm'), aspect.ratio = 1)
# ggplot(statKey, aes(x = meanGen_pNpS, y = meanGen_pi)) +
#   geom_point()
cairo_pdf('./plot/PRJEB34536_key_genes_evolutionary_metrics.pdf', width = 5.80, height = 6.79)
print(sP)
dev.off()


######### evolutionary metrics for each gene along genomes loci#####
load('./data/RData/PRJEB34536_statDF3_all.RData')
load('./data/RData/PRJEB34536_ll.RData')
table(statDF3_all$Species)
keyGenes = c('DDE_Tnp_1', 'DDE_Tnp_IS66', 'Gram_pos_anchor',
             'Methylase_S', 'MFS_1', 'Muc_B2', 'MucBP_2',  'MucBP','Nramp', 'PIN_3',
             'Relaxase', 'RelB', 'Rib', 'YkuD')
keyGenes = c('MucBP_2','MucBP', 'Muc_B2', 'Gram_pos_anchor')
statDF3_all$Key_gene =  statDF3_all$gene_name
statDF3_all$Key_gene[!statDF3_all$Key_gene %in% keyGenes] = 'Other'
library(paletteer)
pal = paletteer_d("ggsci::nrc_npg")[c(1,3,10,5,6,7,8,9,4)] %>% as.character()
cairo_pdf(file =  './plot/PRJEB34536_gene_for_each_species_evolutionary_metrics.pdf', width = 11.14, height = 7.14, onefile = T)
for(sp in unique(statDF3_all$Species)){
  statDF3_all_sub = statDF3_all[statDF3_all$Species == sp, ]
  length(unique(statDF3_all_sub$gene))
  range(statDF3_all_sub$breadth_minCov)
  geneOrd = data.frame(
    Gene = unique(statDF3_all_sub$gene),
    Order = naturalorder(unique(statDF3_all_sub$gene))
  )
  statDF3_all_sub$gene_order = geneOrd$Order[match(statDF3_all_sub$gene, geneOrd$Gene)]
  # ggplot(statDF3_all, aes(x = gene_order, y = pNpS)) +
  #   geom_point()
  statDF3_all_sub = statDF3_all_sub[!is.na(statDF3_all_sub$meanGen_pNpS), ]
  statDF3_all_sub$density = get_density(statDF3_all_sub$gene_order, statDF3_all_sub$meanGen_pNpS, n = 100)
  pp1 = ggplot() + 
    geom_point(statDF3_all_sub[statDF3_all_sub$Key_gene == 'Other', ], 
               mapping = aes(x = gene_order, y = meanGen_pNpS, color = density),
               stroke = 0.5, alpha = 0.5) +
    geom_point(statDF3_all_sub[statDF3_all_sub$Key_gene != 'Other', ], size = 2,
               mapping = aes(x = gene_order, y = meanGen_pNpS, shape = Key_gene,  fill = Key_gene),
               color = 'black', stroke = 0.5, alpha = 0.75) +
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'red', linewidth = 0.6) +
    labs(x = 'Gene', y = 'pN/pS', size = 'pN/pS', fill = 'Gene', color = 'Density', shape = 'Gene', title = sp) +
    scale_y_continuous(expand = expansion(mult = c(0.025, 0.1))) +
    # scale_size_continuous(range = c(0.1, 3)) +
    scale_color_viridis() + 
    scale_fill_manual(values = pal) +
    scale_shape_manual(values = c(21:25)) +
    mytheme + theme(strip.text = element_text(size = 13), aspect.ratio = 1)
  print(pp1)
}
dev.off()
###### plot dN/dS > 1 frequency > 0.5 ####
load('./data/RData/PRJEB34536_statDF3_all.RData')
load('./data/RData/PRJEB34536_ll.RData')
statDF_dNdSDF = l[["dNdS>1 Freq >= 0.5" ]]

if(T){
  keyGenes = c('DDE_Tnp_1', 'DDE_Tnp_IS66', 'Gram_pos_anchor',
               'Methylase_S', 'MFS_1', 'Muc_B2', 'MucBP_2', 'MucBP', 'Nramp', 'PIN_3',
               'Relaxase', 'RelB', 'Rib', 'YkuD')
  nrow(statDF3_all)
  statplot = statDF3_all[statDF3_all$gene_name %in% keyGenes, ]
  statplot = statDF_dNdSDF[statDF_dNdSDF$gene_name %in% keyGenes, ]
  # statplot = statDF3_all[statDF3_all$NSd_bigger_1_freq > 0.5, ]
  nrow(statplot)
  statplot$Species_short = paste0(str_sub(statplot$Species, 1, 1), '. ',
                                  str_split_fixed(statplot$Species, ' ', 2)[, 2])
  statplot$Species_gene = paste0(
    statplot$gene_name,
    ' (', statplot$Species_short, ')'
  )
  statplot = statplot[, c('NSp_bigger_1_freq','NSd_bigger_1_freq','NSdc_bigger_1_freq',
                          'Species_gene', 'gene')]
  x = 'NSd_bigger_1_freq'
  statplot$x = statplot[, x]
  statplot = statplot[order(statplot$x, decreasing = T), ]
  statplot = statplot[!duplicated(statplot$Species_gene), ]
  statplot$Species_gene = factor(statplot$Species_gene, 
                                 levels = statplot$Species_gene[order(statplot$x)])
  dP = ggplot(statplot[statplot$x > 0, ], aes(x = x*100, y = Species_gene)) +
    geom_col(width = 0.8, color = '#386CB0', fill = '#386CB0', alpha = 0.9) +
    scale_x_continuous(expand = c(0,0), limits = c(0,100)) +
    labs(x = 'Fraction of dN/dS > 1, %', y = 'Gene (species)') +
    mytheme3 + theme(plot.margin = unit(c(0.2,1,0.2,0.2), 'cm'))
  print(dP)
}
###### plot pN/pS > 1 frequency > 0.5 ####
load('./data/RData/PRJEB34536_statDF3_all.RData')
load('./data/RData/PRJEB34536_ll.RData')
statDF_pNpSDF = l[["pNpS>1" ]]
if(T){
  # keyGenes = c('DDE_Tnp_1', 'DDE_Tnp_IS66', 'Gram_pos_anchor',
  #              'Methylase_S', 'MFS_1', 'Muc_B2', 'MucBP_2',  'MucBP','Nramp', 'PIN_3',
  #              'Relaxase', 'RelB', 'Rib', 'YkuD')
  keyGenes = c('Gram_pos_anchor', 'Muc_B2', 'MucBP_2',  'MucBP')
  nrow(statDF3_all)
  statplot = statDF3_all[statDF3_all$gene_name %in% keyGenes, ]
  statplot = statDF_pNpSDF[statDF_pNpSDF$gene_name %in% keyGenes, ]
  nrow(statplot)
  statplot$Species_short = paste0(str_sub(statplot$Species, 1, 1), '. ',
                                  str_split_fixed(statplot$Species, ' ', 2)[, 2])
  statplot$Species_gene = paste0(
    statplot$gene_name,
    ' (', statplot$Species_short, ')'
  )
  statplot = statplot[, c('NSp_bigger_1_freq','NSd_bigger_1_freq','NSdc_bigger_1_freq',
                          'Species_gene', 'gene')]
  x = 'NSp_bigger_1_freq'
  statplot$x = statplot[, x]
  statplot = statplot[order(statplot$x, decreasing = T), ]
  statplot = statplot[!duplicated(statplot$Species_gene), ]
  statplot$Species_gene = factor(statplot$Species_gene, 
                                 levels = statplot$Species_gene[order(statplot$x)])
  dPs2 = ggplot(statplot[statplot$x > 0, ], aes(x = x*100, y = Species_gene)) +
    geom_col(width = 0.7, color = '#386CB0', fill = '#386CB0', alpha = 0.9) +
    scale_x_continuous(expand = c(0,0), limits = c(0,100)) +
    labs(x = 'Fraction of pN/pS > 1, %', y = 'Gene (species)') +
    mytheme3 + theme(plot.margin = unit(c(0.2,1,0.2,0.2), 'cm'))
  print(dPs2)
}
###### plot common dN/dS > 1 frequency > 0.5 ####
load('./data/RData/PRJEB34536_statDF3_all.RData')
load('./data/RData/PRJEB34536_ll.RData')
statDF_dNdScDF = l[["common dNdS>1 Freq >= 0.5" ]]
if(T){
  keyGenes = c('DDE_Tnp_1', 'DDE_Tnp_IS66', 'Gram_pos_anchor',
               'Methylase_S', 'MFS_1', 'Muc_B2', 'MucBP_2',  'MucBP','Nramp', 'PIN_3',
               'Relaxase', 'RelB', 'Rib', 'YkuD')
  statplot = statDF_dNdScDF[statDF_dNdScDF$gene_name %in% keyGenes, ]
  nrow(statplot)
  statplot$Species_short = paste0(str_sub(statplot$Species, 1, 1), '. ',
                                  str_split_fixed(statplot$Species, ' ', 2)[, 2])
  statplot$Species_gene = paste0(
    statplot$gene_name,
    ' (', statplot$Species_short, ')'
  )
  statplot = statplot[, c('NSp_bigger_1_freq','NSd_bigger_1_freq','NSdc_bigger_1_freq',
                          'Species_gene', 'gene')]
  x = 'NSdc_bigger_1_freq'
  statplot$x = statplot[, x]
  statplot = statplot[order(statplot$x, decreasing = T), ]
  statplot = statplot[!duplicated(statplot$Species_gene), ]
  statplot$Species_gene = factor(statplot$Species_gene, 
                                 levels = statplot$Species_gene[order(statplot$x)])
  dPs3 = ggplot(statplot[statplot$x > 0, ], aes(x = x*100, y = Species_gene)) +
    geom_col(width = 0.7, color = '#386CB0', fill = '#386CB0', alpha = 0.9) +
    scale_x_continuous(expand = c(0,0), limits = c(0,100)) +
    labs(x = 'Fraction of common dN/dS > 1, %', y = 'Gene (species)') +
    mytheme3 + theme(plot.margin = unit(c(0.2,1.2,0.2,0.2), 'cm'))
  print(dPs3)
}
cairo_pdf('./plot/PRJEB34536_dnds_freq_barplot.pdf', width = 5.21, height = 5.60, onefile = T)
print(dP)
print(dPs2)
print(dPs3)
dev.off()
###### plot mean dN/dS > 1  ####
load('./data/RData/PRJEB34536_statDF3_all.RData')
load('./data/RData/PRJEB34536_ll.RData')
statDF_dNdS_big = l[["dNdS>1"  ]]
statDF_pNpS_big = l[["pNpS>1"]]
statDF_dNdSc_big = l[["common dNdS>1"]]
if(T){
  statDF_big = rbind(statDF_dNdS_big, statDF_pNpS_big, statDF_dNdSc_big)
  nrow(statDF_dNdS_big)
  keyGenes = c('DDE_Tnp_1', 'DDE_Tnp_IS66', 'Gram_pos_anchor',
               'Methylase_S', 'MFS_1', 'Muc_B2', 'MucBP_2',  'MucBP','Nramp', 'PIN_3',
               'Relaxase', 'RelB', 'Rib', 'YkuD')
  statplot = statDF_dNdS_big[statDF_dNdS_big$gene_name %in% keyGenes, ]
  nrow(statplot)
  statplot$Species_short = paste0(str_sub(statplot$Species, 1, 1), '. ',
                                  str_split_fixed(statplot$Species, ' ', 2)[, 2])
  statplot$Species_gene = paste0(
    statplot$gene_name,
    ' (', statplot$Species_short, ')'
  )
  statplot = statplot[, c('NSp_bigger_1_freq', 'meanGen_pNpS',
                          'NSd_bigger_1_freq', 'meanGen_dNdS',
                          'NSdc_bigger_1_freq', 'meanGen_dNdS_common',
                          'Species_gene', 'gene')]
  x = 'meanGen_dNdS'
  statplot$x = statplot[, x]
  statplot = statplot[order(statplot$x, decreasing = T), ]
  statplot = statplot[!duplicated(statplot$Species_gene), ]
  statplot$Species_gene = factor(statplot$Species_gene, 
                                 levels = statplot$Species_gene[order(statplot$x)])
  dP = ggplot(statplot, aes(x = x, y = Species_gene)) +
    geom_col(width = 0.8, color = '#386CB0', fill = '#386CB0', alpha = 0.9) +
    labs(x = 'dN/dS', y = 'Gene (species)') +
    mytheme3
  print(dP)
}
###### plot mean pN/pS > 1  ####
load('./data/RData/PRJEB34536_statDF3_all.RData')
load('./data/RData/PRJEB34536_ll.RData')
statDF_dNdS_big = l[["dNdS>1"  ]]
statDF_pNpS_big = l[["pNpS>1"]]
statDF_dNdSc_big = l[["common dNdS>1"]]
if(T){
  statDF_big = rbind(statDF_dNdS_big, statDF_pNpS_big, statDF_dNdSc_big)
  nrow(statDF_pNpS_big)
  keyGenes = c('DDE_Tnp_1', 'DDE_Tnp_IS66', 'Gram_pos_anchor', 'MucBP',
               'Methylase_S', 'MFS_1', 'Muc_B2', 'MucBP_2', 'Nramp', 'PIN_3',
               'Relaxase', 'RelB', 'Rib', 'YkuD')
  statplot = statDF_pNpS_big[statDF_pNpS_big$gene_name %in% keyGenes, ]
  nrow(statplot)
  statplot$Species_short = paste0(str_sub(statplot$Species, 1, 1), '. ',
                                  str_split_fixed(statplot$Species, ' ', 2)[, 2])
  statplot$Species_gene = paste0(
    statplot$gene_name,
    ' (', statplot$Species_short, ')'
  )
  statplot = statplot[, c('NSp_bigger_1_freq', 'meanGen_pNpS',
                          'NSd_bigger_1_freq', 'meanGen_dNdS',
                          'NSdc_bigger_1_freq', 'meanGen_dNdS_common',
                          'Species_gene', 'gene')]
  x = 'meanGen_pNpS'
  statplot$x = statplot[, x]
  statplot = statplot[order(statplot$x, decreasing = T), ]
  statplot = statplot[!duplicated(statplot$Species_gene), ]
  statplot$Species_gene = factor(statplot$Species_gene, 
                                 levels = statplot$Species_gene[order(statplot$x)])
  statplot$Gene = 'Other'
  statplot$Gene[grep('^Muc', statplot$Species_gene)] = 'Mucin binding gene'
  statplot$Gene[grep('^MFS', statplot$Species_gene)] = 'Major Facilitator Superfamily'
  dP2 = ggplot(statplot, aes(x = x, y = Species_gene), color = 'gray', fill = 'gray') +
    geom_col(width = 0.8,alpha = 0.5) +
    labs(x = 'pN/pS', y = 'Gene (species)') +
    # scale_color_manual(values = c('Mucin binding gene' = '#F87A71',
    #                               'Major Facilitator Superfamily' = '#20C5CA',
    #                               'Other' = 'gray')) +
    # scale_fill_manual(values = c('Mucin binding gene' = '#F87A71',
    #                               'Major Facilitator Superfamily' = '#20C5CA',
    #                              'Other' = 'gray')) +
    mytheme3 + theme(legend.position = 'top', legend.margin = margin(0,0,-5,-200), legend.title = element_blank(),
                     plot.margin = unit(c(0.1,1,0.1,0.4), 'cm'))
  print(dP2)
}
###### plot mean pN/pS > 1  ####
load('./data/RData/PRJEB34536_statDF3_all.RData')
load('./data/RData/PRJEB34536_ll.RData')
statDF_dNdS_big = l[["dNdS>1"  ]]
statDF_pNpS_big = l[["pNpS>1"]]
statDF_dNdSc_big = l[["common dNdS>1"]]
if(T){
  statDF_big = rbind(statDF_dNdS_big, statDF_pNpS_big, statDF_dNdSc_big)
  nrow(statDF_dNdSc_big)
  keyGenes = c('DDE_Tnp_1', 'DDE_Tnp_IS66', 'Gram_pos_anchor', 'MucBP',
               'Methylase_S', 'MFS_1', 'Muc_B2', 'MucBP_2', 'Nramp', 'PIN_3',
               'Relaxase', 'RelB', 'Rib', 'YkuD')
  statplot = statDF_dNdSc_big[statDF_dNdSc_big$gene_name %in% keyGenes, ]
  nrow(statplot)
  statplot$Species_short = paste0(str_sub(statplot$Species, 1, 1), '. ',
                                  str_split_fixed(statplot$Species, ' ', 2)[, 2])
  statplot$Species_gene = paste0(
    statplot$gene_name,
    ' (', statplot$Species_short, ')'
  )
  statplot = statplot[, c('NSp_bigger_1_freq', 'meanGen_pNpS',
                          'NSd_bigger_1_freq', 'meanGen_dNdS',
                          'NSdc_bigger_1_freq', 'meanGen_dNdS_common',
                          'Species_gene', 'gene')]
  x = 'meanGen_dNdS_common'
  statplot$x = statplot[, x]
  statplot = statplot[order(statplot$x, decreasing = T), ]
  statplot = statplot[!duplicated(statplot$Species_gene), ]
  statplot$Species_gene = factor(statplot$Species_gene, 
                                 levels = statplot$Species_gene[order(statplot$x)])
  dP3 = ggplot(statplot, aes(x = x, y = Species_gene)) +
    geom_col(width = 0.8, color = '#386CB0', fill = '#386CB0', alpha = 0.9) +
    labs(x = 'common dN/dS', y = 'Gene (species)') +
    mytheme3 + theme(plot.margin = unit(c(0.2,1,0.2,0.2), 'cm'))
  print(dP3)
}
cairo_pdf('./plot/PRJEB34536_dnds_bigger1_barplot.pdf', width = 5.21, height = 5.60, onefile = T)
print(dP)
print(dP2)
print(dP3)
dev.off()

