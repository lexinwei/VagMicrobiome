source('./script/library_function.R')
load('./data/RData/MTMG_metaData.RData')
load('./data/RData/MTMG_statDF3_all.RData')
load('./data/RData/colorList.RData')
load('./data/RData/MTMG_geneDF.RData')
if(F){
  baseDir = './instrain_res/MTMG/gene'
  statList = list.files(baseDir,pattern = '^inStrain_gene_statistics')
  dataList = setdiff(list.files(baseDir)[endsWith(list.files(baseDir), '.csv')], statList)
  
  geneDF = data.frame()
  for (sp in dataList) {
    gene_df = read.csv(paste0(baseDir, '/', sp))
    gene_df$Species = str_extract(sp, 's__.*') %>% str_remove('\\.fasta.csv') %>% 
      str_remove('s__') %>%
      str_replace_all('_', ' ')
    geneDF = rbind(geneDF, gene_df)
  }
  table(geneDF$Species)
  geneDF_sub = geneDF[geneDF$gene == 'NZ_CP039266.1_1594', ]
  
  if(T){
    geneDF$seq_ID = geneDF$sample
    geneDF$pt_ID.u = metaData$UID[match(geneDF$seq_ID, metaData$seq_ID)]
    geneDF$Ethnicity = metaData$Ethnictiy[match(geneDF$seq_ID, metaData$seq_ID)]
    geneDF$Age = metaData$Age_cat[match(geneDF$seq_ID, metaData$seq_ID)]
    geneDF$pt_ID = metaData$SID[match(geneDF$seq_ID, metaData$seq_ID)]
    table(geneDF$pt_ID)
  }
  save(geneDF, file = './data/RData/MTMG_geneDF.RData')
}
######### evolutionary metrics for each gene (pN/pS)  along genomes loci #####
load('./data/RData/MTMG_geneDF.RData')
table(geneDF$Species)
a  = geneDF[grep('muc', geneDF$gene_description, ignore.case = T), c('gene_name', 'gene_description')]
a = a[!duplicated(a), ]
keyGenes = c('DDE_Tnp_1', 'DDE_Tnp_IS66', 'Gram_pos_anchor',
  'Methylase_S', 'MFS_1', 'Muc_B2', 'MucBP_2',  'MucBP','Nramp', 'PIN_3',
  'Relaxase', 'RelB', 'Rib', 'YkuD')
keyGenes = c('Muc_B2', 'MucBP', 'MucBP_2', 'Gram_pos_anchor')
geneDF$Key_gene =  geneDF$gene_name
geneDF$Key_gene[!geneDF$Key_gene %in% keyGenes] = 'Other'
table(geneDF$Key_gene)
# myPal = c('#E4211C', '#FF7F00', '#4DAF4A', '#377EB8')
# myPal = c('#FBB4AE', '#FED9A6', '#CCEBC5', '#B3CDE3')
myPal = c('#FB8072', '#FDB462', '#B3DE68', '#80B1D3')
names(myPal) = keyGenes
myShape = c(22:25)
names(myShape) = keyGenes
table(geneDF$Ethnicity)
geneDF$Ethnicity[is.na(geneDF$Ethnicity)] = 'Other'
geneDF$Ethnicity %<>% factor(., levels = c('White', 'Asian', 'Latina', 'Black', 'Other'))
cairo_pdf(file =  './plot/MTMG_gene_for_each_species_evolutionary_metrics_pNpS_facet_byEthniticy.pdf', width = 11.14, height = 4, onefile = T)
for(sp in unique(geneDF$Species)){
  geneDF_sub = geneDF[geneDF$Species == sp  & (!is.na(geneDF$pNpS)), ]
  length(unique(geneDF_sub$sample))
  length(unique(geneDF_sub$gene))
  range(geneDF_sub$breadth_minCov)
  geneOrd = data.frame(
    Gene = unique(geneDF_sub$gene),
    Order = naturalorder(unique(geneDF_sub$gene))
  )
  geneDF_sub$gene_order = geneOrd$Order[match(geneDF_sub$gene, geneOrd$Gene)]
  geneDF_sub = geneDF_sub[(!is.na(geneDF_sub$pNpS)), ]
  geneDF_sub$density = get_density(geneDF_sub$gene_order, geneDF_sub$pNpS, n = 100)
  pp1 = ggplot() + 
    geom_point(geneDF_sub[geneDF_sub$Key_gene == 'Other', ], 
               mapping = aes(x = gene_order, y = pNpS, color = density),
               stroke = 0.1, alpha = 0.5) +
    geom_point(geneDF_sub[geneDF_sub$Key_gene != 'Other', ], size = 2,
               mapping = aes(x = gene_order, y = pNpS, shape = Key_gene, fill = Key_gene),
               color = 'black', stroke = 0.5, alpha = 0.9) +
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'red', linewidth = 0.6) +
    labs(x = 'Gene', y = 'pN/pS', size = 'pN/pS', fill = 'Gene',color = 'Density',   shape = 'Gene', title = sp) +
    scale_y_continuous(expand = expansion(mult = c(0.025, 0.1))) +
    # scale_size_continuous(range = c(0.1, 3)) +
    scale_color_viridis() + 
    scale_fill_manual(values = myPal) +
    scale_shape_manual(values = myShape) +
    facet_grid(.~Ethnicity) + mytheme + theme(strip.text = element_text(size = 13), aspect.ratio = 1.5/1)
  print(pp1)
  geneDF_sub_mean = geneDF_sub %>% 
    group_by(Ethnicity, gene_order) %>% summarise(meanpNpS = mean(pNpS, na.rm = T))
  geneDF_sub_mean$Key_gene = geneDF_sub$Key_gene[match(geneDF_sub_mean$gene_order, geneDF_sub$gene_order)]
  geneDF_sub_mean = geneDF_sub_mean[!is.na(geneDF_sub_mean$meanpNpS), ]
  geneDF_sub_mean$density = get_density(geneDF_sub_mean$gene_order, geneDF_sub_mean$meanpNpS, n = 100)
  pp2 = ggplot() + 
    geom_point(geneDF_sub_mean[geneDF_sub_mean$Key_gene == 'Other', ], 
               mapping = aes(x = gene_order, y = meanpNpS, color = density),
               stroke = 0.1, alpha = 0.5) +
    geom_point(geneDF_sub_mean[geneDF_sub_mean$Key_gene != 'Other', ], size = 3,
               mapping = aes(x = gene_order, y = meanpNpS, shape = Key_gene, fill = Key_gene),
               color = 'black', stroke = 0.7, alpha = 0.9) +
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'red', linewidth = 0.6) +
    labs(x = 'Gene', y = 'Avg. pN/pS', size = 'pN/pS', fill = 'Gene', color = 'Density', shape = 'Gene', title = sp) +
    scale_y_continuous(expand = expansion(mult = c(0.025, 0.1))) +
    # scale_size_continuous(range = c(0.1, 3)) +
    scale_color_viridis() + 
    scale_fill_manual(values = myPal) +
    scale_shape_manual(values = myShape) +
    facet_grid(.~Ethnicity) + mytheme + theme(strip.text = element_text(size = 13), aspect.ratio = 1.5/1)
  print(pp2)
}
dev.off()

cairo_pdf(file =  './plot/MTMG_gene_for_each_species_evolutionary_metrics_2.pdf', width = 7.09, height = 3.49, onefile = T)
for(sp in unique(geneDF$Species)){
  geneDF_sub = geneDF[geneDF$Species == sp, ]
  length(unique(geneDF_sub$sample))
  length(unique(geneDF_sub$gene))
  range(geneDF_sub$breadth_minCov)
  geneOrd = data.frame(
    Gene = unique(geneDF_sub$gene),
    Order = naturalorder(unique(geneDF_sub$gene))
  )
  geneDF_sub$gene_order = geneOrd$Order[match(geneDF_sub$gene, geneOrd$Gene)]
  geneDF_sub = geneDF_sub[(!is.na(geneDF_sub$pNpS)), ]
  geneDF_sub$density = get_density(geneDF_sub$gene_order, geneDF_sub$pNpS, n = 100)
  pp1 = ggplot() + 
    geom_point(geneDF_sub[geneDF_sub$Key_gene == 'Other', ], 
               mapping = aes(x = gene_order, y = pNpS, color = density),
               stroke = 0.1, alpha = 0.5) +
    geom_point(geneDF_sub[geneDF_sub$Key_gene != 'Other', ], size = 2,
               mapping = aes(x = gene_order, y = pNpS, shape = Key_gene, fill = Key_gene),
               color = 'black', stroke = 0.5, alpha = 0.9) +
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'red', linewidth = 0.6) +
    labs(x = 'Gene', y = 'pN/pS', size = 'pN/pS', fill = 'Gene',color = 'Density',   shape = 'Gene', title = sp) +
    scale_y_continuous(expand = expansion(mult = c(0.025, 0.1))) +
    # scale_size_continuous(range = c(0.1, 3)) +
    scale_color_viridis() + 
    scale_fill_manual(values = myPal) +
    scale_shape_manual(values = myShape) +
    mytheme + theme(strip.text = element_text(size = 13), aspect.ratio = 1)
  print(pp1)
  geneDF_sub_mean = geneDF_sub %>% 
    group_by(gene_order) %>% summarise(N = n(), meanpNpS = mean(pNpS, na.rm = T))
  geneDF_sub_mean$Key_gene = geneDF_sub$Key_gene[match(geneDF_sub_mean$gene_order, geneDF_sub$gene_order)]
  geneDF_sub_mean = geneDF_sub_mean[!is.na(geneDF_sub_mean$meanpNpS), ]
  geneDF_sub_mean$density = get_density(geneDF_sub_mean$gene_order, geneDF_sub_mean$meanpNpS, n = 100)
  pp2 = ggplot() + 
    geom_point(geneDF_sub_mean[geneDF_sub_mean$Key_gene == 'Other', ], 
               mapping = aes(x = gene_order, y = meanpNpS, color = density),
               stroke = 0.1, alpha = 0.5) +
    geom_point(geneDF_sub_mean[geneDF_sub_mean$Key_gene != 'Other', ], size = 3,
               mapping = aes(x = gene_order, y = meanpNpS, shape = Key_gene, fill = Key_gene),
               color = 'black', stroke = 0.7, alpha = 0.9) +
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'red', linewidth = 0.6) +
    labs(x = 'Gene', y = 'Avg. pN/pS', size = 'pN/pS', fill = 'Gene', color = 'Density', shape = 'Gene', title = sp) +
    scale_y_continuous(expand = expansion(mult = c(0.025, 0.1))) +
    # scale_size_continuous(range = c(0.1, 3)) +
    scale_color_viridis() + 
    scale_fill_manual(values = myPal) +
    scale_shape_manual(values = myShape) +
    mytheme + theme(strip.text = element_text(size = 13), aspect.ratio = 1)
  print(pp2)
}
dev.off()
######### evolutionary metrics for each gene (dN/dS) along genomes loci #####
load('./data/RData/MTMG_geneDF.RData')
table(geneDF$Species)
a  = geneDF[grep('muc', geneDF$gene_description, ignore.case = T), c('gene_name', 'gene_description')]
a = a[!duplicated(a), ]
keyGenes = c('DDE_Tnp_1', 'DDE_Tnp_IS66', 'Gram_pos_anchor',
             'Methylase_S', 'MFS_1', 'Muc_B2', 'MucBP_2',  'MucBP','Nramp', 'PIN_3',
             'Relaxase', 'RelB', 'Rib', 'YkuD')
keyGenes = c('Muc_B2', 'MucBP', 'MucBP_2', 'Gram_pos_anchor')
geneDF$Key_gene =  geneDF$gene_name
geneDF$Key_gene[!geneDF$Key_gene %in% keyGenes] = 'Other'
table(geneDF$Key_gene)
# myPal = c('#E4211C', '#FF7F00', '#4DAF4A', '#377EB8')
# myPal = c('#FBB4AE', '#FED9A6', '#CCEBC5', '#B3CDE3')
myPal = c('#FB8072', '#FDB462', '#B3DE68', '#80B1D3')
names(myPal) = keyGenes
myShape = c(22:25)
names(myShape) = keyGenes
cairo_pdf(file =  './plot/MTMG_gene_for_each_species_evolutionary_metrics_dNdS_facet_byEthniticy.pdf', 
          width = 11.14, height = 4, onefile = T)
for(sp in unique(geneDF$Species)){
  geneDF_sub = geneDF[geneDF$Species == sp, ]
  length(unique(geneDF_sub$sample))
  length(unique(geneDF_sub$gene))
  range(geneDF_sub$breadth_minCov)
  geneOrd = data.frame(
    Gene = unique(geneDF_sub$gene),
    Order = naturalorder(unique(geneDF_sub$gene))
  )
  geneDF_sub$gene_order = geneOrd$Order[match(geneDF_sub$gene, geneOrd$Gene)]
  geneDF_sub = geneDF_sub[(!is.na(geneDF_sub$dNdS)), ]
  geneDF_sub$density = get_density(geneDF_sub$gene_order, geneDF_sub$dNdS, n = 100)
  pp1 = ggplot() + 
    geom_point(geneDF_sub[geneDF_sub$Key_gene == 'Other', ], 
               mapping = aes(x = gene_order, y = dNdS, color = density),
               stroke = 0.1, alpha = 0.5) +
    geom_point(geneDF_sub[geneDF_sub$Key_gene != 'Other', ], size = 2,
               mapping = aes(x = gene_order, y = dNdS, shape = Key_gene, fill = Key_gene),
               color = 'black', stroke = 0.5, alpha = 0.9) +
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'red', linewidth = 0.6) +
    labs(x = 'Gene', y = 'dN/dS', size = 'dN/dS', fill = 'Gene',color = 'Density',   shape = 'Gene', title = sp) +
    scale_y_continuous(expand = expansion(mult = c(0.025, 0.1))) +
    # scale_size_continuous(range = c(0.1, 3)) +
    scale_color_viridis() + 
    scale_fill_manual(values = myPal) +
    scale_shape_manual(values = myShape) +
    facet_grid(.~Ethnicity) + mytheme + theme(strip.text = element_text(size = 13), aspect.ratio = 1.5/1)
  print(pp1)
  geneDF_sub_mean = geneDF_sub %>% 
    group_by(Ethnicity, gene_order) %>% summarise(meandNdS = mean(dNdS, na.rm = T))
  geneDF_sub_mean$Key_gene = geneDF_sub$Key_gene[match(geneDF_sub_mean$gene_order, geneDF_sub$gene_order)]
  geneDF_sub_mean = geneDF_sub_mean[!is.na(geneDF_sub_mean$meandNdS), ]
  geneDF_sub_mean$density = get_density(geneDF_sub_mean$gene_order, geneDF_sub_mean$meandNdS, n = 100)
  pp2 = ggplot() + 
    geom_point(geneDF_sub_mean[geneDF_sub_mean$Key_gene == 'Other', ], 
               mapping = aes(x = gene_order, y = meandNdS, color = density),
               stroke = 0.1, alpha = 0.5) +
    geom_point(geneDF_sub_mean[geneDF_sub_mean$Key_gene != 'Other', ], size = 3,
               mapping = aes(x = gene_order, y = meandNdS, shape = Key_gene, fill = Key_gene),
               color = 'black', stroke = 0.7, alpha = 0.9) +
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'red', linewidth = 0.6) +
    labs(x = 'Gene', y = 'Avg. dN/dS', size = 'dN/dS', fill = 'Gene', color = 'Density', shape = 'Gene', title = sp) +
    scale_y_continuous(expand = expansion(mult = c(0.025, 0.1))) +
    # scale_size_continuous(range = c(0.1, 3)) +
    scale_color_viridis() + 
    scale_fill_manual(values = myPal) +
    scale_shape_manual(values = myShape) +
    facet_grid(.~Ethnicity) + mytheme + theme(strip.text = element_text(size = 13), aspect.ratio = 1.5/1)
  print(pp2)
}
dev.off()

cairo_pdf(file =  './plot/MTMG_gene_for_each_species_evolutionary_metrics_2_dNdS.pdf', width = 7.09, height = 3.49, onefile = T)
for(sp in unique(geneDF$Species)){
  geneDF_sub = geneDF[geneDF$Species == sp, ]
  length(unique(geneDF_sub$sample))
  length(unique(geneDF_sub$gene))
  range(geneDF_sub$breadth_minCov)
  geneOrd = data.frame(
    Gene = unique(geneDF_sub$gene),
    Order = naturalorder(unique(geneDF_sub$gene))
  )
  geneDF_sub$gene_order = geneOrd$Order[match(geneDF_sub$gene, geneOrd$Gene)]
  geneDF_sub = geneDF_sub[(!is.na(geneDF_sub$dNdS)), ]
  geneDF_sub$density = get_density(geneDF_sub$gene_order, geneDF_sub$dNdS, n = 100)
  pp1 = ggplot() + 
    geom_point(geneDF_sub[geneDF_sub$Key_gene == 'Other', ], 
               mapping = aes(x = gene_order, y = dNdS, color = density),
               stroke = 0.1, alpha = 0.5) +
    geom_point(geneDF_sub[geneDF_sub$Key_gene != 'Other', ], size = 2,
               mapping = aes(x = gene_order, y = dNdS, shape = Key_gene, fill = Key_gene),
               color = 'black', stroke = 0.5, alpha = 0.9) +
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'red', linewidth = 0.6) +
    labs(x = 'Gene', y = 'dN/dS', size = 'dN/dS', fill = 'Gene',color = 'Density',   shape = 'Gene', title = sp) +
    scale_y_continuous(expand = expansion(mult = c(0.025, 0.1))) +
    # scale_size_continuous(range = c(0.1, 3)) +
    scale_color_viridis() + 
    scale_fill_manual(values = myPal) +
    scale_shape_manual(values = myShape) +
    mytheme + theme(strip.text = element_text(size = 13), aspect.ratio = 1)
  print(pp1)
  geneDF_sub_mean = geneDF_sub %>% 
    group_by(gene_order) %>% summarise(N = n(), meandNdS = mean(dNdS, na.rm = T))
  geneDF_sub_mean$Key_gene = geneDF_sub$Key_gene[match(geneDF_sub_mean$gene_order, geneDF_sub$gene_order)]
  geneDF_sub_mean = geneDF_sub_mean[!is.na(geneDF_sub_mean$meandNdS), ]
  geneDF_sub_mean$density = get_density(geneDF_sub_mean$gene_order, geneDF_sub_mean$meandNdS, n = 100)
  pp2 = ggplot() + 
    geom_point(geneDF_sub_mean[geneDF_sub_mean$Key_gene == 'Other', ], 
               mapping = aes(x = gene_order, y = meandNdS, color = density),
               stroke = 0.1, alpha = 0.5) +
    geom_point(geneDF_sub_mean[geneDF_sub_mean$Key_gene != 'Other', ], size = 3,
               mapping = aes(x = gene_order, y = meandNdS, shape = Key_gene, fill = Key_gene),
               color = 'black', stroke = 0.7, alpha = 0.9) +
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'red', linewidth = 0.6) +
    labs(x = 'Gene', y = 'Avg. dN/dS', size = 'dN/dS', fill = 'Gene', color = 'Density', shape = 'Gene', title = sp) +
    scale_y_continuous(expand = expansion(mult = c(0.025, 0.1))) +
    # scale_size_continuous(range = c(0.1, 3)) +
    scale_color_viridis() + 
    scale_fill_manual(values = myPal) +
    scale_shape_manual(values = myShape) +
    mytheme + theme(strip.text = element_text(size = 13), aspect.ratio = 1)
  print(pp2)
}
dev.off()

######### evolutionary metrics for each gene (pi)  along genomes loci#####
load('./data/RData/MTMG_geneDF.RData')
table(geneDF$Species)
a  = geneDF[grep('muc', geneDF$gene_description, ignore.case = T), c('gene_name', 'gene_description')]
a = a[!duplicated(a), ]
keyGenes = c('DDE_Tnp_1', 'DDE_Tnp_IS66', 'Gram_pos_anchor',
             'Methylase_S', 'MFS_1', 'Muc_B2', 'MucBP_2',  'MucBP','Nramp', 'PIN_3',
             'Relaxase', 'RelB', 'Rib', 'YkuD')
keyGenes = c('Muc_B2', 'MucBP', 'MucBP_2', 'Gram_pos_anchor')
geneDF$Key_gene =  geneDF$gene_name
geneDF$Key_gene[!geneDF$Key_gene %in% keyGenes] = 'Other'
table(geneDF$Key_gene)
# myPal = c('#E4211C', '#FF7F00', '#4DAF4A', '#377EB8')
# myPal = c('#FBB4AE', '#FED9A6', '#CCEBC5', '#B3CDE3')
myPal = c('#FB8072', '#FDB462', '#B3DE68', '#80B1D3')
names(myPal) = keyGenes
myShape = c(22:25)
names(myShape) = keyGenes
cairo_pdf(file =  './plot/gene_for_each_species_evolutionary_metrics_Pi_facet_byEthniticy.pdf', width = 11.14, height = 4, onefile = T)
for(sp in unique(geneDF$Species)){
  geneDF_sub = geneDF[geneDF$Species == sp, ]
  length(unique(geneDF_sub$sample))
  length(unique(geneDF_sub$gene))
  range(geneDF_sub$breadth_minCov)
  geneOrd = data.frame(
    Gene = unique(geneDF_sub$gene),
    Order = naturalorder(unique(geneDF_sub$gene))
  )
  geneDF_sub$gene_order = geneOrd$Order[match(geneDF_sub$gene, geneOrd$Gene)]
  geneDF_sub = geneDF_sub[(!is.na(geneDF_sub$nucl_diversity)), ]
  geneDF_sub$density = get_density(geneDF_sub$gene_order, geneDF_sub$nucl_diversity, n = 100)
  pp1 = ggplot() + 
    geom_point(geneDF_sub[geneDF_sub$Key_gene == 'Other', ], 
               mapping = aes(x = gene_order, y = nucl_diversity, color = density),
               stroke = 0.1, alpha = 0.5) +
    geom_point(geneDF_sub[geneDF_sub$Key_gene != 'Other', ], size = 2,
               mapping = aes(x = gene_order, y = nucl_diversity, shape = Key_gene, fill = Key_gene),
               color = 'black', stroke = 0.5, alpha = 0.9) +
    # geom_hline(yintercept = 1, linetype = 'dashed', color = 'red', linewidth = 0.3) +
    labs(x = 'Gene', y = expression(pi), fill = 'Gene',color = 'Density',   shape = 'Gene', title = sp) +
    scale_y_continuous(expand = expansion(mult = c(0.025, 0.1))) +
    # scale_size_continuous(range = c(0.1, 3)) +
    scale_color_viridis() + 
    scale_fill_manual(values = myPal) +
    scale_shape_manual(values = myShape) +
    facet_grid(.~Ethnicity) + mytheme + theme(strip.text = element_text(size = 13), aspect.ratio = 1.5/1)
  print(pp1)
  geneDF_sub_mean = geneDF_sub %>% 
    group_by(Ethnicity, gene_order) %>% summarise(meanPi = mean(nucl_diversity, na.rm = T))
  geneDF_sub_mean$Key_gene = geneDF_sub$Key_gene[match(geneDF_sub_mean$gene_order, geneDF_sub$gene_order)]
  geneDF_sub_mean = geneDF_sub_mean[!is.na(geneDF_sub_mean$meanPi), ]
  geneDF_sub_mean$density = get_density(geneDF_sub_mean$gene_order, geneDF_sub_mean$meanPi, n = 100)
  pp2 = ggplot() + 
    geom_point(geneDF_sub_mean[geneDF_sub_mean$Key_gene == 'Other', ], 
               mapping = aes(x = gene_order, y = meanPi, color = density),
               stroke = 0.1, alpha = 0.5) +
    geom_point(geneDF_sub_mean[geneDF_sub_mean$Key_gene != 'Other', ], size = 3,
               mapping = aes(x = gene_order, y = meanPi, shape = Key_gene, fill = Key_gene),
               color = 'black', stroke = 0.7, alpha = 0.9) +
    # geom_hline(yintercept = 1, linetype = 'dashed', color = 'red', linewidth = 0.3) +
    labs(x = 'Gene', y = expression(~ "Avg." ~ pi), size = expression(pi), fill = 'Gene', color = 'Density', shape = 'Gene', title = sp) +
    scale_y_continuous(expand = expansion(mult = c(0.025, 0.1))) +
    # scale_size_continuous(range = c(0.1, 3)) +
    scale_color_viridis() + 
    scale_fill_manual(values = myPal) +
    scale_shape_manual(values = myShape) +
    facet_grid(.~Ethnicity) + mytheme + theme(strip.text = element_text(size = 13), aspect.ratio = 1.5/1)
  print(pp2)
}
dev.off()

cairo_pdf(file =  './plot/gene_for_each_species_evolutionary_metrics_2_Pi.pdf', width = 7.09, height = 3.49, onefile = T)
for(sp in unique(geneDF$Species)){
  geneDF_sub = geneDF[geneDF$Species == sp, ]
  length(unique(geneDF_sub$sample))
  length(unique(geneDF_sub$gene))
  range(geneDF_sub$breadth_minCov)
  geneOrd = data.frame(
    Gene = unique(geneDF_sub$gene),
    Order = naturalorder(unique(geneDF_sub$gene))
  )
  geneDF_sub$gene_order = geneOrd$Order[match(geneDF_sub$gene, geneOrd$Gene)]
  geneDF_sub = geneDF_sub[(!is.na(geneDF_sub$nucl_diversity)), ]
  geneDF_sub$density = get_density(geneDF_sub$gene_order, geneDF_sub$nucl_diversity, n = 100)
  pp1 = ggplot() + 
    geom_point(geneDF_sub[geneDF_sub$Key_gene == 'Other', ], 
               mapping = aes(x = gene_order, y = nucl_diversity, color = density),
               stroke = 0.1, alpha = 0.5) +
    geom_point(geneDF_sub[geneDF_sub$Key_gene != 'Other', ], size = 2,
               mapping = aes(x = gene_order, y = nucl_diversity, shape = Key_gene, fill = Key_gene),
               color = 'black', stroke = 0.5, alpha = 0.9) +
    # geom_hline(yintercept = 1, linetype = 'dashed', color = 'red', linewidth = 0.3) +
    labs(x = 'Gene', y = expression(pi), fill = 'Gene',color = 'Density',   shape = 'Gene', title = sp) +
    scale_y_continuous(expand = expansion(mult = c(0.025, 0.1))) +
    # scale_size_continuous(range = c(0.1, 3)) +
    scale_color_viridis() + 
    scale_fill_manual(values = myPal) +
    scale_shape_manual(values = myShape) +
    mytheme + theme(strip.text = element_text(size = 13), aspect.ratio = 1)
  print(pp1)
  geneDF_sub_mean = geneDF_sub %>% 
    group_by(gene_order) %>% summarise(N = n(), meanPi = mean(nucl_diversity, na.rm = T))
  geneDF_sub_mean$Key_gene = geneDF_sub$Key_gene[match(geneDF_sub_mean$gene_order, geneDF_sub$gene_order)]
  geneDF_sub_mean = geneDF_sub_mean[!is.na(geneDF_sub_mean$meanPi), ]
  geneDF_sub_mean$density = get_density(geneDF_sub_mean$gene_order, geneDF_sub_mean$meanPi, n = 100)
  pp2 = ggplot() + 
    geom_point(geneDF_sub_mean[geneDF_sub_mean$Key_gene == 'Other', ], 
               mapping = aes(x = gene_order, y = meanPi, color = density),
               stroke = 0.1, alpha = 0.5) +
    geom_point(geneDF_sub_mean[geneDF_sub_mean$Key_gene != 'Other', ], size = 3,
               mapping = aes(x = gene_order, y = meanPi, shape = Key_gene, fill = Key_gene),
               color = 'black', stroke = 0.7, alpha = 0.9) +
    # geom_hline(yintercept = 1, linetype = 'dashed', color = 'red', linewidth = 0.3) +
    labs(x = 'Gene', y = expression(~ "Avg." ~ pi), fill = 'Gene', color = 'Density', shape = 'Gene', title = sp) +
    scale_y_continuous(expand = expansion(mult = c(0.025, 0.1))) +
    # scale_size_continuous(range = c(0.1, 3)) +
    scale_color_viridis() + 
    scale_fill_manual(values = myPal) +
    scale_shape_manual(values = myShape) +
    mytheme + theme(strip.text = element_text(size = 13), aspect.ratio = 1)
  print(pp2)
}
dev.off()

######### evolutionary metrics for keyGenes #####
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
        cat(gs, ', P =', signif(corT_Pearson$p.value, 3), 'R =', signif(corT_Pearson$estimate, 3), '\n')
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
load('./data/RData/geneDF.RData')
keyGenes = c('Muc_B2', 'MucBP', 'MucBP_2', 'Gram_pos_anchor')
geneDF_key = geneDF[geneDF$gene_name %in% keyGenes, ]
table(geneDF_key$Species)
cairo_pdf('./plot/keyGenes_pi_vs_pNpS.pdf', width = 13, height = 3.39, onefile = T)
if(T){
  geneDF_key$gene_name %<>% factor(., levels = keyGenes)
  pp = ggplot(geneDF_key, aes(x = nucl_diversity, y = pNpS, fill = gene_name)) +
    geom_point(aes(shape = gene_name), size = 2.5, alpha = 0.75, stroke = 0.2) + 
    geom_smooth(aes(color = gene_name), method = 'lm', se = F) + 
    scale_shape_manual(values = myShape, name = 'Gene') +
    scale_color_manual(values = myPal, name = 'Gene') +
    scale_fill_manual(values = myPal, name = 'Gene') +
    labs(x = expression(pi), y = 'pN/pS') +
    mytheme + theme(aspect.ratio = 1)
  print(pp)
  corFitFun(df = geneDF_key, xval = 'nucl_diversity', yval = 'pNpS', gval =  'gene_name',xLab = 'pi', yLab = 'pN/pS')
  
  pp2 = ggplot(geneDF_key, aes(x = nucl_diversity, y = pNpS, fill = gene_name)) +
    geom_point(aes(shape = gene_name), size = 2.5, alpha = 0.75, stroke = 0.2) + 
    geom_smooth(aes(color = gene_name), method = 'lm', se = F) + 
    scale_shape_manual(values = myShape, name = 'Gene') +
    scale_color_manual(values = myPal, name = 'Gene') +
    scale_fill_manual(values = myPal, name = 'Gene') +
    labs(x = expression(pi), y = 'pN/pS') +
    facet_grid(.~Ethnicity, scales = 'free') +
    mytheme + theme(aspect.ratio = 1, strip.text = element_text(size = 13),
                    axis.text.x = element_text(angle = 45, hjust = 1))
  print(pp2)
  corFitFun(df = geneDF_key, xval = 'nucl_diversity', yval = 'pNpS', zval = 'Ethnicity',
            gval =  'gene_name',xLab = 'pi', yLab = 'pN/pS')
  
  pp1 = ggplot(geneDF_key, aes(x = nucl_diversity, y = dNdS, fill = gene_name)) +
    geom_point(aes(shape = gene_name), size = 2.5, alpha = 0.75, stroke = 0.2) + 
    geom_smooth(aes(color = gene_name), method = 'lm', se = F) + 
    scale_shape_manual(values = myShape, name = 'Gene') +
    scale_color_manual(values = myPal, name = 'Gene') +
    scale_fill_manual(values = myPal, name = 'Gene') +
    labs(x = expression(pi), y = 'dN/dS') +
    mytheme + theme(aspect.ratio = 1)
  print(pp1)
  pp3 = ggplot(geneDF_key, aes(x = nucl_diversity, y = dNdS, fill = gene_name)) +
    geom_point(aes(shape = gene_name), size = 2.5, alpha = 0.75, stroke = 0.2) + 
    geom_smooth(aes(color = gene_name), method = 'lm', se = F) + 
    scale_shape_manual(values = myShape, name = 'Gene') +
    scale_color_manual(values = myPal, name = 'Gene') +
    scale_fill_manual(values = myPal, name = 'Gene') +
    labs(x = expression(pi), y = 'dN/dS') +
    facet_grid(.~Ethnicity, scales = 'free') +
    mytheme + theme(aspect.ratio = 1, strip.text = element_text(size = 13),
                    axis.text.x = element_text(angle = 30, hjust = 1))
  print(pp3)
}
dev.off()

cairo_pdf('./plot/keyGenes_pi_vs_pNpS_for_each_species.pdf', width = 10.80, height = 3.39, onefile = T)
sp = 'Lactobacillus crispatus'
for(sp in unique(geneDF_key$Species)){
  geneDF_key_sub =  geneDF_key[geneDF_key$Species == sp, ]
  geneDF_key_sub = geneDF_key_sub[!is.na(geneDF_key_sub$pNpS), ]
  pp = ggplot(geneDF_key_sub, aes(x = nucl_diversity, y = pNpS, fill = gene_name)) +
    geom_point(aes(shape = gene_name), size = 3, alpha = 0.75, stroke = 0.2) + 
    geom_smooth(aes(color = gene_name), method = 'lm', se = F) + 
    scale_shape_manual(values = myShape, name = 'Gene') +
    scale_color_manual(values = myPal, name = 'Gene') +
    scale_fill_manual(values = myPal, name = 'Gene') +
    labs(x = expression(pi), y = 'pN/pS', title = sp) +
    mytheme + theme(aspect.ratio = 1)
  print(pp)
  corFitFun(df = geneDF_key_sub, xval = 'nucl_diversity', yval = 'pNpS',
            gval =  'gene_name',xLab = 'pi', yLab = 'pN/pS')
  
  pp2 = ggplot(geneDF_key_sub, aes(x = nucl_diversity, y = pNpS, fill = gene_name)) +
    geom_point(aes(shape = gene_name), size = 3, alpha = 0.75, stroke = 0.2) + 
    geom_smooth(aes(color = gene_name), method = 'lm', se = F) + 
    scale_shape_manual(values = myShape, name = 'Gene') +
    scale_color_manual(values = myPal, name = 'Gene') +
    scale_fill_manual(values = myPal, name = 'Gene') +
    labs(x = expression(pi), y = 'pN/pS', title = sp) +
    facet_grid(.~Ethnicity, scales = 'free') +
    mytheme + theme(aspect.ratio = 1, strip.text = element_text(size = 13),
                    axis.text.x = element_text(angle = 30, hjust = 1))
  print(pp2)
  corFitFun(df = geneDF_key_sub, xval = 'nucl_diversity', yval = 'pNpS', zval = 'Ethnicity',
            gval =  'gene_name',xLab = 'pi', yLab = 'pN/pS')
  
  geneDF_key_sub =  geneDF_key[geneDF_key$Species == sp, ]
  geneDF_key_sub = geneDF_key_sub[!is.na(geneDF_key_sub$dNdS), ]
  pp1 = ggplot(geneDF_key_sub, aes(x = nucl_diversity, y = dNdS, fill = gene_name)) +
    geom_point(aes(shape = gene_name), size = 3, alpha = 0.75, stroke = 0.2) + 
    geom_smooth(aes(color = gene_name), method = 'lm', se = F) + 
    scale_shape_manual(values = myShape, name = 'Gene') +
    scale_color_manual(values = myPal, name = 'Gene') +
    scale_fill_manual(values = myPal, name = 'Gene') +
    labs(x = expression(pi), y = 'dN/dS', title = sp) +
    mytheme + theme(aspect.ratio = 1)
  print(pp1)
  pp3 = ggplot(geneDF_key_sub, aes(x = nucl_diversity, y = dNdS, fill = gene_name)) +
    geom_point(aes(shape = gene_name), size = 3, alpha = 0.75, stroke = 0.2) + 
    geom_smooth(aes(color = gene_name), method = 'lm', se = F) + 
    scale_shape_manual(values = myShape, name = 'Gene') +
    scale_color_manual(values = myPal, name = 'Gene') +
    scale_fill_manual(values = myPal, name = 'Gene') +
    labs(x = expression(pi), y = 'dN/dS', title = sp) +
    facet_grid(.~Ethnicity, scales = 'free') +
    mytheme + theme(aspect.ratio = 1, strip.text = element_text(size = 13),
                    axis.text.x = element_text(angle = 30, hjust = 1))
  print(pp3)
  # sp = 'Lactobacillus coleohominis'
  # sp = 'Lactobacillus gasseri'
  # geneDF_key_sub =  geneDF_key[geneDF_key$Species == sp, ]
  # table(geneDF_key_sub$gene)
  # ggplot(geneDF_key_sub, aes(x = nucl_diversity, y = pNpS, fill = paste0(gene, '|', gene_name))) +
  #   geom_point(size = 3, shape = 21,alpha = 0.75, stroke = 0.2) +
  #   geom_smooth(aes(color = paste0(gene, '|', gene_name)), method = 'lm', se = F) +
  #   # scale_shape_manual(values = myShape) +
  #   # scale_color_manual(values = myPal) +
  #   # scale_fill_manual(values = myPal) +
  #   labs(x = expression(pi), y = 'pN/pS', title = sp) +
  #   facet_grid(.~Ethnicity) +
  #   mytheme
}
dev.off()

##### find gene whose pi and dN/dS have trends among variables ####
plotFlag = F
#### single variable comparison lm & lme #####
# keySpecies = c("Fannyhessea vaginae", "Gardnerella vaginalis",
#                # "Lactobacillus coleohominis",
#                "Lactobacillus crispatus", "Lactobacillus gasseri", "Lactobacillus iners", 
#                "Lactobacillus jensenii")
# 0
# statDFsel = statDF3_all[statDF3_all$meanGen_dNdS_common > 1 & !is.na(statDF3_all$meanGen_dNdS_common) &
#                           statDF3_all$Species %in% keySpecies &
#                           statDF3_all$sampleNum > 20, ]
# 1
# keyGenes = c('DDE_Tnp_1', 'DDE_Tnp_IS66', 'Gram_pos_anchor',
#              'Methylase_S', 'MFS_1', 'Muc_B2', 'MucBP_2',  'MucBP','Nramp', 'PIN_3',
#              'Relaxase', 'RelB', 'Rib', 'YkuD')
keykeyGene = c('Gram_pos_anchor', 'Muc_B2', 'MucBP_2',  'MucBP')
statDFsel = statDF3_all[statDF3_all$gene_name %in% keykeyGene & statDF3_all$sampleNum >10, ]
length(unique(statDFsel$gene))
geneDF1 = geneDF[geneDF$gene %in% unique(statDFsel$gene), ]
nrow(geneDF1)
length(unique(statDFsel$gene))
# # 2
# keyGenes = c('DDE_Tnp_1', 'DDE_Tnp_IS66', 'Gram_pos_anchor',
#              'Methylase_S', 'MFS_1', 'Muc_B2', 'MucBP_2',  'MucBP','Nramp', 'PIN_3',
#              'Relaxase', 'RelB', 'Rib', 'YkuD')
# statDFsel = statDF3_all[statDF3_all$gene_name %in% keyGenes, ]
# length(unique(statDFsel$gene))
# geneDF1 = geneDF[geneDF$gene %in% unique(statDFsel$gene), ]
# nrow(geneDF1)
# length(unique(statDFsel$gene))

geneDF1 = geneDF[geneDF$gene %in% unique(statDFsel$gene), ]
nrow(geneDF1)
length(unique(statDFsel$gene))
if(T){
  # ys = c('nucl_diversity', 'dNdS', 'pNpS', 'dNdS_common')
  ys = c('nucl_diversity', 'dNdS', 'pNpS')
  ysNorm = list(nucl_diversity = 'Nucleotide diversity', dNdS = 'dN/dS', pNpS = 'pN/pS', dNdS_common = 'dN/dS (common)')
  xs = c("Ethnicity", "Trimester", "Housing", "Employment", "Term", "Age", "BMI","Marriage","FOB", "Depression")
  comp_DF = data.frame()
  uniRes = data.frame()
  # ge = unique(geneDF$gene)[1]
  # x = xs[1]
  # uniGene = unique(geneDF1$gene)
  uniGene = unique(geneDF1$gene)
  length(uniGene)
  for (gei in 1:length(uniGene)) {
    if(gei %% 100 == 0){
      cat(gei, '/', length(uniGene))
    }
    ge = uniGene[gei]
    tmpDF = geneDF1[geneDF1$gene == ge, ]
    # tmpDF = geneDF1[geneDF1$gene_name == ge, ]
    for(y in ys){
      # cat('  y =', y, '\n')
      for (x in xs) {
        # cat('    x =', x, '\n')
        tmpDF2 = data.frame(
          'y' = tmpDF[, y],
          'x' = tmpDF[, x],
          trimester = tmpDF$Trimester,
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
          lmeTry = try({lmer(as.formula(paste0('y ~ x + (1|pt_ID)')),
                             data = tmpDF2, REML = F)}, silent = F)
          if(class(lmeTry) != 'try-error'){
            lmeUni = lmeTry
            lmeUni.R2 = r.squaredGLMM(lmeUni)[1]
            lmeUni.aov = anova(lmeUni)
            sv = c(ge, y, x, (lmUni.aov['x', c('F value', 'Pr(>F)')]) %>% unlist(), lmUni.summary$r.squared %>% as.numeric(),
                   formula_char, lmeUni.aov$`F value`[1], lmeUni.aov$`Pr(>F)`[1], lmeUni.R2)
          }else{
            sv = c(ge, y, x, (lmUni.aov['x', c('F value', 'Pr(>F)')]) %>% unlist(), lmUni.summary$r.squared %>% as.numeric(),
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
          if(nrow(tmp_comp_DF) > 0){
            tmp_comp_DF$y = y
            tmp_comp_DF$x = x
            tmp_comp_DF$Gene = ge
            tmp_comp_DF$size1 = size$Freq[match(tmp_comp_DF$group1, size$Var1)]
            tmp_comp_DF$size2 = size$Freq[match(tmp_comp_DF$group2, size$Var1)]
            
            tmp_comp_DF$mean1 = mm$y[match(tmp_comp_DF$group1, mm$x)]
            tmp_comp_DF$mean2 = mm$y[match(tmp_comp_DF$group2, mm$x)]
            
            tmp_comp_DF$log2_FC = log2(tmp_comp_DF$mean1/tmp_comp_DF$mean2)
            
            tmp_comp_DF$med1 = mmed$y[match(tmp_comp_DF$group1, mmed$x)]
            tmp_comp_DF$med2 = mmed$y[match(tmp_comp_DF$group2, mmed$x)]
            
            tmp_comp_DF = tmp_comp_DF[, c('Gene', 'y', 'x', 'group1', 'group2', 'size1', 'size2', 'mean1', 'mean2', 'log2_FC', 'med1', 'med2',
                                          'p', 'p.adj', 'p.format', 'p.signif')]
            comp_DF = rbind(comp_DF, tmp_comp_DF)
          }
        }
      }
    }
  }
  colnames(uniRes) = c('gene', 'y', 'x', 'lm.uni.F_value', 'lm.uni.P_value', 'lm.uni.R2',
                       'lme.uni.Model', 'lme.uni.F_value', 'lme.uni.P_value', 'lme.uni.R2')
  uniRes$lm.uni.F_value %<>% as.numeric()
  uniRes$lm.uni.P_value %<>% as.numeric()
  uniRes$lm.uni.R2 %<>% as.numeric()
  uniRes$lme.uni.F_value %<>% as.numeric()
  uniRes$lme.uni.P_value %<>% as.numeric()
  uniRes$lme.uni.R2 %<>% as.numeric()
  uniRes$gene_name = geneDF1$gene_name[match(uniRes$gene, geneDF1$gene)]
  uniRes$gene_description = geneDF1$gene_description[match(uniRes$gene, geneDF1$gene)]
  uniRes$Species = geneDF1$Species[match(uniRes$gene, geneDF1$gene)]
  uniRes$sampleNum = statDFsel$sampleNum[match(uniRes$gene, statDFsel$gene)]
  comp_DF$gene_name = geneDF1$gene_name[match(comp_DF$Gene, geneDF1$gene)]
  comp_DF$gene_description = geneDF1$gene_description[match(comp_DF$Gene, geneDF1$gene)]
  comp_DF$Species = geneDF1$Species[match(comp_DF$Gene, geneDF1$gene)]
  comp_DF$sampleNum = statDFsel$sampleNum[match(comp_DF$Gene, statDFsel$gene)]
  if(plotFlag){
    openxlsx::write.xlsx(comp_DF, file = './data/gene_microdiversity_pairwise_group_comp_v2.xlsx')  
    openxlsx::write.xlsx(uniRes, file = './data/gene_microdiversity_uni_regression_v2.xlsx')
  }
}

##### lme regression fixed model #####
uniRes = readxl::read_xlsx('./data/gene_microdiversity_uni_regression_v2.xlsx')
# nucl div associations to metadata variables
colnames(geneDF1)
xs = c(#'Insulin', 'B6', 
  'Ethnicity', 'Employment', 'Housing',
  'Term', # 'Abortions', 
  #"Marriage",  
  'BMI','Marriage', # 'BMI_char',  
  'Age', 'Trimester',
  'FOB','Depression')  # interest priority
ys = c('nucl_diversity', 'dNdS', 'pNpS', 'dNdS_common')
x_random = c('pt_ID')

# sp = "Lactobacillus_vaginalis (NCBI)" 
# y = 'nucl_diversity'
lmeRes = data.frame()
length(unique(unique(geneDF1$gene)))
for (ge in unique(geneDF1$gene)) {
  cat('  gene =', ge, '\n')
  for (y in ys) {
    cat('  y =', y, '\n')
    dat = geneDF1[geneDF1$gene == ge, c(xs, x_random, y)]
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
        vv_aov_df = data.frame(gene = ge,
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
table(lmeRes$gene)
# merge
res = merge(uniRes, lmeRes, by = c('gene', 'y', 'x'), all = T)
openxlsx::write.xlsx(res, './data/gene_microdiveristy_lme_regression_fixed_v2.xlsx')
res = res[res$sampleNum > 10 & res$y != 'dNdS_common',]
table(res$lme.uni.P_value < 0.05)
table(res$lme.multi.P_value < 0.05)
keyGenes = c('Muc_B2', 'MucBP', 'MucBP_2', 'Gram_pos_anchor')
table(res$gene_name %in% keyGenes)
res_unisig = res[which(res$lme.uni.P_value < 0.05), ]
nrow(res_unisig)
table(res_unisig$gene_name %in% keyGenes)
res_mulisig = res[which(res$lme.multi.P_value < 0.05), ]
nrow(res_mulisig)
table(res_mulisig$gene_name %in% keyGenes)
# res_mulisig = res[which(res$lme.multi.P_value < 0.05 & res$lme.uni.P_value < 0.05),]
# nrow(res_mulisig)
res_mulisig = arrange(res_mulisig, gene, y, x, lme.multi.P_value)
openxlsx::write.xlsx(res_mulisig, './data/gene_microdiveristy_lme_regression_fixed_sig_v2.xlsx')

##### plot boxplot for candidate aka significant species and variable #####
nrow(res_mulisig)
res_mulisig2 = res_mulisig
cand_sp_var = res_mulisig2[, c('gene', 'y', 'x')]

res_unisig2 = res_unisig[res_unisig$gene_name %in% keyGenes, ]
cand_sp_var = res_unisig2[, c('gene', 'y', 'x')]

nrow(cand_sp_var)
plotFlag = T
ysNorm = list(nucl_diversity = 'Nucleotide diversity', dNdS = 'dN/dS', pNpS = 'pN/pS', SNV_per_kbp = 'SNVs/kbp', D_rev = "1 - D'")

if(plotFlag){
  pdf('./plot/gene_microdiveristy_uni_sig_v1.pdf', width = 20, height = 5)
}
for (pr in 1:nrow(cand_sp_var)) {
  sp = cand_sp_var[pr, ][1] %>% unlist
  y = cand_sp_var[pr, ][2] %>% unlist
  x = cand_sp_var[pr, ][3] %>% unlist
  tmpDF = geneDF1[geneDF1$gene == sp, ]
  tmpDF2 = data.frame(
    'fq_ID' = tmpDF$sample,
    'pt_ID' = tmpDF$pt_ID,
    'pt_ID.u' = tmpDF$pt_ID.u,
    'gene' = tmpDF$gene,
    'gene_name' = tmpDF$gene_name,
    'species' = tmpDF$Species,
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
        mytheme3 + theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = 1) +
        stat_compare_means(comparisons=CN, method = "wilcox.test", p.adjust.method = "BH", aes(label=..p.adj..)) +
        labs(title = paste0(x, '|', y, '|', tmpDF$gene_name[1], '|', sp, '|', tmpDF$Species[1]), 
             color = x,
             y = ysNorm[[y]], x = x) + 
        scale_color_brewer(palette = 'Dark2')
      print(p)
    }else{
      p2 = ggplot(na.omit(tmpDF2), aes(x = x, y = y)) +
        geom_smooth(method = 'lm', color = '#969696', size = 0.5, fill = "#D5D5D5") +
        geom_point(aes(color = 'Blue', fill = 'Blue'), shape = 21, size = 2, alpha = 0.5) +
        labs(x = x, y = ysNorm[[y]], title = paste0(x, '|', y, '|', tmpDF$gene_name[1], '|', sp, '|', tmpDF$Species[1])) +
        stat_poly_eq(formula = y ~ x,
                     aes(group=1, label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), 
                     parse = TRUE) +
        mytheme + theme(strip.text = element_text(size = 12), legend.position = 'none', aspect.ratio = 1)
      print(p2)
    }
  }
}
dev.off()
##### again,plot boxplot for candidate aka significant species and variable #####
focus = c(
  'Marriage|dNdS|Muc_B2|NZ_CP018809.1_5|Lactobacillus jensenii',
  'Ethnicity|nucl_diversity|Muc_B2|NZ_CP018809.1_5|Lactobacillus jensenii',
  'Marriage|dNdS|Gram_pos_anchor|NZ_CP039266.1_48|Lactobacillus crispatus', #*
  'Employment|pNpS|MucBP|NZ_CP045664.1_1224|Lactobacillus iners',
  'Ethnicity|nucl_diversity|MucBP|NZ_CP045664.1_439|Lactobacillus iners', #*
  'Ethnicity|nucl_diversity|MucBP|NZ_CP045664.1_6|Lactobacillus iners',
  'Marriage|dNdS|MucBP|NZ_CP045664.1_6|Lactobacillus iners',
  'Marriage|nucl_diversity|Muc_B2|NZ_CP104399.1_1163|Lactobacillus vaginalis',
  'Ethnicity|nucl_diversity|MucBP|NZ_CP104399.1_1335|Lactobacillus vaginalis',
  'Ethnicity|nucl_diversity|Muc_B2|NZ_CP104399.1_26|Lactobacillus vaginalis',
  'Ethnicity|dNdS|Muc_B2|NZ_GG698802.1_114|Lactobacillus coleohominis',
  'Ethnicity|nucl_diversity|Muc_B2|NZ_GG698802.1_114|Lactobacillus coleohominis',
  'Housing|nucl_diversity|Muc_B2|NZ_GG698802.1_114|Lactobacillus coleohominis',
  'Term|nucl_diversity|Muc_B2|NZ_GG698802.1_114|Lactobacillus coleohominis',
  'Term|nucl_diversity|Muc_B2|NZ_GG698802.1_408|Lactobacillus coleohominis',
  'Ethnicity|nucl_diversity|Muc_B2|NZ_GG698802.1_408|Lactobacillus coleohominis',
  'FOB|pNpS|MucBP|NZ_GG698803.1_366|Lactobacillus coleohominis',
  'Ethnicity|nucl_diversity|Muc_B2|NZ_GG698804.1_112|Lactobacillus coleohominis',
  'Housing|nucl_diversity|Muc_B2|NZ_GG698804.1_112|Lactobacillus coleohominis',
  'Term|nucl_diversity|Muc_B2|NZ_GG698804.1_112|Lactobacillus coleohominis',
  'Trimester|pNpS|Muc_B2|NZ_GG698804.1_112|Lactobacillus coleohominis',
  'Housing|pNpS|Muc_B2|NZ_GG698805.1_117|Lactobacillus coleohominis', 
  'Trimester|pNpS|Muc_B2|NZ_GG698805.1_117|Lactobacillus coleohominis', 
'Ethnicity|pNpS|Muc_B2|NZ_GG698805.1_135|Lactobacillus coleohominis',
'Employment|nucl_diversity|Muc_B2|NZ_GG698806.1_18|Lactobacillus coleohominis',
'Ethnicity|nucl_diversity|Muc_B2|NZ_GG698806.1_18|Lactobacillus coleohominis',
'Housing|nucl_diversity|Muc_B2|NZ_GG698806.1_18|Lactobacillus coleohominis',
'FOB|dNdS|MucBP_2|NZ_WBOA01000001.1_609|Lactobacillus gasseri',
'Marriage|dNdS|MucBP_2|NZ_WBOA01000001.1_609|Lactobacillus gasseri',
'Ethnicity|nucl_diversity|MucBP_2|NZ_WBOA01000001.1_609|Lactobacillus gasseri',
'Ethnicity|nucl_diversity|Muc_B2|NZ_WBOA01000001.1_610|Lactobacillus gasseri',
'Trimester|pNpS|Muc_B2|NZ_WBOA01000001.1_610|Lactobacillus gasseri',
'Employment|nucl_diversity|Muc_B2|NZ_WBOA01000002.1_1|Lactobacillus gasseri',
'Employment|nucl_diversity|Muc_B2|NZ_WBOA01000002.1_363|Lactobacillus gasseri',
'Term|pNpS|Muc_B2|NZ_WBOA01000002.1_364|Lactobacillus gasseri',
'Ethnicity|nucl_diversity|Muc_B2|NZ_WBOA01000002.1_365|Lactobacillus gasseri',
'Depression|nucl_diversity|Muc_B2|NZ_WBOA01000004.1_17|Lactobacillus gasseri',
'FOB|pNpS|Muc_B2|NZ_WBOA01000004.1_17|Lactobacillus gasseri',
'Marriage|pNpS|Muc_B2|NZ_WBOA01000004.1_17|Lactobacillus gasseri',
'Depression|nucl_diversity|Muc_B2|NZ_WBOA01000004.1_9|Lactobacillus gasseri',
'Employment|nucl_diversity|Muc_B2|NZ_WBOA01000004.1_9|Lactobacillus gasseri',
'Ethnicity|nucl_diversity|Muc_B2|NZ_WBOA01000004.1_9|Lactobacillus gasseri'
)
focusDF = str_split_fixed(focus, '\\|', 5)
colnames(focusDF) = c('x', 'y', 'gene_name', 'gene', 'species')
cand_sp_var = focusDF[, c('gene', 'y', 'x')]

nrow(cand_sp_var)
plotFlag = T
ysNorm = list(nucl_diversity = 'Nucleotide diversity', dNdS = 'dN/dS', pNpS = 'pN/pS', SNV_per_kbp = 'SNVs/kbp', D_rev = "1 - D'")
if(plotFlag){
  pdf('./plot/gene_microdiveristy_uni_sig_v1_again.pdf', width = 20, height = 5)
}
for (pr in 1:nrow(cand_sp_var)) {
  sp = cand_sp_var[pr, ][1] %>% unlist
  y = cand_sp_var[pr, ][2] %>% unlist
  x = cand_sp_var[pr, ][3] %>% unlist
  tmpDF = geneDF1[geneDF1$gene == sp, ]
  tmpDF2 = data.frame(
    'fq_ID' = tmpDF$sample,
    'pt_ID' = tmpDF$pt_ID,
    'pt_ID.u' = tmpDF$pt_ID.u,
    'gene' = tmpDF$gene,
    'gene_name' = tmpDF$gene_name,
    'species' = tmpDF$Species,
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
        mytheme3 + theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = 1) +
        stat_compare_means(comparisons=CN, method = "wilcox.test", p.adjust.method = "BH", aes(label=..p.adj..)) +
        labs(title = paste0(x, '|', y, '|', tmpDF$gene_name[1], '|', sp, '|', tmpDF$Species[1]), 
             color = x,
             y = ysNorm[[y]], x = x) + 
        scale_color_brewer(palette = 'Dark2')
      print(p)
    }else{
      p2 = ggplot(na.omit(tmpDF2), aes(x = x, y = y)) +
        geom_smooth(method = 'lm', color = '#969696', size = 0.5, fill = "#D5D5D5") +
        geom_point(aes(color = 'Blue', fill = 'Blue'), shape = 21, size = 2, alpha = 0.5) +
        labs(x = x, y = ysNorm[[y]], title = paste0(x, '|', y, '|', tmpDF$gene_name[1], '|', sp, '|', tmpDF$Species[1])) +
        stat_poly_eq(formula = y ~ x,
                     aes(group=1, label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), 
                     parse = TRUE) +
        mytheme + theme(strip.text = element_text(size = 12), legend.position = 'none', aspect.ratio = 1)
      print(p2)
    }
  }
}
dev.off()

##### a little carefully!! ####

relAbuPlot = function(data0, xval, yval, zval = NULL, yLab = 'Value',  titleLab = 'title', jitterSize = 1.8,freeMode = 'fixed',spaceMode = 'fixed',
                      y.position = NULL,step.increase = 0.12, expandMult = NULL,leName = NULL, addJitter = T,aspectRatio = NULL,zvalX = NULL,
                      x.angle = 0, CLD = F, hide.ns = F, legend.position = 'right', if_log10 = F, show_Padj = F){
  data0$x = data0[, xval]
  if(if_log10){
    data0$y = log10(data0[, yval])
  }else{
    data0$y = data0[, yval]
  }
  
  if(endsWith(xval, '2')){xLab = str_remove(xval, '2$')}else{xLab = xval}
  if(!is.null(zval)){data0$z = data0[, zval]}
  if(!is.null(zvalX)){data0$zX = data0[, zvalX]}
  if(x.angle == 0){x.hjust = 0.5}else{x.hjust = 1}
  if(is.null(leName)){
    leName = xLab
  }
  # base ggplot
  if(addJitter){
    pbO = ggboxplot(data0, x = 'x', y = 'y',
                    color = 'x', title = titleLab,
                    add = "jitter", 
                    add.params = list(shape = 21, size = jitterSize, color = 'x', fill = 'x', alpha = 0.7),
                    bxp.errorbar = T, notch = F) +
      scale_color_manual(values = get(paste0('color', xval)), name = leName) +
      scale_fill_manual(values = get(paste0('color', xval)), name = leName) +
      labs(y = yLab,  x = xLab) +
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
        if(!is.null(zvalX)){
          stat.test <- data0 %>%
            group_by(z, zX) %>%
            wilcox_test(y ~ x) %>%
            adjust_pvalue(method = "BH") %>%
            add_significance("p.adj") %>%
            add_xy_position(x = xval,dodge = 0, step.increase = step.increase)
        }else{
          stat.test <- data0 %>%
            group_by(z) %>%
            wilcox_test(y ~ x) %>%
            adjust_pvalue(method = "BH") %>%
            add_significance("p.adj") %>%
            add_xy_position(x = xval,dodge = 0, step.increase = step.increase)
        }
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
    if(!is.null(zvalX)){
      pbO = pbO + facet_grid(zX~z, scales = freeMode, space = spaceMode)
    }else{
      pbO = pbO + facet_grid(.~z, scales = freeMode, space = spaceMode)
    }

  }
  if(!is.null(aspectRatio)){
    pbO = pbO + theme(aspect.ratio  = aspectRatio)
  }
  print(pbO)
}
cairo_pdf('./plot/gene_super_good_results.pdf', width = 6.05,height = 3.55, onefile = T)
if(T){
  for(ge in keyGenes){
    tmpDF = geneDF1[geneDF1$gene_name == ge, ]
    relAbuPlot(data0 = tmpDF, xval = 'Ethnicity', yLab = expression(pi),titleLab = ge,
               yval = 'nucl_diversity', if_log10 = F, hide.ns = T, 
               # y.position = seq(0.29, 0.29 + 0.04*5, 0.04), 
               aspectRatio = 1/1.2,
               expandMult = expansion(mult = c(0.03, 0.1)))
    try({    relAbuPlot(data0 = tmpDF, xval = 'Ethnicity', yLab = "dN/dS",titleLab = ge,
                        yval = 'dNdS', if_log10 = F, hide.ns = T, 
                        # y.position = seq(0.29, 0.29 + 0.04*5, 0.04), 
                        aspectRatio = 1/1.2,
                        expandMult = expansion(mult = c(0.03, 0.1)))})

    try({    relAbuPlot(data0 = tmpDF, xval = 'Ethnicity', yLab = "pN/pS",titleLab = ge,
               yval = 'pNpS', if_log10 = F, hide.ns = T, 
               # y.position = seq(0.29, 0.29 + 0.04*5, 0.04), 
               aspectRatio = 1/1.2,
               expandMult = expansion(mult = c(0.03, 0.1)))})
  }


}
dev.off()
cairo_pdf('./plot/gene_super_good_results_3.pdf', width = 6.05,height = 3.55, onefile = T)
geneDF1 = geneDF
if(T){
  sp = 'Lactobacillus iners'
  tmpDF = geneDF1[geneDF1$Species == sp & geneDF1$gene == 'NZ_CP045664.1_439', ]
  relAbuPlot(data0 = tmpDF, xval = 'Ethnicity', yLab = expression(pi), 
             titleLab = 'Lactobacillus iners | MucBP | NZ_CP045664.1_439',
             yval = 'nucl_diversity',
             if_log10 = F, hide.ns = T, aspectRatio = 1/1.2,
             y.position = c(0.025, 0.028,0.031),
             expandMult = expansion(mult = c(0.03, 0.1)))
  sp = 'Lactobacillus iners'
  tmpDF = geneDF1[geneDF1$Species == sp & geneDF1$gene_name == 'MucBP', ]
  relAbuPlot(data0 = tmpDF, xval = 'Ethnicity', yLab = expression(pi), 
             titleLab = 'Lactobacillus iners | MucBP',
             yval = 'nucl_diversity',
             if_log10 = F, hide.ns = T, aspectRatio = 1/1.2,
             y.position =c(0.025, 0.028,0.031, 0.034)*1.2 + 0.005,
             expandMult = expansion(mult = c(0.03, 0.1)))
  sp = 'Lactobacillus iners'
  tmpDF = geneDF1[geneDF1$Species == sp & geneDF1$gene == 'NZ_CP045664.1_439', ]
  relAbuPlot(data0 = tmpDF, xval = 'Ethnicity', yLab = 'dN/dS',
             titleLab = 'Lactobacillus iners | MucBP | NZ_CP045664.1_439',
             yval = 'dNdS',
             if_log10 = F, hide.ns = T, aspectRatio = 1/1.2,
             # y.position = seq(0.4, 0.4 + 0.05* 2, 0.05)+0.02, 
             expandMult = expansion(mult = c(0.03, 0.1)))
  sp = 'Lactobacillus iners'
  tmpDF = geneDF1[geneDF1$Species == sp & geneDF1$gene_name == 'MucBP', ]
  relAbuPlot(data0 = tmpDF, xval = 'Ethnicity', yLab = 'dN/dS',
             titleLab = 'Lactobacillus iners | MucBP',
             yval =  'dNdS',
             if_log10 = F, hide.ns = T, aspectRatio = 1/1.2,
             y.position = seq(0.4, 0.4 + 0.05* 2, 0.05)+0.02, 
             expandMult = expansion(mult = c(0.03, 0.1)))
  sp = 'Lactobacillus iners'
  tmpDF = geneDF1[geneDF1$Species == sp & geneDF1$gene == 'NZ_CP045664.1_439', ]
  relAbuPlot(data0 = tmpDF, xval = 'Ethnicity', yLab = 'pN/pS',
             titleLab = 'Lactobacillus iners | MucBP | NZ_CP045664.1_439',
             yval = 'pNpS',
             if_log10 = F, hide.ns = T, aspectRatio = 1/1.2,
             y.position = seq(2, 2 + 0.4*3, 0.4),
             expandMult = expansion(mult = c(0.03, 0.1)))
  sp = 'Lactobacillus iners'
  tmpDF = geneDF1[geneDF1$Species == sp & geneDF1$gene_name == 'MucBP', ]
  relAbuPlot(data0 = tmpDF, xval = 'Ethnicity', yLab = 'pN/pS',
             titleLab = 'Lactobacillus iners | MucBP',
             yval =  'pNpS',
             if_log10 = F, hide.ns = T, aspectRatio = 1/1.2,
             y.position = seq(2, 2 + 0.4*3, 0.4),
             expandMult = expansion(mult = c(0.03, 0.1)))
  sp = 'Lactobacillus coleohominis'
  tmpDF = geneDF1[geneDF1$Species == sp & geneDF1$gene_name == 'Muc_B2', ]
  tmpDF$Term2 = tmpDF$Term
  relAbuPlot(data0 = tmpDF, xval = 'Term2', 
             yLab = expression(pi),
             titleLab = 'Lactobacillus coleohominis | MucBP',
             yval = 'nucl_diversity',
             if_log10 = F, hide.ns = T, aspectRatio = 1/1.2,
             # y.position = seq(2, 2 + 0.4*3, 0.4),
             expandMult = expansion(mult = c(0.03, 0.1)))
  sp = 'Lactobacillus coleohominis'
  tmpDF = geneDF1[geneDF1$Species == sp & geneDF1$gene_name == 'Muc_B2', ]
  tmpDF$Term2 = tmpDF$Term
  relAbuPlot(data0 = tmpDF, xval = 'Housing', 
             yLab = expression(pi),
             titleLab = 'Lactobacillus coleohominis | MucBP',
             yval = 'nucl_diversity',
             if_log10 = F, hide.ns = T, aspectRatio = 1/1.2,x.angle = 30,
             y.position = c(0.011,0.012),
             expandMult = expansion(mult = c(0.03, 0.1)))
  sp = 'Lactobacillus gasseri'
  tmpDF = geneDF1[geneDF1$Species == sp & geneDF1$gene_name == 'Muc_B2', ]
  # tmpDF$Term2 = tmpDF$Term
  relAbuPlot(data0 = tmpDF, xval = 'Employment', 
             yLab = expression(pi),
             titleLab = 'Lactobacillus coleohominis | MucBP',
             yval = 'nucl_diversity',
             if_log10 = F, hide.ns = T, aspectRatio = 1/1.2,x.angle = 30,
             y.position = c(0.011),
             expandMult = expansion(mult = c(0.03, 0.1)))
}
dev.off()
#### single variable comparison lm & lme (for gene name) #####
keykeyGene = c('Gram_pos_anchor', 'Muc_B2', 'MucBP_2',  'MucBP')
statDFsel = statDF3_all[statDF3_all$gene_name %in% keykeyGene & statDF3_all$sampleNum >10, ]
length(unique(statDFsel$gene))
geneDF1 = geneDF[geneDF$gene %in% unique(statDFsel$gene), ]
nrow(geneDF1)
geneDF1$spGene = paste0(geneDF1$Species, ' | ', geneDF1$gene_name)
length(unique(geneDF1$spGene))
if(T){
  ys = c('nucl_diversity', 'dNdS', 'pNpS')
  ysNorm = list(nucl_diversity = 'Nucleotide diversity', dNdS = 'dN/dS', pNpS = 'pN/pS', dNdS_common = 'dN/dS (common)')
  xs = c("Ethnicity", "Trimester", "Housing", "Employment", "Term", "Age", "BMI","Marriage","FOB", "Depression")
  comp_DF = data.frame()
  uniRes = data.frame()
  uniGene = unique(geneDF1$spGene)
  length(uniGene)
  for (gei in 1:length(uniGene)) {
    if(gei %% 100 == 0){
      cat(gei, '/', length(uniGene))
    }
    ge = uniGene[gei]
    tmpDF = geneDF1[geneDF1$spGene == ge, ]
    # tmpDF = geneDF1[geneDF1$gene_name == ge, ]
    for(y in ys){
      # cat('  y =', y, '\n')
      for (x in xs) {
        # cat('    x =', x, '\n')
        tmpDF2 = data.frame(
          'y' = tmpDF[, y],
          'x' = tmpDF[, x],
          trimester = tmpDF$Trimester,
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
          lmeTry = try({lmer(as.formula(paste0('y ~ x + (1|pt_ID)')),
                             data = tmpDF2, REML = F)}, silent = F)
          if(class(lmeTry) != 'try-error'){
            lmeUni = lmeTry
            lmeUni.R2 = r.squaredGLMM(lmeUni)[1]
            lmeUni.aov = anova(lmeUni)
            sv = c(ge, y, x, (lmUni.aov['x', c('F value', 'Pr(>F)')]) %>% unlist(), lmUni.summary$r.squared %>% as.numeric(),
                   formula_char, lmeUni.aov$`F value`[1], lmeUni.aov$`Pr(>F)`[1], lmeUni.R2)
          }else{
            sv = c(ge, y, x, (lmUni.aov['x', c('F value', 'Pr(>F)')]) %>% unlist(), lmUni.summary$r.squared %>% as.numeric(),
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
          if(nrow(tmp_comp_DF) > 0){
            tmp_comp_DF$y = y
            tmp_comp_DF$x = x
            tmp_comp_DF$Gene = ge
            tmp_comp_DF$size1 = size$Freq[match(tmp_comp_DF$group1, size$Var1)]
            tmp_comp_DF$size2 = size$Freq[match(tmp_comp_DF$group2, size$Var1)]
            
            tmp_comp_DF$mean1 = mm$y[match(tmp_comp_DF$group1, mm$x)]
            tmp_comp_DF$mean2 = mm$y[match(tmp_comp_DF$group2, mm$x)]
            
            tmp_comp_DF$log2_FC = log2(tmp_comp_DF$mean1/tmp_comp_DF$mean2)
            
            tmp_comp_DF$med1 = mmed$y[match(tmp_comp_DF$group1, mmed$x)]
            tmp_comp_DF$med2 = mmed$y[match(tmp_comp_DF$group2, mmed$x)]
            
            tmp_comp_DF = tmp_comp_DF[, c('Gene', 'y', 'x', 'group1', 'group2', 'size1', 'size2', 'mean1', 'mean2', 'log2_FC', 'med1', 'med2',
                                          'p', 'p.adj', 'p.format', 'p.signif')]
            comp_DF = rbind(comp_DF, tmp_comp_DF)
          }
        }
      }
    }
  }
  colnames(uniRes) = c('gene', 'y', 'x', 'lm.uni.F_value', 'lm.uni.P_value', 'lm.uni.R2',
                       'lme.uni.Model', 'lme.uni.F_value', 'lme.uni.P_value', 'lme.uni.R2')
  uniRes$lm.uni.F_value %<>% as.numeric()
  uniRes$lm.uni.P_value %<>% as.numeric()
  uniRes$lm.uni.R2 %<>% as.numeric()
  uniRes$lme.uni.F_value %<>% as.numeric()
  uniRes$lme.uni.P_value %<>% as.numeric()
  uniRes$lme.uni.R2 %<>% as.numeric()
  if(plotFlag){
    openxlsx::write.xlsx(comp_DF, file = './data/gene_microdiversity_pairwise_group_comp_v2_geneName.xlsx')  
    openxlsx::write.xlsx(uniRes, file = './data/gene_microdiversity_uni_regression_v2_geneName.xlsx')
  }
}

##### lme regression fixed model (for gene name) #####
uniRes = readxl::read_xlsx('./data/gene_microdiversity_uni_regression_v2_geneName.xlsx')
# nucl div associations to metadata variables
colnames(geneDF1)
xs = c(#'Insulin', 'B6', 
  'Ethnicity', 'Employment', 'Housing',
  'Term', # 'Abortions', 
  #"Marriage",  
  'BMI','Marriage', # 'BMI_char',  
  'Age', 'Trimester',
  'FOB','Depression')  # interest priority
ys = c('nucl_diversity', 'dNdS', 'pNpS', 'dNdS_common')
x_random = c('pt_ID')
lmeRes = data.frame()
length(unique(geneDF1$spGene))
for (ge in unique(geneDF1$spGene)) {
  cat('  gene =', ge, '\n')
  for (y in ys) {
    cat('  y =', y, '\n')
    dat = geneDF1[geneDF1$spGene == ge, c(xs, x_random, y)]
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
        vv_aov_df = data.frame(gene = ge,
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
table(lmeRes$gene)
# merge
res = merge(uniRes, lmeRes, by = c('gene', 'y', 'x'), all = T)
openxlsx::write.xlsx(res, './data/gene_microdiveristy_lme_regression_fixed_v2_geneName.xlsx')
table(res$lme.uni.P_value < 0.05)
table(res$lme.multi.P_value < 0.05)
keyGenes = c('Muc_B2', 'MucBP', 'MucBP_2', 'Gram_pos_anchor')
table(res$gene_name %in% keyGenes)
res_unisig = res[which(res$lme.uni.P_value < 0.05), ]
nrow(res_unisig)
table(res_unisig$gene_name %in% keyGenes)
res_mulisig = res[which(res$lme.multi.P_value < 0.05), ]
nrow(res_mulisig)
table(res_mulisig$gene_name %in% keyGenes)
# res_mulisig = res[which(res$lme.multi.P_value < 0.05 & res$lme.uni.P_value < 0.05),]
# nrow(res_mulisig)
res_mulisig = arrange(res_mulisig, gene, y, x, lme.multi.P_value)
openxlsx::write.xlsx(res_mulisig, './data/gene_microdiveristy_lme_regression_fixed_sig_v2_geneName.xlsx')


###### boxplot for cand comparsion (for gene Name) ####
res_unisig2 = res_unisig
cand_sp_var = res_unisig2[, c('gene', 'y', 'x')]

nrow(cand_sp_var)
plotFlag = T
ysNorm = list(nucl_diversity = 'Nucleotide diversity', dNdS = 'dN/dS', pNpS = 'pN/pS', SNV_per_kbp = 'SNVs/kbp', D_rev = "1 - D'")

if(plotFlag){
  pdf('./plot/geneName_microdiveristy_uni_sig_v1.pdf', width = 20, height = 5)
}
for (pr in 1:nrow(cand_sp_var)) {
  sp = cand_sp_var[pr, ][1] %>% unlist
  y = cand_sp_var[pr, ][2] %>% unlist
  x = cand_sp_var[pr, ][3] %>% unlist
  tmpDF = geneDF1[geneDF1$spGene == sp, ]
  tmpDF2 = data.frame(
    'fq_ID' = tmpDF$sample,
    'pt_ID' = tmpDF$pt_ID,
    'pt_ID.u' = tmpDF$pt_ID.u,
    'gene' = tmpDF$gene,
    'gene_name' = tmpDF$gene_name,
    'species' = tmpDF$Species,
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
        mytheme3 + theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = 1) +
        stat_compare_means(comparisons=CN, method = "wilcox.test", p.adjust.method = "BH", aes(label=..p.adj..)) +
        labs(title = paste0(x, '|', y, '|', tmpDF$gene_name[1], '|', sp, '|', tmpDF$Species[1]), 
             color = x,
             y = ysNorm[[y]], x = x) + 
        scale_color_brewer(palette = 'Dark2')
      print(p)
    }else{
      p2 = ggplot(na.omit(tmpDF2), aes(x = x, y = y)) +
        geom_smooth(method = 'lm', color = '#969696', size = 0.5, fill = "#D5D5D5") +
        geom_point(aes(color = 'Blue', fill = 'Blue'), shape = 21, size = 2, alpha = 0.5) +
        labs(x = x, y = ysNorm[[y]], title = paste0(x, '|', y, '|', tmpDF$gene_name[1], '|', sp, '|', tmpDF$Species[1])) +
        stat_poly_eq(formula = y ~ x,
                     aes(group=1, label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), 
                     parse = TRUE) +
        mytheme + theme(strip.text = element_text(size = 12), legend.position = 'none', aspect.ratio = 1)
      print(p2)
    }
  }
}
dev.off()

#### boxplot carefully ####
focus = c(
  'Term|nucl_diversity|Muc_B2|Lactobacillus coleohominis|Muc_B2|Lactobacillus coleohominis',
  'Trimester|pNpS|Muc_B2|Lactobacillus coleohominis|Muc_B2|Lactobacillus coleohominis',
  'Marriage|dNdS|Gram_pos_anchor|Lactobacillus crispatus|Gram_pos_anchor|Lactobacillus crispatus',
  'Trimester|nucl_diversity|Muc_B2|Lactobacillus gasseri|Muc_B2|Lactobacillus gasseri',
  'Ethnicity|pNpS|Muc_B2|Lactobacillus gasseri | Muc_B2|Lactobacillus gasseri',
  'Ethnicity|nucl_diversity|MucBP|Lactobacillus iners |MucBP|Lactobacillus iners'
)
cairo_pdf('./plot/gene_super_good_results_4.pdf', width = 6.05,height = 3.31, onefile = T)
geneDF1 = geneDF
if(T){
  sp = 'Lactobacillus coleohominis'
  tmpDF = geneDF1[geneDF1$Species == sp & geneDF1$gene_name == 'Muc_B2', ]
  tmpDF$Term2 = tmpDF$Term
  relAbuPlot(data0 = tmpDF, xval = 'Term2', 
             yLab = expression(pi),
             titleLab = 'Lactobacillus coleohominis | Muc_B2',
             yval = 'nucl_diversity',
             if_log10 = F, hide.ns = T, aspectRatio = 1/1.2,
             # y.position = seq(2, 2 + 0.4*3, 0.4),
             expandMult = expansion(mult = c(0.03, 0.1)))
  sp = 'Lactobacillus coleohominis'
  tmpDF = geneDF1[geneDF1$Species == sp & geneDF1$gene_name == 'Muc_B2', ]
  relAbuPlot(data0 = tmpDF[!is.na(tmpDF$pNpS), ], xval = 'Trimester', 
             yLab = 'pN/pS',
             titleLab = 'Lactobacillus coleohominis | Muc_B2',
             yval = 'pNpS',CLD = T,
             if_log10 = F, hide.ns = T, aspectRatio = 1/1.2,
             # y.position = c(0.011,0.012),
             expandMult = expansion(mult = c(0.03, 0.5)))
  
  sp = 'Lactobacillus crispatus'
  tmpDF = geneDF1[geneDF1$Species == sp & geneDF1$gene_name == 'Gram_pos_anchor', ]
  relAbuPlot(data0 = tmpDF[!is.na(tmpDF$dNdS) & !is.na(tmpDF$Marriage), ], xval = 'Marriage', 
             yLab = 'dN/dS',
             titleLab = 'Lactobacillus crispatus | Gram_pos_anchor',
             yval = 'dNdS',CLD = F,
             if_log10 = F, hide.ns = T, aspectRatio = 1/1.2,
             y.position = 2.7,
             expandMult = expansion(mult = c(0.03, 0.1)))

  sp = 'Lactobacillus gasseri'
  tmpDF = geneDF1[geneDF1$Species == sp & geneDF1$gene_name == 'Muc_B2', ]
  # tmpDF$Term2 = tmpDF$Term
  relAbuPlot(data0 = tmpDF, xval = 'Trimester', 
             yLab = expression(pi),
             titleLab = 'Lactobacillus gasseri | Muc_B2',
             yval = 'nucl_diversity',
             if_log10 = F, hide.ns = T, aspectRatio = 1/1.2,
             y.position =seq(0.01, 0.01 + 0.002*5, 0.002)-0.002,
             expandMult = expansion(mult = c(0.03, 0.1)))
  sp = 'Lactobacillus gasseri'
  tmpDF = geneDF1[geneDF1$Species == sp & geneDF1$gene_name == 'Muc_B2', ]
  relAbuPlot(data0 = tmpDF, xval = 'Ethnicity', 
             yLab = 'pN/pS',
             titleLab = 'Lactobacillus gasseri | Muc_B2',
             yval = 'pNpS',CLD = T,
             if_log10 = F, hide.ns = T, aspectRatio = 1/1.2,
             y.position =seq(0.01, 0.01 + 0.002*5, 0.002)-0.002,
             expandMult = expansion(mult = c(0.03, 0.2)))

  sp = 'Lactobacillus iners'
  tmpDF = geneDF1[geneDF1$Species == sp & geneDF1$gene == 'NZ_CP045664.1_439', ]
  relAbuPlot(data0 = tmpDF, xval = 'Ethnicity', yLab = expression(pi), 
             titleLab = 'Lactobacillus iners | MucBP | NZ_CP045664.1_439',
             yval = 'nucl_diversity',
             if_log10 = F, hide.ns = T, aspectRatio = 1/1.2,
             y.position = c(0.025, 0.028,0.031),
             expandMult = expansion(mult = c(0.03, 0.1)))
  sp = 'Lactobacillus iners'
  tmpDF = geneDF1[geneDF1$Species == sp & geneDF1$gene_name == 'MucBP', ]
  relAbuPlot(data0 = tmpDF, xval = 'Ethnicity', yLab = expression(pi), 
             titleLab = 'Lactobacillus iners | MucBP',
             yval = 'nucl_diversity',
             if_log10 = F, hide.ns = T, aspectRatio = 1/1.2,
             y.position =c(0.025, 0.028,0.031, 0.034)*1.2 + 0.005,
             expandMult = expansion(mult = c(0.03, 0.1)))
  sp = 'Lactobacillus iners'
  tmpDF = geneDF1[geneDF1$Species == sp & geneDF1$gene == 'NZ_CP045664.1_439', ]
  relAbuPlot(data0 = tmpDF, xval = 'Ethnicity', yLab = 'dN/dS',
             titleLab = 'Lactobacillus iners | MucBP | NZ_CP045664.1_439',
             yval = 'dNdS',
             if_log10 = F, hide.ns = T, aspectRatio = 1/1.2,
             # y.position = seq(0.4, 0.4 + 0.05* 2, 0.05)+0.02, 
             expandMult = expansion(mult = c(0.03, 0.1)))
  sp = 'Lactobacillus iners'
  tmpDF = geneDF1[geneDF1$Species == sp & geneDF1$gene_name == 'MucBP', ]
  relAbuPlot(data0 = tmpDF, xval = 'Ethnicity', yLab = 'dN/dS',
             titleLab = 'Lactobacillus iners | MucBP',
             yval =  'dNdS',
             if_log10 = F, hide.ns = T, aspectRatio = 1/1.2,
             y.position = seq(0.4, 0.4 + 0.05* 2, 0.05)+0.02, 
             expandMult = expansion(mult = c(0.03, 0.1)))
  sp = 'Lactobacillus iners'
  tmpDF = geneDF1[geneDF1$Species == sp & geneDF1$gene == 'NZ_CP045664.1_439', ]
  relAbuPlot(data0 = tmpDF, xval = 'Ethnicity', yLab = 'pN/pS',
             titleLab = 'Lactobacillus iners | MucBP | NZ_CP045664.1_439',
             yval = 'pNpS',
             if_log10 = F, hide.ns = T, aspectRatio = 1/1.2,
             y.position = seq(2, 2 + 0.4*3, 0.4),
             expandMult = expansion(mult = c(0.03, 0.1)))
  sp = 'Lactobacillus iners'
  tmpDF = geneDF1[geneDF1$Species == sp & geneDF1$gene_name == 'MucBP', ]
  relAbuPlot(data0 = tmpDF, xval = 'Ethnicity', yLab = 'pN/pS',
             titleLab = 'Lactobacillus iners | MucBP',
             yval =  'pNpS',
             if_log10 = F, hide.ns = T, aspectRatio = 1/1.2,
             y.position = seq(2, 2 + 0.4*3, 0.4),
             expandMult = expansion(mult = c(0.03, 0.1)))
}
dev.off()
######### comparsion evolutionary metrics for each gene #####
load('./data/RData/geneDF.RData')
keyGenes = c('Muc_B2', 'MucBP', 'MucBP_2', 'Gram_pos_anchor')
myPal = c('#FB8072', '#FDB462', '#B3DE68', '#80B1D3')
names(myPal) = keyGenes
colorGene = myPal
statKey = statDF3_all[statDF3_all$gene_name %in% keyGenes & statDF3_all$sampleNum > 10,]
nrow(statKey)
myShape = c(22:25)
names(myShape) = keyGenes
# for(ge in keyGenes){
if(T){
  tmpDF = geneDF[geneDF$gene_name %in% keyGenes, ]
  tmpDF$Gene = tmpDF$gene_name
  tmpDF$Species %<>% str_replace_all(., 'Lactobacillus', 'L.')
  
  stat.test = tmpDF %>% 
    group_by(Species) %>%
    wilcox_test(nucl_diversity ~ Gene) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj") %>%
    add_xy_position(x = "Gene", group = 'Species', scales = "free")
  stat.test = tmpDF %>% filter(!is.na(pNpS)) %>%
    group_by(Species) %>%
    wilcox_test(pNpS ~ Gene) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj") %>%
    add_xy_position(x = "Gene", group = 'Species', scales = "free")
  stat.test = tmpDF %>% filter(!is.na(dNdS)) %>%
    group_by(Species) %>%
    wilcox_test(dNdS ~ Gene) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj") %>%
    add_xy_position(x = "Gene", group = 'Species', scales = "free")
  pbO = ggboxplot(tmpDF, x = 'Gene', y = 'nucl_diversity',
                  color = 'Gene', 
                  # add = "jitter",
                  add.params = list(shape = 21, size = 1.5, color = 'Gene', fill = 'Gene', alpha = 0.7),
                  bxp.errorbar = T, notch = F) +
    scale_color_manual(values = colorGene) +
    scale_fill_manual(values = colorGene) +
    # stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.05, step.increase = 0.01) +
    labs(y = expression(pi),  x = 'Gene') +
    facet_grid(.~Species, scales = 'free', space = 'free') +
    mytheme + theme(strip.text = element_text(size = 12),
                    # axis.text.x = element_text(angle = 30, hjust = 1),
                    axis.text.x = element_blank(), legend.margin = margin(5,5,-5,5),
                     plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 'cm'), legend.position = 'top')
  cairo_pdf('./plot/compare_gene_Pi_for_each_species.pdf', width = 8.01, height = 2.98)
  print(pbO)
  dev.off()
  pbO = ggboxplot(tmpDF, x = 'Gene', y = 'pNpS',
                  color = 'Gene', 
                  # add = "jitter",
                  add.params = list(shape = 21, size = 1.5, color = 'Gene', fill = 'Gene', alpha = 0.7),
                  bxp.errorbar = T, notch = F) +
    scale_color_manual(values = colorGene) +
    scale_fill_manual(values = colorGene) +
    # stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.05, step.increase = 0.01) +
    labs(y = 'pN/pS',  x = 'Gene') +
    facet_grid(.~Species, scales = 'free', space = 'free') +
    mytheme + theme(strip.text = element_text(size = 12),
                    # axis.text.x = element_text(angle = 30, hjust = 1),
                    axis.text.x = element_blank(), legend.margin = margin(5,5,-5,5),
                    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 'cm'), legend.position = 'top')
  cairo_pdf('./plot/compare_gene_pNpS_for_each_species.pdf', width = 8.01, height = 2.98)
  print(pbO)
  dev.off()
  pbO = ggboxplot(tmpDF, x = 'Gene', y = 'dNdS',
                  color = 'Gene', 
                  # add = "jitter",
                  add.params = list(shape = 21, size = 1.5, color = 'Gene', fill = 'Gene', alpha = 0.7),
                  bxp.errorbar = T, notch = F) +
    scale_color_manual(values = colorGene) +
    scale_fill_manual(values = colorGene) +
    # stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.05, step.increase = 0.01) +
    labs(y = 'dN/dS',  x = 'Gene') +
    facet_grid(.~Species, scales = 'free', space = 'free') +
    mytheme + theme(strip.text = element_text(size = 12),
                    # axis.text.x = element_text(angle = 30, hjust = 1),
                    axis.text.x = element_blank(), legend.margin = margin(5,5,-5,5),
                    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 'cm'), legend.position = 'top')
  cairo_pdf('./plot/compare_gene_dNdS_for_each_species.pdf', width = 8.01, height = 2.98)
  print(pbO)
  dev.off()
}
  # pbO = ggboxplot(tmpDF, x = 'Gene', y = 'nucl_diversity',
  #                 color = 'Gene', add = "jitter", 
  #                 add.params = list(shape = 21, size = 1.5, color = 'Gene', fill = 'Gene', alpha = 0.7),
  #                 bxp.errorbar = T, notch = F) +
  #   scale_color_manual(values = colorGene) +
  #   scale_fill_manual(values = colorGene) +
  #   stat_compare_means() +
  #   labs(y = expression(pi),  x = 'Gene') +
  #   facet_grid(Ethnicity~Species, scales = 'free', space = 'free') +
  #   mytheme + theme(strip.text = element_text(size = 12), axis.text.x = element_text(angle = x.angle, hjust = x.hjust),
  #                    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 'cm'), legend.position = legend.position)
  # print(pbO)
######### comparsion evolutionary metrics for each species #####
load('./data/RData/geneDF.RData')
keyGenes = c('Muc_B2', 'MucBP', 'MucBP_2', 'Gram_pos_anchor')
myPal = c('#FB8072', '#FDB462', '#B3DE68', '#80B1D3')
names(myPal) = keyGenes
colorGene = myPal
statKey = statDF3_all[statDF3_all$gene_name %in% keyGenes & statDF3_all$sampleNum > 10,]
nrow(statKey)
myShape = c(22:25)
names(myShape) = keyGenes
# for(ge in keyGenes){
if(T){
  tmpDF = geneDF[geneDF$gene_name %in% keyGenes, ]
  tmpDF$Gene = tmpDF$gene_name
  tmpDF$Species %<>% str_replace_all(., 'Lactobacillus', 'L.')
  names(colorSpecies) %<>% str_replace_all(., 'Lactobacillus', 'L.')
  stat.test = tmpDF %>%
    group_by(Gene) %>%
    wilcox_test(nucl_diversity ~ Species) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj") %>%
    add_xy_position(x = "Species", group = 'Gene', scales = "free")
  print(stat.test, n = Inf)
  stat.test = tmpDF %>% filter(!is.na(tmpDF$pNpS)) %>%
    group_by(Gene) %>%
    wilcox_test(pNpS ~ Species) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj") %>%
    add_xy_position(x = "Species", group = 'Gene', scales = "free")
  print(stat.test, n = Inf)
  stat.test = tmpDF %>% filter(!is.na(tmpDF$dNdS)) %>%
    group_by(Gene) %>%
    wilcox_test(dNdS ~ Species) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj") %>%
    add_xy_position(x = "Species", group = 'Gene', scales = "free")
  print(stat.test, n = Inf)
  pbO = ggboxplot(tmpDF, x = 'Species', y = 'nucl_diversity',
                  color = 'Species', 
                  # add = "jitter",
                  add.params = list(shape = 21, size = 1.5, color = 'Species', fill = 'Species', alpha = 0.7),
                  bxp.errorbar = T, notch = F) +
    scale_color_manual(values = colorSpecies) +
    scale_fill_manual(values = colorSpecies) +
    # stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.05, step.increase = 0.01) +
    labs(y = expression(pi),  x = 'Gene') +
    facet_grid(.~Gene, scales = 'free', space = 'free') +
    mytheme + theme(strip.text = element_text(size = 12),
                    # axis.text.x = element_text(angle = 30, hjust = 1),
                    axis.text.x = element_blank(), legend.margin = margin(5,5,-5,5),
                    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 'cm'), legend.position = 'top')
  cairo_pdf('./plot/compare_gene_Pi_for_each_genes.pdf', width = 6.65, height = 3.48)
  print(pbO)
  dev.off()
  pbO = ggboxplot(tmpDF, x = 'Species', y = 'pNpS',
                  color = 'Species', 
                  # add = "jitter",
                  add.params = list(shape = 21, size = 1.5, color = 'Species', fill = 'Species', alpha = 0.7),
                  bxp.errorbar = T, notch = F) +
    scale_color_manual(values = colorSpecies) +
    scale_fill_manual(values = colorSpecies) +
    # stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.05, step.increase = 0.01) +
    labs(y = 'pN/pS',  x = 'Species') +
    facet_grid(.~Gene, scales = 'free', space = 'free') +
    mytheme + theme(strip.text = element_text(size = 12),
                    # axis.text.x = element_text(angle = 30, hjust = 1),
                    axis.text.x = element_blank(), legend.margin = margin(5,5,-5,5),
                    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 'cm'), legend.position = 'top')
  cairo_pdf('./plot/compare_gene_pNpS_for_each_genes.pdf', width = 6.65, height = 3.48)
  print(pbO)
  dev.off()
  pbO = ggboxplot(tmpDF, x = 'Species', y = 'dNdS',
                  color = 'Species', 
                  # add = "jitter",
                  add.params = list(shape = 21, size = 1.5, color = 'Species', fill = 'Species', alpha = 0.7),
                  bxp.errorbar = T, notch = F) +
    scale_color_manual(values = colorSpecies) +
    scale_fill_manual(values = colorSpecies) +
    # stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.05, step.increase = 0.01) +
    labs(y = 'dN/dS',  x = 'Species') +
    facet_grid(.~Gene, scales = 'free', space = 'free') +
    mytheme + theme(strip.text = element_text(size = 12),
                    # axis.text.x = element_text(angle = 30, hjust = 1),
                    axis.text.x = element_blank(), legend.margin = margin(5,5,-5,5),
                    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 'cm'), legend.position = 'top')
  cairo_pdf('./plot/compare_gene_dNdS_for_each_genes.pdf', width = 6.65, height = 3.48)
  print(pbO)
  dev.off()
}

