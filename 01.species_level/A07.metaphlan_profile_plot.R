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

###### heatmap #####
cairo_pdf('./plot/abu_heatmap_6.pdf', width = 9.82, height =6.21, onefile = T)
if(T){
  dim(abuMatOther)
  mat = abuMat2
  kingdom = c()
  for (s in row.names(abuMat2)) {
    if(s %in% c('Actinomycetaceae_SGB989', 'Hungateiclostridiaceae_SGB4003' )){
      kk = 'Bacteria'
    }else{
      kk = clade2taxidDF$clade_name[grep(s, clade2taxidDF$clade_name)][1] %>% str_extract(., 'k__\\w+')
    }
    kingdom = c(kingdom, kk)
  }
  kingdom %<>% str_remove('k__') %>% str_replace(., 'Eukaryota', 'Fungi')
  table(kingdom)
  row_ha = rowAnnotation(Kingdom = kingdom, col = list(Kingdom = unlist(colorKingdom)))
  Vagitype = factor(sampleVagitypes$vagitypes_Fettweis2[match(colnames(abuMat2), sampleVagitypes$seq_ID)], levels = names(colorVagitype))
  CST = factor(sampleVagitypes$vagitypes_Susan[match(colnames(abuMat2), sampleVagitypes$seq_ID)], levels = names(colorCST))
  column_ha = HeatmapAnnotation(CST = CST,Vagitype = Vagitype, 
                                Term = swabData$Term2[match(colnames(abuMat2),swabData$seq_ID)], 
                                # Trimester = swabData$trimester[match(colnames(abuMat2),swabData$seq_ID)],
                                Ethnicity = swabData$Ethnicity2[match(colnames(abuMat2),swabData$seq_ID)], 
                                "Shannon diveristy" = anno_barplot(shannonDF0$diversity_shannon[match(colnames(abuMat2), shannonDF0$sample)],
                                                                 gp = gpar(border = NA, fill="#FDB462", lty="blank"),  
                                                                 axis_param = list(gp = gpar(fontsize=8))),
                                annotation_height = c(0.55, 0.55,0.55, 0.55, 1.7), height = unit(3, 'cm'),
                                col = list(Vagitype = colorVagitype, 
                                           CST = colorCST,
                                           Term = colorTerm2,
                                           Trimester = colorTrimester,
                                           Ethnicity = colorEthnicity2),
                                annotation_legend_param = list(
                                  Vagitype = list( 
                                    # title_gp = gpar(fontsize = 16), 
                                    # labels_gp = gpar(fontsize = 8), 
                                    ncol=2)))
  row.names(mat) %<>% str_replace_all(., '_', ' ')
  htColor = colorRampPalette(c('black', '#00EAE7'))(100)
  htColor = colorRampPalette(c('black', 'green'))(100)
  htColor = colorRampPalette(c("black","#007FFF", "cyan",
                               "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(100)
  ht = Heatmap(mat, name = "Relative abundance", column_dend_side = 'top',
          top_annotation = column_ha, cluster_columns = F,
          clustering_method_columns = 'ward.D2',
          heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(3.3, "cm")),
          left_annotation = row_ha, show_column_names = F, 
          col = htColor)
  draw(ht, merge_legend = TRUE, heatmap_legend_side = "bottom",
       annotation_legend_side = "bottom", padding = unit(2, 'cm'))
  # ht2 = Heatmap(mat, name = "Relative abundance", column_dend_side = 'top',
  #              top_annotation = column_ha, cluster_columns = T,
  #              clustering_method_columns = 'ward.D2',
  #              right_annotation = row_ha, show_column_names = F, 
  #              col = colorRampPalette(c('black', '#00EAE7'))(100))
  # draw(ht2, merge_legend = TRUE, heatmap_legend_side = "botetom",
  #      annotation_legend_side = "bottom", padding = unit(2, 'cm'))
}
dev.off()

###### stacked barplot ######
if(T){
  colSums(abuMat3) %>% range()
  rowSums(abuMat3) %>% sort(.,decreasing = T)
  metaphlanDF = melt(abuMat3)
  colnames(metaphlanDF) = c('Species', 'seq_ID', 'relative_abundance')
  range(metaphlanDF$relative_abundance)
  metaphlanDF$patient_id = swabData$SampleID[match(metaphlanDF$seq_ID, swabData$seq_ID)]
  metaphlanDF$sample_id = swabData$pt_ID.u[match(metaphlanDF$seq_ID, swabData$seq_ID)]
  
  # metaphlanDF$patient_id %<>% str_replace_all(., 'KIT', 'NC')
  # orderpt = sort(metaphlanDF$patient_id) %>% unique()
  # metaphlanDF$patient_id %<>% factor(., levels = c(orderpt[orderpt != 'NC'], 'NC'))
  # metaphlanDF$timepoint = str_split_fixed(metaphlanDF$sample_id, '_|-', 2)[, 2]
  # metaphlanDF$timepoint = str_split_fixed(metaphlanDF$timepoint, '_|-', 2)[, 1]
  # metaphlanDF$timepoint %<>% str_replace_all(., 'N', '')
  # metaphlanDF$timepoint %<>% factor(., levels = unique(str_sort(metaphlanDF$timepoint, numeric = TRUE)))
  
  metaphlanDF$Species  %<>% as.character() %>% str_replace_all(., '_', ' ')
  table(metaphlanDF$Species)
  metaphlanDF$Species %<>% factor(., levels = names(colorSpecies))
  metaphlanDF$pt_ID = str_split_fixed(metaphlanDF$sample_id, '_', 2)[, 1]
  
  metaphlanDF$Pregnancy.period = metaData$Pregnancy.period[match(metaphlanDF$patient_id, metaData$SampleID)]
  metaphlanDF$trimester = swabData$trimester[match(metaphlanDF$sample_id, swabData$pt_ID.u)]
  metaphlanDF$week = swabData$Sample_GA[match(metaphlanDF$sample_id, swabData$pt_ID.u)]
  metaphlanDF$week_PP = metaphlanDF$week
  metaphlanDF$week_PP[metaphlanDF$trimester == 'P'] = 'P'
  metaphlanDF$order = f_pad_zero(str_split_fixed(metaphlanDF$sample_id, '_', 2)[, 2])
  metaphlanDF$order[metaphlanDF$trimester == 'P'] = 'P'
  metaphlanDF$Term_char = metaData$Term_char[match(metaphlanDF$patient_id, metaData$SampleID)]
  metaphlanDF$Term_char2 = metaData$Term_char2[match(metaphlanDF$patient_id, metaData$SampleID)]
  metaphlanDF$Ethnicity = metaData$Ethnicity[match(metaphlanDF$patient_id, metaData$SampleID)] %>% as.character()
  metaphlanDF$Ethnicity[metaphlanDF$Ethnicity == 'Pacific Islander'] = 'Pacific'
  metaphlanDF$Ethnicity %<>% factor(., levels = c('White', 'Asian', 'Latina', 'Black', 'Pacific', 'Other'))
  metaphlanDF$Ethnicity2 = metaData$Ethnicity2[match(metaphlanDF$patient_id, metaData$SampleID)]
  # metaphlanDF$xs = paste0(metaphlanDF$patient_id , '_',  metaphlanDF$trimester, '_', metaphlanDF$week)
  # metaphlanDF$xs[grep('KIT',metaphlanDF$pt_ID)] = metaphlanDF$sample_id[grep('KIT',metaphlanDF$pt_ID)]
  # metaphlanDF$xs %<>% str_remove_all(., '_plate|KIT-')
  # metaphlanDF$xs %<>% factor(., levels = unique(metaphlanDF$xs))
  # metaphlanDF = arrange(metaphlanDF, Term_char,Ethnicity,  patient_id)
  metaphlanDF = arrange(metaphlanDF, patient_id)
  metaphlanDF$patient_id_term = paste0(metaphlanDF$patient_id, ' ',
                                       # metaphlanDF$Term_char, ' ',
                                       metaphlanDF$Ethnicity)
  # metaphlanDF$patient_id_term[metaphlanDF$patient_id_term == 'NC NA'] = 'NC'
  orderpt = metaphlanDF$patient_id_term %>% unique()
  metaphlanDF$patient_id_term %<>% factor(., levels = orderpt)
}
save(metaphlanDF, file = './data/RData/metaphlanDF.RData')
###### stacked barplot - overall ######
cairo_pdf('./plot/sp_abu_percentage_barplot_by_metaphlan4_3.pdf', 
          width = 7.48, height = 8.48, onefile = T)
cairo_pdf('./plot/sp_abu_percentage_barplot_by_metaphlan4_4.pdf', 
          width = 8.42, height = 7.7, onefile = T)
if(T){
  # p1 = ggplot(metaphlanDF, aes(x = week_PP, y = relative_abundance, fill = Species)) +
  #   geom_bar(width = 0.95, position = position_fill(reverse = TRUE), stat = 'identity') + 
  #   facet_wrap(~patient_id_term, scales="free_x", nrow = 6, ncol = 6) +
  #   xlab("Sample") +
  #   # ylab('Relative abundance') +
  #   ylab('') +
  #   scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  #   scale_fill_manual(name = 'Species', values = colorSpecies) +
  #   guides(fill = guide_legend(nrow = 4, title.position = 'top')) +
  #   mytheme2 + theme(legend.position = 'bottom', 
  #                    strip.text.x = element_text(size = 13), 
  #                    axis.text.x = element_text(size = 5))
  # print(p1)
  
  metaphlanDF2 = metaphlanDF[metaphlanDF$Species != 'Unclassified',]
  metaphlanDF2$Species %<>% factor(., levels = names(colorSpecies[1:(length(colorSpecies)-1)]))
  p2 = ggplot(metaphlanDF2, aes(x = order, y = relative_abundance, fill = Species)) +
    geom_bar(width = 0.95, position = position_stack(reverse = TRUE), stat = 'identity') + 
    facet_wrap(~patient_id_term, scales="free_x", nrow = 5, ncol =7) +
    # facet_wrap(~patient_id_term, scales="free_x", nrow = 3, ncol = 12) +
    xlab("Sample") +
    ylab('Relative abundance (%)') +
    # ylab('') +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    scale_fill_manual(name = 'Species', values = colorSpecies[1:(length(colorSpecies)-1)]) +
    guides(fill = guide_legend(nrow = 6, title.position = 'top')) +
    mytheme2 + theme(legend.position = 'bottom', legend.box.margin = margin(-10,35,0,-20),
                     legend.margin = margin(0,0,0,0),
                     legend.spacing.x = unit(0.1,'cm'),legend.key.size = unit(0.5, 'cm'),
                     strip.text.x = element_text(size = 14, margin = margin(0.1,0,0.1,0, 'cm')), 
                     axis.text.x = element_text(size = 4), panel.spacing = unit(0.1, 'cm'),
                     axis.ticks.length = unit(0.8,'mm'))
  print(p2)
  
}
dev.off()
###### stacked barplot - facet by trimester ######
# for each person, pick the first sample of each trimester
if(T){
  swabData_sub = swabData %>% filter(Trimester != 'NC') %>% 
    group_by(pt_ID, trimester) %>% arrange(Sample_GA) %>% filter(row_number()==1)
  table(swabData_sub$pt_ID)
  table(swabData_sub$trimester, swabData_sub$Term)
  metaphlanDF_sub = metaphlanDF[metaphlanDF$seq_ID %in% swabData_sub$seq_ID,]
  metaphlanDF_sub$Vagitype = sampleVagitypes$vagitypes_Fettweis[match(metaphlanDF_sub$seq_ID, sampleVagitypes$seq_ID)]
  # abuMat3_sub = t(abuMat3[match(names(colorSpecies), row.names(abuMat3) %>% str_replace_all(., '_', ' ')), colnames(abuMat3) %in% swabData_sub$seq_ID]) %>% as.data.frame()
  # abuMat3_sub$Vagitype = sampleVagitypes$vagitypes_Fettweis[match(row.names(abuMat3_sub), sampleVagitypes$seq_ID)]
  # sampleOrder = abuMat3_sub[do.call(order, c(abuMat3_sub, list(decreasing=TRUE))),] %>% row.names()
  # sampleOrder = abuMat3_sub[order(abuMat3_sub$Unclassified),] %>% row.names()
  metaphlanDF_sub = arrange(metaphlanDF_sub, Vagitype, -relative_abundance)
  sampleOrder = unique(metaphlanDF_sub$seq_ID)
    
  setdiff(names(colorSpecies), row.names(abuMat3) %>% str_replace_all(., '_', ' '))
  setdiff(row.names(abuMat3) %>% str_replace_all(., '_', ' '), names(colorSpecies))
}
cairo_pdf("./plot/sp_abu_percentage_barplot_by_term_and_trimester.pdf", onefile = T,
          width = 9.2, height = 8)
if(T){
  metaphlanDF_sub$seq_ID %<>% factor(., levels = sampleOrder)
  metaphlanDF_sub_TB = metaphlanDF_sub[metaphlanDF_sub$Term_char2 == 'Preterm',]
  metaphlanDF_sub_PTB = metaphlanDF_sub[metaphlanDF_sub$Term_char2 == 'Full-term',]
  # metaphlanDF_sub$seq_ID %<>% as.character()
  if(T){
    p1_TB = ggplot(metaphlanDF_sub_TB, aes(x = seq_ID, y = relative_abundance, fill = Species)) +
      geom_bar(width = 1, position = position_fill(reverse = TRUE), stat = 'identity') + 
      facet_grid(Term_char2~trimester, scales = 'free_x', space="free") +
      xlab("") +
      ylab('Relative abundance (%)') +
      scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_fill_manual(name = 'Species', values = colorSpecies) +
      mytheme2 + theme(legend.position = 'none', 
                       strip.text = element_text(size = 13), 
                       axis.title.x = element_blank(),
                       axis.text.x = element_blank())
    p1_PTB = ggplot(metaphlanDF_sub_PTB, aes(x = seq_ID, y = relative_abundance, fill = Species)) +
      geom_bar(width = 1, position = position_fill(reverse = TRUE), stat = 'identity') + 
      facet_grid(Term_char2~trimester, scales = 'free_x', space="free") +
      xlab("Sample") +
      ylab('Relative abundance (%)') +
      scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_fill_manual(name = 'Species', values = colorSpecies) +
      guides(fill = guide_legend(nrow = 6, title.position = 'top')) +
      mytheme2 + theme(legend.position = 'bottom', 
                       strip.text = element_text(size = 13), 
                       axis.text.x = element_blank())
  
    gt1_TB = ggplot_gtable(ggplot_build(p1_TB))
    facetInd = gt1_TB$layout$l[grep('panel-1-*', gt1_TB$layout$name)]
    gt1_PTB = ggplot_gtable(ggplot_build(p1_PTB))
    gt1_PTB$layout$l[grep('panel-1-*', gt1_PTB$layout$name)]
    gt1_TB$widths[facetInd] = gt1_PTB$widths[facetInd] = (gt1_PTB$widths[facetInd] + gt1_TB$widths[facetInd])/2
  
    grid.arrange(grobs=list(gt1_TB, gt1_PTB), nrow=2,
                 heights=unit(c(6,11.7), c("cm", "cm")))
  }
  if(T){
    metaphlanDF2_sub_TB = metaphlanDF_sub_TB[metaphlanDF_sub_TB$Species != 'Unclassified',]
    metaphlanDF2_sub_TB$Species %<>% factor(., levels = names(colorSpecies[1:(length(colorSpecies)-1)]))
    metaphlanDF2_sub_PTB = metaphlanDF_sub_PTB[metaphlanDF_sub_PTB$Species != 'Unclassified',]
    metaphlanDF2_sub_PTB$Species %<>% factor(., levels = names(colorSpecies[1:(length(colorSpecies)-1)]))
    p2_TB = ggplot(metaphlanDF2_sub_TB, aes(x = seq_ID, y = relative_abundance, fill = Species)) +
      geom_bar(width = 1, position = position_stack(reverse = TRUE), stat = 'identity') + 
      facet_grid(Term_char2~trimester, scales = 'free_x', space="free") +
      xlab("") +
      ylab('Relative abundance (%)') +
      scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_fill_manual(name = 'Species', values = colorSpecies[1:(length(colorSpecies)-1)]) +
      mytheme2 + theme(legend.position = 'none', 
                       strip.text = element_text(size = 13), 
                       axis.title.x = element_blank(),
                       axis.text.x = element_blank())
    p2_PTB = ggplot(metaphlanDF2_sub_PTB, aes(x = seq_ID, y = relative_abundance, fill = Species)) +
      geom_bar(width = 1, position = position_stack(reverse = TRUE), stat = 'identity') + 
      facet_grid(Term_char2 ~ trimester, scales = 'free_x', space="free") +
      xlab("Sample") +
      ylab('Relative abundance (%)') +
      scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_fill_manual(name = 'Species', values = colorSpecies[1:(length(colorSpecies)-1)]) +
      guides(fill = guide_legend(nrow = 6, title.position = 'top')) +
      mytheme2 + theme(legend.position = 'bottom', 
                       strip.text = element_text(size = 13), 
                       axis.text.x = element_blank())
    
    gt2_TB = ggplot_gtable(ggplot_build(p2_TB))
    facetInd = gt2_TB$layout$l[grep('panel-1-*', gt2_TB$layout$name)]
    gt2_PTB = ggplot_gtable(ggplot_build(p2_PTB))
    gt2_PTB$layout$l[grep('panel-1-*', gt2_PTB$layout$name)]
    gt2_TB$widths[facetInd] = gt2_PTB$widths[facetInd] = (gt2_PTB$widths[facetInd] + gt1_TB$widths[facetInd])/2
    grid.arrange(grobs=list(gt2_TB, gt2_PTB), nrow=2,
                 heights=unit(c(6,11.7), c("cm", "cm")))
  }
}
dev.off()
###### stacked barplot - facet by ethnicity vs trimester (T2,T3) t######
if(T){
  swabDataT2 = swabData[swabData$trimester == 'T2',] %>% filter(Trimester != 'NC') %>% 
    group_by(pt_ID) %>% arrange(Sample_GA) %>% filter(row_number()==1)
  swabDataT3 = swabData[swabData$trimester == 'T3',] %>% filter(Trimester != 'NC') %>% 
    group_by(pt_ID) %>% arrange(-Sample_GA) %>% filter(row_number()==1)
  swabData_sub = rbind(swabDataT2, swabDataT3)
  metaphlanDF_sub = metaphlanDF[metaphlanDF$seq_ID %in% swabData_sub$seq_ID,]
  metaphlanDF_sub$Vagitype = sampleVagitypes$vagitypes_Fettweis[match(metaphlanDF_sub$seq_ID, sampleVagitypes$seq_ID)]
  # abuMat3_sub = t(abuMat3[match(names(colorSpecies), row.names(abuMat3) %>% str_replace_all(., '_', ' ')), colnames(abuMat3) %in% swabData_sub$seq_ID]) %>% as.data.frame()
  # abuMat3_sub$Vagitype = sampleVagitypes$vagitypes_Fettweis[match(row.names(abuMat3_sub), sampleVagitypes$seq_ID)]
  # sampleOrder = abuMat3_sub[do.call(order, c(abuMat3_sub, list(decreasing=TRUE))),] %>% row.names()
  # sampleOrder = abuMat3_sub[order(abuMat3_sub$Unclassified),] %>% row.names()
  metaphlanDF_sub = arrange(metaphlanDF_sub, Vagitype, -relative_abundance)
  sampleOrder = unique(metaphlanDF_sub$seq_ID)
  
  setdiff(names(colorSpecies), row.names(abuMat3) %>% str_replace_all(., '_', ' '))
  setdiff(row.names(abuMat3) %>% str_replace_all(., '_', ' '), names(colorSpecies))
}
cairo_pdf("./plot/sp_abu_percentage_barplot_by_ethnicity_and_trimesterT2T3t.pdf", onefile = T,
          width = 7.71, height = 4.02)
if(T){
  metaphlanDF_sub = metaphlanDF_sub[metaphlanDF_sub$Species != 'Unclassified',]
  metaphlanDF_sub$Species %<>% factor(., levels = names(colorSpecies[1:(length(colorSpecies)-1)]))
  metaphlanDF_sub$seq_ID %<>% factor(., levels = sampleOrder)
  # metaphlanDF_sub$seq_ID %<>% as.character()
  if(T){
    metaphlanDF_sub$trimester %<>% as.character()
    metaphlanDF_sub_List = metaphlanDF_sub %>% group_by(trimester) %>% group_split()
    p3_T2 = ggplot(metaphlanDF_sub_List[[1]], aes(x = seq_ID, y = relative_abundance, fill = Species)) +
      geom_bar(width = 1, position = position_stack(reverse = TRUE), stat = 'identity') + 
      facet_grid(trimester~Ethnicity2, scales = 'free_x', space="free") +
      xlab("") +
      ylab('') +
      scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_fill_manual(name = 'Species', values = colorSpecies[names(colorSpecies) != "Unclassified"]) +
      mytheme2 + theme(legend.position = 'none', 
                       strip.text = element_text(size = 13), 
                       axis.title.x = element_blank(),
                       axis.text.x = element_blank())
    p3_T3 = ggplot(metaphlanDF_sub_List[[2]], aes(x = seq_ID, y = relative_abundance, fill = Species)) +
      geom_bar(width = 1, position = position_stack(reverse = TRUE), stat = 'identity') + 
      facet_grid(trimester~Ethnicity2, scales = 'free_x', space="free") +
      xlab("") +
      ylab('') +
      scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_fill_manual(name = 'Species', values = colorSpecies[names(colorSpecies) != "Unclassified"]) +
      guides(fill = guide_legend(nrow = 4, title.position = 'left')) +
      mytheme2 + theme(legend.position = 'bottom', 
                       legend.box.margin = margin(0,0,0,-15),
                       legend.margin =  margin(0,0,0,-15),
                       strip.background.x = element_blank(),
                       strip.text.x = element_blank(),
                       strip.text.y = element_text(size = 13), 
                       axis.title.x = element_blank(),
                       axis.text.x = element_blank())
    gt3_T2 = ggplot_gtable(ggplot_build(p3_T2))
    facetInd = gt3_T2$layout$l[grep('panel-1-*', gt3_T2$layout$name)]
    gt3_T3 = ggplot_gtable(ggplot_build(p3_T3))
    gt3_T3$layout$l[grep('panel-1-*', gt3_T3$layout$name)]
    gt3_T2$widths[facetInd] =  gt3_T3$widths[facetInd] = (gt3_T2$widths[facetInd] + gt3_T3$widths[facetInd])/2
    grid.arrange(grobs=list(gt3_T2, gt3_T3), nrow=2,
                 heights=unit(c(4,6.5), c("cm", "cm")))
  }
}
dev.off()

###### stacked barplot - facet by ethnicity vs trimester (T2,T3) ######
if(T){
  swabDataT2 = swabData[swabData$trimester == 'T2',] %>% filter(Trimester != 'NC') %>% 
    group_by(pt_ID) %>% arrange(Sample_GA) %>% filter(row_number()==1)
  swabDataT3 = swabData[swabData$trimester == 'T3',] %>% filter(Trimester != 'NC') %>% 
    group_by(pt_ID) %>% arrange(-Sample_GA) %>% filter(row_number()==1)
  swabData_sub = rbind(swabDataT2, swabDataT3)
  metaphlanDF_sub = metaphlanDF[metaphlanDF$seq_ID %in% swabData_sub$seq_ID,]
  metaphlanDF_sub$Vagitype = sampleVagitypes$vagitypes_Fettweis[match(metaphlanDF_sub$seq_ID, sampleVagitypes$seq_ID)]
  # abuMat3_sub = t(abuMat3[match(names(colorSpecies), row.names(abuMat3) %>% str_replace_all(., '_', ' ')), colnames(abuMat3) %in% swabData_sub$seq_ID]) %>% as.data.frame()
  # abuMat3_sub$Vagitype = sampleVagitypes$vagitypes_Fettweis[match(row.names(abuMat3_sub), sampleVagitypes$seq_ID)]
  # sampleOrder = abuMat3_sub[do.call(order, c(abuMat3_sub, list(decreasing=TRUE))),] %>% row.names()
  # sampleOrder = abuMat3_sub[order(abuMat3_sub$Unclassified),] %>% row.names()
  metaphlanDF_sub = arrange(metaphlanDF_sub, Vagitype, -relative_abundance)
  sampleOrder = unique(metaphlanDF_sub$seq_ID)
  
  setdiff(names(colorSpecies), row.names(abuMat3) %>% str_replace_all(., '_', ' '))
  setdiff(row.names(abuMat3) %>% str_replace_all(., '_', ' '), names(colorSpecies))
}
cairo_pdf("./plot/sp_abu_percentage_barplot_by_ethnicity_and_trimesterT2T3.pdf", onefile = T,
          width = 4.19, height = 6.36)
if(T){
  metaphlanDF_sub = metaphlanDF_sub[metaphlanDF_sub$Species != 'Unclassified',]
  metaphlanDF_sub$Species %<>% factor(., levels = names(colorSpecies[1:(length(colorSpecies)-1)]))
  metaphlanDF_sub$seq_ID %<>% factor(., levels = sampleOrder)
  # metaphlanDF_sub$seq_ID %<>% as.character()
  if(T){
    metaphlanDF_sub_List = metaphlanDF_sub %>% group_by(Ethnicity2) %>% group_split()
    p3_white = ggplot(metaphlanDF_sub_List[[1]], aes(x = seq_ID, y = relative_abundance, fill = Species)) +
      geom_bar(width = 1, position = position_stack(reverse = TRUE), stat = 'identity') + 
      facet_grid(Ethnicity2~trimester, scales = 'free_x', space="free") +
      xlab("") +
      ylab('') +
      scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_fill_manual(name = 'Species', values = colorSpecies[names(colorSpecies) != "Unclassified"]) +
      mytheme2 + theme(legend.position = 'none', 
                       strip.text = element_text(size = 13), 
                       axis.title.x = element_blank(),
                       axis.text.x = element_blank())
    p3_asian = ggplot(metaphlanDF_sub_List[[2]], aes(x = seq_ID, y = relative_abundance, fill = Species)) +
      geom_bar(width = 1, position = position_stack(reverse = TRUE), stat = 'identity') + 
      facet_grid(Ethnicity2~trimester, scales = 'free_x', space="free") +
      xlab("") +
      ylab('') +
      scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_fill_manual(name = 'Species', values = colorSpecies[names(colorSpecies) != "Unclassified"]) +
      mytheme2 + theme(legend.position = 'none', 
                       strip.background.x = element_blank(),
                       strip.text.x = element_blank(),
                       strip.text.y = element_text(size = 13),
                       axis.title.x = element_blank(),
                       axis.text.x = element_blank())
    p3_latino = ggplot(metaphlanDF_sub_List[[3]], aes(x = seq_ID, y = relative_abundance, fill = Species)) +
      geom_bar(width = 1, position = position_stack(reverse = TRUE), stat = 'identity') + 
      facet_grid(Ethnicity2~trimester, scales = 'free_x', space="free") +
      xlab("") +
      ylab('') +
      scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0,1)) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_fill_manual(name = 'Species', values = colorSpecies[names(colorSpecies) != "Unclassified"]) +
      mytheme2 + theme(legend.position = 'none', 
                       strip.background.x = element_blank(),
                       strip.text.x = element_blank(),
                       strip.text.y = element_text(size = 13),
                       axis.title.x = element_blank(),
                       axis.text.x = element_blank())
    p3_black = ggplot(metaphlanDF_sub_List[[4]], aes(x = seq_ID, y = relative_abundance, fill = Species)) +
      geom_bar(width = 1, position = position_stack(reverse = TRUE), stat = 'identity') + 
      facet_grid(Ethnicity2~trimester, scales = 'free_x', space="free") +
      xlab("") +
      ylab('') +
      scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_fill_manual(name = 'Species', values = colorSpecies[names(colorSpecies) != "Unclassified"]) +
      mytheme2 + theme(legend.position = 'none', 
                       strip.background.x = element_blank(),
                       strip.text.x = element_blank(),
                       strip.text.y = element_text(size = 13), 
                       axis.title.x = element_blank(),
                       axis.text.x = element_blank())
    p3_others = ggplot(metaphlanDF_sub_List[[5]], aes(x = seq_ID, y = relative_abundance, fill = Species)) +
      geom_bar(width = 1, position = position_stack(reverse = TRUE), stat = 'identity') + 
      facet_grid(Ethnicity2~trimester, scales = 'free_x', space="free") +
      xlab("") +
      ylab('') +
      scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0,1)) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_fill_manual(name = 'Species', values = colorSpecies[names(colorSpecies) != "Unclassified"]) +
      guides(fill = guide_legend(ncol= 3, title.position = 'top')) +
      mytheme2 + theme(legend.position = 'none', 
                       strip.background.x = element_blank(),
                       strip.text.x = element_blank(),
                       strip.text.y = element_text(size = 13), 
                       axis.title.x = element_blank(),
                       axis.text.x = element_blank())
    gt3_white = ggplot_gtable(ggplot_build(p3_white))
    facetInd = gt3_white$layout$l[grep('panel-1-*', gt3_white$layout$name)]
    gt3_asian = ggplot_gtable(ggplot_build(p3_asian))
    gt3_asian$layout$l[grep('panel-1-*', gt3_asian$layout$name)]
    gt3_latino = ggplot_gtable(ggplot_build(p3_latino))
    gt3_latino$layout$l[grep('panel-1-*', gt3_latino$layout$name)]
    gt3_black = ggplot_gtable(ggplot_build(p3_black))
    gt3_black$layout$l[grep('panel-1-*', gt3_black$layout$name)]
    gt3_other = ggplot_gtable(ggplot_build(p3_others))
    gt3_other$layout$l[grep('panel-1-*', gt3_other$layout$name)]
    gt3_white$widths[facetInd] =  gt3_asian$widths[facetInd] = gt3_latino$widths[facetInd] = gt3_black$widths[facetInd] = gt3_other$widths[facetInd]=
      (gt3_white$widths[facetInd] + gt3_asian$widths[facetInd] + gt3_latino$widths[facetInd] + gt3_black$widths[facetInd] + gt3_other$widths[facetInd])/2
    grid.arrange(grobs=list(gt3_white, gt3_asian, gt3_latino, gt3_black,gt3_other), nrow=5,
                 heights=unit(c(6, 5,5,5,5)-2, "cm"))
    p3_othersLegend = ggplot(metaphlanDF_sub_List[[5]], aes(x = seq_ID, y = relative_abundance, fill = Species)) +
      geom_bar(width = 1, position = position_stack(reverse = TRUE), stat = 'identity') + 
      facet_grid(Ethnicity2~trimester, scales = 'free_x', space="free") +
      xlab("") +
      ylab('') +
      scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0,1)) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_fill_manual(name = 'Species', values = colorSpecies[names(colorSpecies) != "Unclassified"]) +
      guides(fill = guide_legend(nrow= 4, title.position = 'left')) +
      mytheme2 + theme(legend.position = 'bottom', 
                       legend.title = element_text(size = 18),
                       legend.text = element_text(size = 18),
                       legend.spacing.x = unit(0.2, 'cm'),
                       strip.background.x = element_blank(),
                       strip.text.x = element_blank(),
                       strip.text.y = element_text(size = 13), 
                       axis.title.x = element_blank())

  }
}
dev.off()
cairo_pdf(file = './plot/sp_abu_percentage_barplot_legend.pdf',width = 16.44, height = 3.29)
print(p3_othersLegend)
dev.off()
###### stacked barplot - facet by ethnicity vs trimester ######
if(T){
  swabData_sub = swabData %>% filter(Trimester != 'NC') %>% 
    group_by(pt_ID, trimester) %>% arrange(Sample_GA) %>% filter(row_number()==1)
  metaphlanDF_sub = metaphlanDF[metaphlanDF$seq_ID %in% swabData_sub$seq_ID,]
  metaphlanDF_sub$Vagitype = sampleVagitypes$vagitypes_Fettweis[match(metaphlanDF_sub$seq_ID, sampleVagitypes$seq_ID)]
  # abuMat3_sub = t(abuMat3[match(names(colorSpecies), row.names(abuMat3) %>% str_replace_all(., '_', ' ')), colnames(abuMat3) %in% swabData_sub$seq_ID]) %>% as.data.frame()
  # abuMat3_sub$Vagitype = sampleVagitypes$vagitypes_Fettweis[match(row.names(abuMat3_sub), sampleVagitypes$seq_ID)]
  # sampleOrder = abuMat3_sub[do.call(order, c(abuMat3_sub, list(decreasing=TRUE))),] %>% row.names()
  # sampleOrder = abuMat3_sub[order(abuMat3_sub$Unclassified),] %>% row.names()
  metaphlanDF_sub = arrange(metaphlanDF_sub, Vagitype, -relative_abundance)
  sampleOrder = unique(metaphlanDF_sub$seq_ID)
  
  setdiff(names(colorSpecies), row.names(abuMat3) %>% str_replace_all(., '_', ' '))
  setdiff(row.names(abuMat3) %>% str_replace_all(., '_', ' '), names(colorSpecies))
}
cairo_pdf("./plot/sp_abu_percentage_barplot_by_ethnicity_and_trimester.pdf", onefile = T,
          width = 9.2, height = 9)
if(T){
  metaphlanDF_sub = metaphlanDF_sub[metaphlanDF_sub$Species != 'Unclassified',]
  metaphlanDF_sub$Species %<>% factor(., levels = names(colorSpecies[1:(length(colorSpecies)-1)]))
  metaphlanDF_sub$seq_ID %<>% factor(., levels = sampleOrder)
  # metaphlanDF_sub$seq_ID %<>% as.character()
  if(T){
    metaphlanDF_sub_List = metaphlanDF_sub %>% group_by(Ethnicity) %>% group_split()
    p3_white = ggplot(metaphlanDF_sub_List[[1]], aes(x = seq_ID, y = relative_abundance, fill = Species)) +
      geom_bar(width = 1, position = position_stack(reverse = TRUE), stat = 'identity') + 
      facet_grid(Ethnicity~trimester, scales = 'free_x', space="free") +
      xlab("") +
      ylab('Relative abundance (%)') +
      scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_fill_manual(name = 'Species', values = colorSpecies[names(colorSpecies) != "Unclassified"]) +
      mytheme2 + theme(legend.position = 'none', 
                       strip.text = element_text(size = 13), 
                       axis.title.x = element_blank(),
                       axis.text.x = element_blank())
    p3_latino = ggplot(metaphlanDF_sub_List[[3]], aes(x = seq_ID, y = relative_abundance, fill = Species)) +
      geom_bar(width = 1, position = position_stack(reverse = TRUE), stat = 'identity') + 
      facet_grid(Ethnicity~trimester, scales = 'free_x', space="free") +
      xlab("") +
      ylab('Relative abundance (%)') +
      scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0,1)) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_fill_manual(name = 'Species', values = colorSpecies[names(colorSpecies) != "Unclassified"]) +
      mytheme2 + theme(legend.position = 'none', 
                       strip.text = element_text(size = 13), 
                       axis.title.x = element_blank(),
                       axis.text.x = element_blank())
    p3_black = ggplot(metaphlanDF_sub_List[[4]], aes(x = seq_ID, y = relative_abundance, fill = Species)) +
      geom_bar(width = 1, position = position_stack(reverse = TRUE), stat = 'identity') + 
      facet_grid(Ethnicity~trimester, scales = 'free_x', space="free") +
      xlab("") +
      ylab('Relative abundance (%)') +
      scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_fill_manual(name = 'Species', values = colorSpecies[names(colorSpecies) != "Unclassified"]) +
      guides(fill = guide_legend(nrow = 6, title.position = 'top')) +
      mytheme2 + theme(legend.position = 'bottom', 
                       strip.text = element_text(size = 13), 
                       axis.title.x = element_blank(),
                       axis.text.x = element_blank())
    gt3_white = ggplot_gtable(ggplot_build(p3_white))
    facetInd = gt3_white$layout$l[grep('panel-1-*', gt3_white$layout$name)]
    gt3_latino = ggplot_gtable(ggplot_build(p3_latino))
    gt3_latino$layout$l[grep('panel-1-*', gt3_latino$layout$name)]
    gt3_black = ggplot_gtable(ggplot_build(p3_black))
    gt3_black$layout$l[grep('panel-1-*', gt3_black$layout$name)]
    gt3_white$widths[facetInd] = gt3_latino$widths[facetInd] = gt3_black$widths[facetInd] =
      (gt3_white$widths[facetInd] + gt3_latino$widths[facetInd] + gt3_black$widths[facetInd])/2
    grid.arrange(grobs=list(gt3_white, gt3_latino, gt3_black), nrow=3,
                 heights=unit(c(6, 6,11), c("cm", "cm")))
  }
}
dev.off()
###### stacked barplot - facet by ethnicity vs group ######
# group = early and late
# collected early (the first sample before 23 weeks of gestation) and 
# late (last sample collected after 32 weeks of gestation) during pregnancy.
if(T){
  swabData_frist = swabData %>% filter(Trimester != 'NC' & trimester != 'P' & Sample_GA <= 23) %>% 
    group_by(pt_ID) %>% arrange(Sample_GA) %>% filter(row_number()==1)
  swabData_last = swabData %>% filter(Trimester != 'NC' & trimester != 'P' & Sample_GA >= 32) %>% 
    group_by(pt_ID) %>% arrange(Sample_GA) %>% filter(row_number()==n())
  swabData_post = swabData %>% filter(trimester == 'P') %>% 
    group_by(pt_ID) %>% arrange(Sample_GA) %>% filter(row_number()==1)
  swabData_frist$Group = 'Early'
  swabData_last$Group = 'Late'
  swabData_post$Group = 'Postpartum'
  swabData_sub = rbind(swabData_frist, swabData_last, swabData_post)
  table(swabData_sub$pt_ID)
  table(swabData_sub$Group, swabData_sub$Ethnicity)

  metaphlanDF_sub = metaphlanDF[metaphlanDF$seq_ID %in% swabData_sub$seq_ID,]
  metaphlanDF_sub$Group = swabData_sub$Group[match(metaphlanDF_sub$seq_ID, swabData_sub$seq_ID)]
  metaphlanDF_sub$Vagitype = sampleVagitypes$vagitypes_Fettweis[match(metaphlanDF_sub$seq_ID, sampleVagitypes$seq_ID)]
  # abuMat3_sub = t(abuMat3[match(names(colorSpecies), row.names(abuMat3) %>% str_replace_all(., '_', ' ')), colnames(abuMat3) %in% swabData_sub$seq_ID]) %>% as.data.frame()
  # abuMat3_sub$Vagitype = sampleVagitypes$vagitypes_Fettweis[match(row.names(abuMat3_sub), sampleVagitypes$seq_ID)]
  # sampleOrder = abuMat3_sub[do.call(order, c(abuMat3_sub, list(decreasing=TRUE))),] %>% row.names()
  # sampleOrder = abuMat3_sub[order(abuMat3_sub$Unclassified),] %>% row.names()
  metaphlanDF_sub = arrange(metaphlanDF_sub, Vagitype, -relative_abundance)
  sampleOrder = unique(metaphlanDF_sub$seq_ID)
  
  setdiff(names(colorSpecies), row.names(abuMat3) %>% str_replace_all(., '_', ' '))
  setdiff(row.names(abuMat3) %>% str_replace_all(., '_', ' '), names(colorSpecies))
}
cairo_pdf("./plot/sp_abu_percentage_barplot_by_ethnicity_and_group_2.pdf", onefile = T,
          width = 9.2, height = 15)
if(T){
  metaphlanDF_sub = metaphlanDF_sub[metaphlanDF_sub$Species != 'Unclassified',]
  metaphlanDF_sub$Species %<>% factor(., levels = names(colorSpecies[1:(length(colorSpecies)-1)]))
  metaphlanDF_sub$seq_ID %<>% factor(., levels = sampleOrder)
  table(metaphlanDF_sub$Ethnicity)
  if(T){
    metaphlanDF_sub_List = metaphlanDF_sub %>% group_by(Ethnicity) %>% group_split()
    p3_white = ggplot(metaphlanDF_sub_List[[1]], aes(x = seq_ID, y = relative_abundance, fill = Species)) +
      geom_bar(width = 1, position = position_stack(reverse = TRUE), stat = 'identity') + 
      facet_grid(Ethnicity~Group, scales = 'free_x', space="free") +
      xlab("") +
      ylab('Relative abundance (%)') +
      scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_fill_manual(name = 'Species', values = colorSpecies) +
      mytheme2 + theme(legend.position = 'none', 
                       strip.text = element_text(size = 13), 
                       axis.title.x = element_blank(),
                       axis.text.x = element_blank())
    p3_asian = ggplot(metaphlanDF_sub_List[[2]], aes(x = seq_ID, y = relative_abundance, fill = Species)) +
      geom_bar(width = 1, position = position_stack(reverse = TRUE), stat = 'identity') + 
      facet_grid(Ethnicity~Group, scales = 'free_x', space="free") +
      xlab("") +
      ylab('Relative abundance (%)') +
      scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_fill_manual(name = 'Species', values = colorSpecies) +
      mytheme2 + theme(legend.position = 'none', 
                       strip.text = element_text(size = 13), 
                       axis.title.x = element_blank(),
                       axis.text.x = element_blank())
    p3_latino = ggplot(metaphlanDF_sub_List[[3]], aes(x = seq_ID, y = relative_abundance, fill = Species)) +
      geom_bar(width = 1, position = position_stack(reverse = TRUE), stat = 'identity') + 
      facet_grid(Ethnicity~Group, scales = 'free_x', space="free") +
      xlab("") +
      ylab('Relative abundance (%)') +
      scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0,1)) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_fill_manual(name = 'Species', values = colorSpecies) +
      mytheme2 + theme(legend.position = 'none', 
                       strip.text = element_text(size = 13), 
                       axis.title.x = element_blank(),
                       axis.text.x = element_blank())
    p3_black = ggplot(metaphlanDF_sub_List[[4]], aes(x = seq_ID, y = relative_abundance, fill = Species)) +
      geom_bar(width = 1, position = position_stack(reverse = TRUE), stat = 'identity') + 
      facet_grid(Ethnicity~Group, scales = 'free_x', space="free") +
      xlab("") +
      ylab('Relative abundance (%)') +
      scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_fill_manual(name = 'Species', values = colorSpecies[names(colorSpecies) != "Unclassified"]) +
      guides(fill = guide_legend(nrow = 6, title.position = 'top')) +
      mytheme2 + theme(legend.position = 'none', 
                       strip.text = element_text(size = 13), 
                       axis.title.x = element_blank(),
                       axis.text.x = element_blank())
    p3_Other = ggplot(rbind(metaphlanDF_sub_List[[5]], metaphlanDF_sub_List[[6]]), aes(x = seq_ID, y = relative_abundance, fill = Species)) +
      geom_bar(width = 1, position = position_stack(reverse = TRUE), stat = 'identity') + 
      facet_grid(Ethnicity2~Group, scales = 'free_x', space="free") +
      xlab("") +
      ylab('Relative abundance (%)') +
      scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0,1)) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_fill_manual(name = 'Species', values = colorSpecies[names(colorSpecies) != "Unclassified"]) +
      guides(fill = guide_legend(nrow = 6, title.position = 'top')) +
      mytheme2 + theme(legend.position = 'bottom', 
                       strip.text = element_text(size = 13), 
                       axis.title.x = element_blank(),
                       axis.text.x = element_blank())
    gt3_white = ggplot_gtable(ggplot_build(p3_white))
    facetInd = gt3_white$layout$l[grep('panel-1-*', gt3_white$layout$name)]
    gt3_asian = ggplot_gtable(ggplot_build(p3_asian))
    gt3_asian$layout$l[grep('panel-1-*', gt3_asian$layout$name)]
    gt3_latino = ggplot_gtable(ggplot_build(p3_latino))
    gt3_latino$layout$l[grep('panel-1-*', gt3_latino$layout$name)]
    gt3_black = ggplot_gtable(ggplot_build(p3_black))
    gt3_black$layout$l[grep('panel-1-*', gt3_black$layout$name)]
    gt3_Other = ggplot_gtable(ggplot_build(p3_Other))
    gt3_Other$layout$l[grep('panel-1-*', gt3_Other$layout$name)]
    gt3_white$widths[facetInd] = gt3_asian$widths[facetInd] = gt3_latino$widths[facetInd] = gt3_black$widths[facetInd] = gt3_Other$widths[facetInd] =
      (gt3_white$widths[facetInd] +  gt3_asian$widths[facetInd] + gt3_latino$widths[facetInd] + gt3_black$widths[facetInd])/2
    grid.arrange(grobs=list(gt3_white, gt3_asian, gt3_latino, gt3_black, gt3_Other), nrow = 5,
                 heights=unit(c(6, 6, 6, 6, 11), c("cm", "cm")))
  }
}
dev.off()
table(metaData$Term_char, metaData$Ethnicity)
