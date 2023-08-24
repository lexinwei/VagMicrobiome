source('./script/library_function.R')
load('./data/RData/swabData.RData')
load('./data/RData/metaData.RData')
load('./data/RData/colorList.RData')
SGB2sp = read.csv('./data/SGB_species_to_taxid.csv')
plotFlag = T
plotFlag = F
list.dirs('./strainphlan_results')
#### get distmat of alignment ####
if(T){
  subDir = 'without_NCBI_references'
  # subDir = 'with_NCBI_reference'
  # subDir = 'with_metaphlan4_bins'
  if(plotFlag){pdf(paste0('./plot/kimura_matrix_', subDir, '_heatmap2.pdf'), width = 25, height = 25)}
  dismatList = list()
  strainDir = paste0("./strainphlan_results/", subDir)
  spList = (list.files(strainDir) %>% str_split_fixed(., '\\.', 2))[,1] %>% unique() %>% sort()
  # rmSp = c('Lactobacillus_brevis', 'Lactobacillus_rhamnosus')
  # fqID_without_swabID = read.delim2('./data/fqID_without_swabID.txt')[, 1]
  for(sp in spList){
    spName = SGB2sp$species_name[SGB2sp$SGB_name == sp] %>% str_remove(., 's__') %>% str_replace_all(., '_', ' ')
  
    if(length(spName) == 0){
      spName = sp
      key = sp
    }else{
      key = paste0(spName, ' (', sp, ')')
    }
    if(key == 'SGB4003'){
      key = 'GGB3012 SGB4003'
    }else if(key == 'SGB989'){
      key = 'GGB753 SGB989'
    }
    print(key)
    dismat0 = read.delim2(paste0(strainDir, '/', sp, '.dismat'), skip = 8, sep = '\t', header = F)
    dismat0 = dismat0[, colSums(is.na(dismat0)) == 0]
    colnames(dismat0) = (dismat0[, ncol(dismat0)] %>% str_split_fixed(., ' ', 2))[, 1] 
    dismat0 = dismat0[, 1:(ncol(dismat0)-1)]
    refName = colnames(dismat0)[!startsWith(colnames(dismat0), 'Ming_nova_VS_')]
    dismat0 = apply(dismat0, 2, as.numeric)
    row.names(dismat0) = colnames(dismat0) %<>% str_remove_all(., 'Ming_nova_VS_|_bwa.sorted')
    keepSample = intersect(colnames(dismat0), c(refName, swabData$seq_ID[swabData$Trimester != 'NC']))
    keepSampleName = swabData$pt_ID.u[match(keepSample, swabData$seq_ID)]
    keepSampleName[is.na(keepSampleName)] = paste0("Ref_", keepSample[is.na(keepSampleName)])
    dismat1 = dismat0[match(keepSample, row.names(dismat0)), match(keepSample, colnames(dismat0))]
    row.names(dismat1) = colnames(dismat1) = keepSampleName
    # row.names(dismat0)[is.na(row.names(dismat0))] = 'Reference'
    dismat1[is.na(dismat1)] = t(dismat1)[is.na(dismat1)]
    dismat1 = dismat1/max(dismat1, na.rm = T)
    dismat1[is.nan(dismat1)] = NA
    dismat2 = dismat1
    naInd = which(colSums(is.na(dismat2)) > 0)
    while (length(naInd) > 0) {
      tmpMax = sort(colSums(is.na(dismat2)), decreasing = T)[1]
      dismat2 = dismat2[row.names(dismat2) != names(tmpMax), colnames(dismat2) != names(tmpMax)]
      if(!is.matrix(dismat2)){
        break
      }
      naInd = which(colSums(is.na(dismat2)) > 0)
    }
    
    if(length(refName) >= (ncol(dismat1)-1)){
      cat(spName, '--', sp, 'all are reference, skip...\n')
      next
    }
    if(!is.matrix(dismat2)){
      next
    }
    pS = sum(table(str_extract(colnames(dismat2), 'SF\\d+')) >= 2)
    if(pS < 2 | sum(dismat2) == 0){
      next
    }
    dismatList[[key]] = dismat2
    if(plotFlag){
      # ph = try({
        annoColumns = data.frame(Participant = str_extract(colnames(dismat2), 'SF\\d+'))
        row.names(annoColumns) = colnames(dismat2)
        annoColor = list(Participant = colorParticipant[names(colorParticipant) %in% annoColumns$Participant])
        annoRows = data.frame(Participant = str_extract(row.names(dismat2), 'SF\\d+'))
        row.names(annoRows) = row.names(dismat2)
        ph = pheatmap(1-dismat2, scale = 'none', border_color = NA,
                 color = brewer.pal(n=10, name="PuOr"), main = key, 
                 cellwidth = 10, cellheight = 10, annotation_row = annoRows,
                 annotation_col = annoColumns, annotation_colors = annoColor)
      # })
      # if(class(ph) != 'try-error'){
        print(ph)
      # }
    }
  }
  save(dismatList, file = paste0('./data/RData/', subDir, '.RData'))
  if(plotFlag){dev.off()}
}

##### within and between participant distance ####
load('./data/RData/without_NCBI_references.RData')
load('./data/RData/without_NCBI_references2.RData')
spN_old = names(dismatList) 
spN = str_split_fixed(spN_old, '\\(', 2)[, 1] %>% str_trim()
spN[spN %in% spN[duplicated(spN)]] = spN_old[spN %in% spN[duplicated(spN)]]
spN[spN == "GGB3012 SGB4003"] = "Hungateiclostridiaceae SGB4003"
spN[spN == "GGB753 SGB989"] = "Actinomycetaceae SGB989"
names(dismatList) = spN

# save(dismatList, file = './data/RData/without_NCBI_references.RData')
# combine each samples
comb_dist = data.frame()
for (j in names(dismatList)) {
  cat(j, '\n')
  dismat = dismatList[[j]]
  dt = melt(dismat)
  dt = dt[dt$X1 != dt$X2,]
  dt$Y1 = str_split_fixed(dt$X1, '_', 2)[, 1]
  dt$Y2 = str_split_fixed(dt$X2, '_', 2)[, 1]
  dt$group[dt$Y1 == dt$Y2] = 'intra-participant'
  dt$group[dt$Y1 != dt$Y2] = 'inter-participant'
  print(table(dt$group))
  if(length(table(dt$group)) == 2){
    dt$species = paste0(j, ': ', length(dt$X1 %>% unique()))
    comb_dist = rbind(comb_dist, dt)
  }
}

agg = aggregate(comb_dist[comb_dist$group == 'inter-participant', c("value", "species")], 
                list(comb_dist$species[comb_dist$group == 'inter-participant']), function(x){median(as.numeric(x), na.rm = T)})
comb_dist$species %<>% factor(., levels = agg$Group.1[order(agg$value)])
library(ggforestplot)
stat.test <- comb_dist %>%
  group_by(species) %>%
  wilcox_test(value ~ group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "x", dodge = 0)
# stat.test$p.adj %<>% signif(., 2)
p = ggplot(comb_dist, aes(x = value*100, y  = species, fill = group)) + 
  geom_stripes(odd = "#F5F5F5", even = "white") +
  geom_boxplot(outlier.colour = NA, size = 0.25, width = 0.6) + mytheme + 
  theme(panel.grid= element_blank(), panel.border = element_rect(size = 0.6), 
        legend.position = 'top', legend.title = element_blank(), legend.direction = 'vertical',
        legend.margin = margin(1,1,-2,1,'mm')) +
  labs(x = 'SNP dissimilarity (%)', y = '', fill ='Group') +
  scale_fill_manual(values = c('#FDC086', '#B3CDE3')) +
  scale_x_continuous(expand = expansion(mult = c(0.05,0.1)), breaks = seq(0,100,25), labels = seq(0,100,25)) +
  annotate('text', x = 110, y = stat.test$species, label = stat.test$p.adj.signif)
if(plotFlag){cairo_pdf('./plot/distance_of_inter_and_intra_participant.pdf', width = 5.71, height = 7)}
print(p)
if(plotFlag){dev.off()}

#### Make the dendrogram based on Kimura matrix ####
load('./data/RData/without_NCBI_references.RData')
load('./data/RData/without_NCBI_references2.RData')
if(plotFlag){pdf('./plot/dendrogram_of_Kimura_2.pdf', width = 6, height = 5)}
clusterDF = data.frame()
for (l in names(dismatList)) {
  mat = dismatList[[l]]
  row.names(mat) = swabData$SampleID.u[match(row.names(mat), swabData$pt_ID.u)]
  colnames(mat) = swabData$SampleID.u[match(colnames(mat), swabData$pt_ID.u)]
  mat[is.nan(mat)] = NA
  colSums(is.nan(mat))
  # mat = mat[order(colnames(mat)), order(colnames(mat))]
  # colnames(mat) = row.names(mat) = paste(swabData$pt_ID, swabData$trimester, swabData$Sample_GA, sep = '_')[match(row.names(mat), swabData$pt_ID.u)]
  # hcT = try({
    hc_complete <- hclust(as.dist(mat), method = "complete")
  # })
  # if(class(hcT) == 'try-error'){
  #   cat(l, 'can not do hclust\n')
  #   next
  # }
  hc_complete = dendsort(hc_complete)
  library(dendextend)
  clist = cutree(hc_complete, h = 0.7, order_clusters_as_data = F)
  clusterDF = rbind(clusterDF, data.frame(species = l,
                                samples = names(clist),
                                cluster = clist))
  k = table(clist) %>% length()
  dend = as.dendrogram(hc_complete)
  if(plotFlag){
    par(mfrow = c(1, 1), las=1)
    dend %>% 
      set("branches_k_color", k = k, value = colorStrain[1:k]) %>%
      set("labels_col", colorStrain[1:k], k = k) %>%
      set("labels_cex", 0.2) %>%
      plot(xlim = c(1, length(clist)), main = l, ylab = 'Dissimilarity')
    abline(h = 0.7, lty = 2)
  }
}
if(plotFlag){dev.off()}
save(clusterDF, file = './data/RData/clusterDF.RData')
length(unique(clusterDF$species))
###### careful plot using complex heatmap ####
load('./data/RData/without_NCBI_references.RData')
# load('./data/RData/without_NCBI_references2.RData')
plotFlag = T
if(plotFlag){pdf('./plot/heatmap_with_dend_of_Kimura_2.pdf', width = 20, height = 20)}
for (l in names(dismatList)){
  mat = dismatList[[l]]
  row.names(mat) = swabData$SampleID.u[match(row.names(mat), swabData$pt_ID.u)]
  colnames(mat) = swabData$SampleID.u[match(colnames(mat), swabData$pt_ID.u)]
  colSums(is.nan(mat))
  # mat = mat/100
  # mat = mat[order(colnames(mat)), order(colnames(mat))]
  # colnames(mat) = row.names(mat) = paste(swabData$pt_ID, swabData$trimester, swabData$Sample_GA, sep = '_')[match(row.names(mat), swabData$pt_ID.u)]
  hc_complete <- hclust(as.dist(mat), method = "complete")
  hc_complete = dendsort(hc_complete)
  library(dendextend)
  clist = cutree(hc_complete, h = 0.3, order_clusters_as_data = F)
  dend = as.dendrogram(hc_complete)
  k = table(clist) %>% length()
  dend = color_branches(dend, k = k, h= 0.3,col = colorStrain[1:k])
  
  pS = str_extract(colnames(mat), 'P\\d+')
  colorParticipant2 = colorParticipant
  names(colorParticipant2) = metaData$SampleID[match(names(colorParticipant), metaData$pt_ID)]
  topAnno = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = colorStrain[1:k]),
                                               labels = paste0(1:k, ''), 
                                               labels_gp = gpar(col = "black", fontsize = 10)),
                              Participant = pS,
                              border = T, annotation_name_side = 'left',
                              col = list(Participant = colorParticipant2[names(colorParticipant2) %in% pS]))
  leftAnno = rowAnnotation(foo = anno_block(gp = gpar(fill = colorStrain[1:k]),
                                            labels = paste0(1:k, ''), 
                                            labels_gp = gpar(col = "black", fontsize = 10)),
                           Participant = pS,show_annotation_name = F,
                           col = list(Participant = colorParticipant2[names(colorParticipant2) %in% pS]),
                           show_legend = c(TRUE, FALSE))
  topAnno = re_size(topAnno, annotation_height =unit(c(2,2)*0.035*nrow(mat)+3, "mm"))
  leftAnno = re_size(leftAnno, annotation_width =unit(c(2,2)*0.035*nrow(mat)+3, "mm"))
  if(k < 2){
    hp = Heatmap(1-mat, name = "Similarity", cluster_columns = dendsort(dend), cluster_rows = dendsort(dend),
                 col = brewer.pal(n  =10, name="PuOr"), 
                 # row_split = k, column_split = k,
                 height = nrow(mat)*unit(2, "mm"), width = ncol(mat)*unit(2, "mm"), 
                 row_names_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 5),
                 top_annotation = topAnno, left_annotation = leftAnno,
                 column_title = l)
  }else{
    hp = Heatmap(1-mat, name = "Similarity", cluster_columns = dendsort(dend), cluster_rows = dendsort(dend),
                 col = brewer.pal(n  =10, name="PuOr"), 
                 row_split = k, column_split = k,
                 height = nrow(mat)*unit(2, "mm"), width = ncol(mat)*unit(2, "mm"), 
                 row_names_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 5),
                 top_annotation = topAnno, left_annotation = leftAnno,
                 column_title = l)
  }
  draw(hp, merge_legend = TRUE)
}
dev.off()
#### view strain replacement by species ####
load('./data/RData/clusterDF.RData')
load('./data/RData/colorList.RData')
if(T){ # prepare data
  clusterDF$pt_ID = swabData$SampleID[match(clusterDF$samples, swabData$SampleID.u)]
  unique(clusterDF$pt_ID)
  # for each species,remove person with only on sample present
  clusterDF2 = data.frame()
  for(sp in unique(clusterDF$species)){
    clusterDF_sub = clusterDF[clusterDF$species == sp,]
    tmpN = clusterDF_sub %>% group_by(pt_ID) %>% summarise(N = n())
    clusterDF_sub2 =  clusterDF_sub[clusterDF_sub$pt_ID %in% tmpN$pt_ID[tmpN$N > 1], ]
    clusterDF2 = rbind(clusterDF2, clusterDF_sub2)
  }
  unique(clusterDF2$pt_ID)
  clusterDF2$cluster %<>% as.character()
  # for each species, each person add the additional sample even though we can't call strain from them
  clusterDF3 = data.frame()
  metaData3 = data.frame()
  for(sp in unique(clusterDF$species)){
    clusterDF_sub = clusterDF2[clusterDF2$species == sp,]
    swb = swabData[swabData$SampleID %in% clusterDF_sub$pt_ID, ]
    diffS = setdiff(swb$SampleID.u, clusterDF_sub$samples)
    if(length(diffS) > 0){
      addDF = data.frame(
        species = sp,
        samples = diffS,
        cluster = '0',
        # time_point = str_split_fixed(diffS, '-', 2)[, 2]
        pt_ID = str_split_fixed(diffS, '-', 2)[, 1]
      )
      clusterDF_sub = rbind(clusterDF_sub, addDF)
    }
    metaData_sub = metaData[metaData$SampleID %in% clusterDF_sub$pt_ID, ]
    metaData_sub$species = sp
    clusterDF3 = rbind(clusterDF3, clusterDF_sub)
    metaData3 = rbind(metaData3, metaData_sub)
  }
  unique(clusterDF3$pt_ID)
  clusterDF3$time = swabData$Sample_GA[match(clusterDF3$samples, swabData$SampleID.u)]
  clusterDF3$Term = metaData$Term[match(clusterDF3$pt_ID, metaData$SampleID)]
  clusterDF3$Term_char = metaData$Term_char[match(clusterDF3$pt_ID, metaData$SampleID)]
  clusterDF3$Term_char2 = metaData$Term_char2[match(clusterDF3$pt_ID, metaData$SampleID)]
  clusterDF3$Ethnicity = metaData$Ethnicity[match(clusterDF3$pt_ID, metaData$SampleID)]
  clusterDF3$Ethnicity2 = metaData$Ethnicity2[match(clusterDF3$pt_ID, metaData$SampleID)]
  clusterDF3$trimester = swabData$trimester[match(clusterDF3$samples, swabData$SampleID.u)]
  save(clusterDF3, file = './data/RData/clusterDF3.RData')
  clusterDF_T = clusterDF3[clusterDF3$trimester != 'P', ]
  clusterDF_P = clusterDF3[clusterDF3$trimester == 'P', ]
  # clusterDF_T$pt_ID %<>% factor(., levels = metaData_sub$pt_ID[order(metaData_sub$Term)])
  cairo_pdf('./plot/strain_transition_facet_by_term.pdf', width = 15, height = 40)
  tp = ggplot(data = clusterDF_T) +
    geom_point(data = clusterDF_T, mapping = aes(x = time, y = pt_ID, fill = cluster), alpha = 0.65, size = 3, shape = 21, stroke = 0.2) +
    geom_point(data = clusterDF_P, mapping = aes(x = 45, y = pt_ID, fill = cluster), alpha = 0.65, size = 3, shape = 21, stroke = 0.2) +
    geom_point(data = metaData3, mapping = aes(x = Term, y = pt_ID, color = Term_char2), shape = ')', size = 4) +
    geom_vline(xintercept = 39, linetype = 'dashed', color = 'red', linewidth = 0.2) +
    geom_vline(xintercept = c(14,26), linetype = 'dashed', color = 'black', linewidth = 0.2) +
    scale_fill_manual(values = colorStrain, name = 'Strain', guide = guide_legend(order = 1)) +
    scale_color_manual(values = colorTerm2, name = 'Term', guide = guide_legend(order = 2)) +
    scale_shape_manual(values = c(2, 24), name = 'CST changed', guide = guide_legend(order = 3)) +
    ylab('Participants') + xlab('Gestational Weeks') +
    facet_grid(species~Term_char2, scales = 'free_y', space = 'free') +
    mytheme
  print(tp)
  dev.off()
  cairo_pdf('./plot/strain_transition_facet_by_ethnicity.pdf', width = 15, height = 50)
  tp = ggplot(data = clusterDF_T) +
    geom_point(data = clusterDF_T, mapping = aes(x = time, y = pt_ID, fill = cluster), alpha = 0.65, size = 3, shape = 21, stroke = 0.2) +
    geom_point(data = clusterDF_P, mapping = aes(x = 45, y = pt_ID, fill = cluster), alpha = 0.65, size = 3, shape = 21, stroke = 0.2) +
    geom_point(data = metaData3, mapping = aes(x = Term, y = pt_ID, color = Term_char2), shape = ')', size = 4) +
    geom_vline(xintercept = 39, linetype = 'dashed', color = 'red', linewidth = 0.2) +
    geom_vline(xintercept = c(14,26), linetype = 'dashed', color = 'black', linewidth = 0.2) +
    scale_fill_manual(values = colorStrain, name = 'Strain', guide = guide_legend(order = 1)) +
    scale_color_manual(values = colorTerm2, name = 'Term', guide = guide_legend(order = 2)) +
    scale_shape_manual(values = c(2, 24), name = 'CST changed', guide = guide_legend(order = 3)) +
    ylab('Participants') + xlab('Gestational Weeks') +
    facet_grid(species~Ethnicity2, scales = 'free_y', space = 'free') +
    mytheme
  print(tp)
  dev.off()
}
###### select species ######
# for term
unique(clusterDF3$species) %>% sort()
selSp = c(
  "Aerococcus christensenii",
  "Fannyhessea vaginae (SGB990)",
  "Fannyhessea vaginae (SGB991)",
  "Gardnerella vaginalis (SGB17302)",  
  "Gardnerella vaginalis (SGB21500)",
  "GGB3012 SGB4003",
  "Lactobacillus iners",
  "Megasphaera genomosp type 1" 
)

for(sp in selSp){
  clusterDF_sub = clusterDF3[clusterDF3$species == sp & clusterDF3$trimester != 'P', ]
  clusterDF_subP = clusterDF3[clusterDF3$species == sp & clusterDF3$trimester == 'P', ]
  metaData_sub = metaData[metaData$pt_ID %in% clusterDF_sub$pt_ID, ]
  clusterDF_sub$pt_ID %<>% factor(., levels = metaData_sub$pt_ID[order(metaData_sub$Term)])
  pp = ggplot(data = clusterDF_sub) +
    geom_point(data = clusterDF_sub, mapping = aes(x = time, y = pt_ID, fill = cluster), alpha = 0.65, size = 3, shape = 21, stroke = 0.2) +
    geom_point(data = clusterDF_subP, mapping = aes(x = 45, y = pt_ID, fill = cluster), alpha = 0.65, size = 3, shape = 21, stroke = 0.2) +
    geom_point(data = metaData_sub, mapping = aes(x = Term, y = pt_ID, color = Term_char2), shape = ')', size = 4) +
    geom_vline(xintercept = 39, linetype = 'dashed', color = 'red', linewidth = 0.2) +
    geom_vline(xintercept = c(14,26), linetype = 'dashed', color = 'black', linewidth = 0.2) +
    scale_fill_manual(values = colorStrain, name = 'Strain', guide = guide_legend(order = 1)) +
    scale_color_manual(values = colorTerm2, name = 'Term', guide = guide_legend(order = 2)) +
    scale_shape_manual(values = c(2, 24), name = 'CST changed', guide = guide_legend(order = 3)) +
    ylab('Participants') + xlab('Gestational Weeks') +
    facet_grid(Term_char2~species, scales = 'free_y', space = 'free') + 
    mytheme + theme(strip.text = element_text(size = 12.5), strip.text.y.right = element_text(angle = 90))
  cairo_pdf(paste0('./plot/strain_selected_for_term_', sp,'.pdf'), width = 6.78, height = 0.2*length(unique(clusterDF_sub$pt_ID))+1, onefile = T)
  print(pp)
  dev.off()
}

# for ethnicity
selSp = c(
  "Anaerococcus tetradius",
  "Fannyhessea vaginae (SGB991)",
  "Gardnerella vaginalis (SGB17301)",
  "Gardnerella vaginalis (SGB17302)",  
  "Gardnerella vaginalis (SGB17307)",  
  "Hungateiclostridiaceae SGB4003",
  "Lactobacillus iners",
  "Megasphaera genomosp type 1",
  'Lactobacillus crispatus'
)
for(sp in selSp){
  clusterDF_sub = clusterDF3[clusterDF3$species == sp & clusterDF3$trimester != 'P', ]
  clusterDF_subP = clusterDF3[clusterDF3$species == sp & clusterDF3$trimester == 'P', ]
  metaData_sub = metaData[metaData$SampleID %in% clusterDF_sub$pt_ID, ]
  clusterDF_sub$pt_ID %<>% factor(., levels = metaData_sub$SampleID[order(metaData_sub$Term)])
  pp = ggplot(data = clusterDF_sub) +
    geom_point(data = clusterDF_sub, mapping = aes(x = time, y = pt_ID, fill = cluster), alpha = 0.65, size = 3, shape = 21, stroke = 0.2) +
    geom_point(data = clusterDF_subP, mapping = aes(x = 45, y = pt_ID, fill = cluster), alpha = 0.65, size = 3, shape = 21, stroke = 0.2) +
    geom_point(data = metaData_sub, mapping = aes(x = Term, y = SampleID, color = Term_char2), shape = ')', size = 4) +
    geom_vline(xintercept = 39, linetype = 'dashed', color = 'red', linewidth = 0.2) +
    geom_vline(xintercept = c(14,26), linetype = 'dashed', color = 'black', linewidth = 0.2) +
    scale_fill_manual(values = colorStrain, name = 'Strain', guide = guide_legend(order = 1)) +
    scale_color_manual(values = colorTerm2, name = 'Term', guide = guide_legend(order = 2)) +
    scale_shape_manual(values = c(2, 24), name = 'CST changed', guide = guide_legend(order = 3)) +
    ylab('Participants') + xlab('Gestational Weeks') +
    facet_grid(Ethnicity2~species, scales = 'free_y', space = 'free') + 
    mytheme + theme(strip.text = element_text(size = 12.5), strip.text.y.right = element_text(angle = 90))
  cairo_pdf(paste0('./plot/strain_selected_for_ethnicity_', sp,'.pdf'), width = 6.78, height = 0.2*length(unique(clusterDF_sub$pt_ID))+1, onefile = T)
  print(pp)
  dev.off()
}
# for subspecies
selSp = c(
  "Fannyhessea vaginae (SGB990)",
  "Fannyhessea vaginae (SGB991)",
  "Gardnerella vaginalis (SGB17301)",
  "Gardnerella vaginalis (SGB17301)",
  "Gardnerella vaginalis (SGB17302)",  
  "Gardnerella vaginalis (SGB17305)",
  "Gardnerella vaginalis (SGB17307)",   
  "Gardnerella vaginalis (SGB21500)",
  "Gardnerella vaginalis (SGB7097)"
)
for(sp in selSp){
  clusterDF_sub = clusterDF3[clusterDF3$species == sp & clusterDF3$trimester != 'P', ]
  clusterDF_subP = clusterDF3[clusterDF3$species == sp & clusterDF3$trimester == 'P', ]
  metaData_sub = metaData[metaData$pt_ID %in% clusterDF_sub$pt_ID, ]
  clusterDF_sub$pt_ID %<>% factor(., levels = metaData_sub$pt_ID[order(metaData_sub$Term)])
  pp = ggplot(data = clusterDF_sub) +
    geom_point(data = clusterDF_sub, mapping = aes(x = time, y = pt_ID, fill = cluster), alpha = 0.65, size = 3, shape = 21, stroke = 0.2) +
    geom_point(data = clusterDF_subP, mapping = aes(x = 45, y = pt_ID, fill = cluster), alpha = 0.65, size = 3, shape = 21, stroke = 0.2) +
    geom_point(data = metaData_sub, mapping = aes(x = Term, y = pt_ID, color = Term_char2), shape = ')', size = 4) +
    geom_vline(xintercept = 39, linetype = 'dashed', color = 'red', linewidth = 0.2) +
    geom_vline(xintercept = c(14,26), linetype = 'dashed', color = 'black', linewidth = 0.2) +
    scale_fill_manual(values = colorStrain, name = 'Strain', guide = guide_legend(order = 1)) +
    scale_color_manual(values = colorTerm2, name = 'Term', guide = guide_legend(order = 2)) +
    scale_shape_manual(values = c(2, 24), name = 'CST changed', guide = guide_legend(order = 3)) +
    ylab('Participants') + xlab('Gestational Weeks') + labs(title = sp) +
    # facet_grid(Ethnicity2~species, scales = 'free_y', space = 'free') + 
    mytheme + theme(strip.text = element_text(size = 12.5), strip.text.y.right = element_text(angle = 90))
  cairo_pdf(paste0('./plot/strain_selected_for_subspecies_', sp,'.pdf'), width = 6.78, height = 0.2*length(unique(clusterDF_sub$pt_ID))+1, onefile = T)
  print(pp)
  dev.off()
}
#### view strain replacement by heatmap only replacement ####
clusterDF = data.frame()
for (l in names(dismatList)) {
  mat = dismatList[[l]]
  mat = mat[!grepl('KIT', colnames(mat)), !grepl('KIT', colnames(mat))]
  mat = mat[order(colnames(mat)), order(colnames(mat))]
  # colnames(mat) = row.names(mat) = paste(swabData$pt_ID, swabData$trimester, swabData$Sample_GA, sep = '_')[match(row.names(mat), swabData$pt_ID)]
  hc_complete <- hclust(as.dist(mat), method = "complete")
  hc_complete = dendsort(hc_complete)
  clist = cutree(hc_complete, h = 0.7,order_clusters_as_data = F)
  clusterDF = rbind(clusterDF, data.frame(species = l,
                                          samples = names(clist),
                                          cluster = clist))
}
clusterDF = clusterDF[(!grepl('kit', clusterDF$samples, ignore.case = T)) & clusterDF$samples != 'Reference', ]
clusterDF$time_point = str_split_fixed(clusterDF$samples, '_', 2)[, 2]
clusterDF$time_point = factor(clusterDF$time_point, 
                              levels = table(clusterDF$time_point) %>% names %>% gtools::mixedsort())
clusterDF$person = str_split_fixed(clusterDF$samples, '_', 2)[, 1]
sp = unique(clusterDF$species)[1]
clusterDF$week = swabData$Sample_GA[match(clusterDF$samples, swabData$pt_ID.u)]
clusterDF$Trimester = swabData$trimester[match(clusterDF$samples, swabData$pt_ID.u)]
clusterDF$Trimester %<>% as.character()
clusterDF$Trimester %<>% as.factor()
clusterDF_Repl = data.frame()
for(sp in unique(clusterDF$species)){
  subDF = clusterDF[clusterDF$species == sp, ]
  for(ps in unique(subDF$person)){
    ssubDF = subDF[subDF$person == ps, ]
    if(length(table(ssubDF$cluster[ssubDF$cluster != 0])) > 1){ 
      adduS = swabData$pt_ID.u[(swabData$pt_ID == ps) & (!swabData$pt_ID.u %in% ssubDF$samples)]
      if(length(adduS)>0){
        ssubDF = rbind(ssubDF, data.frame(
          species = sp,
          samples = adduS,
          cluster = 0,
          time_point = str_split_fixed(adduS, '_', 2)[, 2],
          person = str_split_fixed(adduS, '_', 2)[, 1],
          week = swabData$Sample_GA[match(adduS, swabData$pt_ID.u)],
          Trimester = swabData$trimester[match(adduS, swabData$pt_ID.u)]
        ))
      }
      ssubDF = ssubDF[ssubDF$Trimester != 'P', ]
      clusterDF_Repl = rbind(clusterDF_Repl, ssubDF)
    }
  }
}
clusterDF_Repl$d = paste0(clusterDF_Repl$Trimester, '\n(', clusterDF_Repl$week, ')')
clusterDF_Repl = clusterDF_Repl[clusterDF_Repl$species %in% c('Atopobium_vaginae', 
                                                              'Aerococcus_christensenii', 
                                                              'Gardnerella_vaginalis', 
                                                              'Lactobacillus_crispatus', 
                                                              'Lactobacillus_iners'), ]
clusterDF_Repl$cluster %<>% as.factor()
clusterDF_Repl$person %<>% factor(., levels = rev(sort(unique(clusterDF_Repl$person))))
if(plotFlag){pdf('./plot/strain_replacement_point_new_2.pdf', width = 12, height = 12)}
p = ggplot(clusterDF_Repl, aes(x = week, y = person, species, fill = cluster), 
       color = 'black') + 
  geom_point(shape = 21, size = 4,stroke = 0.2) +
  geom_vline(xintercept = c(14, 26), linetype= "dashed", 
             color = "blue", size= 0.3) +
  scale_fill_manual(values = clusterCol %>% alpha(.,1.7)) + 
  labs(x = 'Gestational Weeks', fill = 'Strain', y = 'Person') +
  facet_grid(rows = vars(species), scales="free_y") +
  mytheme +
  theme(strip.text.y = element_text(size = 12))
print(p)
dev.off()

if(plotFlag){pdf('./plot/strain_replacement_heatmap_label_only_no_replace_3_1.pdf', width = 12, height = 25)}
if(plotFlag){
  subMat = dcast(clusterDF_Repl, person + species ~ time_point, value.var = 'cluster')
  subMat = subMat[order(subMat$species, subMat$person), ]
  subMat2 = dcast(clusterDF_Repl, person + species ~ time_point, value.var = 'd')
  subMat2 = subMat2[order(subMat2$species, subMat2$person), ]
  mat = subMat[, 3:ncol(subMat)] %>% as.matrix()
  row.names(mat)=paste0(subMat$species, '_',subMat$person )
  mat2 = subMat2[, 3:ncol(subMat2)] %>% as.matrix()
  row.names(mat2)=paste0(subMat2$species, '_',subMat2$person )
  uniC = unique(c(mat) %>% na.omit()) %>% gtools::mixedsort()
  colors = clusterCol[match(uniC, names(clusterCol))]
  load('./data/col_vector.RData')
  
  # haCol = brewer.pal(11, 'Paired') %>% sample(11)
  haCol = col_vector[match(unique(subMat$species), names(col_vector))]
  names(haCol) = unique(subMat$species)
  ha = rowAnnotation(labels = anno_text(subMat$person, which = "row"),
                     'Ethnicity' = metaData$Ethnicity[match(subMat$person, metaData$pt_ID)],
                     'Term' = metaData$Term_char[match(subMat$person, metaData$pt_ID)],
                     'Species' = subMat$species, 
                     col = list('Species' = haCol,
                                 'Ethnicity' = c('White' = '#96A87F' %>% lighten(., factor = 1.3), 
                                                 'Asian' = '#CDA245' %>% alpha(., alpha = 0.6), 
                                                 'Latina' = '#6B988B' %>% lighten(., factor = 1.3), 
                                                 'Black' = '#847879' %>% lighten(., factor = 1.3), 
                                                 'Pacific Islander' = '#C17E96' %>% lighten(., factor = 1.3), 
                                                 'Others' = '#B3B3B3' %>% lighten(., factor = 1.2)),
                                 'Term' = c('PTB' = '#2E3338' %>%  lighten(., factor = 1.9), 
                                            'ETB' = '#467283' %>%  lighten(., factor = 1.4),
                                            'TB' = '#C9DAD6' %>%  lighten(., factor = 1)), 
                                  'Replacement' = c('No' = '#D1D1E7', 'Yes' = '#8B87BE')))
  ph = Heatmap(mat, name = paste0("Strain cluster"), cluster_rows = F, cluster_columns = F,
               column_names_side = "top", row_names_side = "left", column_title_side = 'top', 
               row_names_centered = F, column_names_centered = T, show_row_names = F,
               rect_gp = gpar(col = "white", lwd = 2), na_col = brewer.pal(11, 'RdBu')[6],
               column_names_rot = T, col = colors, 
               width = ncol(mat)*unit(7, "mm"), 
               height = nrow(mat)*unit(7, "mm"),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if(!is.na(mat2[i, j])){
                   grid.text(sprintf(mat2[i, j]), x, y, gp = gpar(fontsize = 5, col = 'black', alpha = 0.7))
                 }
               },
               row_title = "Participant", column_title = "Sampling order", 
               left_annotation = ha
  )
  draw(ph, heatmap_legend_side = "right", annotation_legend_side = "bottom", merge_legend = TRUE)
  
  # sort by person
  subMat = subMat[order(subMat$person, subMat$species), ]
  subMat2 = subMat2[order(subMat2$person, subMat2$species), ]
  mat = subMat[, 3:ncol(subMat)] %>% as.matrix()
  row.names(mat)=paste0(subMat$person, '_',subMat$species)
  mat2 = subMat2[, 3:ncol(subMat2)] %>% as.matrix()
  row.names(mat2)=paste0(subMat2$person, '_',subMat2$species)
  uniC = unique(c(mat) %>% na.omit()) %>% gtools::mixedsort()
  colors = clusterCol[match(uniC, names(clusterCol))]
  ha = rowAnnotation(labels = anno_text(subMat2$person, which = "row"),
                     'Ethnicity' = metaData$Ethnicity[match(subMat$person, metaData$pt_ID)],
                     'Term' = metaData$Term_char[match(subMat$person, metaData$pt_ID)],
                     'Species' = subMat$species, col = list('Species' = haCol,
                                                            'Ethnicity' = c('White' = '#96A87F' %>% lighten(., factor = 1.3), 
                                                                            'Asian' = '#CDA245' %>% alpha(., alpha = 0.6), 
                                                                            'Latina' = '#6B988B' %>% lighten(., factor = 1.3), 
                                                                            'Black' = '#847879' %>% lighten(., factor = 1.3), 
                                                                            'Pacific Islander' = '#C17E96' %>% lighten(., factor = 1.3), 
                                                                            'Others' = '#B3B3B3' %>% lighten(., factor = 1.2)),
                                                            'Term' = c('PTB' = '#2E3338' %>%  lighten(., factor = 1.9), 
                                                                       'ETB' = '#467283' %>%  lighten(., factor = 1.4),
                                                                       'TB' = '#C9DAD6' %>%  lighten(., factor = 1)))
                     )
  ph = Heatmap(mat, name = paste0("Strain cluster"), cluster_rows = F, cluster_columns = F,
               column_names_side = "top", row_names_side = "left", column_title_side = 'top', 
               row_names_centered = F, column_names_centered = T,show_row_names = F,
               rect_gp = gpar(col = "white", lwd = 2), na_col = brewer.pal(11, 'RdBu')[6],
               column_names_rot = T, col = colors, 
               width = ncol(mat)*unit(7, "mm"), 
               height = nrow(mat)*unit(7, "mm"),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if(!is.na(mat2[i, j])){
                   grid.text(sprintf(mat2[i, j]), x, y, gp = gpar(fontsize = 5, col = 'black', alpha = 0.7))
                 }
               },
               row_title = "Participant", column_title = "Sampling order", 
               left_annotation = ha
  )
  draw(ph, heatmap_legend_side = "right", annotation_legend_side = "bottom", merge_legend = TRUE)
}
if(plotFlag){dev.off()}

