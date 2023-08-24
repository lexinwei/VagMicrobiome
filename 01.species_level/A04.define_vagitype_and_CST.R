source('./script/library_function.R')
load('./data/RData/abuList.RData')
load('./data/RData/swabData.RData')
load('./data/RData/metaData.RData')
load('./data/RData/colorList.RData')
load('./data/RData/sampleVagitypes.RData')
if(T){ # prepare abundance matrix
  abuS = abuList[['s']]
  # discard NC samples
  removeSamples = swabData$seq_ID[grepl('^KIT', swabData$pt_ID)]
  length(removeSamples)
  abuS1 = abuS[, !colnames(abuS) %in% removeSamples]
  range(rowSums(abuS1[, 4:ncol(abuS1)])) 
  # discard unclassified
  abuS2 = abuS1[abuS1$taxa != 'unclassified', ]
  abuS_mat = abuS2[, 4:ncol(abuS2)]
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
#### vagitypes defined by Fettweis et al, Nature Medicine paper ####
if(T){
  sampleVtype = c()
  for(i in 4:ncol(abuS2)){
    maxInd = order(abuS2[, i], decreasing = T)[1]
    vType = abuS2$taxa[maxInd]
    if(abuS2[maxInd, i] > 30 & vType != 'unclassified'){
      sampleVtype = c(sampleVtype, vType)
    }else{
      sampleVtype = c(sampleVtype, 'None dominant')
    }
  }
  sampleVagitypes = data.frame(seq_ID = colnames(abuS2[, 4:ncol(abuS2)]), vagitypes_Fettweis = sampleVtype)
  sampleVagitypes %<>% add_column(pt_ID.u = swabData$pt_ID.u[match(sampleVagitypes$seq_ID, swabData$seq_ID)],
                                  .before = 'seq_ID')
  sampleVagitypes$vagitypes_Fettweis %<>% as.character()
  sampleVagitypes$vagitypes_Fettweis %<>% str_replace_all(., '_', ' ')
  table(sampleVagitypes$vagitypes_Fettweis) %>% sort()
  sampleVagitypes$vagitypes_Fettweis2 = sampleVagitypes$vagitypes_Fettweis %>% as.character()
  sampleVagitypes$vagitypes_Fettweis2[sampleVagitypes$vagitypes_Fettweis %in% names(table(sampleVagitypes$vagitypes_Fettweis)[table(sampleVagitypes$vagitypes_Fettweis)<10])] = 'Other'
  sampleVagitypes$vagitypes_Fettweis %<>% factor(., levels = c(names(colorVagitype)[1:7], 
                                                               'Lactobacillus vaginalis', 'Megasphaera genomosp type 1',
                                                               'Lacticaseibacillus rhamnosus', 'Mobiluncus curtisii',
                                                               'Enterococcus faecalis', 'Streptococcus anginosus', 
                                                               'Streptococcus agalactiae', 'Candida albicans', 
                                                               'Lawsonella SGB3665', 'Hungateiclostridiaceae SGB4003', 'Actinomycetaceae SGB989', 'None dominant'))
  sampleVagitypes$vagitypes_Fettweis2 %<>% str_replace_all(., '_', ' ')
  sampleVagitypes$vagitypes_Fettweis2 %<>% factor(., levels = names(colorVagitype))
  mytable = table(sampleVagitypes$vagitypes_Fettweis2)
  lbls <- paste(names(mytable), "\n", mytable, sep="")
  pie(mytable, labels = lbls, col = colorVagitype[names(mytable)],
      main = "Pie Chart of Vagitypes")

}

###### Susan's PNAS paper ####
library(phyloseq)
if(T){
  dim(abuMat_T)
  bray_curtis_dist <- vegan::vegdist(abuMat_T, method = "bray")
  pcoa <- cmdscale(bray_curtis_dist, k = (nrow(abuMat) - 1), eig = TRUE)
  points = as.data.frame(pcoa$points)
  eig = pcoa$eig
  print(eig[1:50])
  print(tail(eig))
  h_sub5 <- hist(eig[6:length(eig)], 100)
  plot(h_sub5$mids, h_sub5$count, log="y", type='h', lwd=10, lend=2)
  pamPCoA = function(x, k) {
    list(cluster = pam(x[,1:2], k, cluster.only = TRUE))
  }
  NDIM = 7
  gs = clusGap(points[, 1:NDIM], FUN = pamPCoA, K.max = 12, B = 50)
  plot_clusgap(gs) + scale_x_continuous(breaks=c(seq(0, 12, 2)))
  
  K <- 5 # 3,7
  x <- points[, 1:NDIM]
  clust <- as.factor(pam(x, k=K, cluster.only=T))
  table(clust)
  table(clust, sampleVagitypes$vagitypes_Fettweis[match(names(clust), sampleVagitypes$seq_ID)])
  table(clust, sampleVagitypes$vagitypes_Fettweis2[match(names(clust), sampleVagitypes$seq_ID)])
  vagitypes_Susan = as.character(clust)
  names(vagitypes_Susan) = names(clust)
  vagitypes_Susan[clust == '1'] = 'I'
  vagitypes_Susan[clust == '5'] = 'II'
  vagitypes_Susan[clust == '3'] = 'III'
  vagitypes_Susan[clust == '4'] = 'IV'
  vagitypes_Susan[clust == '2'] = 'V'
  sampleVagitypes$vagitypes_Susan = vagitypes_Susan[match(sampleVagitypes$seq_ID,names(clust))]
  # get Median Bray-Curtis dissimilarity
  bcDissimilarity = as.matrix(bray_curtis_dist)
  bcdmp_mean = bcdmp_median = c()
  ptList = unique(swabData$pt_ID[!grepl('^KIT', swabData$pt_ID)])
  for(p in ptList){
    ptL = swabData$seq_ID[swabData$pt_ID == p & swabData$trimester != "P"]
    bcSub = bcDissimilarity[match(ptL, row.names(bcDissimilarity)), match(ptL, colnames(bcDissimilarity))]
    bcdmp_mean = c(bcdmp_mean, bcSub[upper.tri(bcSub)] %>% mean())
    bcdmp_median = c(bcdmp_median, bcSub[upper.tri(bcSub)] %>% median())
  }
  bcDissimilarityMedianPerson = data.frame(pt_ID = ptList,
                                           bcdmp_mean = bcdmp_mean,
                                           bcdmp_median = bcdmp_median)
  save(bcDissimilarityMedianPerson, file = './data/RData/bcDissimilarityMedianPerson.RData')
  sampleVagitypes$vagitypes_Susan %<>% factor()
  save(sampleVagitypes, file = './data/RData/sampleVagitypes.RData')
  

}

###### plot CST dynamics ####
load('./data/RData/sampleVagitypes.RData')
load('./data/RData/colorList.RData')
if(T){ # prepare data
  sampleVagitypes$CST = sampleVagitypes$vagitypes_Susan
  sampleVagitypes$Vagitypes = sampleVagitypes$vagitypes_Fettweis2
  sampleVagitypes$time = swabData$Sample_GA[match(sampleVagitypes$pt_ID.u, swabData$pt_ID.u)]
  sampleVagitypes$pt_ID = swabData$pt_ID[match(sampleVagitypes$pt_ID.u, swabData$pt_ID.u)]
  sampleVagitypes$Term = metaData$Term[match(sampleVagitypes$pt_ID, metaData$pt_ID)]
  sampleVagitypes$Term_char = metaData$Term_char[match(sampleVagitypes$pt_ID, metaData$pt_ID)]
  sampleVagitypes$Term_char2 = metaData$Term_char2[match(sampleVagitypes$pt_ID, metaData$pt_ID)]
  sampleVagitypes$Ethnicity = metaData$Ethnicity[match(sampleVagitypes$pt_ID, metaData$pt_ID)]
  sampleVagitypes$Ethnicity2 = metaData$Ethnicity2[match(sampleVagitypes$pt_ID, metaData$pt_ID)]
  sampleVagitypes$trimester = swabData$trimester[match(sampleVagitypes$pt_ID.u, swabData$pt_ID.u)]
  ddtt = sampleVagitypes[sampleVagitypes$trimester != "P", ]
  ddtt = ddtt[!is.na(ddtt$pt_ID),]
  range(ddtt$bcMedian)
  unique(ddtt$pt_ID)
  ddtt = ddtt[naturalorder(ddtt$pt_ID.u), ]
  CST_chs = Vagitype_chs = CST_first = Vagitype_first = c()
  for (p in unique(ddtt$pt_ID)){
    CST_first = c(CST_first, as.character(ddtt$CST)[ddtt$pt_ID == p][1])
    Vagitype_first = c(Vagitype_first, as.character(ddtt$Vagitypes)[ddtt$pt_ID == p][1])
    if(length(unique(ddtt$CST[ddtt$pt_ID == p])) > 1){CST_chs = c(CST_chs, '1')}else{CST_chs = c(CST_chs, '0')}
    if(length(unique(ddtt$Vagitypes[ddtt$pt_ID == p])) > 1){Vagitype_chs = c(Vagitype_chs, '1')}else{Vagitype_chs = c(Vagitype_chs, '0')}
  }
  ddttChs = data.frame(pt_ID = unique(ddtt$pt_ID), 
                       CST_changes = CST_chs, Vagitype_changes = Vagitype_chs,
                       CST_first = CST_first, 
                       Vagitype_first = Vagitype_first)
  ddttChs$Ethnicity = metaData$Ethnicity2[match(ddttChs$pt_ID, metaData$pt_ID)]
  ddttChs$Term = metaData$Term_char2[match(ddttChs$pt_ID, metaData$pt_ID)]
  table(ddttChs$CST_changes, ddttChs$Ethnicity)
  chisq.test(table(ddttChs$CST_changes, ddttChs$Ethnicity))
  
  table(ddttChs$Vagitype_changes, ddttChs$Ethnicity)
  chisq.test(table(ddttChs$Vagitype_changes, ddttChs$Ethnicity))
  
  table(ddttChs$CST_changes, ddttChs$Term)
  chisq.test(table(ddttChs$CST_changes, ddttChs$Term))
  
  table(ddttChs$Vagitype_changes, ddttChs$Term)
  chisq.test(table(ddttChs$Vagitype_changes, ddttChs$Term))
}
save(sampleVagitypes, file = './data/RData/sampleVagitypes.RData')
if(T){ # plot order by the time of the last samples
  ddtt$pt_ID = factor(ddtt$pt_ID, levels = unique(ddtt$pt_ID[order(ddtt$Term, decreasing = F)]))
  p1 = ggplot() +
    geom_point(data = ddtt, mapping = aes(x = time, y = pt_ID, fill = CST), alpha = 0.65, size = 3, shape = 21, stroke = 0.1) +
    geom_point(data = metaData, mapping = aes(x = Term, y = pt_ID, color = Term_char2), shape = ')', size = 4) +
    geom_point(data = ddttChs, mapping = aes(x = 42, y = pt_ID, shape = CST_changes), size = 3,fill = '#BF5B16') +
    geom_vline(xintercept = 39, linetype = 'dashed', color = 'red', size = 0.3) +
    geom_vline(xintercept = c(14,26), linetype = 'dashed', color = 'black', size = 0.3) +
    scale_fill_manual(values = colorCST, name = 'CST', guide = guide_legend(order = 1)) +
    scale_color_manual(values = colorTerm2, name = 'Term', guide = guide_legend(order = 2)) +
    scale_shape_manual(values = c(2, 24), name = 'CST changed', guide = guide_legend(order = 3)) +
    ylab('Participants') + xlab('Gestational Weeks') +
    mytheme + theme(aspect.ratio = 1.5/1)
  p2 = ggplot() +
    geom_point(data = ddtt, mapping = aes(x = time, y = pt_ID, fill = Vagitypes), alpha = 0.65, size = 3, shape = 21, stroke = 0.1) +
    geom_point(data = metaData, mapping = aes(x = Term, y = pt_ID, color = Term_char2), shape = ')', size = 4) +
    geom_point(data = ddttChs, mapping = aes(x = 42, y = pt_ID, shape = Vagitype_changes), size = 3,fill = '#BF5B16') +
    geom_vline(xintercept = 39, linetype = 'dashed', color = 'red', size = 0.3) +
    geom_vline(xintercept = c(14,26), linetype = 'dashed', color = 'black', size = 0.3) +
    scale_fill_manual(values = colorVagitype, name = 'Vagitype', guide = guide_legend(order = 1)) +
    scale_color_manual(values = colorTerm2, name = 'Term', guide = guide_legend(order = 2)) +
    scale_shape_manual(values = c(2, 24), name = 'Vagitype changed', guide = guide_legend(order = 3)) +
    ylab('Participants') + xlab('Gestational Weeks') +
    mytheme + theme(aspect.ratio = 1.5/1)
  cairo_pdf(file = './plot/CST_during_pergnant_order_by_term_2.pdf', width = 7.86, height = 7.4, onefile = T)
  print(p1)
  print(p2)
  dev.off()
}
if(T){ # plot order by the time of the median bray curtis dissimilarity
  load('data/RData/bcDissimilarityMedianPerson.RData')
  ddtt$pt_ID = factor(ddtt$pt_ID, levels = unique(bcDissimilarityMedianPerson$pt_ID[order(bcDissimilarityMedianPerson$bcdmp_median, decreasing = T)]))
  ddttChs$Vagitype_first = factor(ddttChs$Vagitype_first, levels = names(colorVagitype))
  metaData$Vagitype_first = ddttChs$Vagitype_first[match(metaData$pt_ID, ddttChs$pt_ID)]
  ddtt$Vagitype_first = ddttChs$Vagitype_first[match(ddtt$pt_ID, ddttChs$pt_ID)]
  p1 = ggplot() +
    geom_point(data = ddtt, mapping = aes(x = time, y = pt_ID, fill = CST), alpha = 0.65, size = 3, shape = 21, stroke = 0.1) +
    geom_point(data = metaData, mapping = aes(x = Term, y = pt_ID, color = Term_char2), shape = ')', size = 4) +
    geom_point(data = ddttChs, mapping = aes(x = 42, y = pt_ID, shape = CST_changes), size = 3,fill = '#BF5B16') +
    geom_vline(xintercept = 39, linetype = 'dashed', color = 'red', size = 0.3) +
    geom_vline(xintercept = c(14,26), linetype = 'dashed', color = 'black', size = 0.3) +
    scale_fill_manual(values = colorCST, name = 'CST', guide = guide_legend(order = 1)) +
    scale_color_manual(values = colorTerm2, name = 'Term', guide = guide_legend(order = 2)) +
    scale_shape_manual(values = c(2, 24), name = 'CST changed', guide = guide_legend(order = 3)) +
    ylab('Participants') +
    xlab('Gestational Weeks') +
    mytheme + 
    theme(strip.text.y = element_blank(), legend.position = 'left',
          plot.margin = unit(c(5, 0, 5, 5), 'mm'))
  p2 = ggplot() +
    geom_point(data = ddtt, mapping = aes(x = time, y = pt_ID, fill = Vagitypes), alpha = 0.65, size = 3, shape = 21, stroke = 0.1) +
    geom_point(data = metaData, mapping = aes(x = Term, y = pt_ID, color = Term_char2), shape = ')', size = 4) +
    geom_point(data = ddttChs, mapping = aes(x = 42, y = pt_ID, shape = Vagitype_changes), size = 3,fill = '#BF5B16') +
    geom_vline(xintercept = 39, linetype = 'dashed', color = 'red', size = 0.3) +
    geom_vline(xintercept = c(14,26), linetype = 'dashed', color = 'black', size = 0.3) +
    scale_fill_manual(values = colorVagitype, name = 'Vagitype', guide = guide_legend(order = 1)) +
    scale_color_manual(values = colorTerm2, name = 'Term', guide = guide_legend(order = 2)) +
    scale_shape_manual(values = c(2, 24), name = 'Vagitype changed', guide = guide_legend(order = 3)) +
    ylab('Participants') +
    xlab('Gestational Weeks') +
    mytheme +
    theme(strip.text.y = element_blank(), legend.position = 'left',
          plot.margin = unit(c(5, 0, 5, 5), 'mm'))
  ddttChs$bcdmp_median = bcDissimilarityMedianPerson$bcdmp_median[match(ddttChs$pt_ID, ddttChs$pt_ID)]
  ddttChs$pt_ID = factor(ddttChs$pt_ID, levels = unique(ddttChs$pt_ID[order(ddttChs$bcdmp_median, decreasing = T)]))
  sy = ggplot(ddttChs,aes(x=1,y=pt_ID,fill=bcdmp_median))+
    geom_tile() + ylab('y') +
    # scale_y_discrete(expand=c(0,0))+
    scale_x_continuous(expand=c(0,0), breaks = 1, label = '1')+
    scale_fill_gradient(low = "#132B43", high = "#56B1F7", breaks=seq(0, 1, 0.25), 
                        limits = c(0, 1), name = 'Within-individual median\nBC dissimilarity') +
    mytheme + theme(legend.position="right", strip.text.y = element_blank(), 
                    axis.title.y=element_blank(), axis.text.y=element_blank(),
                    axis.ticks.y=element_blank(), axis.text.x=element_text(color="white"),
                    axis.title.x=element_text(color="white"), 
                    plot.margin = unit(c(5, 5, 5, 0), 'mm'),
                    axis.ticks.x=element_line(color="white"), panel.grid = element_blank())
  #Define layout for the plots (2 rows, 2 columns)
  layt<-grid.layout(nrow=1,ncol=2,heights=c(7/8, 7/8),widths=c(5.1/8, 1.9/8),default.units=c('null','null'))
  #View the layout of plots
  grid.show.layout(layt)
  #Draw plots one by one in their positions
  cairo_pdf(file = './plot/CST_during_pergnant_order_by_dissimilarity_2.pdf',
            width = 10.26, height = 7.95, onefile = T)
  grid.newpage()
  pushViewport(viewport(layout=layt))
  print(p1,vp=viewport(layout.pos.row=1,layout.pos.col=1))
  print(sy,vp=viewport(layout.pos.row=1,layout.pos.col=2))
  grid.newpage()
  pushViewport(viewport(layout=layt))
  print(p2,vp=viewport(layout.pos.row=1,layout.pos.col=1))
  print(sy,vp=viewport(layout.pos.row=1,layout.pos.col=2))
  dev.off()
}
if(T){ # plot order by the time of the median bray curtis dissimilarity add heatmap
  bcDissimilarityMedianPerson
  bcDissimilarityMedianPerson$SampleID =metaData$SampleID[match(bcDissimilarityMedianPerson$pt_ID, metaData$pt_ID)]
  ddtt$SampleID = metaData$SampleID[match(ddtt$pt_ID, metaData$pt_ID)]
  ddtt$SampleID = factor(ddtt$SampleID, levels = unique(bcDissimilarityMedianPerson$SampleID[order(bcDissimilarityMedianPerson$bcdmp_median, decreasing = T)]))
  ddttChs$SampleID =  metaData$SampleID[match(ddttChs$pt_ID, metaData$pt_ID)]
  ddttChs$Vagitype_first = factor(ddttChs$Vagitype_first, levels = names(colorVagitype))
  metaData$Vagitype_first = ddttChs$Vagitype_first[match(metaData$pt_ID, ddttChs$pt_ID)]
  ddtt$Vagitype_first = ddttChs$Vagitype_first[match(ddtt$pt_ID, ddttChs$pt_ID)]
  # ddttChs$Vagitype_first = factor(ddttChs$CST_first, levels = names(colorCST))
  # metaData$Vagitype_first = ddttChs$Vagitype_first[match(metaData$pt_ID, ddttChs$pt_ID)]
  # ddtt$Vagitype_first = ddttChs$Vagitype_first[match(ddtt$pt_ID, ddttChs$pt_ID)]
  ddttChs$Ethnicity = metaData$Ethnicity2[match(ddttChs$pt_ID, metaData$pt_ID)]
  ddttChs$Term = metaData$Term_char[match(ddttChs$pt_ID, metaData$pt_ID)]
  cs = table(ddttChs$CST_first, ddttChs$Term)
  chisq_test(cs)
  cs_pre = apply(as.matrix(cs), 2, function(x){round(x/sum(x), 3)})
  write.xlsx(as.table(cs_pre), file = './table/CST_first_vs_term.xlsx')
  
  
  cs = table(ddttChs$Vagitype_first, ddttChs$Ethnicity)
  chisq_test(cs)
  cs_pre = apply(as.matrix(cs), 2, function(x){round(x/sum(x), 2)})
  write.xlsx(as.table(cs_pre), file = './table/CST_first_vs_ethnicity.xlsx')
  p1 = ggplot() +
    geom_point(data = ddtt, mapping = aes(x = time, y = SampleID, fill = CST), alpha = 0.65, size = 3, shape = 21, stroke = 0.1) +
    geom_point(data = metaData, mapping = aes(x = Term, y = SampleID, color = Term_char2), shape = ')', size = 4) +
    geom_point(data = ddttChs, mapping = aes(x = 42, y = SampleID, shape = CST_changes), size = 3,fill = '#BF5B16') +
    geom_vline(xintercept = 39, linetype = 'dashed', color = 'red', size = 0.3) +
    geom_vline(xintercept = c(14,26), linetype = 'dashed', color = 'black', size = 0.3) +
    scale_fill_manual(values = colorCST, name = 'CST', guide = guide_legend(order = 1)) +
    scale_color_manual(values = colorTerm2, name = 'Term', guide = guide_legend(order = 2)) +
    scale_shape_manual(values = c(2, 24), name = 'CST changed', guide = guide_legend(order = 3)) +
    ylab('Participants') +
    xlab('Gestational Weeks') +
    mytheme + 
    theme(strip.text.y = element_blank(), legend.position = 'left', 
          plot.margin = unit(c(5, 0, 5, 5), 'mm'),
          axis.title = element_text(size = 16)) +
    facet_grid(Vagitype_first ~ ., scales = 'free_y', space="free")
  p2 = ggplot() +
    geom_point(data = ddtt, mapping = aes(x = time, y = SampleID, fill = Vagitypes), alpha = 0.65, size = 3, shape = 21, stroke = 0.1) +
    geom_point(data = metaData, mapping = aes(x = Term, y = SampleID, color = Term_char2), shape = ')', size = 4) +
    geom_point(data = ddttChs, mapping = aes(x = 42, y = SampleID, shape = Vagitype_changes), size = 3,fill = '#BF5B16') +
    geom_vline(xintercept = 39, linetype = 'dashed', color = 'red', size = 0.3) +
    geom_vline(xintercept = c(14,26), linetype = 'dashed', color = 'black', size = 0.3) +
    scale_fill_manual(values = colorVagitype, name = 'Vagitype', guide = guide_legend(order = 1)) +
    scale_color_manual(values = colorTerm2, name = 'Term', guide = guide_legend(order = 2)) +
    scale_shape_manual(values = c(2, 24), name = 'Vagitype changed', guide = guide_legend(order = 3)) +
    ylab('Participants') +
    xlab('Gestational Weeks') +
    mytheme +
    theme(strip.text.y = element_blank(), legend.position = 'left',
          plot.margin = unit(c(5, 0, 5, 5), 'mm'), axis.title = element_text(size = 16)) +
    facet_grid(Vagitype_first~., scales = 'free_y', space="free")
  ddttChs$bcdmp_median = bcDissimilarityMedianPerson$bcdmp_median[match(ddttChs$SampleID, bcDissimilarityMedianPerson$SampleID)]
  ddttChs$SampleID = factor(ddttChs$SampleID, levels = unique(ddttChs$SampleID[order(ddttChs$bcdmp_median, decreasing = T)]))
  sy = ggplot(ddttChs,aes(x=1,y=pt_ID,fill=bcdmp_median))+
    geom_tile() + ylab('y') +
    # scale_y_discrete(expand=c(0,0))+
    scale_x_continuous(expand=c(0,0), breaks = 1, label = '1')+
    scale_fill_gradient(low = "#132B43", high = "#56B1F7", breaks=seq(0, 1, 0.25), 
                        limits = c(0, 1), name = 'Within-individual median\nBC dissimilarity',
                        guide = guide_colorbar(title.position = 'top',direction = 'horizontal',title.hjust = 0.5, barwidth = unit(5,'cm'))) +
    mytheme + theme(legend.position="right", strip.text.y = element_blank(), 
          axis.title.y=element_blank(), axis.text.y=element_blank(),
          axis.ticks.y=element_blank(), axis.text.x=element_text(color="white"),
          axis.title.x=element_text(color="white", size = 16),
          plot.margin = unit(c(5, 5, 5, 0), 'mm'),
          axis.ticks.x=element_line(color="white"), panel.grid = element_blank()) + 
    facet_grid(Vagitype_first~., scales = 'free_y', space="free")

  #Draw plots one by one in their positions
  cairo_pdf(file = './plot/CST_during_pergnant_order_by_dissimilarity_add_bc_3.pdf',
            width = 9.47, height = 7, onefile = T)
  #Define layout for the plots (2 rows, 2 columns)
  layt<-grid.layout(nrow=1,ncol=2,heights=c(7/8, 7/8),widths=c(4.8/7, 2/7),default.units=c('null','null'))
  #View the layout of plots
  grid.show.layout(layt)
  grid.newpage()
  pushViewport(viewport(layout=layt))
  print(p1,vp=viewport(layout.pos.row=1,layout.pos.col=1))
  print(sy,vp=viewport(layout.pos.row=1,layout.pos.col=2))
  dev.off()
  cairo_pdf(file = './plot/Vagitype_during_pergnant_order_by_dissimilarity_add_bc_2.pdf',
            width = 10.47, height = 7, onefile = T)
  layt<-grid.layout(nrow=1,ncol=2,heights=c(7/8, 7/8),widths=c(6/8, 2.2/8),default.units=c('null','null'))
  #View the layout of plots
  grid.show.layout(layt)
  grid.newpage()
  pushViewport(viewport(layout=layt))
  print(p2,vp=viewport(layout.pos.row=1,layout.pos.col=1))
  print(sy,vp=viewport(layout.pos.row=1,layout.pos.col=2))
  dev.off()
}

table(ddttChs$Vagitype_first)
ggboxplot(ddttChs, x = 'Vagitype_first', y = 'bcdmp_median',
            ylab = 'bcdmp_median', color = 'Vagitype_first', 
            add = "jitter", 
            add.params = list(shape = 21, size = 1.8, color = 'Vagitype_first', fill = 'Vagitype_first', alpha = 0.7),
            bxp.errorbar = T, notch = F) +
  mytheme + theme(axis.text.x = element_text(angle = 45, hjust = 1))

####### statistics for vagitypes ######
load('data/RData/sampleVagitypes.RData')
statDF = sampleVagitypes
# statDF = sampleVagitypes[naturalorder(sampleVagitypes$pt_ID.u), ]
# statDF = statDF[!duplicated(statDF$pt_ID), ]

if(T){
  cairo_pdf('./plot/vagitype_cst_stats_table.pdf',  width = 4.98, height = 6.22)
  # relationship between vagitype and CST
  mytable1 = table(statDF$vagitypes_Fettweis, statDF$vagitypes_Susan)
  grid.table(mytable1)
  dev.off()
  cairo_pdf('./plot/vagitype_cst_stats_table_1.pdf',  width = 7, height = 6.22)
  mytable1_pre = apply(as.matrix(mytable1), 2, function(x){paste0(x, ' (', round(x/sum(x), 3), ')')})
  row.names(mytable1_pre) = row.names(mytable1)
  colnames(mytable1_pre) = colnames(mytable1)
  grid.table(mytable1_pre)
  dev.off()
  cairo_pdf('./plot/vagitype_cst_stats_table_2.pdf',  width = 7, height = 6.22)
  mytable1_pre2 = apply(as.matrix(mytable1), 2, function(x){round(x/sum(x), 3)})
  row.names(mytable1_pre2) = row.names(mytable1)
  colnames(mytable1_pre2) = colnames(mytable1)
  grid.table(mytable1_pre2)
  dev.off()
  l = list(as.matrix(mytable1), as.table(mytable1_pre2), as.table(mytable1_pre))
  write.xlsx(l, file = './table/CST_vs_Vagitype.xlsx')
  
  mytable2 = table(statDF$vagitypes_Fettweis2, statDF$vagitypes_Susan)
  chisq.test(mytable1)
  chisq.test(mytable2)
  # relationship between term and vagitype/CST
  mytable3 = table(statDF$vagitypes_Fettweis, statDF$Term_char)
  round(mytable3/colSums(mytable3), 2)
  mytable4 = table(statDF$vagitypes_Fettweis2, statDF$Term_char)
  mytable5 = table(statDF$vagitypes_Susan, statDF$Term_char)
  chisq.test(mytable3)
  chisq.test(mytable4)
  
  # relationship between ethnicity and vagitype/CST
  statDF$Ethnicity2 = statDF$Ethnicity %>% as.character()
  statDF$Ethnicity2[statDF$Ethnicity2  == "Pacific Islander"] = 'Other'
  statDF$Ethnicity2 %<>% factor(., levels = names(colorEthnicity2))
  mytable6 = table(statDF$vagitypes_Fettweis, statDF$Ethnicity2)
  chisq_test(mytable6)
  
  mytable6_pre = apply(as.matrix(mytable6), 2, function(x){round(x/sum(x), 3)})
  write.xlsx(as.table(mytable6_pre), file = './table/Ethnicity_vs_Vagitype.xlsx')
  mytable7 = table(statDF$CST, statDF$Ethnicity2)
  chisq_test(mytable7)
  mytable7_pre = apply(as.matrix(mytable7), 2, function(x){round(x/sum(x), 3)})
  write.xlsx(as.table(mytable7_pre), file = './table/CST_vs_Vagitype.xlsx')
  
  # relationship between term and vagitype/CST
  mytable6 = table(statDF$CST, statDF$Term_char)
  chisq_test(mytable6)
  mytable6_pre = apply(as.matrix(mytable6), 2, function(x){round(x/sum(x), 3)})
  write.xlsx(as.table(mytable6_pre), file = './table/Term_vs_Vagitype.xlsx')
  mytable7 = table(statDF$CST, statDF$Term_char)
  chisq_test(mytable7)
  mytable7_pre = apply(as.matrix(mytable7), 2, function(x){round(x/sum(x), 3)})
  write.xlsx(as.table(mytable7_pre), file = './table/Term_vs_cst.xlsx')

  
  
  cairo_pdf('./plot/CST_propotion.pdf', width = 8.15, height = 6.84, onefile = T)
  if(T){
    CSTvsEthnicity = table(statDF$Ethnicity2, statDF$vagitypes_Susan)
    CSTvsEthnicity2 = apply(CSTvsEthnicity, 2, function(x){x*100/sum(x)})
    CSTvsEthnicity = melt(CSTvsEthnicity)
    CSTvsEthnicity2 = melt(CSTvsEthnicity2)
    colnames(CSTvsEthnicity) = c('Ethnicity', 'CST', 'Count')
    colnames(CSTvsEthnicity2) = c('Ethnicity', 'CST', 'Frequency')
    CSTvsEthnicity = merge(CSTvsEthnicity, CSTvsEthnicity2)
    ethP = ggplot(CSTvsEthnicity, aes(x = "", y = Frequency, fill = Ethnicity)) +
      geom_col(color = "black", width = 0.5) +
      geom_text(aes(label = ifelse(round(Count) != 0, round(Count), '')), 
                position = position_stack(vjust = 0.5), color = 'black') +
      coord_polar(theta = "y") +
      scale_fill_manual(values = colorEthnicity2) +
      facet_grid(.~ CST) +theme_void() +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            panel.grid = element_blank(),
            legend.title = element_text(size = 13),
            legend.text = element_text(size = 13),
            strip.text = element_text(size = 13))
    print(ethP)
  }
  if(T){
    CSTvsTerm = table(statDF$Term_char, statDF$vagitypes_Susan)
    CSTvsTerm2 = apply(CSTvsTerm, 2, function(x){x*100/sum(x)})
    CSTvsTerm = melt(CSTvsTerm)
    CSTvsTerm2 = melt(CSTvsTerm2)
    colnames(CSTvsTerm) = c('Term', 'CST', 'Count')
    colnames(CSTvsTerm2) = c('Term', 'CST', 'Frequency')
    CSTvsTerm = merge(CSTvsTerm, CSTvsTerm2)
    terP = ggplot(CSTvsTerm, aes(x = "", y = Frequency, fill = Term)) +
      geom_col(color = "black", width = 0.5) +
      geom_text(aes(label = ifelse(round(Count) != 0, round(Count), '')), 
                position = position_stack(vjust = 0.5), color = 'white') +
      coord_polar(theta = "y") +
      scale_fill_manual(values = colorTerm) +
      facet_grid(.~ CST) +theme_void() +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            panel.grid = element_blank(),
            legend.title = element_text(size = 13),
            legend.text = element_text(size = 13),
            strip.text = element_text(size = 13))
    print(terP)
  }
  dev.off()
  cairo_pdf('./plot/VAG_propotion.pdf', width = 8.15, height = 7.31, onefile = T)
  if(T){
    VAGvsEthnicity = table(statDF$Ethnicity2, statDF$vagitypes_Fettweis2)
    VAGvsEthnicity2 = apply(VAGvsEthnicity, 2, function(x){x*100/sum(x)})
    VAGvsEthnicity = melt(VAGvsEthnicity)
    VAGvsEthnicity2 = melt(VAGvsEthnicity2)
    colnames(VAGvsEthnicity) = c('Ethnicity', 'VAG', 'Count')
    colnames(VAGvsEthnicity2) = c('Ethnicity', 'VAG', 'Frequency')
    VAGvsEthnicity = merge(VAGvsEthnicity, VAGvsEthnicity2)
    ethP = ggplot(VAGvsEthnicity, aes(x = "", y = Frequency, fill = Ethnicity)) +
      geom_col(color = "black", width = 0.5) +
      geom_text(aes(label = ifelse(round(Count) != 0, round(Count), '')), 
                position = position_stack(vjust = 0.5), color = 'black') +
      coord_polar(theta = "y") +
      scale_fill_manual(values = colorEthnicity2) +
      facet_wrap(.~ VAG) +theme_void() +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            panel.grid = element_blank(),
            legend.title = element_text(size = 13),
            legend.text = element_text(size = 13),
            strip.text = element_text(size = 13))
    print(ethP)
  }
  if(T){
    VAGvsTerm = table(statDF$Term_char, statDF$vagitypes_Fettweis2)
    VAGvsTerm2 = apply(VAGvsTerm, 2, function(x){x*100/sum(x)})
    VAGvsTerm = melt(VAGvsTerm)
    VAGvsTerm2 = melt(VAGvsTerm2)
    colnames(VAGvsTerm) = c('Term', 'VAG', 'Count')
    colnames(VAGvsTerm2) = c('Term', 'VAG', 'Frequency')
    VAGvsTerm = merge(VAGvsTerm, VAGvsTerm2)
    terP = ggplot(VAGvsTerm, aes(x = "", y = Frequency, fill = Term)) +
      geom_col(color = "black", width = 0.5) +
      geom_text(aes(label = ifelse(round(Count) != 0, round(Count), '')), 
                position = position_stack(vjust = 0.5), color = 'white') +
      coord_polar(theta = "y") +
      scale_fill_manual(values = colorTerm) +
      facet_wrap(.~ VAG) +theme_void() +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            panel.grid = element_blank(),
            legend.title = element_text(size = 13),
            legend.text = element_text(size = 13),
            strip.text = element_text(size = 13))
    print(terP)
  }
  dev.off()
}
###### estimate the transition of CST/vagitype #####
# the distribution of space between two sequential samples
if(T){
  spaceDays = c()
  for(p in metaData$pt_ID){
    subDF = swabData[swabData$pt_ID == p & swabData$trimester != 'P', ]
    subDF = subDF[naturalorder(subDF$pt_ID.u), ]
    spaceDays = c(spaceDays, (subDF$Sample_GA[2:nrow(subDF)] - subDF$Sample_GA[1:(nrow(subDF)-1)])*7)
    if(55.3 %in% spaceDays){print(p)}
  }
  range(spaceDays)
  dp = ggplot(data.frame(x = spaceDays), aes(x = x)) + 
    geom_histogram(binwidth = 1, fill = '#FDF4D6', color = 'black') +
    scale_x_continuous(breaks = seq(0,55,5)) +
    labs(x = 'Interval days of sequential samples', y = 'Count') +
    mytheme
  cairo_pdf('./plot/Interval days of sequential samples.pdf', width = 4.14, height = 4.04)
  print(dp)
  dev.off()
}
# set window of sequential samples and filter out samples
window = c(11, 17) # biweekly, 11, 17 about two weeks; 
# window = c(4, 10) # weekly, 4, 10, ~ about 1 weeks
# window = c(3.5, 55.3) # all
space = 'biweekly'
# space = 'weekly'
# space = 'all'
if(T){
  sequentialDF = data.frame()
  for(p in metaData$pt_ID){
    subDF = swabData[swabData$pt_ID == p & swabData$trimester != 'P', ]
    subDF = subDF[naturalorder(subDF$pt_ID.u), ]
    spaceDays = (subDF$Sample_GA[2:nrow(subDF)] - subDF$Sample_GA[1:(nrow(subDF)-1)])*7
    preInd = which(spaceDays >= window[1] & spaceDays <= window[2])
    nextInd = preInd + 1
    sequentialDF = rbind(sequentialDF, data.frame(
      preSample = subDF$pt_ID.u[preInd],
      nextSample = subDF$pt_ID.u[nextInd]
    ))
  }
  nrow(sequentialDF)
}

cairo_pdf(paste0("./plot/cst_markov_chain_", space, ".pdf"), 
          width=4.41, height=3.72, onefile = T)
if(T){ # transition rate matrix and Markov Chain
  using_Col = 'vagitypes_Susan' # vagitypes_Fettweis     vagitypes_Fettweis2 vagitypes_Susan CST               Vagitypes
  sequentialDF$preState = sampleVagitypes[match(sequentialDF$preSample, sampleVagitypes$pt_ID.u), using_Col]
  sequentialDF$nextState = sampleVagitypes[match(sequentialDF$nextSample, sampleVagitypes$pt_ID.u), using_Col]
  nstates = length(levels(sampleVagitypes[, using_Col]))
  CSTs = levels(sampleVagitypes[, using_Col])
  ttab <- table(sequentialDF$preState, sequentialDF$nextState) # prevstate=row, curstate=col
  trans <- matrix(ttab, nrow=nstates)
  trans <- trans/rowSums(trans)  # Normalize row sums to 1
  CSTtrans <- trans
  colnames(CSTtrans) <- CSTs
  rownames(CSTtrans) <- CSTs
  t_persist <- -1/log(diag(CSTtrans))
  table_CSTtrans = as.data.frame(CSTtrans)
  row.names(table_CSTtrans) = row.names(CSTtrans)
  l = list(ttab, table_CSTtrans)
  write.xlsx(l, file = './table/CST_biweekly_trans.xlsx', row.name =T)
  sgrid.table(CSTtrans %>% round(., 3)) # Paper
  t_persist # Paper
  mcPreg <- new("markovchain", states=CSTs,
                transitionMatrix = trans, name="PregCST")
  mcPreg
}
if(T){ # Set up igraph of the markov chain
  library(markovchain)
  library(igraph)
  netMC <- markovchain:::.getNet(mcPreg, round = TRUE)
  wts <- E(netMC)$weight/100
  
  edgel <- get.edgelist(netMC)
  elcat <- paste(edgel[,1], edgel[,2])
  elrev <- paste(edgel[,2], edgel[,1])
  edge.curved <- sapply(elcat, function(x) x %in% elrev)
  
  samdf_def <- sampleVagitypes[!is.na(sampleVagitypes$Term_char2),] # Only those definitely assigned, i.e. not marginal
  
  premat <- table(samdf_def[, using_Col], samdf_def$Term_char2)
  rownames(premat) <- markovchain::states(mcPreg)
  colnames(premat) <- c("Preterm", "Full-term")
  premat
  premat <- premat/rowSums(premat)
  grid.newpage()
  grid.table(premat %>% round(., 3))
  CSTColors <- brewer.pal(6,"Paired")[c(1,3,2,5,4,6)] # Length 6 for consistency with pre-revision CST+ coloration
  names(CSTColors) <- CSTs
  vert.CSTclrs <- CSTColors
}

if(T){
  library(igraph)
  default.par <- par(no.readonly = TRUE)
  # Define color scale
  # Plotting function for markov chain
  plotMC <- function(object, ...) {
    netMC <- markovchain:::.getNet(object, round = TRUE)
    plot.igraph(x = netMC, ...)  
  }
  # Color bar for the markov chain visualization, gradient in strength of preterm association
  color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title=NULL) {
    scale = (length(lut)-1)/(max-min)
    
    #    dev.new(width=1.75, height=5)
    # cur.par <- par(no.sreadonly=T)
    par(mar=c(0,7,1,7)+0.1, oma=c(0,0,0,0)+0.1, mgp=c(2,0.2,0))
    # par(ps = 10, cex = 0.8)
    # par(tcl=-0.2, cex.axis=0.8, cex.lab = 0.8)
    par(tcl=-0.2)
    plot(c(min,max), c(0,10), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
    axis(1, seq(0, 1, 0.25))
    for (i in 1:(length(lut)-1)) {
      x = (i-1)/scale + min
      rect(x,0,x+1/scale,10, col=lut[i], border=NA)
    }
  }

  pal <- colorRampPalette(c("grey50", "maroon", "magenta2"))(101)
  library(viridis)
  pal <- magma(101)
  vert.clrs <- sapply(states(mcPreg), function(x) pal[1+round(100*premat[x,"Preterm"])])
  vert.sz <- 4 + 2*sapply(states(mcPreg), 
                          function(x) length(unique(samdf_def[samdf_def[, using_Col]==x,"pt_ID"])))
  vert.sz <- vert.sz * 0.85
  vert.font.clrs <- c("white", "white", "white", "white", "white")
  # E(netMC) to see edge list, have to define loop angles individually by the # in edge list, not vertex
  edge.loop.angle = c(0, 0, 0, 0, 3.14, 3.14, 0, 0, 0, 3.14, 3.14, 0.60, 0, 0, 0, 0.45, 0,  0, 0, 0, 0)-0.45
  
  edge.loop.angle = rep(0, length(E(netMC)))
  edgesMat = ends(netMC,es = E(netMC))
  edgesList = paste0(edgesMat[, 1], ' -> ', edgesMat[, 2])
  edge.loop.angle[edgesList == 'I -> I'] = -0.45
  edge.loop.angle[edgesList == 'II -> II'] = -0.45
  edge.loop.angle[edgesList == 'III -> III'] = 3.25
  edge.loop.angle[edgesList == 'IV -> IV'] = 0.45

  layout <- matrix(c(0.6,0.95, 0.43,1, 0.3,0.66, 0.55,0.3, 0.75,0.65), nrow=5, ncol=2, byrow=T)
  
  # Colored by association with preterm birth
  layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE), heights=c(1,10))
  color.bar(pal, min=0, max=1, nticks=6, title="Fraction preterm")
  par(mar=c(0,1,1,1)+0.1)
  edge.arrow.size=0.6
  edge.arrow.width=1.2
  edge.width = (15*wts + 0.1)*0.6
  edge.labels <- as.character(round(E(netMC)$weight/100, 3))
  # edge.labels[edge.labels<0.4] <- NA  # labels only for self-loops
  par(family = "serif")
  plotMC(mcPreg, 
         edge.arrow.size=edge.arrow.size, edge.arrow.width = edge.arrow.width,
         edge.label = edge.labels, edge.label.cex=1.2, edge.label.color="black",
         # FIX EDGE LABELS FOR PUBLICATION IN POST-PROCESSING
         edge.width=edge.width, edge.curved=edge.curved, edge.color = 'lightgray',
         vertex.color=vert.clrs, vertex.size=(vert.sz),
         vertex.label.font = 2, vertex.label.cex = 1,
         vertex.label.color = vert.font.clrs, vertex.frame.color = 'black', 
         layout=layout, edge.loop.angle = edge.loop.angle)
  #dev.off()
  par(default.par)
}
dev.off()
# *********************************
# vagitypes transition
# *********************************

cairo_pdf(paste0("./plot/vagitype_markov_chain_", space, ".pdf"), 
          width=11, height=4.28, onefile = T)
if(T){ # transition rate matrix and Markov Chain
  sampleVagitypes$vagitypes_Fettweis3 = sampleVagitypes$vagitypes_Fettweis2 %>% as.character()
  table(sampleVagitypes$vagitypes_Fettweis3)
  sampleVagitypes$vagitypes_Fettweis3 %<>% str_replace_all(., 'Lactobacillus', 'L.') %>% str_replace_all(., 'Fannyhessea', 'F.') %>% str_replace_all(., 'Gardnerella', 'G.') %>% str_replace_all(., 'Gardnerella', 'G.')
  sampleVagitypes$vagitypes_Fettweis3 %<>% factor(., levels = names(colorVagitypeShort))
  using_Col = 'vagitypes_Fettweis3' # vagitypes_Fettweis     vagitypes_Fettweis2 vagitypes_Susan CST               Vagitypes
  sequentialDF$preState = sampleVagitypes[match(sequentialDF$preSample, sampleVagitypes$pt_ID.u), using_Col]
  sequentialDF$nextState = sampleVagitypes[match(sequentialDF$nextSample, sampleVagitypes$pt_ID.u), using_Col]
  nstates = length(levels(sampleVagitypes[, using_Col]))
  CSTs = levels(sampleVagitypes[, using_Col])
  ttab <- table(preState = sequentialDF$preState, nextState = sequentialDF$nextState) # prevstate=row, curstate=col
  trans <- matrix(ttab, nrow=nstates)
  trans <- trans/rowSums(trans)  # Normalize row sums to 1
  CSTtrans <- trans
  colnames(CSTtrans) <- CSTs
  rownames(CSTtrans) <- CSTs
  t_persist <- -1/log(diag(CSTtrans))
  grid.table(CSTtrans %>% round(., 3) )# Paper
  t_persist # Paper
  mcPreg <- new("markovchain", states=CSTs,
                transitionMatrix = trans, name="PregCST")
  mcPreg
  
  table_CSTtrans = as.data.frame(CSTtrans)
  row.names(table_CSTtrans) = row.names(CSTtrans)
  l = list(ttab, table_CSTtrans)
  write.xlsx(l, file = './table/Vagitype_biweekly_trans.xlsx', row.name =T)

}
if(T){ # Set up igraph of the markov chain
  netMC <- markovchain:::.getNet(mcPreg, round = TRUE)
  wts <- E(netMC)$weight/100
  
  edgel <- get.edgelist(netMC)
  elcat <- paste(edgel[,1], edgel[,2])
  elrev <- paste(edgel[,2], edgel[,1])
  edge.curved <- sapply(elcat, function(x) x %in% elrev)

  samdf_def <- sampleVagitypes[!is.na(sampleVagitypes$Term_char2),] # Only those definitely assigned, i.e. not marginal
  premat <- table(samdf_def[, using_Col], samdf_def$Term_char2)
  rownames(premat) <- markovchain::states(mcPreg)
  colnames(premat) <- c("Preterm", "Full-term")
  premat
  premat <- premat/rowSums(premat)
  grid.newpage()
  grid.table(premat %>% round(., 3) )
  CSTColors <- colorVagitype # Length 6 for consistency with pre-revision CST+ coloration
  names(CSTColors) <- CSTs
  vert.CSTclrs <- CSTColors
}

if(T){
  library(igraph)
  default.par <- par(no.readonly = TRUE)
  # Define color scale
  # Plotting function for markov chain
  plotMC <- function(object, ...) {
    netMC <- markovchain:::.getNet(object, round = TRUE)
    plot.igraph(x = netMC, ...)  
  }
  # Color bar for the markov chain visualization, gradient in strength of preterm association
  color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title=NULL) {
    scale = (length(lut)-1)/(max-min)
    
    #    dev.new(width=1.75, height=5)
    # cur.par <- par(no.sreadonly=T)
    par(mar=c(0,27,1,27)+0.1, oma=c(0,0,0,0)+0.1, mgp=c(2,0.2,0))
    # par(ps = 10, cex = 0.8)
    # par(tcl=-0.2, cex.axis=0.8, cex.lab = 0.8)
    par(tcl=-0.2)
    plot(c(min,max), c(0,10), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
    axis(1, seq(0, 1, 0.25))
    for (i in 1:(length(lut)-1)) {
      x = (i-1)/scale + min
      rect(x,0,x+1/scale,10, col=lut[i], border=NA)
    }
  }
  
  pal <- colorRampPalette(c("grey50", "maroon", "magenta2"))(101)
  library(viridis)
  pal <- magma(101)
  vert.clrs <- sapply(states(mcPreg), function(x) pal[1+round(100*premat[x,"Preterm"])])
  vert.sz <- 4 + 2*sapply(states(mcPreg), 
                          function(x) length(unique(samdf_def[samdf_def[, using_Col]==x,"pt_ID"])))
  vert.sz <- vert.sz * 0.85
  vert.font.clrs <- c("white", "white", "white", "white", "white")
  # E(netMC) to see edge list, have to define loop angles individually by the # in edge list, not vertex
  edge.loop.angle = rep(0, length(E(netMC)))
  edgesMat = ends(netMC,es = E(netMC))
  edgesList = paste0(edgesMat[, 1], ' -> ', edgesMat[, 2])
  edge.loop.angle[str_count(edgesList, 'F. vaginae') == 2] = 3.25
  edge.loop.angle[str_count(edgesList, 'L. jensenii') == 2] = 3.25
  edge.loop.angle[str_count(edgesList, 'L. gasseri') == 2] = 3.25
  edge.loop.angle[str_count(edgesList, 'L. coleohominis') == 2] = 2.20
  edge.loop.angle[str_count(edgesList, 'L. iners') == 2] = -0.45
  edge.loop.angle[str_count(edgesList, 'None dominant') == 2] = 0.45
  edge.loop.angle[str_count(edgesList, 'Other') == 2] = 0.45
  layouts <- grep("^layout_", ls("package:igraph"), value=TRUE)[-1] 
  layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama|tree", layouts)]
  layout = "layout_in_circle"
  l <- do.call(layout, list(netMC)) 
  # Remove layouts that do not apply to our graph.
  
  # Colored by association with preterm birth
  layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE), heights=c(1,10))
  color.bar(pal, min=0, max=1, nticks=6, title="Fraction preterm")
  par(mar=c(1,1,1,1)+0.1)
  edge.arrow.size=0.6
  edge.arrow.width=1.2
  edge.width = (15*wts + 0.1)*0.6
  edge.labels <- as.character(round(E(netMC)$weight/100, 3))
  edge.labels[edge.labels<0.25] <- NA  # labels only for self-loops
  
  par(family = "serif")
  plotMC(mcPreg, edge.arrow.size=edge.arrow.size, edge.arrow.width = edge.arrow.width,
         edge.label = edge.labels, edge.label.cex=1.2, edge.label.color="black",
         # FIX EDGE LABELS FOR PUBLICATION IN POST-PROCESSING
         edge.width=edge.width, edge.curved=edge.curved, edge.color = 'lightgray',
         vertex.color=vert.clrs, vertex.size=(vert.sz),
         vertex.label.font = 4, vertex.label.cex = 1,
         vertex.label.color = 'red', vertex.frame.color = NA, 
         edge.loop.angle = edge.loop.angle,
         layout=l)
  #dev.off()
  par(default.par)
}
dev.off()

##### VALENCIA ######
# load the reference centroids csv
load('./data/RData/taxa2taxidDF.RData')
refCSV = read.csv('./reference_code/VALENCIA/CST_centroids_012920.csv', check.names = F)
refTaxa = data.frame(ref_taxa = colnames(refCSV)[2:ncol(refCSV)])
refTaxa$ref_taxa2[!grepl('^[dkpcofgs]_', refTaxa$ref_taxa)] = paste0('s__', refTaxa$ref_taxa[!grepl('^[dkpcofgs]_', refTaxa$ref_taxa)])
for(l in c('^d_', '^k_', '^p_', '^c_', '^o_', '^f_', '^g_', '^s_')){
  refTaxa$ref_taxa2[grepl(l, refTaxa$ref_taxa)] = refTaxa$ref_taxa[grepl(l, refTaxa$ref_taxa)] %>% str_replace(., l, paste0(str_sub(l, 2, 3), '_'))
}
refTaxa$ref_taxa2[refTaxa$ref_taxa2 == 's__Atopobium_vaginae'] = 's__Fannyhessea_vaginae'
refTaxa$ref_taxa2[refTaxa$ref_taxa2 == 'd__Eukaryota'] = 'k__Eukaryota'
refTaxa$ref_taxa2[refTaxa$ref_taxa2 == 'g__Escherichia.Shigella'] = 'g__Escherichia'
refTaxa$ref_taxa2[refTaxa$ref_taxa2 == 's__Lactobacillus_fermentum'] = 's__Limosilactobacillus_fermentum'
refTaxa$ref_taxa2[refTaxa$ref_taxa2 == 'g__Burkholderia.Paraburkholderia'] = 'o__Burkholderiales'
refTaxa$ref_taxa2[refTaxa$ref_taxa2 == 's__Atopobium_parvulum'] = 's__Lancefieldella_parvula'
# refTaxa$ref_taxa2[refTaxa$ref_taxa2 == 'f__Bacteroidales_S24.7_group'] = 'f__Bacteroidales_unclassified'
# refTaxa$ref_taxa2[refTaxa$ref_taxa2 == 'g__Acidaminococcus'] = 'f__Acidaminococcaceae'

refTaxa$my_taxa = taxa2taxidDF$taxa[match(refTaxa$ref_taxa2, taxa2taxidDF$taxa)]


