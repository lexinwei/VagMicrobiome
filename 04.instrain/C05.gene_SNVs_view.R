source('./script/library_function.R')
library(trackViewer)
load('./data/RData/MTMG_metaData.RData')
load('./data/RData/MTMG_statDF3_all.RData')
load('./data/RData/colorList.RData')
load('./data/RData/MTMG_geneDF.RData')
##### pick genes ####
if(T){
  # L. crisptus, Muc_B2, NZ_CP039266.1_1594
  
  # L. crisptus, Gram_pos_anchor, NZ_CP039266.1_48
  
  # L. iners, Gram_pos_anchor, NZ_CP045664.1_1
  
  # L. jensenii, Muc_B2, NZ_CP018809.1_5
}
lolPlot = function(., type = 'AA', keep_mut_on_domain = F){ # or type = 'NT'
  if(type == 'AA'){
    wholeWidth = geneLength/3
  }else{
    wholeWidth = geneLength
  }
  snvDF = snvDF[snvDF$sample %in% geneDF$sample[geneDF$gene == geneID], ]
  GeneDirection = geneDF$direction[geneDF$gene == geneID][1]
  table(snvDF$class)
  # a = snvDF[, c('position', 'mutation', 'class')]
  # a = a[!duplicated(a), ]
  # table(duplicated(a$mutation))
  classDF = snvDF %>% group_by(mutation, class) %>% 
    summarise(class_counts = n())
  classMat0 = reshape2::dcast(classDF, mutation ~ class)
  classMat0[is.na(classMat0)] = 0
  classMat = data.frame(mutation = classMat0[, 1], classMat0[, -1]/rowSums(classMat0[, -1]))
  # b = classDF[classDF$mutation %in% names(table(classDF$mutation))[table(classDF$mutation) > 2], ]
  classDFtop = classDF %>% group_by(mutation) %>% 
    arrange(desc(class_counts), .by_group = TRUE) %>% top_n(1, class_counts)
  classDFtop$class2 = classDFtop$class %>% str_replace_all(., 'con_|pop_', '')
  
  snvDF$position = snvDF$position + 1
  snvDF$NT_pos_rel = str_extract(snvDF$mutation, '\\d+') %>% as.numeric()
  snvDF$seq_ID = str_extract(snvDF$sample, '\\d+_\\d+')
  snvCountDF = snvDF[, c('mutation', 'sample')] %>% group_by(mutation) %>% 
    summarise(Counts = n())
  snvCountDF$Freq = snvCountDF$Counts/length(unique(snvDF$sample))*100
  snvCountDF$NT_pos = snvDF$position[match(snvCountDF$mutation, snvDF$mutation)]
  snvCountDF$NT_pos_rel = snvCountDF$NT_pos - scaffoldStart + 1
  snvCountDF$NT_mut = str_replace(snvCountDF$mutation, '\\d+', as.character(snvCountDF$NT_pos_rel))
  if(GeneDirection == '-1'){
    snvCountDF$NT_pos_rel = geneLength - snvCountDF$NT_pos_rel + 1
  }
  snvCountDF$AA_pos = ceil(snvCountDF$NT_pos_rel/3)
  snvCountDF$AA_mut = str_replace(snvCountDF$mutation, '\\d+', as.character(snvCountDF$AA_pos))
  table(duplicated(snvCountDF$AA_mut))
  snvCountDF = merge(snvCountDF, classMat0, by.x = 'mutation', by.y = 'mutation', all.x = T)
  snvCountDF$SNV.all = rowSums(snvCountDF[, colnames(snvCountDF) %in% c('con_SNV', 'pop_SNV', 'SNV')])
  if(type == 'AA'){
    snvCountDF$X = snvCountDF$AA_pos
  }else{
    snvCountDF$X = snvCountDF$NT_pos_rel
  }
  snvCountDF = snvCountDF[!snvCountDF$NT_pos_rel %in% excludePos, ]
  fileN = ''
  if(keep_mut_on_domain){
    keepV = c()
    for(d in 1:length(domainStart)){
      keepV = c(keepV, seq(domainStart[d], domainStart[d] + domainLength[d]))
    }
    snvCountDF = snvCountDF[snvCountDF$X %in% keepV, ]
    fileN = '_on_domain'
  }
  range(snvCountDF$NT_pos_rel)
  snvCountDF$mutation_type = snvDF$mutation_type[match(snvCountDF$mutation, snvDF$mutation)]
  snvCountDF$class = classDFtop$class2[match(snvCountDF$mutation, classDFtop$mutation)]
  snvCountDF = arrange(snvCountDF, X)
  
  features <- GRanges("chr1", IRanges(c(1, domainStart), 
                                      height=c(0.02, domainHeight),
                                      width = c(wholeWidth, domainLength-1), 
                                      fill = c('gray', domainColor),
                                      color = c('gray', domainColor),
                                      names = c(geneName, domainName)))
  keepInd = which(snvCountDF$Counts > 1)
  if(type == 'AA'){
    SNP = snvCountDF$AA_pos[keepInd]
    SNP_name = snvCountDF$AA_mut[keepInd]
  }else{
    SNP = snvCountDF$NT_pos_rel[keepInd]
    SNP_name = snvCountDF$NT_mut[keepInd]
  }
  SNP_name[startsWith(SNP_name, 'S:')] = ''
  SNP_name %<>% str_remove_all(., 'N:')
  SNP.gr <- GRanges("chr1", IRanges(SNP, width=1, names = SNP_name))
  # SNP.gr$color = ifelse(snvCountDF$mutation_type[keepInd] == 'N', fillSyno[1], fillSyno[2])
  SNP.gr$color = ifelse(snvCountDF$class[keepInd] == 'SNS', fillClass[1], fillClass[2])
  SNP.gr$alpha = 0.5
  SNP.gr$dashline.col = ifelse(snvCountDF$mutation_type[keepInd] == 'N', colorSyno[1], NA)
  SNP.gr$border = ifelse(snvCountDF$mutation_type[keepInd] == 'N', colorSyno[1], colorSyno[2])
  SNP.gr$cex = ifelse(snvCountDF$mutation_type[keepInd] == 'N', 0.7, 0.3)
  # SNP.gr$value1 = snvCountDF$pretermFreq[keepInd]
  # SNP.gr$value2 = 1 - snvCountDF$pretermFreq[keepInd]
  SNP.gr$label.parameter.rot = 60
  SNP.gr$score = snvCountDF$Freq[keepInd]
  # label.parameter.gp.SNS = list(gpar(col = colorClass[1])) # , fontsize = 7
  # label.parameter.gp.SNV = list(gpar(col = colorClass[2]))
  # SNP.gr$label.parameter.gp = ifelse(snvCountDF$class[keepInd] == 'SNS', label.parameter.gp.SNS, label.parameter.gp.SNV)
  SNP.gr$label.parameter.gp = gpar(fontsize = 5)
  SNP.gr$lwd = ifelse(snvCountDF$mutation_type[keepInd] == 'N', 1, 0.5)
  
  SNP2.gr = SNP.gr
  SNP2.gr$cex = ifelse(snvCountDF$mutation_type[keepInd] == 'N', 1, 0.3)
  SNP2.gr$border = ifelse(snvCountDF$mutation_type[keepInd] == 'N', 'black', 'gray')
  SNP2.gr$value1 = snvCountDF$SNS[keepInd]*100/snvCountDF$Counts[keepInd]
  SNP2.gr$value2 = snvCountDF$SNV.all[keepInd]*100/snvCountDF$Counts[keepInd]
  # SNP2.gr$value2 = classMat$SNV[match(snvCountDF$mutation[keepInd], classMat$mutation)]
  # SNP2.gr$value3 = classMat$con_SNV[match(snvCountDF$mutation[keepInd], classMat$mutation)]
  # SNP2.gr$value4 = classMat$pop_SNV[match(snvCountDF$mutation[keepInd], classMat$mutation)]
  SNP2.gr$color <- rep(list(colorClass[1:2]), length(SNP))
  SNP2.gr$alpha = 1
  SNP2.gr$lwd =0.2
  
  legend = list(
    list(labels=c('Non-synonymous', 'Synonymous'), cex = 1, fill= 'white', col = colorSyno), # fill = fillSyno
    list(labels=c('SNS', 'SNV'), cex = 1, fill = fillClass, col = 'white')
    # list(labels=c('SNS', 'SNV', 'con_SNV', 'pop_SNV'), fill = colorClass)
  )
  cairo_pdf(paste0('./plot/MTMG_lol_mutation_', species, '_', geneID, '_', type, fileN, '.pdf'), width = 9.65, height = 6.43)
  lolliplot(list(A = SNP2.gr, B = SNP.gr), 
            list(x = features, y = features), 
            ylab = c('', 'Frequency of mutation'), 
            type=c("pie", "circle"), legend = legend, label_on_feature = F,
            xaxis = xAxis, jitter = c("node", "label"),
            xaxis.gp = gpar(lex = 1, lineheight = 0.7, fontsize=13, lwd=1.5),
            yaxis = seq(0, 100, 25), 
            yaxis.gp = gpar(lex = 1, lineheight = 0.7, fontsize=13, lwd=1.5),
            ranges = GRanges("chr1", IRanges(1, wholeWidth)))
  dev.off()
}

lolPlot2 = function(., type = 'AA', keep_mut_on_domain = F, pdfHei = 3.44){ # or type = 'NT'
  if(type == 'AA'){
    wholeWidth = geneLength/3
  }else{
    wholeWidth = geneLength
  }
  snvDF = snvDF[snvDF$sample %in% geneDF$sample[geneDF$gene == geneID], ]
  GeneDirection = geneDF$direction[geneDF$gene == geneID][1]
  table(snvDF$class)
  snvDF$class2 = str_replace_all(snvDF$class, 'con_|pop_', '')
  table(snvDF$class2)
  
  snvDF$position = snvDF$position + 1
  snvDF$NT_pos_rel = str_extract(snvDF$mutation, '\\d+') %>% as.numeric()
  snvDF$seq_ID = str_extract(snvDF$sample, '\\d+_\\d+')
  snvCountDF = snvDF[, c('mutation', 'sample', 'class2')] %>% group_by(mutation, class2) %>% 
    summarise(Counts = n())
  snvCountDF$Freq = snvCountDF$Counts/length(unique(snvDF$sample))*100
  snvCountDF$NT_pos = snvDF$position[match(snvCountDF$mutation, snvDF$mutation)]
  snvCountDF$NT_pos_rel = snvCountDF$NT_pos - scaffoldStart + 1
  snvCountDF$NT_mut = str_replace(snvCountDF$mutation, '\\d+', as.character(snvCountDF$NT_pos_rel))
  if(GeneDirection == '-1'){
    snvCountDF$NT_pos_rel = geneLength - snvCountDF$NT_pos_rel + 1
  }
  snvCountDF$AA_pos = ceil(snvCountDF$NT_pos_rel/3)
  snvCountDF$AA_mut = str_replace(snvCountDF$mutation, '\\d+', as.character(snvCountDF$AA_pos))
  table(duplicated(snvCountDF$AA_mut))
  if(type == 'AA'){
    snvCountDF$X = snvCountDF$AA_pos
  }else{
    snvCountDF$X = snvCountDF$NT_pos_rel
  }
  snvCountDF = snvCountDF[!snvCountDF$NT_pos_rel %in% excludePos, ]
  fileN = ''
  if(keep_mut_on_domain){
    keepV = c()
    for(d in 1:length(domainStart)){
      keepV = c(keepV, seq(domainStart[d], domainStart[d] + domainLength[d]))
    }
    snvCountDF = snvCountDF[snvCountDF$X %in% keepV, ]
    fileN = '_on_domain'
  }
  range(snvCountDF$NT_pos_rel)
  snvCountDF$mutation_type = snvDF$mutation_type[match(snvCountDF$mutation, snvDF$mutation)]
  snvCountDF$class = snvCountDF$class2
  
  snvCountDF$Score = convScore$Score[match(snvCountDF$AA_pos, convScore$Pos)]
  range(snvCountDF$Score)
  snvCountDF = arrange(snvCountDF, X)
  
  features <- GRanges("chr1", IRanges(c(1, domainStart), 
                                      height=c(0.02, domainHeight)/2,
                                      width = c(wholeWidth, domainLength), 
                                      fill = c('gray', domainColor),
                                      color = c('gray', domainColor),
                                      names = c(geneName, domainName)))
  keepInd = which(snvCountDF$Counts > 1)
  if(type == 'AA'){
    SNP = snvCountDF$AA_pos[keepInd]
    SNP_name = snvCountDF$AA_mut[keepInd]
  }else{
    SNP = snvCountDF$NT_pos_rel[keepInd]
    SNP_name = snvCountDF$NT_mut[keepInd]
  }
  SNP_name[startsWith(SNP_name, 'S:')] = ''
  SNP_name %<>% str_remove_all(., 'N:')
  SNP.gr = GRanges("chr1", IRanges(SNP, width = 1, names = SNP_name))
  # SNP.gr$color = ifelse(snvCountDF$mutation_type[keepInd] == 'N', fillSyno[1], fillSyno[2])
  SNP.gr$color = ifelse(snvCountDF$class[keepInd] == 'SNS', fillClass[1], fillClass[2])
  SNP.gr$alpha = 0.5
  SNP.gr$dashline.col = NA
  SNP.gr$dashline.col[snvCountDF$class[keepInd] == 'SNS'] = colorSyno[1]
  SNP.gr$dashline.col[snvCountDF$class[keepInd] == 'SNV'] = colorSyno[1]
  SNP.gr$dashline.col[snvCountDF$mutation_type[keepInd] == 'S'] = NA
  SNP.gr$border = 'gray'
  SNP.gr$border[snvCountDF$class[keepInd] == 'SNS'] = colorSyno[1]
  SNP.gr$border[snvCountDF$class[keepInd] == 'SNV'] = colorSyno[1]
  SNP.gr$border[snvCountDF$mutation_type[keepInd] == 'S'] = 'gray'
  # SNP.gr$cex = snvCountDF$Score[keepInd]
  SNP.gr$cex = ifelse(snvCountDF$mutation_type[keepInd] == 'N', 0.7, 0.3)
  # SNP.gr$value1 = snvCountDF$pretermFreq[keepInd]
  # SNP.gr$value2 = 1 - snvCountDF$pretermFreq[keepInd]
  SNP.gr$label.parameter.rot = 60
  SNP.gr$score = snvCountDF$Freq[keepInd]
  # label.parameter.gp.SNS = list(gpar(col = colorClass[1])) # , fontsize = 7
  # label.parameter.gp.SNV = list(gpar(col = colorClass[2]))
  # SNP.gr$label.parameter.gp = ifelse(snvCountDF$class[keepInd] == 'SNS', label.parameter.gp.SNS, label.parameter.gp.SNV)
  SNP.gr$label.parameter.gp = gpar(fontsize = 5)
  SNP.gr$lwd = ifelse(snvCountDF$mutation_type[keepInd] == 'N', 1, 0.5)
  SNP.gr$SNPsideID = 'top'
  SNP.gr$SNPsideID[snvCountDF$class[keepInd] == 'SNV'] = "bottom"
  
  legend = list(
    list(labels=c('Non-synonymous', 'Synonymous'), cex = 1, fill= 'white', col = colorSyno), # fill = fillSyno
    list(labels=c('SNS', 'SNV'), cex = 1, fill = fillClass, col = 'white')
    # list(labels=c('SNS', 'SNV', 'con_SNV', 'pop_SNV'), fill = colorClass)
  )
  cairo_pdf(paste0('./plot/MTMG_bilol_mutation_', species, '_', geneID, '_', type, fileN, '.pdf'), 
            width = 9.65, height = pdfHei, onefile = T)
  par(mar=c(0.5, 0,0,.5))
  lolliplot(SNP.gr, features,
            ylab = 'Frequency of mutation\n', 
            type= "circle", legend = legend, label_on_feature = F,
            xaxis = xAxis[c(1,length(xAxis))],
            jitter = c("node", "label"),
            # xaxis.gp = gpar(lex = 1, lineheight = 0.7, fontsize=13, lwd=1.5),
            xaxis.gp = gpar(lex = 0, lineheight = 0, fontsize=13, lwd=0),
            yaxis = seq(0, 100, 25), 
            yaxis.gp = gpar(lex = 1, lineheight = 0.7, fontsize=13, lwd=1.5),
            ranges = GRanges("chr1", IRanges(1, wholeWidth)))
  
  filterDF = snvCountDF[snvCountDF$mutation_type == 'N', c('AA_mut', 'Counts')] %>% group_by(AA_mut) %>% 
    summarise(Counts = sum(Counts))
  filterDF$AA_pos = str_extract(filterDF$AA_mut, '\\d+') %>% as.numeric()
  countDF = filterDF[, c('AA_pos', 'Counts')] %>% group_by(AA_pos) %>% 
    summarise(Counts = max(Counts))
  countDF$Freq = countDF$Counts / length(unique(snvDF$sample))
  range(countDF$Freq)
  convScore$Counts = countDF$Counts[match(convScore$Pos, countDF$AA_pos)]
  convScore$Counts[is.na(convScore$Counts)] = 0
  hist(convScore$Counts)
  convScore$Freq = countDF$Freq[match(convScore$Pos, countDF$AA_pos)]
  convScore$Freq[is.na(convScore$Freq)] = 0
  
  hist(convScore$Freq)
  # b = snvCountDF[snvCountDF$Counts > 1 & snvCountDF$AA_pos %in% c(seq(1476, 1476 + 79), seq(1597, 1597 + 79), seq(1718, 1718 + 79),
  #                                                                  seq(1839, 1839 + 79), seq(1960, 1960 + 79), seq(2081, 2081 + 79)), ]
  # c = snvDF[snvDF$mutation %in% b$mutation, ]
  # c$SampleID = swabData$SampleID.u[match(c$seq_ID, swabData$seq_ID)]
  # pheatmap(table(c$mutation, c$SampleID))
  # convScore$Group = 'WT'
  # convScore$Group[convScore$Freq > 0] = 'Mut'
  # table(convScore$Group)
  # convScore$Group %<>% factor(., levels = c('WT', 'Mut'))
  # table(convScore$Group)
  convScore$Group = 'WT'
  convScore$Group[convScore$Counts > 1] = 'Mut'
  table(convScore$Group)
  convScore$Group %<>% factor(., levels = c('WT', 'Mut'))
  table(convScore$Group)
  
  p1 = ggplot(convScore, aes(x = Pos, y = Score)) +
    geom_line(color = 'gray', linewidth = 0.3) +
    geom_point(shape = 21, color = '#377EB8', fill = '#377EB8', alpha = 0.5) +
    labs(x = 'Position', y = 'Conservation score') +
    scale_x_continuous(expand = expansion(mult = c(0,0))) +
    mytheme + theme(aspect.ratio = 1/5)
  p2 = ggplot(convScore, aes(x = Pos, y = Freq)) +
    geom_line(color = 'gray', linewidth = 0.3) +
    geom_point(shape = 21, color = '#377EB8', fill = '#377EB8', alpha = 0.5) +
    labs(x = 'Position', y = 'Frequency of mutation') +
    scale_x_continuous(expand = expansion(mult = c(0,0))) +
    mytheme + theme(aspect.ratio = 1/5)
  print(ggarrange(p1, p2, nrow = 2, ncol = 1))
  p3 = ggplot(convScore, aes(x = Freq, y = Score)) + 
    geom_point(shape = 21, color = '#377EB8', fill = '#377EB8', alpha = 0.5) +
    labs(x = 'Count of mutation', y = 'Conservation score') +
    mytheme + theme(aspect.ratio = 1)
  snvCountDF$Score = convScore$Score[match(snvCountDF$AA_pos, convScore$Pos)]
  p4 = ggplot(snvCountDF, aes(x = Freq, y = Score)) + 
    geom_point(shape = 21, aes(color = class, fill = class), alpha = 0.5) +
    labs(x = 'Frequency of mutation', y = 'Conservation score') +
    mytheme + theme(aspect.ratio = 1)
  print(p3 | p4)

  p5 = ggboxplot(convScore,x = 'Group', y = 'Score', color = 'Group',add = 'jitter', add.params = list(size = 0.3)) +
    labs(x = 'Group', y = 'Conservation score') +
    stat_compare_means(comparisons = list(c('WT', 'Mut')), 
                       method = 'wilcox.test') +
    scale_color_manual(values = c('#7FC97F', '#BEAED4')) + 
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    mytheme3 + theme(aspect.ratio = 1.5/1, legend.position = 'top',
                     legend.margin = margin(0,0,-10,0))
  print(p5)
  ht = Heatmap( t(convScore[, 'Score'] ), name = 'Conservation score', 
                col = colorRampPalette(brewer.pal(8, "PuBuGn") %>% rev())(100),
                heatmap_legend_param = list(legend_width = unit(3.5, "cm"), direction = "horizontal", title_position = 'topcenter'),
                cluster_columns = F, cluster_rows = F, height = unit(0.6, 'cm'))
  draw(ht,heatmap_legend_side = 'bottom')
  dev.off()
}
##### L. crisptus, Muc_B2, NZ_CP039266.1_1594 #####
if(T){ # NT
  library(trackViewer)
  species = 'Lactobacillus_crispatus'
  geneID = 'NZ_CP039266.1_1594'
  GeneDirection = geneDF$direction[geneDF$gene == geneID][1]
  geneName = 'Muc_B2'
  snvDF = read.csv(paste0('./instrain_res/MTMG/gene/s__', species, '.fasta_', geneID, '.csv'))
  
  domainStart = c(2, 32, 171) *3
  domainLength = c(44, 115, 61) *3
  domainHeight = c(0.07, 0.05, 0.07)
  domainColor = c("#F99300", "#A2D1E1", "#F99300")
  domainName = c('Mub B2-like domain', 
                 'PBP2_NikA_DppA_OppA_like', # 	The substrate-binding domain of an ABC-type nickel/oligopeptide-like import system contains the type 2 periplasmic binding fold
                 'Mub B2-like domain')
  geneLength = 1239
  scaffoldStart = 1644128
  
  fillSyno = c('red', 'gray')
  colorSyno = c('#E5C494', 'gray')
  fillClass = c('#8EC80F', '#7BC5FB')
  colorClass = c('#8EC80F', '#7BC5FB', '#1F78B4', 'red')

  # excludePos = seq(337, 418, 1)
  excludePos = c()
  xAxis = c(1, seq(100, 1100, 100), geneLength)
  lolPlot(type = 'NT')
  lolPlot(type = 'NT', keep_mut_on_domain = T)
}
if(T){ # AA
  library(trackViewer)
  species = 'Lactobacillus_crispatus'
  geneID = 'NZ_CP039266.1_1594'
  GeneDirection = geneDF$direction[geneDF$gene == geneID][1]
  geneName = 'Muc_B2'
  snvDF = read.csv(paste0('./instrain_res/MTMG/gene/s__', species, '.fasta_', geneID, '.csv'))
  
  domainStart = c(2, 32, 171)
  domainLength = c(44, 115, 61)
  domainHeight = c(0.07, 0.05, 0.07)
  domainColor = c("#F99300", "#A2D1E1", "#F99300")
  domainName = c('Mub B2-like domain', 
                 'PBP2_NikA_DppA_OppA_like', # 	The substrate-binding domain of an ABC-type nickel/oligopeptide-like import system contains the type 2 periplasmic binding fold
                 'Mub B2-like domain')
  geneLength = 1239
  scaffoldStart = 1644128
  
  fillSyno = c('red', 'gray')
  colorSyno = c('#E5C494', 'gray')
  fillClass = c('#8EC80F', '#7BC5FB')
  colorClass = c('#8EC80F', '#7BC5FB', '#1F78B4', 'red')
  
  # excludePos = seq(337, 418, 1)
  excludePos = c()
  xAxis = c(1, seq(50, 350, 50), geneLength/3)
  # lolPlot(type = 'AA')
  # lolPlot(type = 'AA', keep_mut_on_domain = T)
  convScore = read.table('./colabfold/1594_conserv.txt', header = T)
  lolPlot2(type = 'AA', keep_mut_on_domain = F, pdfHei = 6.43)
}
##### L. crisptus, Gram_pos_anchor, NZ_CP039266.1_48 #####
if(T){ # NT
  library(trackViewer)
  species = 'Lactobacillus_crispatus'
  geneID = 'NZ_CP039266.1_48'
  GeneDirection = geneDF$direction[geneDF$gene == geneID][1]
  geneName = 'Gram_pos_anchor'
  snvDF = read.csv(paste0('./instrain_res/s__', species, '.fasta_', geneID, '.csv'))
  
  domainStart = c(52)*3
  domainLength = c(92)*3
  domainHeight = c(0.06)
  domainColor = c("#D49119")
  
  geneLength = 612
  scaffoldStart = 48780
  
  fillSyno = c('red', 'gray')
  colorSyno = c('#E5C494', 'gray')
  fillClass = c('#8EC80F', '#7BC5FB')
  colorClass = c('#8EC80F', '#7BC5FB', '#1F78B4', 'red')
  domainName = c('PRK09418')
  domainName = c("Bifunctional 2',3'-cyclic-nucleotide 2'-phosphodiesterase/3'-nucleotidase")
  excludePos = c(seq(49079, 49105, 1), seq(49215, 49248, 1)) - scaffoldStart
  excludePos =c()
  xAxis = c(1, seq(100,500,100), geneLength)
  lolPlot(type = 'NT')
  lolPlot(type = 'NT', keep_mut_on_domain = T)
}
if(T){ # AA
  library(trackViewer)
  species = 'Lactobacillus_crispatus'
  geneID = 'NZ_CP039266.1_48'
  GeneDirection = geneDF$direction[geneDF$gene == geneID][1]
  geneName = 'Gram_pos_anchor'
  snvDF = read.csv(paste0('./instrain_res/MTMG/gene/s__', species, '.fasta_', geneID, '.csv'))
  
  domainStart = c(52)
  domainLength = c(92)
  domainHeight = c(0.06)
  domainColor = c("#D49119")
  
  geneLength = 612
  scaffoldStart = 48780
  
  fillSyno = c('red', 'gray')
  colorSyno = c('#E5C494', 'gray')
  fillClass = c('#8EC80F', '#7BC5FB')
  colorClass = c('#8EC80F', '#7BC5FB', '#1F78B4', 'red')
  domainName = c('PRK09418')
  domainName = c("Bifunctional 2',3'-cyclic-nucleotide 2'-phosphodiesterase/3'-nucleotidase")
  # excludePos = c(seq(49079, 49105, 1), seq(49215, 49248, 1)) - scaffoldStart
  excludePos =c()
  xAxis = c(1, seq(25,175,25), geneLength/3)
  # lolPlot(type = 'AA')
  # lolPlot(type = 'AA', keep_mut_on_domain = T)
  convScore = read.table('./colabfold/48_conserv.txt', header = T)
  lolPlot2(type = 'AA', keep_mut_on_domain = F, pdfHei = 6.43)
}

##### L. iners, Gram_pos_anchor, NZ_CP045664.1_1 #####
if(T){ # NT
  library(trackViewer)
  species = 'Lactobacillus_iners'
  geneID = 'NZ_CP045664.1_1'
  GeneDirection = geneDF$direction[geneDF$gene == geneID][1]
  geneName = 'Gram_pos_anchor'
  snvDF = read.csv(paste0('./instrain_res/MTMG/gene/s__', species, '.fasta_', geneID, '.csv'))
  
  domainStart = c(49)*3
  domainLength = c(92)*3
  domainHeight = c(0.06)
  domainColor = c("#F3685E")
  
  geneLength = 276
  scaffoldStart = 3
  
  fillSyno = c('red', 'gray')
  colorSyno = c('#E5C494', 'gray')
  fillClass = c('#8EC80F', '#7BC5FB')
  colorClass = c('#8EC80F', '#7BC5FB', '#1F78B4', 'red')
  domainName = c("Gram_pos_anchor super family")
  domainName = c('Gram_pos_anchor')
  domainName = c("LPXTG cell wall anchor motif")
  excludePos =c()
  xAxis = c(1, seq(25,250,25), geneLength)
  lolPlot(type = 'NT')
  lolPlot(type = 'NT', keep_mut_on_domain = T)
}
if(T){ # AA
  library(trackViewer)
  species = 'Lactobacillus_iners'
  geneID = 'NZ_CP045664.1_1'
  GeneDirection = geneDF$direction[geneDF$gene == geneID][1]
  geneName = 'Gram_pos_anchor'
  snvDF = read.csv(paste0('./instrain_res/MTMG/gene/s__', species, '.fasta_', geneID, '.csv'))
  
  domainStart = c(49)
  domainLength = c(92)
  domainHeight = c(0.06)
  domainColor = c("#F3685E")
  
  geneLength = 276
  scaffoldStart = 3
  
  fillSyno = c('red', 'gray')
  colorSyno = c('#E5C494', 'gray')
  fillClass = c('#8EC80F', '#7BC5FB')
  colorClass = c('#8EC80F', '#7BC5FB', '#1F78B4', 'red')
  domainName = c("Gram_pos_anchor super family")
  domainName = c('Gram_pos_anchor')
  domainName = c("LPXTG cell wall anchor motif")
  excludePos =c()
  xAxis = c(1, seq(10,80,10), geneLength/3)
  lolPlot(type = 'AA')
  lolPlot(type = 'AA', keep_mut_on_domain = T)
  convScore = read.table('./colabfold/1_conserv.txt', header = T)
  lolPlot2(type = 'AA', keep_mut_on_domain = F, pdfHei = 6.43)
}
##### L. jensenii, Muc_B2, NZ_CP018809.1_5 #####
if(T){ # NT
  library(trackViewer)
  species = 'Lactobacillus_jensenii'
  geneID = 'NZ_CP018809.1_5'
  GeneDirection = geneDF$direction[geneDF$gene == geneID][1]
  geneName = 'Muc_B2'
  snvDF = read.csv(paste0('./instrain_res/MTMG/gene/s__', species, '.fasta_', geneID, '.csv'))
  
  domainStart = c(12, 865, 967, 1476, 1597, 1718, 1839, 1960, 2081, 2198, 2292, 2394)*3
  domainLength = c(42, 96, 98, rep(79, 6), 70, 72, 44)*3
  domainName = c("YSIRK_signal", 'RibLong', 'YPDG_rpt', rep('Mub_B2', 7), 'Rib', 'Gram_pos_anchor')
  domainName = c("Gram-positive signal peptide, YSIRK family", 'Long Rib domain', 'YPDG domain', 
                 rep('Mub B2-like domain', 7), 'Rib/alpha-like repeat', 'LPXTG cell wall anchor motif')
  domainHeight = rep(0.06, 12)
  domainColor = c('#F781BF', '#BC80BD', '#33A02B', rep("#F99300", 7), '#6A3D9A', '#B15928')

  geneLength = 7314
  scaffoldStart = 4465
  
  fillSyno = c('red', 'gray')
  colorSyno = c('#E5C494', 'gray')
  fillClass = c('#8EC80F', '#7BC5FB')
  colorClass = c('#8EC80F', '#7BC5FB', '#1F78B4', 'red')

  excludePos =c()
  xAxis = c(1, seq(500,7000,500), geneLength)
  lolPlot(type = 'NT')
  lolPlot(type = 'NT', keep_mut_on_domain = T)
}
if(T){ # AA
  library(trackViewer)
  species = 'Lactobacillus_jensenii'
  geneID = 'NZ_CP018809.1_5'
  GeneDirection = geneDF$direction[geneDF$gene == geneID][1]
  geneName = 'Muc_B2'
  snvDF = read.csv(paste0('./instrain_res/MTMG/gene/s__', species, '.fasta_', geneID, '.csv'))
  
  domainStart = c(12, 865, 967, 1476, 1597, 1718, 1839, 1960, 2081, 2198, 2292, 2394)
  domainLength = c(42, 96, 98, rep(79, 6), 70, 72, 44)
  domainName = c("YSIRK_signal", 'RibLong', 'YPDG_rpt', rep('Mub_B2', 7), 'Rib', 'Gram_pos_anchor')
  domainName = c("Gram-positive signal peptide, YSIRK family", 'Long Rib domain', 'YPDG domain', 
                 rep('Mub B2-like domain', 7), 'Rib/alpha-like repeat', 'LPXTG cell wall anchor motif')
  domainHeight = rep(0.06, 12)
  domainColor = c('#F781BF', '#BC80BD', '#33A02B', rep("#F99300", 7), '#6A3D9A', '#B15928')
  
  geneLength = 7314
  scaffoldStart = 4465
  
  fillSyno = c('red', 'gray')
  colorSyno = c('#E5C494', 'gray')
  fillClass = c('#8EC80F', '#7BC5FB')
  colorClass = c('#8EC80F', '#7BC5FB', '#1F78B4', 'red')
  
  excludePos =c()
  xAxis = c(1, seq(200, 2200, 200), geneLength/3)
  # lolPlot(type = 'AA')
  # lolPlot(type = 'AA', keep_mut_on_domain = T)
  # 
  convScore = read.table('./colabfold/5_conserv.txt', header = T)
  lolPlot2(type = 'AA', keep_mut_on_domain = F, pdfHei = 6.43)
}
##### trackViewer L. iners, MucBP, NZ_CP045664.1_439 #####
if(F){
  snvDF = read.csv('./instrain_res/MTMG/gene/s__Lactobacillus_iners.fasta_NZ_CP045664.1_439.csv')
  snvDF$gene_pos <- str_extract(snvDF$mutation, '\\d+') %>% as.numeric()
  snvDF$seq_ID = str_extract(snvDF$sample, '\\d+_\\d+')
  snvCountDF = read.csv('./instrain_res/MTMG/gene/s__Lactobacillus_iners.fasta_NZ_CP045664.1_439_stat.csv')
  snvCountDF$NT_pos = snvDF$position[match(snvCountDF$mutation, snvDF$mutation)]
  snvCountDF$gene_pos <- str_extract(snvCountDF$mutation, '\\d+') %>% as.numeric()
  range(snvCountDF$gene_pos)
  snvCountDF$mutation_type = snvDF$mutation_type[match(snvCountDF$mutation, snvDF$mutation)]
  library(trackViewer)
  features <- GRanges("chr1", IRanges(c(908, 1016, 1124, 1232, 1340,1388)*3, 
                                      width=c(72, 72, 72, 72, 72, 112)*3,
                                      names= c('MucBP', 'MucBP', 'MucBP', 'MucBP', 'MucBP', 'Trypan_PARP')))
  features$fill = c("#FF8833", "#FF8833", "#FF8833", "#FF8833", "#FF8833", "#DFA32D")
  keepInd = which(snvCountDF$Counts > 50)
  SNP = snvCountDF$gene_pos[keepInd]
  SNP_name = snvCountDF$mutation[keepInd]
  length(SNP)
  sample.gr <- GRanges("chr1", IRanges(SNP, width=1, names=SNP_name))
  sample.gr$color = ifelse(snvCountDF$mutation_type[keepInd] == 'N', 'red', 'blue')
  sample.gr$node.label <- snvCountDF$Counts[keepInd]
  sample.gr$score <- snvCountDF$Counts[keepInd]
  sample.gr$value1 = snvCountDF$pretermFreq[keepInd]
  sample.gr$value2 = 1 - snvCountDF$pretermFreq[keepInd]
  legend <- list(labels=c('Non-synonymous', 'Synonymous'), 
                 col = 'black', fill=c('red', 'blue'))
  lolliplot(sample.gr, features, legend = legend,
            ranges = GRanges("chr1", IRanges(1, 4554)))
}

##### trackViewer Lactobacillus_crispatus, DDE_Tnp_1_3, NZ_CP039266.1_1431 #####
if(F){
  snvDF = read.csv('./instrain_res/MTMG/gene/s__Lactobacillus_crispatus.fasta_NZ_CP039266.1_1431.csv')
  snvDF$gene_pos <- str_extract(snvDF$mutation, '\\d+') %>% as.numeric()
  snvDF$seq_ID = str_extract(snvDF$sample, '\\d+_\\d+')
  snvCountDF = read.csv('./instrain_res/MTMG/gene/s__Lactobacillus_crispatus.fasta_NZ_CP039266.1_1431_stat.csv')
  snvCountDF$gene_pos <- str_extract(snvCountDF$mutation, '\\d+') %>% as.numeric()
  snvCountDF$NT_pos = snvDF$position[match(snvCountDF$mutation, snvDF$mutation)]
  snvCountDF$NT_pos_rel = snvCountDF$NT_pos
  snvCountDF$mutation_type = snvDF$mutation_type[match(snvCountDF$mutation, snvDF$mutation)]
  library(trackViewer)
  
  features <- GRanges("chr1", IRanges(c(1)*3, 
                                      width=c(108)*3,
                                      names= c('Transposase')))
  features$fill = c("#FF8833")
  keepInd = which(snvCountDF$Counts > 10)
  SNP = snvCountDF$gene_pos[keepInd]
  SNP_name = snvCountDF$mutation[keepInd]
  length(SNP)
  sample.gr <- GRanges("chr1", IRanges(SNP, width=1, names=SNP_name))
  sample.gr$color = ifelse(snvCountDF$mutation_type[keepInd] == 'N', 'red', 'blue')
  sample.gr$node.label <- snvCountDF$Counts[keepInd]
  sample.gr$score <- snvCountDF$Counts[keepInd]
  sample.gr$value1 = snvCountDF$pretermFreq[keepInd]
  sample.gr$value2 = 1 - snvCountDF$pretermFreq[keepInd]
  legend <- list(labels=c('Non-synonymous', 'Synonymous'), 
                 col= 'black', 
                 fill=c('red', 'blue'))
  lolliplot(sample.gr, features, legend = legend,
            ranges = GRanges("chr1", IRanges(1,324)))
}
