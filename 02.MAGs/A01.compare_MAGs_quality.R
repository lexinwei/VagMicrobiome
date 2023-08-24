##### assembly individually, bin across samples for each person ####
##### + assembly coReads by person, bin across samples for each person ####
load('./data/Rdata/colorList.RData')
load('./data/RData/swabData.RData')
source('./script/library_function.R')
library(ggExtra)
#### quality plot #####
if(T){
  statsDir = './MAGs/stats_cross/'
  magsQual = data.frame()
  for(f in list.files(statsDir)){
    tmpDF = read.delim(paste0(statsDir, f), sep = '\t')
    tmpDF$bin_name = paste0(str_remove(f, 'metawrap_50_10_bins.stats'), tmpDF$bin, '.fa')
    magsQual = rbind(magsQual, tmpDF)
  }
  range(magsQual$completeness)
  range(magsQual$contamination)
  save(magsQual, file = './data/RData/magsQual.RData')
  drepBins = read.csv('~/workspace/vag/MAGs/Widb.csv') # /share/home/jianglab/weixin/workspace/vag/drep/cross_Person_ani99_2/data_tables/Widb.csv
  drepBins$sample = str_extract(drepBins$genome, '\\d_\\d+')
  drepBins$pt_ID.u = swabData$pt_ID.u[match(drepBins$sample, swabData$seq_ID)]
  drepBins2 = read.csv('~/workspace/vag/MAGs/Widb_ani95.csv') # /share/home/jianglab/weixin/workspace/vag/drep/cross_Person_ani99_2_ani95/data_tables/Widb.csv
  drepBins2$sample = str_extract(drepBins2$genome, '\\d_\\d+')
  drepBins2$pt_ID.u = swabData$pt_ID.u[match(drepBins2$sample, swabData$seq_ID)]
  drepBins2$name = paste0('Bin.', 1:nrow(drepBins2))
  bin2name = drepBins2[, c('genome', 'name')]
  save(bin2name, file = './data/RData/bin2name.RData')
  intersect(drepBins2$genome, drepBins$genome)
  
  magsQual = magsQual[magsQual$bin_name %in% drepBins$genome, ]
  magsQual$Quality = 'Filter out'
  range(magsQual$completeness)
  range(magsQual$contamination)
  range(drepBins$completeness)
  range(drepBins$contamination)
  magsQual$Quality[magsQual$completeness < 50 & magsQual$contamination < 10] = 'Low'
  magsQual$Quality[magsQual$completeness >= 50 & magsQual$contamination < 10] = 'Medium'
  magsQual$Quality[magsQual$completeness > 90 & magsQual$contamination < 5] = 'High'
  table(magsQual$Quality)
  magsQual$Quality %<>% factor(., levels = names(colorQuality))
  magsQual$score = drepBins$score[match(magsQual$bin_name, drepBins$genome)]
  magsQual = magsQual[!startsWith(magsQual$bin_name, 'NC'), ]
  piris <- ggplot(magsQual, aes(completeness, contamination)) +
    geom_point(aes(size = score, fill = Quality, color = Quality), stroke = 0.2, shape = 21, alpha = 0.5) +
    labs(x = 'Completeness (%)', y = 'Contamination (%)', size = 'dRep score') +
    scale_size_continuous(range = c(1, 2.5)) +
    scale_color_manual(values = colorQuality) +
    scale_fill_manual(values = colorQuality) +
    mytheme + theme(legend.position = 'left')
  cairo_pdf('./plot/MAGs_quality2.pdf', width = 4.91, height = 3.34)
  ggMarginal(piris, fill = '#F1F1F1')
  dev.off()
}
nrow(magsQual)
save(magsQual, file = './data/RData/magsQual.RData')
# /share/home/jianglab/weixin/workspace/vag/drep/cross_Person_ani99_2_ani95 
##### species-level representitive genome #####
ani95DF = read.delim('./MAGs/SRG/Cdb.csv', sep = ',')
srgDF = read.delim('./MAGs/SRG/Wdb.csv', sep = ',')
# lineage
lineageAr = read.table('./MAGs/gtdbtk.ar53.summary.tsv', sep = '\t', header = T)
lineageBac = read.table('./MAGs/gtdbtk.bac120.summary.tsv', sep = '\t', header = T)
lineage = rbind(lineageBac, lineageAr)
bin2rank = data.frame('ID' = lineage$user_genome, str_split_fixed(lineage$classification, ';', 7))
bin2rank = apply(bin2rank, 2, FUN = function(x){str_replace_all(x, '^[dpcofgs]__', '')}) %>% as.data.frame()
colnames(bin2rank) = c('ID', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
table(bin2rank$phylum) %>% length()
save(bin2rank, file = './data/RData/bin2rank.RData')
ani95DF$SRG = 'No'
ani95DF$SRG[ani95DF$genome %in% srgDF$genome] = 'Yes'
ani95DF$Phylum = bin2rank$phylum[match(ani95DF$genome, paste0(bin2rank$ID, '.fa'))]
ani95DF$Species = bin2rank$species[match(ani95DF$genome, paste0(bin2rank$ID, '.fa'))]
ani95DF$Genus = bin2rank$genus[match(ani95DF$genome, paste0(bin2rank$ID, '.fa'))]
table(ani95DF$SRG)

srgDF$SRG_type = 'uSRG'
for(b in srgDF$genome){
  cc = ani95DF$secondary_cluster[ani95DF$genome == b]
  tmp = ani95DF[ani95DF$secondary_cluster == cc, ]
  srgDF$N_of_constructed_genomes[srgDF$genome == b] = nrow(tmp)
  if(sum(tmp$Species != "") > 0){
    srgDF$SRG_type[srgDF$genome == b] = 'kSRG'
  }
}
table(srgDF$N_of_constructed_genomes)
sum(table(srgDF$N_of_constructed_genomes))
table(srgDF$N_of_constructed_genomes, srgDF$SRG_type) %>% colSums()

srgDF$Phylum = ani95DF$Phylum[match(srgDF$genome, ani95DF$genome)]
setdiff(names(table(ani95DF$Phylum)), names(table(srgDF$Phylum )))
table(srgDF$Phylum) %>% length()
table(srgDF$SRG_type)
Phylum2SRGTab = table(srgDF$Phylum, srgDF$SRG_type)
Phylum2SRG = data.frame(
  Phylum = row.names(Phylum2SRGTab),
  kSRG = Phylum2SRGTab[, 1],
  uSRG = Phylum2SRGTab[, 2]
)
sum(Phylum2SRG$kSRG)
sum(Phylum2SRG$uSRG)
sum(Phylum2SRG$uSRG)/sum(Phylum2SRG$SRG)
Phylum2SRG$SRG = rowSums(Phylum2SRG[, -1])
ani95DF =ani95DF[!startsWith(ani95DF$genome, 'NC'), ]
Phylum2SRG$total = table(ani95DF$Phylum)[match(Phylum2SRG$Phylum, names(table(ani95DF$Phylum)))]
sum(Phylum2SRG$total )
Phylum2SRG$Phylum %<>% factor(., levels = c(names(colorPhylumMAGs), 'Thermoproteota'))
blankTheme = mytheme + theme(axis.title.y = element_blank(),
                             panel.border = element_blank(), panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(), 
                             axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), 
                             legend.position = 'none')
p1 <- ggplot(Phylum2SRG, aes(Phylum, uSRG/SRG, fill = Phylum)) + 
  geom_col(color = 'black', width = 0.8, size = 0.3, alpha = 0.7) + 
  geom_text(aes(label = signif(100*(uSRG/SRG), 3)), vjust = -0.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.3), labels = scales::percent) +
  scale_fill_manual(values = c(colorPhylumMAGs, 'Thermoproteota' = '#8DA0CB')) +
  xlab('Percent of uSGB in each Phylum (uSRG/SRG, %)') +
  blankTheme + theme(axis.text.x = element_blank(), plot.margin = unit(c(3,3,0,3), 'mm'))
p2 <- ggplot(Phylum2SRG,aes(Phylum, SRG, fill = Phylum)) + 
  geom_col(color = 'black', width = 0.8, size = 0.3, alpha = 0.7) + 
  geom_text(aes(label = SRG), vjust = -0.5) +
  scale_y_continuous(expand = expansion(mult = c(0.01,0.08)), limits = c(0, 70)) +
  scale_fill_manual(values = c(colorPhylumMAGs, 'Thermoproteota' = '#8DA0CB')) +
  xlab('Number of SRG in each Phylum') +
  blankTheme + theme(axis.text.x = element_blank(), plot.margin = unit(c(0,3,0,3), 'mm'))
p3<- ggplot(Phylum2SRG,aes(Phylum, total, fill = Phylum)) + 
  geom_col(color = 'black', width = 0.8, size = 0.3, alpha = 0.7) + 
  geom_text(aes(label = total), vjust = -0.5) +
  scale_y_continuous(expand = expansion(mult = c(0.01,0.05)), limits = c(0, 300)) +
  scale_fill_manual(values = c(colorPhylumMAGs, 'Thermoproteota' = '#8DA0CB')) +
  xlab('Number of MAG in each Phylum') +
  blankTheme + theme(axis.text.x = element_blank(), plot.margin = unit(c(0,3,0,3), 'mm'))

p4<- ggplot(Phylum2SRG,aes(Phylum, 0.000005, color = Phylum)) + geom_point(size = 4) +
  scale_y_continuous(limits = c(0, 0.00001), expand = c(0, 0)) +
  scale_color_manual(values = c(colorPhylumMAGs, 'Thermoproteota' = '#8DA0CB')) +
  blankTheme + theme(axis.title = element_blank(), axis.text.x = element_text(angle=90, size = 13, hjust = 1, vjust = 0.5), 
                     plot.margin = unit(c(2,3,0,3), 'mm'))

pdf('./plot/SRG2.pdf', width = 4.04, height = 5.29)
egg::ggarrange(p1, p2, p3, p4, heights = c(0.46, 0.46, 0.46,0.08))
dev.off()


