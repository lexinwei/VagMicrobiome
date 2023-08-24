source('~/workspace/library_function.R')
load('./data/metaData.RData')
load('./data/swabData.RData')

##### load ref/mags sample-wise profile #####
selCol = c("genome", "coverage", "breadth", "nucl_diversity", "length",
           "true_scaffolds", "detected_scaffolds", "coverage_median",
           "coverage_std", "coverage_SEM", "breadth_minCov",
           "breadth_expected", "nucl_diversity_rarefied",
           "conANI_reference", "popANI_reference", "SNS_count",
           "SNV_count", "reads_mean_PID", "reads_unfiltered_reads")
##### ref
refDF = data.frame()
for(ss in metaData$pt_ID){
  tmpDat = read.delim2(paste0('./instrain_res/ref_coReads_profile/genome_info/', ss, '.IS_genome_info.tsv'))
  tmpDat = tmpDat[, selCol]
  tmpDat$breadth_minCov %<>% as.numeric()
  refDF_tmp = tmpDat[tmpDat$breadth_minCov >= 0.1, ]
  if(nrow(refDF_tmp) > 0){
    refDF_tmp$sample = ss
    refDF = rbind(refDF, refDF_tmp)
  }
}
refDF$pt_ID = refDF$sample
refDF$nucl_diversity %<>% as.numeric()
refDF$SNV_count %<>% as.numeric()
refDF$pt_ID %<>% factor(., levels = mixedsort(refDF$pt_ID %>% unique))
refDF$species = str_replace(refDF$genome, '.fna', '')
table(refDF$species)
refDF$type = 'Reference'

#### mags 
magDF = data.frame()
for(ss in metaData$pt_ID){
  tmpDat = read.delim2(paste0('./instrain_res/mags_coReads_profile/genome_info/', ss, '.IS_genome_info.tsv'))
  tmpDat = tmpDat[, selCol]
  tmpDat$breadth_minCov %<>% as.numeric()
  magDF_tmp = tmpDat[tmpDat$breadth_minCov >= 0.1, ]
  if(nrow(magDF_tmp) > 0){
    magDF_tmp$sample = ss
    magDF = rbind(magDF, magDF_tmp)
  }
}
magDF$pt_ID = magDF$sample
magDF$nucl_diversity %<>% as.numeric()
magDF$SNV_count %<>% as.numeric()
magDF$pt_ID %<>% factor(., levels = mixedsort(magDF$pt_ID %>% unique))
magDF$genome %<>% str_replace(., '.fna', '')
spName = read.xlsx('./tax_of_mags_co.xlsx')
magDF$species = spName$lca_species[match(magDF$genome, spName$Bin)]
magDF$species = spName$species[match(magDF$genome, spName$Bin)]
magDF = magDF[!is.na(magDF$species), ]
magDF$type = 'MAG'

##### plot nucleotide diversity/coverage/breath ####
impSP = intersect(magDF$species, refDF$species)
allDF = rbind(magDF[magDF$species %in% impSP, ], refDF[refDF$species %in% impSP, ])
allDF$length %<>% as.numeric()
allDF$SNPs_per_bp = allDF$SNV_count/allDF$length
alllVarDF = allDF[, c('species', 'pt_ID', 'coverage', 'breadth', 'breadth_minCov', 'SNPs_per_bp', 'nucl_diversity', 'type')]
alllVarDF = melt(alllVarDF, id.vars = c('species','pt_ID', 'type'))
alllVarDF$value %<>% as.numeric()
if(F){
  pdf('./plot/ref_and_mags_genome_wise_nucl_div_coReads.pdf', width = 13, height = 4)
  for (sp in unique(allDF$species)) {
    tmpDF = allDF[allDF$species == sp, ]
    p = ggplot(tmpDF, aes(x = pt_ID, y = nucl_diversity)) +
      geom_line(aes(group=interaction(pt_ID, type), color = pt_ID, linetype = type), size = 0.3) + 
      geom_point(aes(fill = pt_ID, shape = type), color = 'black', stroke = 0.4, alpha = 0.5) +
      scale_shape_manual(values = c(21, 24)) +
      mytheme + theme(axis.text.x = element_text(angle = 45, size = 5, hjust = 1)) +
      labs(title = paste0(sp), y = paste0('Nucleotide diversity'), x = '')
    print(p)
  }
  dev.off()
  
  pdf('./plot/ref_and_mags_genome_info_coReads.pdf', width = 15, height = 6)
  for (sp in unique(alllVarDF$species)) {
    tmpDF = alllVarDF[alllVarDF$species == sp, ]
    p = ggplot(tmpDF, aes(x = pt_ID, y = value)) +
      geom_line(aes(group=interaction(pt_ID, type), color = pt_ID, linetype = type), size = 0.3) + 
      geom_point(aes(fill = pt_ID, shape = type), color = 'black', stroke = 0.4, alpha = 0.5) +
      scale_shape_manual(values = c(21, 24)) +
      mytheme2 + theme(axis.text.x = element_text(angle = 45, size = 5, hjust = 1),
                       axis.text.y = element_text(size = 10, hjust = 1),
                       strip.text = element_text(size = 10), 
                       panel.grid = element_line(size = 0.3), panel.background = element_rect(fill = '#EEEEEE')) +
      labs(title = paste0(sp), x = 'Sample', y = '') +
      facet_wrap(.~variable, nrow = 5, ncol = 1, scales = 'free_y', strip.position = 'right')
    print(p)
  }
  dev.off()
  
}






























