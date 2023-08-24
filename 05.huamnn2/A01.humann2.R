source('./script/library_function.R')
load('./data/RData/swabData.RData')
load('./data/RData/metaData.RData')
load('./data/RData/sampleVagitypes.RData')
load('./data/RData/colorList.RData')
##### pathabundance ####
if(T){
  datDF = read.delim("./humann/humann_pathabundance_relab.tsv", sep = '\t', check.names = F)
  colnames(datDF) %<>% str_remove('Ming_nova_VS_') %>% str_remove('_bwa.*$')
  datDF = datDF[, c(1, which(colnames(datDF) %in% swabData$seq_ID[!is.na(swabData$trimester)]))]
  datDF = datDF[rowSums(datDF[, -1]) > 0, ]
  dim(datDF)
}
if(T){ # not consider the taxonomy of pathway
  datDF_pw0 = datDF[which(!grepl('s__', datDF[, 1]) & !grepl('unclassified$', datDF[, 1])), ]
  datDF_pw1 = data.frame(Pathway = datDF_pw0[, 1], datDF_pw0[,-1] * 1e6, 
                        check.names = F)
  colSums(datDF_pw1[, 2:5])
  range(datDF_pw1[, -1])
  datDF_pw2 = datDF_pw1[rowSums(datDF_pw1[, -1] > 0.1) > 0.05 * ncol(datDF_pw1[, -1]),]
  datDF_pw3 = data.frame(Pathway = datDF_pw2[, 1], log10(datDF_pw2[, -1] + 0.001), check.names = F)
  annoDF = data.frame(row.names = colnames(datDF_pw3[, -1]))
  annoDF$pt_ID.u = swabData$pt_ID.u[match(row.names(annoDF), swabData$seq_ID)]
  annoDF$pt_ID = swabData$pt_ID[match(row.names(annoDF), swabData$seq_ID)]
  annoDF$Ethnicity2 = metaData$Ethnicity2[match(annoDF$pt_ID, metaData$pt_ID)]
  annoDF$Term2 = metaData$Term_char2[match(annoDF$pt_ID, metaData$pt_ID)]
  annoDF$Term = metaData$Term_char[match(annoDF$pt_ID, metaData$pt_ID)]
  annoDF$CST = sampleVagitypes$vagitypes_Susan[match(annoDF$pt_ID.u, sampleVagitypes$pt_ID.u)]
  annoDF$Vagitype = sampleVagitypes$vagitypes_Fettweis2[match(annoDF$pt_ID.u, sampleVagitypes$pt_ID.u)]
  annoCol = annoDF[, c('Term2', 'Term', 'Ethnicity2', 'CST', 'Vagitype')]
  pheatmap(datDF_pw3[3:nrow(datDF_pw3), -1], scale = 'none', 
           show_rownames = F, show_colnames = F, annotation_col = annoCol,
           annotation_colors = list(Term2 = colorTerm2, Term = colorTerm, Ethnicity2 = colorEthnicity2,
                                    CST = colorCST, Vagitype = colorVagitype))
}
order(rowSums(datDF_pw0[, -1]), decreasing = T)[1:3]
datDF_pw0[154, 1]

##### try to pick one pathway #####
# consider the taxonomy of pathway
if(T){ 

  pw = "PWY-5686: UMP biosynthesis I"
  datDF_pwX0 = datDF[grepl(pw, datDF[, 1]) & datDF[, 1] != pw, ]
  range(datDF_pwX0[, -1])
  datDF_pwX1 = data.frame(Pathway = datDF_pwX0[, 1], datDF_pwX0[, -1], check.names = F)
  range(datDF_pwX1[, -1])
  datDF_pwX1_DF = melt(datDF_pwX1)
  colnames(datDF_pwX1_DF)[2:3] = c('Sample', 'Value')
  datDF_pwX1_DF$Taxonomy = str_split_fixed(datDF_pwX1_DF$Pathway, '\\|', 2)[, 2]
  datDF_pwX1_DF$Pathway = str_split_fixed(datDF_pwX1_DF$Pathway, '\\|', 2)[, 1]
  datDF_pwX1_DF$Genus = str_split_fixed(datDF_pwX1_DF$Taxonomy, '\\.', 2)[, 1]
  datDF_pwX1_DF$Species = str_split_fixed(datDF_pwX1_DF$Taxonomy, '\\.', 2)[, 2]
  datDF_pwX1_DF$Species[datDF_pwX1_DF$Genus == 'unclassified'] = 'unclassified'
  datDF_pwX1_DF$pt_ID = swabData$pt_ID[match(datDF_pwX1_DF$Sample, swabData$seq_ID)]
  datDF_pwX1_DF$Term = metaData$Term_char2[match(datDF_pwX1_DF$pt_ID, metaData$pt_ID)]
  datDF_pwX1_DF$Ethnicity = metaData$Ethnicity2[match(datDF_pwX1_DF$pt_ID, metaData$pt_ID)]
  range(datDF_pwX1_DF$Value)
  ggplot(datDF_pwX1_DF, aes(x = Sample, y = Value, fill = Genus)) +
    geom_bar(width = 0.95, position = position_stack(reverse = TRUE), stat = 'identity') + 
    facet_wrap(~Term, scales="free_x", nrow = 6, ncol = 6) +
    xlab("Sample") +
    ylab('Relative abundance') +
    scale_y_continuous(expand = c(0, 0)) +
    # scale_fill_manual(name = 'Genus') +
    guides(fill = guide_legend(nrow = 4, title.position = 'top')) +
    mytheme2 + theme(legend.position = 'bottom', strip.text.x = element_text(size = 13), 
                     axis.text.x = element_text(size = 5))
  ggplot(datDF_pwX1_DF, aes(x = Sample, y = Value, fill = Species)) +
    geom_bar(width = 0.95, position = position_stack(reverse = TRUE), stat = 'identity') + 
    facet_wrap(~Term, scales="free_x", nrow = 6, ncol = 6) +
    xlab("Sample") +
    ylab('Relative abundance') +
    scale_y_continuous(expand = c(0, 0)) +
    # scale_fill_manual(name = 'Genus') +
    guides(fill = guide_legend(nrow = 4, title.position = 'top')) +
    mytheme2 + theme(legend.position = 'bottom', strip.text.x = element_text(size = 13), 
                     axis.text.x = element_text(size = 5))
}
datDF_tx = datDF[which(((grepl('g__', datDF[, 1]) | grepl('unclassified$', datDF[, 1])) &
                          !startsWith(datDF[, 1], 'UNINTEGRATED')) |
                         datDF[, 1] %in% c('UNMAPPED', 'UNINTEGRATED')), ]
colSums(datDF_tx[, 2:5])
grepl('\\|unclassified$', datDF[, 1])
grepl('^UNINTEGRATED\\|', datDF[, 1])
datDF[, 1] == 'UNMAPPED'
datDF[, 1] == 'UNINTEGRATED'

colSums(datDFs[, 2:5])
datDFs[1:5, 1:5]

##### gene family ####
if(T){
  filePath = "./humann/renamed/humann_genefamilies_relab_uniref90_rename.tsv"
  library(vroom)
  datDF = vroom(filePath)
  class(datDF)
  datDF = datDF[, seq(1,734,2)]
  colnames(datDF) %<>% str_remove('Ming_nova_VS_') %>% str_remove('_bwa.*$')
  datDF = datDF[, match(swabData$seq_ID[!is.na(swabData$trimester)], colnames(datDF))]
  datDF = datDF[rowSums(datDF[, -1]) > 0, ]
  dim(datDF)

}

