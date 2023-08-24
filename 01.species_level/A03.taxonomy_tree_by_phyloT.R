source('./script/library_function.R')
load('./data/RData/abuList.RData')
load('./data/RData/swabData.RData')
load('./data/colorList.RData')
##### get taxid list and prepare annotation #####
if(T){
  abuS = abuList[['s']]
  dim(abuS)
  colnames(abuS)
  range(rowSums(abuS[, 4:ncol(abuS)]))
  # discard NC samples
  removeSamples = swabData$seq_ID[grepl('^KIT', swabData$pt_ID)]
  length(removeSamples)
  abuS1 = abuS[, !colnames(abuS) %in% removeSamples]
  range(rowSums(abuS1[, 4:ncol(abuS1)])) # good! no species exist only in NC samples
  # keep species with taxid
  abuS2 = abuS1[!abuS1$taxid %in% c('', '-1'), ]
  abuS2 = abuS1[!abuS1$taxid %in% c('-1'), ]

  range(colSums(abuS2[,4:ncol(abuS2)]))
  mean(colSums(abuS2[,4:ncol(abuS2)]))
  sd(colSums(abuS2[,4:ncol(abuS2)]))
  if(F){
    newColName = swabData$SampleID.u[match(colnames(abuS2)[4:354], swabData$seq_ID)]
    newColName[is.na(newColName)] = paste0('NA', '.', order(newColName[is.na(newColName)] ))
    colnames(abuS2) = c(colnames(abuS2)[1:3], newColName)
    write.xlsx(abuS2, file = './table/abuS2.xlsx')
  }

  # prepare the taxid and annotation table
  taxDF = data.frame(abuS2[, 1:3], 
             total_relative_abundance = rowSums(abuS2[, 4:ncol(abuS2)]),
             mean_relative_abundance = apply(abuS2[, 4:ncol(abuS2)], 1, mean),
             counts = apply(abuS2[, 4:ncol(abuS2)], 1, function(x){sum(x > 0)}))
  range(taxDF$counts)
  # output taxid list for phyloT
  cat(taxDF$taxid, sep = '\n', file = './data/taxid_list_for_phyloT.txt')
  taxDF$Kingdom = str_split_fixed(taxDF$clade_name, '\\|', 8)[, 1] %>% str_remove('k__')
  taxDF$Phylum = str_split_fixed(taxDF$clade_name, '\\|', 8)[, 2] %>% str_remove('p__')
  taxDF$Class = str_split_fixed(taxDF$clade_name, '\\|', 8)[, 3] %>% str_remove('c__')
  taxDF$Order = str_split_fixed(taxDF$clade_name, '\\|', 8)[, 4] %>% str_remove('o__')
  taxDF$Family = str_split_fixed(taxDF$clade_name, '\\|', 8)[, 5] %>% str_remove('f__')
  taxDF$Genus = str_split_fixed(taxDF$clade_name, '\\|', 8)[, 6] %>% str_remove('g__')
  taxDF$Species = str_split_fixed(taxDF$clade_name, '\\|', 8)[, 7] %>% str_remove('s__')
  taxDF$taxa %<>% str_replace_all(., '_', ' ')
  table(taxDF$Kingdom) %>% sort(., decreasing = T)
  table(taxDF$Phylum)%>% sort(., decreasing = T)
  table(taxDF$Class)%>% sort(., decreasing = T)
  table(taxDF$Order)%>% sort(., decreasing = T)
  table(taxDF$Family)%>% sort(., decreasing = T)
  table(taxDF$Genus)%>% sort(., decreasing = T)
  table(taxDF$Species)%>% sort(., decreasing = T)
  
  table(taxDF$Kingdom, taxDF$Phylum)[, order(colSums(table(taxDF$Kingdom, taxDF$Phylum)), decreasing = T)]
  
  table(taxDF$Kingdom) %>% sort(., decreasing = T) %>% length()
  table(taxDF$Phylum)%>% sort(., decreasing = T) %>% length()
  table(taxDF$Class)%>% sort(., decreasing = T) %>% length()
  table(taxDF$Order)%>% sort(., decreasing = T) %>% length()
  table(taxDF$Family)%>% sort(., decreasing = T) %>% length()
  table(taxDF$Genus)%>% sort(., decreasing = T) %>% length()
  table(taxDF$Species)%>% sort(., decreasing = T)  %>% length()
  

  if(F){
    abuS = abuList[['s']]
    dim(abuS)
    colnames(abuS)
    range(rowSums(abuS[, 4:ncol(abuS)]))
    # discard NC samples
    removeSamples = swabData$seq_ID[grepl('^KIT', swabData$pt_ID)]
    length(removeSamples)
    abuS1 = abuS[, !colnames(abuS) %in% removeSamples]
    range(rowSums(abuS1[, 4:ncol(abuS1)])) # good! no species exist only in NC samples
    # keep species with taxid
    abuS2 = abuS1[!abuS1$taxid %in% c('-1'), ]
    # prepare the taxid and annotation table
    taxDF = data.frame(abuS2[, 1:3], 
                       total_relative_abundance = rowSums(abuS2[, 4:ncol(abuS2)]),
                       mean_relative_abundance = apply(abuS2[, 4:ncol(abuS2)], 1, mean),
                       counts = apply(abuS2[, 4:ncol(abuS2)], 1, function(x){sum(x > 0)}))
    taxDF$freq = taxDF$counts/351
    range(taxDF$total_relative_abundance)
  }
    
}
if(T){
  trPath = './data/phyloT_generated_tree_1674558291_newick.txt'
  tree = read.tree(trPath)
  setdiff(tree$tip.label, taxDF$taxid)
  setdiff(taxDF$taxid, tree$tip.label)
  tree$tip.label[tree$tip.label == 2559073] = 1050843
  tree$tip.label = taxDF$taxa[match(tree$tip.label, taxDF$taxid)]
  branchPhylum = split(taxDF$taxa, taxDF$Phylum)
  annoTreeDF = data.frame(
    ID = tree$tip.label,
    Kingdom = taxDF$Kingdom[match(tree$tip.label, taxDF$taxa)],
    Phylum = taxDF$Phylum[match(tree$tip.label, taxDF$taxa)],
    abuTotal = taxDF$total_relative_abundance[match(tree$tip.label, taxDF$taxa)],
    abuMean = taxDF$mean_relative_abundance[match(tree$tip.label, taxDF$taxa)],
    Counts = taxDF$counts[match(tree$tip.label, taxDF$taxa)],
    Frequency = taxDF$counts[match(tree$tip.label, taxDF$taxa)]/length(4:ncol(abuS2)),
    node = 1:length(tree$tip.label)
  )
  annoTreeDF$Kingdom[annoTreeDF$Kingdom == 'Eukaryota'] = 'Fungi'
  tree = groupOTU(tree, branchPhylum, "Phylum")
  p = ggtree(tree, aes(color = Phylum), size = 0.2, # branch line width
         branch.length=T, layout="circular") + 
    scale_color_manual(values = c('black', unlist(colorPhylum)), guide= 'none') +
    geom_fruit(geom=geom_tile, data = annoTreeDF,
               mapping = aes(y = ID, fill=Kingdom),
               offset = 0.27, width = 8, alpha = 0.3) +
    scale_fill_manual(values = c(unlist(colorKingdom)),
                      guide=guide_legend(keywidth = 0.5, nrow = 2,title.position = 'top',
                                         keyheight=0.5, order = 2)) + 
    new_scale_fill() +
    geom_tiplab(color = 'black', size = 0.87, nudge_x = 0.1) +
    geom_fruit(geom = geom_tile, mapping = aes(fill=Phylum),
               offset = 0.27, width= 0.9, alpha = 1) +
    scale_fill_manual(values = unlist(colorPhylum),
                      guide=guide_legend(keywidth = 0.5, keyheight = 0.5, nrow = 3, order = 3, title.position = 'top')) +
    new_scale_fill() +
    geom_fruit(data = annoTreeDF, geom = geom_tile, 
               mapping = aes(y = ID, fill = scale(abuTotal)),
               offset = 0.065, width= 0.9, alpha = 1) +
    scale_fill_distiller( type = "seq",name = 'Normalized total\nabundance', palette = 'OrRd', direction = 1,
                         guide=guide_colourbar(title.position = 'top', title.hjust = 1, order = 1)) +
    new_scale_fill() +
    geom_fruit(data = annoTreeDF, geom=geom_bar,
               mapping = aes(y = ID, x = Frequency), fill = '#F7AD00', width= 0.8,
               pwidth = 0.3, offset = 0.04,
               orientation = "y", stat = "identity",
               axis.params = list(axis = "x", text.size = 1, vjust = 1, nbreak = 4),
               grid.params = list(size = 0.08, color = '#FDCF99')
    ) + theme(plot.margin = unit(c(-10,10,10,-10), 'mm'), legend.box = 'vertical',
              legend.title = element_text(size = 16), legend.text = element_text(size = 16),
              legend.margin=margin(0,-10,0,-20), legend.box.margin = margin(-10,-10,-10,-20),
              legend.box.just = 'left',legend.spacing = unit(0.3,'cm'),
              legend.box.spacing = unit(1,'cm'),
              legend.key.width = unit(0.9, 'cm'),
              legend.key.height = unit(0.5, 'cm'),
              legend.background = element_rect(),
              # legend.key = element_blank(), # removes the border
              # legend.key.size = unit(0.8, 'cm'), # sets overall area/size of the legend
              legend.position = 'bottom')
  cairo_pdf('./plot/species_taxonomy_tree.pdf', width = 10.36, height = 7.34)
  print(p)
  dev.off()
}
