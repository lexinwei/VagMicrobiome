source('./script/library_function.R')
load('./data/RData/abuList.RData')
load('./data/RData/swabData.RData')
load('./data/RData/colorList.RData')
##### get taxid list and prepare annotation #####
if(T){
  # mags Quality by checkM ran in metawrap
  load('./data/RData/magsQual.RData')
  # taxonmy by GTDB
  lineageBac = read.table('./MAGs/gtdbtk.bac120.summary.tsv', sep = '\t', header = T) # /share/home/jianglab/weixin/workspace/vag/gtdb/byPerson_ani99/gtdbtk.bac120.summary.tsv
  lineage = lineageBac
  bin2rank = data.frame('ID' = lineage$user_genome, 'lineage' = lineage$classification, str_split_fixed(lineage$classification, ';', 7))
  bin2rank = apply(bin2rank, 2, FUN = function(x){str_replace_all(x, '^[dpcofgs]__', '')}) %>% as.data.frame()
  colnames(bin2rank) = c('ID', 'lineage','kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
  table(bin2rank$kingdom) %>% sort(., decreasing = T)
  table(bin2rank$phylum) %>% sort(., decreasing = T)
  table(bin2rank$class) %>% sort(., decreasing = T)
  table(bin2rank$order) %>% sort(., decreasing = T)
  table(bin2rank$family) %>% sort(., decreasing = T)
  table(bin2rank$genus) %>% sort(., decreasing = T)
  table(bin2rank$species) %>% sort(., decreasing = T)
  save(bin2rank, file = './data/RData/bin2rank.Rdata')
  # abundance by coverM
  magAbu0 = data.frame()
  for(f in list.files('./MAGs/coverm/to_MAGs_139_2/',pattern = 'tsv$')){
    tmpf = read.table(paste0('./MAGs/coverm/to_MAGs_139_2/', f), sep = '\t', head = T)
    if(nrow(magAbu0) == 0){
      magAbu0 = tmpf
    }else{
      magAbu0 = merge(magAbu0, tmpf, all = T)
    }
  }
  table(rowSums(magAbu0[, -1]) == 0)
  colnames(magAbu0) %<>% str_extract(., '\\d+_\\d+')
  colnames(magAbu0)[1] = 'Genome'
  magAbu = magAbu0[, !(colnames(magAbu0) %in% swabData$seq_ID[startsWith(swabData$pt_ID.u, 'KIT')])]
  row.names(magAbu) = magAbu0$Genome
  save(magAbu, file = './data/RData/magAbu.RData')
  statAbu = data.frame(ID = row.names(magAbu),
                       total_relative_abundance = rowSums(magAbu[, -1]),
                       mean_relative_abundance = apply(magAbu[, -1], 1, function(x){mean(x)}),
                       counts = apply(magAbu[, -1], 1, function(x){sum(x>0)}))
  statAbu$freq = statAbu$counts/ncol(magAbu[, -1])
  # save(statAbu, file = './data/RData/statAbu.RData')
  if(F){ # for instrain part
    magAbuDF = melt(magAbu0)
    colnames(magAbuDF) = c('Genome', 'seq_ID', 'relative_abundance')
    magAbuDF$pt_ID.u = swabData$pt_ID.u[match(magAbuDF$seq_ID, swabData$seq_ID)]
    save(magAbuDF, file = './data/RData/magAbuDF.RData')
  }
}
if(T){
  trPath = './MAGs/gtdbtk.bac120.iqtree.treefile'
  tree = read.tree(trPath)
  length(tree$tip.label)
  setdiff(tree$tip.label, bin2rank$ID)
  setdiff(tree$tip.label, statAbu$ID)
  magsQual$bin_name = magsQual$bin_name %>% str_remove_all(., '\\.fa')
  setdiff(tree$tip.label, magsQual$bin_name)
  
  taxDF = bin2rank[match(tree$tip.label, bin2rank$ID), ]
  table(taxDF$phylum) %>% sort(., decreasing = T)
  table(taxDF$family) %>% sort(., decreasing = T)
  
  branchPhylum = split(taxDF$ID, taxDF$phylum)
  
  annoTreeDF = data.frame(
    ID = tree$tip.label,
    Kingdom = taxDF$kingdom[match(tree$tip.label, taxDF$ID)],
    Phylum = taxDF$phylum[match(tree$tip.label, taxDF$ID)],
    abuTotal = statAbu$total_relative_abundance[match(tree$tip.label, statAbu$ID)],
    abuMean = statAbu$mean_relative_abundance[match(tree$tip.label, statAbu$ID)],
    Counts = statAbu$counts[match(tree$tip.label, statAbu$ID)],
    Frequency = statAbu$freq[match(tree$tip.label, statAbu$ID)],
    node = 1:length(tree$tip.label)
  )
  tree = groupOTU(tree, branchPhylum, "Phylum")
  p = ggtree(tree, aes(color = Phylum), size = 0.2, # branch line width
             branch.length=T, layout="circular") + 
    scale_color_manual(values = c('black', unlist(colorPhylumMAGs)), guide= 'none') +
    # geom_tiplab(color = 'black', size = 0.87, nudge_x = 0.1) +
    geom_tippoint(mapping=aes(color=Phylum, fill = Phylum), size=1.2, show.legend=FALSE, shape = 21, alpha = 0.7) +
    geom_fruit(geom = geom_tile, mapping = aes(fill=Phylum),
               offset = 0.13, width = 0.15, alpha = 1) +
    scale_fill_manual(values = unlist(colorPhylumMAGs),
                      guide=guide_legend(keywidth = 0.5, keyheight = 0.5, order = 2)) +
    new_scale_fill() +
    geom_fruit(data = annoTreeDF, geom = geom_tile, 
               mapping = aes(y = ID, fill = scale(abuTotal)),
               offset = 0.06, width= 0.15, alpha = 1) +
    scale_fill_distiller(name = 'Normalized total\nabundance', palette = 'OrRd', direction = 1) +
    new_scale_fill() +
    geom_fruit(data = annoTreeDF, geom=geom_bar,
               mapping = aes(y = ID, x = Frequency), fill = '#F7AD00', width= 0.8,
               pwidth = 0.3, offset = 0.04,
               orientation = "y", stat = "identity",
               axis.params = list(axis = "x", text.size = 1, vjust = 1, nbreak = 4),
               grid.params = list(size = 0.08, color = '#FDCF99')
    ) + theme(legend.key.width = unit(0.35, 'cm'), plot.margin = unit(c(-10,10,-10,-10), 'mm'),
              legend.margin=margin(0,-10,0,-20), legend.box.margin = margin(-10,-10,-10,-20))
  cairo_pdf('./plot/MAGs_taxonomy_tree2.pdf', width = 5.02, height = 3.76)
  print(p)
  dev.off()
}
annoTreeDF = merge(annoTreeDF, taxDF)
