source('~/workspace/library_function.R')
spName = read.xlsx('./tax_of_mags_co.xlsx')
spName$lca_species_manual = spName$lca_species
spName$lca_species_manual[!is.na(spName$species)] = spName$species[!is.na(spName$species)]
spName$lca_species_manual[spName$lca_species_manual == 'Life'] = 'cellular organisms'
spName$lca_species_manual %<>% str_replace_all(., '_', ' ')


spName$lca_level_manual = spName$lca_level
spName$lca_level_manual[!is.na(spName$species)] = 'species'
spName$lca_level_manual[spName$lca_level_manual == 'strain'] = 'species' 
spName$lca_level_manual[spName$lca_level_manual == 'superkingdom'] = 'kingdom'
spName$lca_level_manual %<>% capitalize()
spName$lca_level_manual = factor(spName$lca_level_manual, levels = c('Life', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'))
table(spName$lca_level_manual)
name2taxid = read.delim2('./data/name2taxid.txt', sep = '\t')
name2taxid$name %<>% str_remove_all(., "[rdkcofgs]__")
name2taxid$name %<>% str_replace_all(., "_", ' ')
spName$lca_last_taxid = name2taxid$taxid[match(spName$lca_species_manual, name2taxid$name)]
spName$lca_species_manual %<>% str_replace_all(., '_', ' ')

spName$lca_last_taxid
spName$lca_last_taxid[spName$lca_species_manual == 'cellular organisms'] = 131567
spName$lca_last_taxid[spName$lca_species_manual == 'Firmicutes'] = 1239
spName$lca_last_taxid[spName$lca_species_manual == 'Actinomycetia'] = 1760
spName$lca_last_taxid[spName$lca_species_manual == 'Limosilactobacillus reuteri 1063'] = 1273150
spName$lca_last_taxid[spName$lca_species_manual == 'Atopobium minutum'] = 1381
spName$lca_last_taxid[spName$lca_species_manual == 'Actinomyces sp. zg-332'] = 2708340
spName$lca_last_taxid[spName$lca_species_manual == 'Lactobacillus rhamnosus'] = 47715
spName$lca_last_taxid[spName$lca_species_manual == 'Lancefieldella parvulum DSM 20469'] = 521095
spName$lca_last_taxid[spName$lca_species_manual == 'Berryella intestinalis'] = 1531429
spName$lca_last_taxid[spName$lca_species_manual == 'Bifidobacterium scardovii JCM 12489 = DSM 13734'] = 1150461
spName$lca_last_taxid[spName$lca_species_manual == 'Mageeibacillus indolicus '] = 884684
spName$lca_species_manual[spName$lca_species_manual == 'Mageeibacillus indolicus '] = 'Mageeibacillus indolicus'
spName$uni_species = paste0(spName$lca_species_manual, ' (', spName$Bin, ') ')
save(spName, file = './data/spName.RData')

str_c(spName$lca_last_taxid %>% unique, collapse = ',')
mags_linkage = read.delim2('./data/lca_last_taixd_linkage.txt', header = F)
mags_linkage = str_split_fixed(mags_linkage$V3, ';', 7) %>% as.data.frame()
mags_linkage$V1[mags_linkage$V1 == ''] = 'Unclassified' 
mags_linkage[mags_linkage == ''] = NA
colnames(mags_linkage) = c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
cNUm = table(mags_linkage %>% as.matrix())
un = c()
n = c()
ftn = c()
f = c()
t = c()
for (nr in 1:nrow(mags_linkage)) {
  itms = mags_linkage[nr, ] %>% unlist()
  itms = itms[!is.na(itms)]
  if(length(itms) == 1){
    n = c(n, itms[1])
    un = c(un, itms[1])
    f = c(f, itms[1])
    t = c(t, itms[1])
  }else{
    for(im in 1:(length(itms)-1)){
      f = c(f, str_c(itms[1:im], collapse = '/'))
      t = c(t, paste0(str_c(itms[1:im], collapse = '/'), '/', itms[im + 1]))
      n = c(n, itms[im + 1])
      ftn = c(ftn, itms[im + 1])
      if(im == (length(itms)-1)){
        un = c(un, itms[im+1])
      }
    }
  }
}
edgeDF = data.frame(
  from = f,
  to = t,
  short = n
)
table(t)
edgeDF = edgeDF[!duplicated(edgeDF), ]
verticesDF = data.frame(
  name = t,
  nodeName = n,
  level = names(n)
)

verticesDF$size = cNUm[match(verticesDF$nodeName, names(cNUm))]
verticesDF$size2 = table(un)[match(verticesDF$nodeName, names(table(un)))]
verticesDF$size2[!is.na(verticesDF$size2)] = paste0(' (', verticesDF$size2[!is.na(verticesDF$size2)], ')')
verticesDF$size2[is.na(verticesDF$size2)] = ''
verticesDF = verticesDF[!duplicated(verticesDF), ]
verticesDF$level[verticesDF$name == 'Unclassified'] = 'Unclassified'
library(ggraph)
library(influential)
igT = graph_from_data_frame(edgeDF, vertices = verticesDF)
# plot using ggraph
cccolor = brewer.pal(9, 'Set3')[c(1,3:9)]
names(cccolor) = c('Kingdom', 'Phylum','Class', 'Order', 'Family', 'Genus', 'Species', 'Unclassified')
pdf('./plot/mags_taxonomy_tree.pdf', width = 19, height = 12)
ggraph(igT, layout = 'tree', circular = FALSE) + 
  geom_edge_diagonal(color = 'grey', width = 0.3) +
  geom_node_point(aes(color=level %>% factor(., levels = c('Kingdom', 'Phylum','Class', 'Order', 'Family', 'Genus', 'Species', 'Unclassified')),
                      fill = level %>% factor(., levels = c('Kingdom', 'Phylum','Class', 'Order', 'Family', 'Genus', 'Species', 'Unclassified'))), 
                      size = 6.7) +
  geom_node_text(
    aes(label=paste0(size)), size = 4.5
  ) +
  geom_node_text(
    aes(label=paste0(nodeName, size2)), size = 4.5,
    vjust = 0.5,
    hjust = 0,
    nudge_y = 0.08,
    nudge_x = 0
  ) +
  scale_color_manual(values=cccolor %>% alpha(.,0.9)) +
  scale_fill_manual(values=cccolor %>% alpha(.,0.7)) +
  theme_void() + theme(plot.margin = unit(c(0,11,0,0), 'cm'),
                       # legend.key.height= unit(2, 'cm'),
                       # legend.key.width= unit(2, 'cm'),
                       legend.position = c(1.2, 0.5),
                       legend.title = element_text(size=13),
                       legend.text = element_text(size=13)) +
  labs(color = 'Rank', fill = 'Rank') +
  scale_x_reverse() +
  scale_y_reverse() +
  coord_flip(clip = 'off')

dev.off()

{ #### leaf background classification rank ####
  mycol = c('#969696', '#8C96C6', '#A6761D', '#FB8072', '#7FC97F', '#BC80BD', '#FDB462', '#80B1D3') %>% alpha(., 0.5)
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  legend("topleft", legend =paste0(table(spName$lca_level_manual) %>% names(), 
                           ' (', table(spName$lca_level_manual), ')'), 
         pch=16, pt.cex=3, cex=1.5, bty='n',
         col = mycol)
  mtext("Rank", at=0.03, cex=1.7)
  table(spName$lca_level_manual)
  cv = rep(NA, nrow(spName))
  for (l in 1:(table(spName$lca_level_manual) %>% length())) {
    loc = (table(spName$lca_level_manual) %>% names())[l]
    ind = which(spName$lca_level_manual== loc)
    cv[ind] = mycol[l]
  }
  strip = rbind(data.frame(
    "NODE_ID" = spName$lca_species_manual,
    'Range' = 'range',
    "COLOR" = cv,
    "LABEL" = spName$lca_level_manual
  )
  )

  strip = strip[!duplicated(strip), ]
  strip$NODE_ID %<>% str_replace_all(., ' ', '_')

  sink("./data/anno_range_rank.txt")
  cat("TREE_COLORS\nSEPARATOR TAB\nDATA\n")
  write.table(strip, file = "./data/anno_range_rank.txt", 
              sep = '\t', quote = F, col.names = F, row.names = F, append = T)
  sink()
}
# 
# Couldn't find ID Lactobacillus_brevis in the tree
# Couldn't find ID Lactobacillus_rhamnosus in the tree
# Couldn't find ID Mageeibacillus_indolicus_ in the tree



