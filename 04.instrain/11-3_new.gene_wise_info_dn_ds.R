source('~/workspace/library_function.R')
load('./data/metaData.RData')
load('./data/swabData.RData')
load('./data/spName.RData')  # MAGs annotation table

if(T){ # GO term of genes
  pfam2go = read.delim2('data/pfam2go.txt', skip = 5, sep = ';', row.names = NULL)
  pfam2go = cbind(pfam2go$X., pfam2go$row.names %>% str_split_fixed(., '>', 2)) %>% data.frame()
  pfam2go = cbind(pfam2go$X2 %>% str_split_fixed(., ' ', 2), pfam2go$X3 %>% str_remove_all(., ' GO:')) %>% data.frame()
  pfam2go$X3 %<>% str_trim()
  pfam2go$X2 %<>% str_trim()
  pfam2go$X1 %<>% str_replace_all(., 'Pfam:', '')
  table(pfam2go$X3)
  colnames(pfam2go) = c('pfam_accession', 'gene_name', 'GO_term')
  pfam2go$gene_name %<>% str_trim()
  pfam2go = pfam2go[order(pfam2go$pfam_accession, pfam2go$gene_name, pfam2go$GO_term), ]
}

if(F){ # parse hmmscan(pfam) tblout results
  ref_gene_anno = read.delim2('./instrain_res/ref_profile/ref_tblout copy.txt', sep = '\t', header = F)
  target_name = c()
  accession = c()
  query_name = c()
  E_value = c()
  score = c()
  description_of_target = c()
  counter = 1
  for(strg in ref_gene_anno$V1){
    if(counter %% 1000 == 0){
      cat(counter, '/', length(ref_gene_anno$V1), '\n')
    }
    counter = counter + 1
    spres = str_split(strg, '\\s+')[[1]]
    target_name = c(target_name, spres[1])
    accession = c(accession, spres[2])
    query_name = c(query_name, spres[3])
    E_value = c(E_value, spres[5])
    score = c(score, spres[6])
    description_of_target = c(description_of_target, paste0(spres[19:length(spres)], collapse = ' '))
  }
  ref_tbl_table = data.frame(
    query_name = query_name,
    target_accession = accession,
    target_name = target_name,
    E_value = E_value %>% as.numeric(),
    score = score %>% as.numeric(),
    description_of_target = description_of_target
  )
  range(ref_tbl_table$E_value)
  write.csv(ref_tbl_table, './Routput/ref_gene_anno_tbl_table.csv', row.names = F)
  # pick unique target for each query
  ref_tbl_table2 = arrange(ref_tbl_table, desc(score), E_value)
  ref_tbl_table3 = ref_tbl_table2[which(!duplicated(ref_tbl_table2$query_name)), ]
  write.csv(ref_tbl_table3, './Routput/ref_gene_anno_tbl_table_picked.csv', row.names = F)
}

#================== calculate genome wide dnds / pnps ===================#
#### ref ####
if(F){
  # read stb file
  refSTB = read.delim2('./instrain_res/reference_genomes.stb', sep = '\t', header = F, col.names = c('scaffold', 'genome'))
  refDF_gene = data.frame()
  for(ss in swabData$fqID){
    print(ss)
    geneDat = read.delim2(paste0('./instrain_res/ref_profile/gene_info/', ss, '_bwa.IS_gene_info.tsv'))
    geneDat$breadth_minCov %<>% as.numeric()
    geneDat$SNS_N_count %<>% as.numeric()
    geneDat$SNS_S_count %<>% as.numeric()
    geneDat$SNV_S_count %<>% as.numeric()
    geneDat$SNV_N_count %<>% as.numeric()
    geneDat = geneDat[which(geneDat$breadth_minCov >= 0.5), ]
    geneDat$genome = refSTB$genome[match(geneDat$scaffold, refSTB$scaffold)] %>% str_remove_all(., '.fna')
    dt = fread(paste0('./instrain_res/ref_profile/gene_info/', ss, '_bwa_genes_SNP_count.csv.gz'), select = c('gene', 'S_sites', 'N_sites'))
    dt = dt[!duplicated(dt), ]
    geneDat$S_sites = dt$S_sites[match(geneDat$gene, dt$gene)]
    geneDat$N_sites = dt$N_sites[match(geneDat$gene, dt$gene)]
    if(nrow(geneDat) > 0){
      geneDat$dnds = NA
      geneDat$pnps = NA
      geneDat$dn_ds_common = NA
      dn_ds = (geneDat$SNS_N_count / geneDat$N_sites) / (geneDat$SNS_S_count / geneDat$S_sites)
      pn_ps = (geneDat$SNV_N_count / geneDat$N_sites) / (geneDat$SNV_S_count / geneDat$S_sites)
      dn_ds_common = ((geneDat$SNV_N_count + geneDat$SNS_N_count)/geneDat$N_sites) / ((geneDat$SNV_S_count + geneDat$SNS_S_count) / geneDat$S_sites)    
      geneDat$dnds = dn_ds
      geneDat$pnps = pn_ps
      geneDat$dn_ds_common = dn_ds_common
      geneDat$sample = ss
      refDF_gene = rbind(refDF_gene, geneDat)
    }
  }
  refDF_gene$pt_ID = swabData$pt_ID[match(refDF_gene$sample, swabData$fqID)]
  refDF_gene$pt_ID.u = swabData$pt_ID.u[match(refDF_gene$sample, swabData$fqID)]
  refDF_gene$nucl_diversity %<>% as.numeric()
  refDF_gene$SNV_count %<>% as.numeric()
  refDF_gene$species = str_replace(refDF_gene$genome, '.fna', '')
  table(refDF_gene$species)
  refDF_gene$type = 'NCBI'
  refDF_gene$trimester = swabData$trimester[match(refDF_gene$sample, swabData$fqID)]
  refDF_gene$week = swabData$Sample_GA[match(refDF_gene$sample, swabData$fqID)]
  refDF_gene$gene_length %<>% as.numeric()
  refDF_gene$SNV_count %<>% as.numeric()
  refDF_gene$SNPs_per_bp = refDF_gene$SNV_count/(refDF_gene$gene_length*refDF_gene$breadth_minCov)
  save(refDF_gene, file = './data/refDF_gene.RData')
  write.csv(refDF_gene, file = './data/refDF_gene.csv', row.names = F)
}

#### filter out genes with significantly different coverage with other genes on the same genome #####
if(F){
  load('./data/refDF.RData')
  load('./data/refDF_gene.RData')
  refDF_gene$coverage %<>% as.numeric()
  ref_gene_anno = read.csv('./Routput/ref_gene_anno_tbl_table_picked.csv', header = T)
  ref_gene_anno$GO = pfam2go$GO_term[match(str_remove_all(ref_gene_anno$target_accession, '\\.\\d+'), 
                                           pfam2go$pfam_accession)]
  table(is.na(ref_gene_anno$GO))
  table(is.na(ref_gene_anno$target_name))
  table(duplicated(ref_gene_anno$query_name))
  table(duplicated(ref_gene_anno$target_name))

  # remove kit samples
  refDF_gene1 = refDF_gene[!grepl('KIT', refDF_gene$pt_ID.u), ]
  length(unique(refDF_gene1$gene))  # 60220个基因
  nrow(refDF_gene1) # 1362113 
  # add gene annotation
  refDF_gene1$gene_name =  ref_gene_anno$target_name[match(refDF_gene1$gene, ref_gene_anno$query_name)]
  refDF_gene1$description_of_target =  ref_gene_anno$description_of_target[match(refDF_gene1$gene, ref_gene_anno$query_name)]
  refDF_gene1$GO_term = ref_gene_anno$GO[match(refDF_gene1$gene, ref_gene_anno$query_name)]
  table(is.na(refDF_gene1$gene_name))
  length(unique(na.omit(refDF_gene1$gene_name)))
  length(unique(na.omit(refDF_gene1$GO_term)))
  # for each genome, each sample, remove gene beyond 95% quantile and below 55 quantile
  spList = refDF_gene1$genome %>% unique()
  counter = 1
  refDF_gene2 = data.frame()
  for(spp in spList){
    cat('sp =', spp, '  ', counter, '/', length(spList), '\n')
    counter = counter + 1
    refDF_gene_sub = refDF_gene1[refDF_gene1$genome == spp, ]
    counter2 = 1
    tmp_sampleList = unique(refDF_gene_sub$pt_ID.u)
    for(ss in tmp_sampleList){
      refDF_gene_sub_sub = refDF_gene_sub[refDF_gene_sub$pt_ID.u == ss, ]
      if(counter2 %% 100 == 0){
        cat(counter2, '/', length(tmp_sampleList), '\n')
      }
      counter2 = counter2 + 1
      # densityplot(refDF_gene_sub_sub$coverage)
      # shapiro.test(refDF_gene_sub_sub$coverage)
      up95 = quantile(refDF_gene_sub_sub$coverage, probs = 0.95)
      low05 = quantile(refDF_gene_sub_sub$coverage, probs = 0.05)
      refDF_gene_sub_sub_keep = refDF_gene_sub_sub[refDF_gene_sub_sub$coverage < up95 & 
                                                     refDF_gene_sub_sub$coverage > low05, ]
      refDF_gene2 = rbind(refDF_gene2, refDF_gene_sub_sub_keep)
    }
  }
  refDF_gene2$gene_name =  ref_gene_anno$target_name[match(refDF_gene2$gene, ref_gene_anno$query_name)]
  refDF_gene2$description_of_target =  ref_gene_anno$description_of_target[match(refDF_gene2$gene, ref_gene_anno$query_name)]
  refDF_gene2$GO_term = ref_gene_anno$GO[match(refDF_gene2$gene, ref_gene_anno$query_name)]
  length(unique(refDF_gene2$gene)) # 58990
  nrow(refDF_gene2) # 1220463
  length(unique(refDF_gene2$gene)) # 46280
  length(unique(na.omit(refDF_gene2$gene_name))) # 4693
  length(unique(na.omit(refDF_gene2$GO_term))) # 749
  
  
  # remove gene without annotation
  refDF_gene3 = refDF_gene2[!is.na(refDF_gene2$gene_name), ]
  length(unique(refDF_gene3$gene)) # 46280
  length(unique(refDF_gene3$gene_name)) # 4693
  length(unique(na.omit(refDF_gene3$GO_term))) # 749
  nrow(refDF_gene3)
  
  
  densityplot(table(refDF_gene3$gene_name) %>% as.numeric())
  range(table(refDF_gene3$gene_name) %>% sort(., decreasing = T))
  refDF_gene3$GO_term = ref_gene_anno$GO[match(refDF_gene3$gene, ref_gene_anno$query_name)]
  table(refDF_gene3$GO_term %>% is.na())
  save(refDF_gene3, file = 'data/refDF_gene3.RData')
}


##### cat dN/dS>1 or dN/dS<1 of all species, all sample ####
load('data/refDF_gene3.RData')
pdf('./plot/gene_with_dnds_over_1_GO_term.pdf', width = 8.09, height = 5.81)
if(T){
  refDF_cat = refDF_gene3
  refDF_cat$Ethnicity = metaData$Ethnicity[match(refDF_cat$pt_ID, metaData$pt_ID)]
  refDF_cat$Term_char = metaData$Term_char[match(refDF_cat$pt_ID, metaData$pt_ID)]
  refDF_cat$GO_term[is.na(refDF_cat$GO_term)] = 'No GO term'
  
  all_freq = data.frame((refDF_cat$gene_name %>% table() %>% sort(., decreasing = T)))
  genep = ggplot(all_freq[all_freq$Freq >= 3000, ], aes(x = ., y = Freq)) +
    geom_col(fill = '#74ADD1', alpha = 1, width = 0.8) +  mytheme + labs(y = 'Frequency', x = 'Proteins') +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 37000)) + labs(title = ' Proteins with freq>3000 in community of population') +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    geom_text(aes(label=Freq), vjust=-0.5, size = 4)
  print(genep)
  
  sel_pos = refDF_cat[which(refDF_cat$dn_ds_common > 1 & refDF_cat$dn_ds_common != Inf), ]
  length(unique(sel_pos$gene)) # 14470
  length(unique(sel_pos$gene_name)) # 2880
  length(unique(sel_pos$GO_term)) # 544
  nrow(sel_pos)
  sel_pos2 = sel_pos[, c('pt_ID', 'pt_ID.u', 'genome', 'gene_name', 'GO_term')]
  sel_neg = refDF_cat[which(refDF_cat$dn_ds_common < 1 & refDF_cat$dn_ds_common != Inf), ]
  length(unique(sel_neg$gene)) # 44250 
  length(unique(sel_neg$gene_name)) # 4626
  length(unique(sel_neg$GO_term)) # 748
  nrow(sel_neg)
  sel_neg2 = sel_neg[, c('pt_ID', 'pt_ID.u', 'genome', 'gene_name', 'GO_term')]
  # genes with dN/dS >1 and freq>200 in community and population
  pos_freq = data.frame((sel_pos2$gene_name %>% table() %>% sort()))
  pos_freq$all_freq = all_freq$Freq[match(pos_freq$., all_freq$.)]
  pos_freq$ratio  =pos_freq$Freq/pos_freq$all_freq
  genep = ggplot(pos_freq[pos_freq$Freq >= 200, ], aes(x = Freq, y = .)) +
    geom_col(fill = '#74ADD1', alpha = 1, width = 0.8) +  mytheme + labs(x = 'Frequency', y = 'Proteins') +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 1899)) + labs(title = 'Proteins with dN/dS>1 & freq>200 in community and population') +
    geom_text(aes(label=Freq), hjust=-0.2, size = 4.5)
  print(genep)
  
  pos_freq$ratio = pos_freq$Freq/pos_freq$all_freq
  tmp = pos_freq[pos_freq$all_freq >=200, ]
  tmp = pos_freq[pos_freq$ratio >= 0.25 & pos_freq$all_freq >=50, ]
  tmp = tmp[order(tmp$ratio, decreasing = F),]
  tmp$. = factor(tmp$., levels = tmp$.)
  tmp$txt = paste0('(', tmp$Freq, '/', tmp$all_freq, ')')
  genep = ggplot(tmp, aes(x = ratio, y = .)) +
    geom_col(fill = '#74ADD1', alpha = 1, width = 0.8) +  mytheme + labs(x = 'Frequency', y = 'Proteins') +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 0.6)) + labs(title = 'Proteins with dN/dS>1 & freq>50 & pos_ratio>0.25 in community and population') + 
    geom_text(aes(label=txt), hjust=-0.2, size = 4.5)
  print(genep)
  sel_pos2.1 = sel_pos2
  sel_pos2.1 = sel_pos2[sel_pos2$gene_name %in% tmp$., ]
  
  table(duplicated(sel_pos2.1))
  sel_pos3 = ddply(sel_pos2.1,.(pt_ID, pt_ID.u, genome, gene_name, GO_term), nrow)
  sel_pos3 = sel_pos2.1
  sel_pos4 = sel_pos3[sel_pos3$V1 > 2, ]
  sel_pos4 = sel_pos3
  sel_pos5 = rbind(sel_pos4[, c(1,2)] %>% setNames(., c('x', 'y')) %>% unique(), 
                   sel_pos4[, c(2,3)] %>% setNames(., c('x', 'y')) %>% unique(),  
                   sel_pos4[, c(3,4)] %>% setNames(., c('x', 'y')), 
                   sel_pos4[, c(4, 5)] %>% setNames(., c('x', 'y'))
                   # sel_pos4[, c(5,6)] %>% setNames(., c('x', 'y'))
  )
  sel_pos6 = ddply(sel_pos5,.(x,y), nrow)
  # sankey tree to show person-sample-species-genes #####
  # A connection data frame is a list of flows with intensity for each flow
  links <- data.frame(
    source = sel_pos6$x, 
    target =sel_pos6$y,
    value = sel_pos6$V1
  )
  # From these flows we need to create a node data frame: it lists every entities involved in the flow
  nodes <- data.frame(
    name=c(as.character(links$source), 
           as.character(links$target)) %>% unique()
  )
  
  # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
  links$IDsource <- match(links$source, nodes$name)-1 
  links$IDtarget <- match(links$target, nodes$name)-1
  
  # Make the Network
  p <- sankeyNetwork(Links = links, Nodes = nodes,
                     Source = "IDsource", Target = "IDtarget", 
                     Value = "value", NodeID = "name", fontSize = 13, fontFamily = "sans-serif", 
                     sinksRight=FALSE, colourScale = )
  p
  
  
  
  
  
  
  
  
  # genes with dN/dS < 1 and freq>200 in community and population
  neg_freq = data.frame((sel_neg2$gene_name %>% table() %>% sort()))
  genep = ggplot(neg_freq[neg_freq$Freq >= 2000, ], aes(x = Freq, y = .)) +
    geom_col(fill = '#74ADD1', alpha = 1, width = 0.8) +  mytheme + labs(x = 'Frequency', y = 'Proteins') +
    scale_x_continuous(expand = c(0, 0)) + labs(title = 'Genes with dN/dS<1 & freq>2000 in community and population')
  print(genep) 
  
  only_pos_gene = setdiff(pos_freq$., neg_freq$.) # 23
  genep = ggplot(pos_freq[pos_freq$. %in% only_pos_gene, ], aes(x = Freq, y = .)) +
    geom_col(fill = '#74ADD1', alpha = 1, width = 0.8) +  mytheme + labs(x = 'Frequency', y = 'Gene') +
    scale_x_continuous(expand = c(0, 0)) + labs(title = 'Genes with dN/dS>1 & freq>200 in community and population')
  print(genep)
  
  only_neg_gene = setdiff(neg_freq$., pos_freq$.) # 1769
  genep = ggplot(neg_freq[neg_freq$. %in% only_neg_gene, ], aes(x = Freq, y = .)) +
    geom_col(fill = '#74ADD1', alpha = 1, width = 0.8) +  mytheme + labs(x = 'Frequency', y = 'Gene') +
    scale_x_continuous(expand = c(0, 0)) + labs(title = 'Genes with dN/dS>1 & freq>200 in community and population')
  print(genep)
  
  
  
  
  

  
  
  setdiff(sel_neg2$gene_name, sel_pos2$gene_name) %>% length()

  
  
  df = data.frame((sel_neg2$gene_name %>% table() %>% sort()))
  # genes with dN/dS <1 and freq>200 in community and population
  genep = ggplot(df[df$Freq >= 2000, ], aes(x = Freq, y = .)) +
    geom_col(fill = '#74ADD1', alpha = 1, width = 0.8) +  mytheme + labs(x = 'Frequency', y = 'Gene') +
    scale_x_continuous(expand = c(0, 0)) + labs(title = 'Genes with dN/dS<1 & freq>200 in community and population')
  print(genep)
  only_pos_term = setdiff(sel_pos$GO_term, sel_neg$GO_term)
  only_neg_term = setdiff(sel_neg$GO_term, sel_pos$GO_term)
  only_pos_gene = setdiff(sel_pos$gene, sel_neg$gene)
  only_neg_gene = setdiff(sel_neg$gene, sel_pos$gene)
  only_pos_gene_name = setdiff(sel_pos$gene_name, sel_neg$gene_name)
  only_neg_gene_name = setdiff(sel_neg$gene_name, sel_pos$gene_name)
  
  pDF1 = (table(sel_pos$GO_term) %>% sort(., decreasing = T)) %>% as.data.frame
  sel_pos2 = sel_pos[order(sel_pos$dn_ds_common, decreasing = T)[1:100], 
            c('gene_name', 'dn_ds_common', 'description_of_target', 'GO_term')]
  sel_pos2$term = metaData$Term_char[match(sel_pos2$pt_ID, metaData$pt_ID)]
  table(sel_pos2$GO_term, sel_pos2$term)
  pDF2 = pDF1[pDF1$Freq >= 100, ]
  pDF2 = pDF1[order(pDF1$Freq, decreasing = T)[1:20] , ]
  pDF2$Var1 = factor(pDF2$Var1, levels = pDF2$Var1[order(pDF2$Freq)])
  pDF2[1, ]
  genep = ggplot(pDF2[2:nrow(pDF2), ], aes(x = Freq, y = Var1)) +
    geom_col(fill = '#74ADD1', alpha = 1, width = 0.8) +  mytheme + labs(x = 'No. of gene', y = 'GO term') +
    scale_x_continuous(expand = c(0, 0)) + labs(title = 'Top 20 GO term for genes with dN/dS>1')
  print(genep)
}
dev.off()

##### cat gene only have dN/dS<1 of all species ####
refDF_cat = refDF_gene4
refDF_cat$Ethnicity = metaData$Ethnicity[match(refDF_cat$pt_ID, metaData$pt_ID)]
refDF_cat$Term_char = metaData$Term_char[match(refDF_cat$pt_ID, metaData$pt_ID)]
genesSmall = refDF_cat$gene_name[which(refDF_cat$dn_ds_common < 1 & refDF_cat$dn_ds_common != Inf)] %>% unique
genesLarge = refDF_cat$gene_name[which(refDF_cat$dn_ds_common > 1 & refDF_cat$dn_ds_common != Inf)] %>% unique
length(intersect(genesSmall, genesLarge))
geneSuperS = setdiff(genesSmall, genesLarge)
geneSuperSDF = ref_gene_anno[ref_gene_anno$target_name %in% geneSuperS, ]
(table(geneSuperSDF$GO) %>% sort(., decreasing = T))[1:10]
geneSuperL = setdiff(genesLarge, genesSmall)
geneSuperLDF = refDF_cat[which(refDF_cat$gene_name %in% geneSuperL), ]
table(geneSuperLDF$GO)

refDF_cat$GO_term[is.na(refDF_cat$GO_term)] = 'No GO term'
pDF1 = (table(refDF_cat$GO_term) %>% sort(., decreasing = T)) %>% as.data.frame
pDF2 = (table(refDF_cat$gene_name) %>% sort(., decreasing = T)) %>% as.data.frame
refDF_cat2 = refDF_cat[order(refDF_cat$dn_ds_common, decreasing = F)[1:100], 
                       c('gene_name', 'dn_ds_common', 'description_of_target', 'GO_term')]
refDF_cat$term = metaData$Term_char[match(refDF_cat$pt_ID, metaData$pt_ID)]
table(refDF_cat$GO_term, refDF_cat$term)
pDF2 = pDF1[pDF1$Freq >= 100, ]
pDF2$Var1 = factor(pDF2$Var1, levels = pDF2$Var1[order(pDF2$Freq)])
pDF2[1, ]
pdf('./plot/gene_with_dnds_over_1_GO_term.pdf', width = 10, height = 10)
genep = ggplot(pDF2[2:nrow(pDF2), ], aes(x = Freq, y = Var1)) +
  geom_col() +  mytheme + labs(x = 'No. of genes', y = 'GO term') +
  scale_x_continuous(expand = c(0, 0))
print(genep)
dev.off()


##### cat dN/dS of jensenii and gasserii ####
sp = 'Lactobacillus_jensenii'
sp = 'Lactobacillus_gasseri'
refDF_cat = refDF_gene4[refDF_gene4$genome == sp, ]
refDF_cat = refDF_cat[which(refDF_cat$dn_ds_common > 1 & refDF_cat$dn_ds_common != Inf), ]
refDF_cat$GO_term[is.na(refDF_cat$GO_term)] = 'No GO term'
pDF1 = (table(refDF_cat$GO_term) %>% sort(., decreasing = T)) %>% as.data.frame
pDF2 = (table(refDF_cat$gene_name) %>% sort(., decreasing = T)) %>% as.data.frame

refDF_cat$term = metaData$Ethnicity[match(refDF_cat$)]
table(refDF_cat$GO_term, refDF_cat$)

ggplot(pDF1[pDF1$Freq >= 50, ], aes(x = Freq, y = Var1)) +
  geom_col() +  mytheme2 + labs(x = 'No. of genes', y = 'GO term') 

write.csv(refDF_cat, file = paste0('./data/gene_dNdS_', sp, '.csv'), row.names = F)


### 1 find the gene with very high or low nucleotide diveristy #####
load('./data/refDF_gene3.RData')
load('./data/refDF.RData')
# for each genome, each sample, remove gene beyond 95% quantile and below 55 quantile
spList = refDF_gene3$genome %>% unique()
counter = 1
refDF_nc_H = data.frame()
refDF_nc_L = data.frame()
for(spp in spList){
  cat('sp =', spp, '  ', counter, '/', length(spList), '\n')
  counter = counter + 1
  refDF_gene_sub = refDF_gene3[refDF_gene3$genome == spp, ]
  counter2 = 1
  tmp_sampleList = unique(refDF_gene_sub$pt_ID.u)
  for(ss in tmp_sampleList){
    refDF_gene_sub_sub = refDF_gene_sub[refDF_gene_sub$pt_ID.u == ss, ]
    if(counter2 %% 100 == 0){
      cat(counter2, '/', length(tmp_sampleList), '\n')
    }
    counter2 = counter2 + 1
    # densityplot(refDF_gene_sub_sub$nucl_diversity)
    # shapiro.test(refDF_gene_sub_sub$nucl_diversity)
    up95 = quantile(refDF_gene_sub_sub$nucl_diversity, probs = 0.95)
    low05 = quantile(refDF_gene_sub_sub$nucl_diversity, probs = 0.05)
    refDF_gene_sub_sub_keep_H = refDF_gene_sub_sub[refDF_gene_sub_sub$nucl_diversity > up95, ]
    refDF_gene_sub_sub_keep_L = refDF_gene_sub_sub[refDF_gene_sub_sub$nucl_diversity < low05, ]
    refDF_nc_H = rbind(refDF_nc_H, refDF_gene_sub_sub_keep_H)
    refDF_nc_L = rbind(refDF_nc_L, refDF_gene_sub_sub_keep_L)
  }
}

if(F){
  counter = 1
  comped_gene = c()
  wcTest.P = c()
  speciess = c()
  mean1 = c()
  mean2 = c()
  sp_in_upto3_samples = table(refDF_gene3$species)[which(table(refDF_gene3$species) >= 3)] %>% names()
  for(spp in sp_in_upto3_samples){
    cat('sp =', spp, '  ', counter, '/', length(sp_in_upto3_samples), '\n')
    counter = counter + 1
    refDF_gene3_sub = refDF_gene3[refDF_gene3$genome == spp, ]
    counter2 = 1
    for(ggs in unique(refDF_gene3_sub$gene)){
      if(counter2 %% 100 == 0){
        cat(counter2, '/', length(unique(refDF_gene3_sub$gene)), '\n')
      }
      counter2 = counter2 + 1
      cov1 = refDF_gene3_sub$nucl_diversity[refDF_gene3_sub$gene == ggs] %>% as.numeric()
      cov2 = refDF$nucl_diversity[(refDF$species == spp) & (refDF$pt_ID.u %in% refDF_gene3_sub$pt_ID.u)]
      len1 = length(cov1)
      len2 = length(cov2)
      wcTest = wilcox.test(cov1, cov2)
      wcTest.P = c(wcTest.P, wcTest$p.value)
      comped_gene = c(comped_gene, ggs)
      speciess  = c(speciess, spp)
      mean1 = c(mean1, mean(cov1))
      mean2 = c(mean2, mean(cov2))
    }
  }
  gene_nucdiv_wilcoxTest = data.frame(
    species = speciess,
    gene = comped_gene,
    len1 = len1,
    len2 = len2,
    mean1 = mean1,
    mean2 = mean2,
    P_value = wcTest.P
  )
  gene_nucdiv_wilcoxTest$mad = gene_nucdiv_wilcoxTest$mean1 - gene_nucdiv_wilcoxTest$mean2
  range(gene_nucdiv_wilcoxTest$mad)
  table(gene_nucdiv_wilcoxTest$mad < 0)
  
  gene_nucdiv_wilcoxTest2 = merge(gene_nucdiv_wilcoxTest, ref_gene_anno, by.x = 'gene', by.y = 'query_name')
  write.csv(gene_nucdiv_wilcoxTest2, file = './Routput/gene_nucdiv_wilcoxTest_genomeAsBGM.csv', row.names = F)
}


### 2 find the gene with very high or low nucleotide diveristy ######
gene_nucdiv_wilcoxTest = read.csv('./Routput/gene_nucdiv_wilcoxTest_genomeAsBGM.csv')
gene_nucdiv_wilcoxTest$adj.P_value = p.adjust(gene_nucdiv_wilcoxTest$P_value, method = 'BH')
refDF_gene3_nuc_sig = refDF_gene3[refDF_gene3$gene %in% gene_nucdiv_wilcoxTest$gene[gene_nucdiv_wilcoxTest$adj.P_value <= 0.05], ]
nrow(refDF_gene3_nuc_sig)
refDF_gene3_nuc_sigUp = refDF_gene3[refDF_gene3$gene %in% gene_nucdiv_wilcoxTest$gene[gene_nucdiv_wilcoxTest$adj.P_value <= 0.05 & gene_nucdiv_wilcoxTest$mad > 0], ]
refDF_gene3_nuc_sigDown = refDF_gene3[refDF_gene3$gene %in% gene_nucdiv_wilcoxTest$gene[gene_nucdiv_wilcoxTest$adj.P_value <= 0.05 & gene_nucdiv_wilcoxTest$mad < 0], ]

refDF_cat = refDF_gene3_nuc_sigDown
refDF_cat$Ethnicity = metaData$Ethnicity[match(refDF_cat$pt_ID, metaData$pt_ID)]
refDF_cat$Term_char = metaData$Term_char[match(refDF_cat$pt_ID, metaData$pt_ID)]
refDF_cat$GO_term[is.na(refDF_cat$GO_term)] = 'No GO term'
all_freq = data.frame((refDF_gene3_nuc_sig$gene_name %>% table() %>% sort(., decreasing = T)))
genep = ggplot(all_freq[all_freq$Freq >= 1000, ], aes(x = ., y = Freq)) +
  geom_col(fill = '#74ADD1', alpha = 1, width = 0.8) +  mytheme + labs(y = 'Frequency', x = 'Proteins') +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 8500)) + labs(title = ' Proteins with freq>1000 in community of population') +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  geom_text(aes(label=Freq), vjust=-0.5, size = 4)
print(genep)

sel_pos = refDF_gene3_nuc_sigUp
nrow(sel_pos)
length(unique(sel_pos$gene)) # 1106
length(unique(sel_pos$gene_name)) # 668
length(unique(sel_pos$GO_term)) # 147
nrow(sel_pos)

sel_pos2 = sel_pos[, c('pt_ID', 'pt_ID.u', 'genome', 'gene_name', 'GO_term')]

sel_neg = refDF_gene3_nuc_sigDown
nrow(sel_neg)
length(unique(sel_neg$gene)) # 5002 
length(unique(sel_neg$gene_name)) # 1958
length(unique(sel_neg$GO_term)) # 400

sel_neg2 = sel_neg[, c('pt_ID', 'pt_ID.u', 'genome', 'gene_name', 'GO_term')]
# genes with dN/dS >1 and freq>200 in community and population
pos_freq = data.frame((sel_pos2$gene_name %>% table() %>% sort()))
pos_freq$all_freq = all_freq$Freq[match(pos_freq$., all_freq$.)]
pos_freq$ratio  =pos_freq$Freq/pos_freq$all_freq
genep = ggplot(pos_freq[pos_freq$Freq >= 200, ], aes(x = Freq, y = .)) +
  geom_col(fill = '#74ADD1', alpha = 1, width = 0.8) +  mytheme + labs(x = 'Frequency', y = 'Proteins') +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1800)) + labs(title = 'Proteins with above average nucl_div & freq>200 in community and population') +
  geom_text(aes(label=Freq), hjust=-0.2, size = 4.5)
print(genep)

pos_freq$ratio = pos_freq$Freq/pos_freq$all_freq
tmp = pos_freq[pos_freq$ratio >= 0.80 & pos_freq$all_freq >=100, ]
tmp = tmp[order(tmp$ratio, decreasing = F),]
tmp$. = factor(tmp$., levels = tmp$.)
tmp$txt = paste0('(', tmp$Freq, '/', tmp$all_freq, ')')
genep = ggplot(tmp, aes(x = ratio, y = .)) +
  geom_col(fill = '#74ADD1', alpha = 1, width = 0.8) +  mytheme + labs(x = 'Frequency', y = 'Proteins') +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1.2)) + labs(title = 'Proteins with above average nucl_div & freq>100 & pos_ratio>0.8 in community and population') + 
  geom_text(aes(label=txt), hjust=-0.2, size = 4.5)
print(genep)
sel_pos2.1 = sel_pos2
sel_pos2.1 = sel_pos2[sel_pos2$gene_name %in% tmp$., ]
sel_pos2.1$GO_term[is.na(sel_pos2.1$GO_term)] = 'No GO term'
table(duplicated(sel_pos2.1))
sel_pos3 = ddply(sel_pos2.1,.(pt_ID, pt_ID.u, genome, gene_name, GO_term), nrow)
# sel_pos3 = sel_pos2.1
sel_pos4 = sel_pos3[sel_pos3$V1 >3, ]
# sel_pos4 = sel_pos3
sel_pos5 = rbind(sel_pos4[, c(1,2)] %>% setNames(., c('x', 'y')) %>% unique(), 
                 sel_pos4[, c(2,3)] %>% setNames(., c('x', 'y')) %>% unique(),  
                 sel_pos4[, c(3,4)] %>% setNames(., c('x', 'y')), 
                 sel_pos4[, c(4, 5)] %>% setNames(., c('x', 'y'))
                 # sel_pos4[, c(5,6)] %>% setNames(., c('x', 'y'))
)
sel_pos6 = ddply(sel_pos5,.(x,y), nrow)
# sankey tree to show person-sample-species-genes #####
# A connection data frame is a list of flows with intensity for each flow
links <- data.frame(
  source = sel_pos6$x, 
  target =sel_pos6$y,
  value = sel_pos6$V1
)
# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

# Make the Network
library(networkD3)
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget", 
                   Value = "value", NodeID = "name", fontSize = 13, fontFamily = "sans-serif", 
                   sinksRight=FALSE, colourScale = )
p


# 
# 
# gene_nucdiv_wilcoxTest2_sig = gene_nucdiv_wilcoxTest2[gene_nucdiv_wilcoxTest2$adj.P_value <= 0.01, ]
# # top 10 low microdiversity
# gene_nucdiv_wilcoxTest2_sig = gene_nucdiv_wilcoxTest2_sig[order(gene_nucdiv_wilcoxTest2_sig$adj.P_value), ]
# 
# write.csv(gene_nucdiv_wilcoxTest2_sig, file = './Routput/gene_nucdiv_wilcoxTest_genomeAsBGM_sig.csv', row.names = F)
# gene_nucdiv_wilcoxTest2_sig[which(gene_nucdiv_wilcoxTest2_sig$mad < 0), ][1:20, ]
# 
# 
# table(gene_nucdiv_wilcoxTest2_sig$mad > 0)
# sigUP = gene_nucdiv_wilcoxTest2_sig[gene_nucdiv_wilcoxTest2_sig$mad > 0, ]
# sigDOWN = gene_nucdiv_wilcoxTest2_sig[gene_nucdiv_wilcoxTest2_sig$mad < 0, ]
# sigTMP = sigUP
# sigTMP = sigDOWN
# table(sigTMP$species)
# ggg =  table(sigTMP$target_name)[table(sigTMP$target_name) >= 2]
# sort(ggg, decreasing = T)[1:100]
# vsDF= dcast(sigTMP[, c('target_name', 'species')], target_name~species)
# vsMat = as.matrix(vsDF[, -1])
# row.names(vsMat) = vsDF$target_name
# vsMat1 = vsMat[rowSums(vsMat) >= 3, ]
# pheatmap(vsMat1, scale = 'none', cluster_rows = F, cluster_cols = F)
# 
# 
# metaData
# 
# # get top 5 of each species
# t = 1
# upTop = data.frame()
# for(spp in gene_nucdiv_wilcoxTest2_sig$species %>% unique()){
#   sssf = gene_nucdiv_wilcoxTest2_sig[gene_nucdiv_wilcoxTest2_sig$species == spp & gene_nucdiv_wilcoxTest2_sig$mad > 0, ]
#   upTop = rbind(upTop, sssf[order(sssf$P_value)[1:t],])
# }
# upTop %<>% na.omit()
# table(upTop$target_name) %>% sort(., decreasing = T)
# write.csv(upTop, './data/top1_nucl_div_up_each_species.csv', row.names = F)
# 
# downTop = data.frame()
# for(spp in gene_nucdiv_wilcoxTest2_sig$species %>% unique()){
#   sssf = gene_nucdiv_wilcoxTest2_sig[gene_nucdiv_wilcoxTest2_sig$species == spp & gene_nucdiv_wilcoxTest2_sig$mad < 0, ]
#   downTop = rbind(downTop, sssf[order(sssf$P_value)[1:t],])
# }
# downTop %<>% na.omit()
# table(downTop$target_name) %>% sort(., decreasing = T)
# intersect(downTop$target_name, upTop$target_name)
# write.csv(downTop, './data/top1_nucl_div_down_each_species.csv', row.names = F)
# 
# 
# 
# genes with dN/dS >1 and freq>200 in community and population
neg_freq = data.frame((sel_neg2$gene_name %>% table() %>% sort()))
neg_freq$all_freq = all_freq$Freq[match(neg_freq$., all_freq$.)]
neg_freq$ratio  =neg_freq$Freq/neg_freq$all_freq
genep = ggplot(neg_freq[neg_freq$Freq >= 600, ], aes(x = Freq, y = .)) +
  geom_col(fill = '#74ADD1', alpha = 1, width = 0.8) +  mytheme + labs(x = 'Frequency', y = 'Proteins') +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 8500)) + labs(title = 'Proteins with below average nucl_div & freq>600 in community and population') +
  geom_text(aes(label=Freq), hjust=-0.2, size = 4.5)
print(genep)

neg_freq$ratio = neg_freq$Freq/neg_freq$all_freq
tmp = neg_freq[neg_freq$ratio >= 0.80 & neg_freq$all_freq >=600, ]
dim(tmp)
tmp = tmp[order(tmp$ratio, decreasing = F), ]
tmp$. = factor(tmp$., levels = tmp$.)
tmp$txt = paste0('(', tmp$Freq, '/', tmp$all_freq, ')')
genep = ggplot(tmp, aes(x = ratio, y = .)) +
  geom_col(fill = '#74ADD1', alpha = 1, width = 0.8) +  mytheme + labs(x = 'Frequency', y = 'Proteins') +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1.2)) + labs(title = 'Proteins with above average nucl_div & freq>600 & neg_ratio>0.8 in community and population') + 
  geom_text(aes(label=txt), hjust=-0.2, size = 4.5)
print(genep)
sel_neg2.1 = sel_neg2
sel_neg2.1 = sel_neg2[sel_neg2$gene_name %in% tmp$., ]
sel_neg2.1$GO_term[is.na(sel_neg2.1$GO_term)] = 'No GO term'
table(duplicated(sel_neg2.1))
sel_neg3 = ddply(sel_neg2.1,.(pt_ID, pt_ID.u, genome, gene_name, GO_term), nrow)
# sel_neg3 = sel_neg2.1
sel_neg4 = sel_neg3[sel_neg3$V1 > 15 & sel_neg3$gene_name != 'ABC_tran', ]
# sel_neg4 = sel_neg3
sel_neg5 = rbind(sel_neg4[, c(1,2)] %>% setNames(., c('x', 'y')) %>% unique(), 
                 sel_neg4[, c(2,3)] %>% setNames(., c('x', 'y')) %>% unique(),  
                 sel_neg4[, c(3,4)] %>% setNames(., c('x', 'y')), 
                 sel_neg4[, c(4, 5)] %>% setNames(., c('x', 'y'))
                 # sel_neg4[, c(5,6)] %>% setNames(., c('x', 'y'))
)
sel_neg6 = ddply(sel_neg5,.(x,y), nrow)
# sankey tree to show person-sample-species-genes #####
# A connection data frame is a list of flows with intensity for each flow
links <- data.frame(
  source = sel_neg6$x, 
  target =sel_neg6$y,
  value = sel_neg6$V1
)
# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

# Make the Network
library(networkD3)
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget", 
                   Value = "value", NodeID = "name", fontSize = 13, fontFamily = "sans-serif", 
                   sinksRight=FALSE, colourScale = )
p
