source('~/workspace/library_function.R')
load('./data/metaData.RData')
load('./data/swabData.RData')
load('./data/spName.RData')
# parse hmmscan(pfam) tblout results
if(F){
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
  dt = dt[!duplicated(dt),]
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

#### filter out genes with significantly different coverage with other genes on the same genome #####
load('./data/refDF_gene.RData')
load('./data/refDF.RData')
refDF_gene1 = refDF_gene[!grepl('KIT', refDF_gene$pt_ID.u), ]
refDF_gene1$gene_name =  ref_gene_anno$target_name[match(refDF_gene1$gene, ref_gene_anno$query_name)]
refDF_gene1 = refDF_gene1[!is.na(refDF_gene1$gene_name), ]
sp_in_upto3_samples = table(refDF$species)[which(table(refDF$species) >= 3)] %>% names()
refDF_gene2 = refDF_gene1[refDF_gene1$genome %in% sp_in_upto3_samples, ]
gene_in_upto3_samples = table(refDF_gene2$gene)[which(table(refDF_gene2$gene) >= 3)] %>% names()
refDF_gene3 = refDF_gene2[refDF_gene2$gene %in% gene_in_upto3_samples, ]

counter = 1
comped_gene = c()
wcTest.P = c()
speciess = c()
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
    cov1 = refDF_gene3_sub$coverage[refDF_gene3_sub$gene == ggs] %>% as.numeric()
    cov2 = refDF_gene3_sub$coverage[refDF_gene3_sub$gene != ggs] %>% as.numeric()
    wcTest = wilcox.test(cov1, cov2)
    wcTest.P = c(wcTest.P, wcTest$p.value)
    comped_gene = c(comped_gene, ggs)
    speciess  = c(speciess, spp)
  }
}
gene_coverage_wilcoxTest = data.frame(
  species  =speciess,
  gene = comped_gene,
  P_value = wcTest.P
)
write.csv(gene_coverage_wilcoxTest, file = './Routput/gene_coverage_wilcoxTest.csv', row.names = F)

### find the gene with very high or low nucleotide diveristy
gene_coverage_wilcoxTest = read.csv('./Routput/gene_coverage_wilcoxTest.csv')
gene_coverage_wilcoxTest$adj.P_value = p.adjust(gene_coverage_wilcoxTest$P_value, method = 'BH')
refDF_gene4 = refDF_gene3[refDF_gene3$gene %in% gene_coverage_wilcoxTest$gene[gene_coverage_wilcoxTest$adj.P_value >= 0.1], ]

counter = 1
comped_gene = c()
wcTest.P = c()
speciess = c()
med1 = c()
med2 = c()
for(spp in sp_in_upto3_samples){
  cat('sp =', spp, '  ', counter, '/', length(sp_in_upto3_samples), '\n')
  counter = counter + 1
  refDF_gene4_sub = refDF_gene4[refDF_gene4$genome == spp, ]
  counter2 = 1
  for(ggs in unique(refDF_gene4_sub$gene)){
    if(counter2 %% 100 == 0){
      cat(counter2, '/', length(unique(refDF_gene4_sub$gene)), '\n')
    }
    counter2 = counter2 + 1
    cov1 = refDF_gene4_sub$nucl_diversity[refDF_gene4_sub$gene == ggs] %>% as.numeric()
    cov2 = refDF_gene4_sub$nucl_diversity[refDF_gene4_sub$gene != ggs] %>% as.numeric()
    wcTest = wilcox.test(cov1, cov2)
    wcTest.P = c(wcTest.P, wcTest$p.value)
    comped_gene = c(comped_gene, ggs)
    speciess  = c(speciess, spp)
    med1 = c(med1, mean(cov1))
    med2 = c(med2, mean(cov2))
  }
}
gene_nucdiv_wilcoxTest = data.frame(
  species  =speciess,
  gene = comped_gene,
  med1 = med1,
  med2 = med2,
  P_value = wcTest.P
)
gene_nucdiv_wilcoxTest$mad = gene_nucdiv_wilcoxTest$med1 - gene_nucdiv_wilcoxTest$med2
range(gene_nucdiv_wilcoxTest$mad )
table(gene_nucdiv_wilcoxTest$mad < 0)

write.csv(gene_nucdiv_wilcoxTest, file = './Routput/gene_nucdiv_wilcoxTest.csv', row.names = F)
gene_nucdiv_wilcoxTest$adj.P_value = p.adjust(gene_nucdiv_wilcoxTest$P_value, method = 'BH')
refDF_gene4_nuc_sig = refDF_gene4[refDF_gene4$gene %in% gene_nucdiv_wilcoxTest$gene[gene_nucdiv_wilcoxTest$adj.P_value <= 0.05], ]

gene_nucdiv_wilcoxTest2 = merge(gene_nucdiv_wilcoxTest, ref_gene_anno, by.x = 'gene', by.y = 'query_name')
gene_nucdiv_wilcoxTest2_sig = gene_nucdiv_wilcoxTest2[gene_nucdiv_wilcoxTest2$adj.P_value <= 0.05, ]
# top 10 low microdiversity
gene_nucdiv_wilcoxTest2_sig = gene_nucdiv_wilcoxTest2_sig[order(gene_nucdiv_wilcoxTest2_sig$adj.P_value), ]
gene_nucdiv_wilcoxTest2_sig[which(gene_nucdiv_wilcoxTest2_sig$mad < 0), ][1:20, ]
write.csv(gene_nucdiv_wilcoxTest2_sig, file = './Routput/gene_nucdiv_wilcoxTest_sig.csv', row.names = F)










#### associate gene with variables in species wise ####
ref_gene_anno = read.csv('./Routput/ref_gene_anno_tbl_table_picked.csv', header = T)
ys = c('nucl_diversity', 'dnds', 'pnps', 'dn_ds_common')
ysNorm = list(nucl_diversity = 'Nucleotide diversity', dnds = 'dN/dS', pnps = 'pN/pS', dn_ds_common = 'dN/dS (common)')
xs = c('trimester', 'Age_char', 'BMI_char', 'Ethnicity', 'Language', 'Marriage', 'FOB', 'Employed', 'Housing', 
       'Term_char', 'Gravida_char', 'Abortions_char', 'Living_children_char', 
       'Abnormal_pap', 'Depression', 'PTSD', 'PPD_pos', 'DM2', 'Asthma', 'PCN_allergy', 'Sulfa_allergy', 'NKDA',
       'PNV', 'Fe', 'Insulin', 'B6', 'Progesterone',
       'IPV_hx', 'TED_hx', 'Tobacco_hx', 'EtOH_hx', 
       'GBS', 'PCN', 'Antibiotic',
       'Induction', 'Baby_gender', 'Baby_weight_char', 
       # repeated numeric 
       'Age', 'BMI', 'Term', 'Gravida', 'Abortions', 'Living_children', 'Baby_weight')
comp_DF = data.frame()
uniRes = data.frame()
counter = 1
for(sps in sp_in_upto3_samples){
  cat('sp =', sps, '  ', counter, '/', length(sp_in_upto3_samples), '\n')
  counter = counter + 1
  tmpSub = refDF_gene4_nuc_sig[refDF_gene4_nuc_sig$genome == sps, ]
  tmpSub$gene_name = ref_gene_anno$target_name[match(tmpSub$gene, ref_gene_anno$query_name)]
  tmpSub2 = tmpSub[!is.na(tmpSub$gene_name), ]
  # length(unique(tmpSub2$gene))
  # length(unique(tmpSub2$gene_name))
  counter2 = 1
  for(gg in unique(tmpSub2$gene)){
    cat('gene = ', gg, '  ', counter2, '/', length(unique(tmpSub2$gene)), '\n')
    counter2 = counter2 + 1
    tmp = tmpSub2[tmpSub2$gene == gg, ]
    tmp2 = merge(tmp, metaData, by = 'pt_ID', all.x = T)
      for (y in ys) {
        # cat('  y =', y, '\n')
        for (x in xs) {
          # cat('    x =', x, '\n')
          tmpDF2 = data.frame(
            'y' = tmp2[, y],
            'x' = tmp2[, x],
            trimester = tmp2$trimester,
            pt_ID = tmp2$pt_ID
          )
          tmpDF2 = na.omit(tmpDF2)
          tmpDF2 = tmpDF2[!is.infinite(tmpDF2$y), ]
          if(length(na.omit(unique(tmpDF2$x))) >= 2){
            # lm
            lmUni = lm(y ~ x, tmpDF2)
            lmUni.summary = summary(lmUni)
            lmUni.aov = anova(lmUni)
            # lmer
            formula_char = paste0(y, ' ~ (Intercept) + ', x, ' + (1|Subject)')
            lmeTry = try({lmer(as.formula(paste0('y ~ x + trimester + (1|pt_ID)')),
                               data = tmpDF2, REML = F)}, silent = F)
            if(class(lmeTry) != 'try-error'){
              lmeUni = lmeTry
              lmeUni.R2 = r.squaredGLMM(lmeUni)[1]
              lmeUni.aov = anova(lmeUni)
              sv = c(sp, gg, y, x, (lmUni.aov['x', c('F value', 'Pr(>F)')]) %>% unlist(), lmUni.summary$r.squared %>% as.numeric(),
                     formula_char, lmeUni.aov$`F value`[1], lmeUni.aov$`Pr(>F)`[1], lmeUni.R2)
            }else{
              sv = c(sp, gg, y, x, (lmUni.aov['x', c('F value', 'Pr(>F)')]) %>% unlist(), lmUni.summary$r.squared %>% as.numeric(),
                     formula_char, NA, NA, NA)
            }
            uniRes = rbind(uniRes, sv)
          }
          if(class(tmpDF2$x) != 'numeric' & 
             length(na.omit(unique(tmpDF2$x))) >= 2 & length(na.omit(unique(tmpDF2$y))) >= 2 &
             nrow(tmpDF2) != length(unique(tmpDF2$x))){
            mm = aggregate(y ~ x, tmpDF2, mean)
            mmed = aggregate(y ~ x, tmpDF2, median)
            size = table(tmpDF2$x) %>% as.data.frame()
            # do pairwise comparison
            tmp_comp_DF = compare_means(y ~ x, tmpDF2, p.adjust.method = 'BH') %>% as.data.frame()
            tmp_comp_DF$y = y
            tmp_comp_DF$x = x
            tmp_comp_DF$species = sps
            tmp_comp_DF$gene = gg
            tmp_comp_DF$size1 = size$Freq[match(tmp_comp_DF$group1, size$Var1)]
            tmp_comp_DF$size2 = size$Freq[match(tmp_comp_DF$group2, size$Var1)]
            
            tmp_comp_DF$mean1 = mm$y[match(tmp_comp_DF$group1, mm$x)]
            tmp_comp_DF$mean2 = mm$y[match(tmp_comp_DF$group2, mm$x)]
            
            tmp_comp_DF$log2_FC = log2(tmp_comp_DF$mean1/tmp_comp_DF$mean2)
            
            tmp_comp_DF$med1 = mmed$y[match(tmp_comp_DF$group1, mmed$x)]
            tmp_comp_DF$med2 = mmed$y[match(tmp_comp_DF$group2, mmed$x)]
            
            tmp_comp_DF = tmp_comp_DF[, c('species', 'gene', 'y', 'x', 'group1', 'group2', 'size1', 'size2', 'mean1', 'mean2', 'log2_FC', 'med1', 'med2',
                                          'p', 'p.adj', 'p.format', 'p.signif')]
            comp_DF = rbind(comp_DF, tmp_comp_DF)
            
          }
        }
      }
  }
}

colnames(uniRes) = c('species', 'gene', 'y', 'x', 'lm.uni.F_value', 'lm.uni.P_value', 'lm.uni.R2',
                     'lme.uni.Model', 'lme.uni.F_value', 'lme.uni.P_value', 'lme.uni.R2')
uniRes$lm.uni.F_value %<>% as.numeric()
uniRes$lm.uni.P_value %<>% as.numeric()
uniRes$lm.uni.R2 %<>% as.numeric()
uniRes$lme.uni.F_value %<>% as.numeric()
uniRes$lme.uni.P_value %<>% as.numeric()
uniRes$lme.uni.R2 %<>% as.numeric()
if(plotFlag){
  xlsx::write.xlsx(comp_DF, file = './Routput/gene_wise_nucl_div_dnds_metadata_pairwise_comp_2.xlsx')  
  xlsx::write.xlsx(uniRes, file = './Routput/gene_wise_nucl_div_dnds_metadata_uni_model_2.xlsx')
}

unique(refDF_gene$genome)


#### mags ####
# read stb file
magSTB = read.delim2('./instrain_res/co_drep_mags.stb', sep = '\t', header = F, col.names = c('scaffold', 'genome'))
selCol = c("genome", "coverage", "breadth", "nucl_diversity", "length",
           "true_scaffolds", "detected_scaffolds", "coverage_median",
           "coverage_std", "coverage_SEM", "breadth_minCov",
           "breadth_expected", "nucl_diversity_rarefied",
           "conANI_reference", "popANI_reference", "SNS_count",
           "SNV_count", "reads_mean_PID", "reads_unfiltered_reads")
magDF = data.frame()
for(ss in swabData$fqID){
  print(ss)
  geneDat = read.delim2(paste0('./instrain_res/mags_profile/gene_info/', ss, '_bwa.IS_gene_info.tsv'))
  dt = fread(paste0('./instrain_res/mags_profile/gene_info/', ss, '_bwa_genes_SNP_count.csv.gz'), select = c('gene', 'S_sites', 'N_sites'))
  dt = dt[!duplicated(dt),]
  geneDat$S_sites = dt$S_sites[match(geneDat$gene, dt$gene)]
  geneDat$N_sites = dt$N_sites[match(geneDat$gene, dt$gene)]
  geneDat$SNS_S_count %<>% as.numeric()
  geneDat$SNS_N_count %<>% as.numeric()
  geneDat$SNV_S_count %<>% as.numeric()
  geneDat$SNV_N_count %<>% as.numeric()
  geneDat$breadth_minCov %<>% as.numeric()
  geneDat = read.delim2(paste0('./instrain_res/mags_profile/genome_info/', ss, '_bwa.IS_genome_info.tsv'))
  geneDat$breadth_minCov %<>% as.numeric()
  geneDat = geneDat[geneDat$breadth_minCov >= 0.1, selCol]
  if(nrow(geneDat) > 0){
    geneDat$dnds = NA
    geneDat$pnps = NA
    geneDat$dn_ds_common = NA
    for (gg in geneDat$genome) {
      scaffold_list = magSTB$scaffold[magSTB$genome == gg]
      geneDat_sub = geneDat[(geneDat$scaffold %in% scaffold_list) & geneDat$breadth_minCov >= 0.5, ]
      if(nrow(geneDat_sub) > 0){
        dn_ds = (sum(geneDat_sub$SNS_N_count, na.rm = T)/sum(geneDat_sub$N_sites, na.rm = T)) / (sum(geneDat_sub$SNS_S_count, na.rm = T)/sum(geneDat_sub$S_sites, na.rm = T))
        pn_ps = (sum(geneDat_sub$SNV_N_count, na.rm = T)/sum(geneDat_sub$N_sites, na.rm = T)) / (sum(geneDat_sub$SNV_S_count, na.rm = T)/sum(geneDat_sub$S_sites, na.rm = T))
        dn_ds_common = (sum(c(geneDat_sub$SNV_N_count, geneDat_sub$SNS_N_count), na.rm = T)/sum(geneDat_sub$N_sites, na.rm = T)) / (sum(c(geneDat_sub$SNV_S_count, geneDat_sub$SNS_S_count), na.rm = T)/sum(geneDat_sub$S_sites, na.rm = T))  
        geneDat$dnds[geneDat$genome == gg] = dn_ds
        geneDat$pnps[geneDat$genome == gg] = pn_ps
        geneDat$dn_ds_common[geneDat$genome == gg] = dn_ds_common
      }
    }
    geneDat$sample = ss
    magDF = rbind(magDF, geneDat)
  }
}
magDF$pt_ID = swabData$pt_ID[match(magDF$sample, swabData$fqID)]
magDF$pt_ID.u = swabData$pt_ID.u[match(magDF$sample, swabData$fqID)]
magDF$nucl_diversity %<>% as.numeric()
magDF$SNV_count %<>% as.numeric()
magDF$genome %<>% str_replace(., '.fna', '')
magDF$species = spName$lca_species_manual[match(magDF$genome, spName$Bin)]
magDF$type = 'MAG'
magDF$trimester = swabData$trimester[match(magDF$sample, swabData$fqID)]
magDF$week = swabData$Sample_GA[match(magDF$sample, swabData$fqID)]
magDF$length %<>% as.numeric()
magDF$SNV_count %<>% as.numeric()
magDF$SNPs_per_bp = magDF$SNV_count/(magDF$length*magDF$breadth_minCov)
save(magDF, file = './data/magDF.RData')

##### ref plot coverage/breath/nucleotide diversity/dnds/pnps ####
reflVarDF = refDF[, c('species', 'pt_ID.u', 'pt_ID', 'coverage', 'breadth', 'breadth_minCov', 'SNPs_per_bp', 'nucl_diversity', 'dnds', 'pnps', 'dn_ds_common', 'type')]
reflVarDF = melt(reflVarDF, id.vars = c('species', 'pt_ID.u', 'pt_ID', 'type'))
reflVarDF$value %<>% as.numeric()
reflVarDF$pt_ID.u %<>% factor(., levels = mixedsort(reflVarDF$pt_ID.u %>% unique))
if(F){
  pdf('./plot/ref_genome_info_with_dn_ds.pdf', width = 15, height = 10)
  for (sp in unique(reflVarDF$species)) {
    tmpDF = reflVarDF[reflVarDF$species == sp, ]
    p = ggplot(tmpDF, aes(x = pt_ID.u, y = value)) +
      geom_line(aes(group=pt_ID, color = pt_ID), linetype = 'dashed', size = 0.3) + 
      geom_point(aes(fill = pt_ID), shape = 24, color = 'black', stroke = 0.4, alpha = 0.5) +
      mytheme2 + theme(axis.text.x = element_text(angle = 45, size = 5, hjust = 1),
                       axis.text.y = element_text(size = 10, hjust = 1),
                       strip.text = element_text(size = 10), 
                       panel.grid = element_line(size = 0.3), panel.background = element_rect(fill = '#EEEEEE')) +
      labs(title = paste0(sp), x = 'Sample', y = '') +
      facet_wrap(.~variable, nrow = 8, ncol = 1, scales = 'free_y', strip.position = 'right')
    print(p)
  }
  dev.off()
}

##### mag plot coverage/breath/nucleotide diversity/dnds/pnps ####

maglVarDF = magDF[, c('species', 'pt_ID.u', 'pt_ID', 'coverage', 'breadth', 'breadth_minCov', 'SNPs_per_bp', 'nucl_diversity', 'dnds', 'pnps', 'dn_ds_common', 'type')]
maglVarDF = melt(maglVarDF, id.vars = c('species', 'pt_ID.u', 'pt_ID', 'type'))
maglVarDF$value %<>% as.numeric()
maglVarDF$pt_ID.u %<>% factor(., levels = mixedsort(maglVarDF$pt_ID.u %>% unique))
if(F){
  pdf('./plot/mag_genome_info_with_dn_ds.pdf', width = 15, height = 10)
  for (sp in unique(maglVarDF$species)) {
    tmpDF = maglVarDF[maglVarDF$species == sp, ]
    p = ggplot(tmpDF, aes(x = pt_ID.u, y = value)) +
      geom_line(aes(group=pt_ID, color = pt_ID), size = 0.3) + 
      geom_point(aes(fill = pt_ID), shape = 21, color = 'black', stroke = 0.4, alpha = 0.5) +
      mytheme2 + theme(axis.text.x = element_text(angle = 45, size = 5, hjust = 1),
                       axis.text.y = element_text(size = 10, hjust = 1),
                       strip.text = element_text(size = 10), 
                       panel.grid = element_line(size = 0.3), panel.background = element_rect(fill = '#EEEEEE')) +
      labs(title = paste0(sp), x = 'Sample', y = '') +
      facet_wrap(.~variable, nrow = 8, ncol = 1, scales = 'free_y', strip.position = 'right')
    print(p)
  }
  dev.off()
}

##### plot coverage/breath/nucleotide diversity/dnds/pnps ####
magDF$species %<>% str_replace_all(., ' ', '_')
impSP = intersect(magDF$species, refDF$species)
allDF = rbind(magDF[magDF$species %in% imƒpSP, ], refDF[refDF$species %in% impSP, ])
allDF$coverage_log10 = log10(allDF$coverage %>% as.numeric())
alllVarDF = allDF[, c('species', 'pt_ID.u', 'pt_ID', 'coverage_log10', 'breadth', 'breadth_minCov', 'SNPs_per_bp', 'nucl_diversity', 'dnds', 'pnps', 'dn_ds_common', 'type', 'trimester', 'week')]
alllVarDF = melt(alllVarDF, id.vars = c('species', 'pt_ID.u', 'pt_ID', 'type', 'trimester', 'week'))
alllVarDF$value %<>% as.numeric()
alllVarDF = alllVarDF[!grepl('KIT', alllVarDF$pt_ID.u), ]
alllVarDF$person_trimester = paste0(
  alllVarDF$pt_ID, '_',
  alllVarDF$trimester, '_',
  alllVarDF$week
)
if(plotFlag){
  pdf('./plot/ref_and_mags_genome_info_with_dn_ds_new.pdf', width = 15, height = 10)
  labDF = data.frame(
    label = c('Coverage (log10-transformed)', 'Breadth', 'Breadth (coverage>=5)', 'SNPs/bp', 'Nucleotide diversity', 'dN/dS', 'pN/pS', 'dN/dS (common)'),
    variable = paste0(table(alllVarDF$variable) %>% names())
  )
  
  library(randomcoloR)
  n <- 35
  # palette <- distinctColorPalette(n)
  pie(rep(1, n), col=palette)
  mycc = palette
  names(mycc) = unique(alllVarDF$pt_ID) %>% sort()
  labDF$variable %<>% factor(., labDF$variable)
  for (sp in unique(alllVarDF$species)) {
    tmpDF = alllVarDF[alllVarDF$species == sp, ]
    tmpDF = tmpDF[order(tmpDF$pt_ID, tmpDF$week), ]
    tmpDF$person_trimester %<>% factor(., levels = unique(tmpDF$person_trimester))
    tmpDF$variable %<>% factor(., labDF$variable)
    p = ggplot(tmpDF, aes(x = person_trimester, y = value)) +
      geom_line(aes(group=interaction(pt_ID, type), color = pt_ID, linetype = type), size = 0.3) + 
      geom_point(aes(fill = pt_ID, shape = type), color = 'black', stroke = 0.4, alpha = 0.5) +
      scale_shape_manual(values = c(21, 24)) +
      scale_color_manual(values = mycc[match(unique(tmpDF$pt_ID) %>% sort(), names(mycc))] %>% darken(., 1.2) %>% alpha(., 0.7)) +
      scale_fill_manual(values =  mycc[match(unique(tmpDF$pt_ID) %>% sort(), names(mycc))]) +
      scale_y_continuous(expand = expansion(mult = c(0.2, .35))) +
      mytheme2 + theme(axis.text.x = element_text(angle = 45, size = 5, hjust = 1),
                       axis.text.y = element_text(size = 10, hjust = 1),
                       plot.title = element_text(hjust = 0.5),
                       strip.background = element_blank(),
                       # strip.text = element_text(size = 10, hjust = 0),
                       strip.text = element_blank(),
                       panel.grid = element_line(size = 0.3), panel.background = element_rect(fill = '#EEEEEE')) +
      labs(title = paste0(sp), x = 'Sample', y = '', fill = 'Participant', color = 'Participant',
           shape = 'Reference type', linetype = 'Reference type') +
      facet_wrap(.~variable, nrow = 8, ncol = 1, scales = 'free_y', strip.position = 'top') +
      geom_text(
        data    = labDF,
        mapping = aes(x =  -Inf+1, y = -Inf, label = label),
        hjust   = 0,
        vjust   = -8.2, size = 4
      )
    print(p)
  }
  dev.off()
}





