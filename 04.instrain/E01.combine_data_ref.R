source('~/workspace/library_function.R')
#### prepare the metadata for MTMG ####
if(T){
  sampleData0 = read.table('./public_cohort/!!!susan_PNAS_PRJNA288562/PRJNA288562.tsv', sep = '\t', header = T)
  table(sampleData0$scientific_name)
  nrow(sampleData0)
  sampleData0$subjectID = paste0('Subject_', str_split_fixed(sampleData0$sample_title, ' ', 3)[, 2])
  sampleData0$gestational_day = str_extract(sampleData0$sample_title, '\\d+$') %>% as.numeric()
  sampleData0$body_site = str_split_fixed(sampleData0$sample_title, ' ', 3)[, 3] %>% str_remove(., ' collected on gestational day \\d+')
  sampleData0$body_site[sampleData0$body_site == 'distal gut specimen'] = 'Distal_Gut'
  sampleData0$body_site[sampleData0$body_site == 'saliva'] = 'Saliva'
  sampleData0$body_site[sampleData0$body_site == 'vaginal swab'] = 'Vaginal_Swab'
  table(sampleData0$body_site)
  sampleData0$uid = paste0(sampleData0$subjectID, '_',  sampleData0$body_site, '_', sampleData0$gestational_day)
  
  metaData0 = read.xlsx('./public_cohort/!!!susan_PNAS_PRJNA288562/pnas.1502875112.sd02.xlsx')
  metaData0$uid = paste0(metaData0$Subject_ID, '_', metaData0$BodySite, '_', metaData0$GestationalDayOfCollection)
  metaData0$run_accession = sampleData0$run_accession[match(metaData0$uid, sampleData0$uid)]
  length(metaData0$run_accession[!is.na(metaData0$run_accession)] %>% unique())
  table(metaData0$Race, metaData0$Subject_ID)
  table(metaData0$Ethnicity)
  
  metaData1 = metaData0[!is.na(metaData0$run_accession), ]
  table(metaData1$BodySite)
  metaData2 = metaData1[metaData1$BodySite == 'Vaginal_Swab', ]
  metaData3 = metaData2[!endsWith(metaData2$`#SampleID`, '.rs'), ]
  table(metaData3$)
  table(metaData3$PretermBirth)
  table(metaData3$Race)
  table(metaData3$Ethnicity)
  table(metaData3$Ethnicity, metaData3$Subject_ID)
  table(metaData3$Race, metaData3$Subject_ID)
  table(metaData3$Race, metaData3$Ethnicity)
  table(metaData3$Subject_ID)
  metaData3$Subject_ID[metaData3$Race == "Caucasian" & metaData3$Ethnicity == 'hispanic'] %>% unique() %>% length()
  metaData3$Subject_ID[metaData3$Race == "Caucasian" & metaData3$Ethnicity == 'non-hispanic'] %>% unique() %>% length()
}

#### combine profile of genome and gene levels ####
#================== calculate genome wide dnds / pnps ===================#
# read stb file
if(T){
  refSTB = read.delim2('./instrain_res/NCBI_ref/reference_genomes.stb', sep = '\t', header = F, col.names = c('scaffold', 'genome'))
  refDF = data.frame()
  for(ss in metaData3$run_accession[c(1:3, 5:length(metaData3$run_accession))]){
    print(ss)
    geneDat = read.delim2(paste0('./instrain_res/PRJNA288562/gene_info/', ss, '.IS_gene_info.tsv'))
    dt = fread(paste0('./instrain_res/PRJNA288562/genes_SNP_count/', ss, '_genes_SNP_count.csv.gz'), 
               select = c('gene', 'S_sites', 'N_sites'))
    dt = dt[!duplicated(dt), ]
    geneDat$S_sites = dt$S_sites[match(geneDat$gene, dt$gene)]
    geneDat$N_sites = dt$N_sites[match(geneDat$gene, dt$gene)]
    geneDat$SNS_S_count %<>% as.numeric()
    geneDat$SNS_N_count %<>% as.numeric()
    geneDat$SNV_S_count %<>% as.numeric()
    geneDat$SNV_N_count %<>% as.numeric()
    geneDat$breadth_minCov %<>% as.numeric()
    genomeDat = read.delim2(paste0('./instrain_res/PRJNA288562/genome_info/', ss, '.IS_genome_info.tsv'))
    genomeDat$breadth_minCov %<>% as.numeric()
    genomeDat = genomeDat[genomeDat$breadth_minCov >= 0.1, ] # cutoff
    if(nrow(genomeDat) > 0){
      for (gg in genomeDat$genome) {
        scaffold_list = refSTB$scaffold[refSTB$genome == gg]
        geneDat_sub0 = geneDat[(geneDat$scaffold %in% scaffold_list), ] 
        geneDat_sub = geneDat[(geneDat$scaffold %in% scaffold_list) & geneDat$breadth_minCov >= 0.5, ] # cutoff
        if(nrow(geneDat_sub) > 0){
          gene_number_before_filter = nrow(geneDat_sub0)
          SNS_N_count_before_filter = sum(geneDat_sub0$SNS_N_count, na.rm = T)
          SNS_S_count_before_filter = sum(geneDat_sub0$SNS_S_count, na.rm = T)
          SNV_N_count_before_filter = sum(geneDat_sub0$SNV_N_count, na.rm = T)
          SNV_S_count_before_filter = sum(geneDat_sub0$SNV_S_count, na.rm = T)
          SNP_N_count_before_filter = sum(c(geneDat_sub0$SNV_N_count, geneDat_sub0$SNS_N_count), na.rm = T)
          SNP_S_count_before_filter = sum(c(geneDat_sub0$SNV_S_count, geneDat_sub0$SNS_S_count), na.rm = T)
          N_sites_before_filter = sum(geneDat_sub0$N_sites, na.rm = T)
          S_sites_before_filter = sum(geneDat_sub0$S_sites, na.rm = T)
          dn_ds_before_filter = (sum(geneDat_sub0$SNS_N_count, na.rm = T)/sum(geneDat_sub0$N_sites, na.rm = T)) / (sum(geneDat_sub0$SNS_S_count, na.rm = T)/sum(geneDat_sub0$S_sites, na.rm = T))
          pn_ps_before_filter = (sum(geneDat_sub0$SNV_N_count, na.rm = T)/sum(geneDat_sub0$N_sites, na.rm = T)) / (sum(geneDat_sub0$SNV_S_count, na.rm = T)/sum(geneDat_sub0$S_sites, na.rm = T))
          dn_ds_common_before_filter = (sum(c(geneDat_sub0$SNV_N_count, geneDat_sub0$SNS_N_count), na.rm = T)/sum(geneDat_sub0$N_sites, na.rm = T)) / (sum(c(geneDat_sub0$SNV_S_count, geneDat_sub0$SNS_S_count), na.rm = T)/sum(geneDat_sub0$S_sites, na.rm = T))  
          
          genomeDat$gene_number_before_filter[genomeDat$genome == gg] = gene_number_before_filter
          genomeDat$SNS_N_count_before_filter[genomeDat$genome == gg] = SNS_N_count_before_filter
          genomeDat$SNS_S_count_before_filter[genomeDat$genome == gg] = SNS_S_count_before_filter
          genomeDat$SNV_N_count_before_filter[genomeDat$genome == gg] = SNV_N_count_before_filter
          genomeDat$SNV_S_count_before_filter[genomeDat$genome == gg] = SNV_S_count_before_filter
          genomeDat$SNP_N_count_before_filter[genomeDat$genome == gg] = SNP_N_count_before_filter
          genomeDat$SNP_S_count_before_filter[genomeDat$genome == gg] = SNP_S_count_before_filter
          genomeDat$N_sites_before_filter[genomeDat$genome == gg] = N_sites_before_filter
          genomeDat$S_sites_before_filter[genomeDat$genome == gg] = S_sites_before_filter
          genomeDat$dnds_before_filter[genomeDat$genome == gg] = dn_ds_before_filter
          genomeDat$pnps_before_filter[genomeDat$genome == gg] = pn_ps_before_filter
          genomeDat$dn_ds_common_before_filter[genomeDat$genome == gg] = dn_ds_common_before_filter
          
          gene_number = nrow(geneDat_sub)
          SNS_N_count = sum(geneDat_sub$SNS_N_count, na.rm = T)
          SNS_S_count = sum(geneDat_sub$SNS_S_count, na.rm = T)
          SNV_N_count = sum(geneDat_sub$SNV_N_count, na.rm = T)
          SNV_S_count = sum(geneDat_sub$SNV_S_count, na.rm = T)
          SNP_N_count = sum(c(geneDat_sub$SNV_N_count, geneDat_sub$SNS_N_count), na.rm = T)
          SNP_S_count = sum(c(geneDat_sub$SNV_S_count, geneDat_sub$SNS_S_count), na.rm = T)
          N_sites = sum(geneDat_sub$N_sites, na.rm = T)
          S_sites = sum(geneDat_sub$S_sites, na.rm = T)
          dn_ds = (sum(geneDat_sub$SNS_N_count, na.rm = T)/sum(geneDat_sub$N_sites, na.rm = T)) / (sum(geneDat_sub$SNS_S_count, na.rm = T)/sum(geneDat_sub$S_sites, na.rm = T))
          pn_ps = (sum(geneDat_sub$SNV_N_count, na.rm = T)/sum(geneDat_sub$N_sites, na.rm = T)) / (sum(geneDat_sub$SNV_S_count, na.rm = T)/sum(geneDat_sub$S_sites, na.rm = T))
          dn_ds_common = (sum(c(geneDat_sub$SNV_N_count, geneDat_sub$SNS_N_count), na.rm = T)/sum(geneDat_sub$N_sites, na.rm = T)) / (sum(c(geneDat_sub$SNV_S_count, geneDat_sub$SNS_S_count), na.rm = T)/sum(geneDat_sub$S_sites, na.rm = T))  
          genomeDat$gene_number[genomeDat$genome == gg] = gene_number
          genomeDat$SNS_N_count[genomeDat$genome == gg] = SNS_N_count
          genomeDat$SNS_S_count[genomeDat$genome == gg] = SNS_S_count
          genomeDat$SNV_N_count[genomeDat$genome == gg] = SNV_N_count
          genomeDat$SNV_S_count[genomeDat$genome == gg] = SNV_S_count
          genomeDat$SNP_N_count[genomeDat$genome == gg] = SNP_N_count
          genomeDat$SNP_S_count[genomeDat$genome == gg] = SNP_S_count
          genomeDat$N_sites[genomeDat$genome == gg] = N_sites
          genomeDat$S_sites[genomeDat$genome == gg] = S_sites
          genomeDat$dnds[genomeDat$genome == gg] = dn_ds
          genomeDat$pnps[genomeDat$genome == gg] = pn_ps
          genomeDat$dn_ds_common[genomeDat$genome == gg] = dn_ds_common
        }
      }
      genomeDat$sample = ss
      if(!"linked_SNV_count" %in% colnames(genomeDat)){genomeDat %<>% add_column(., .after = "iRep_GC_corrected", linked_SNV_count = NA)}
      if(!"SNV_distance_mean" %in% colnames(genomeDat)){genomeDat %<>% add_column(., .after = "iRep_GC_corrected", SNV_distance_mean = NA)}
      if(!"r2_mean" %in% colnames(genomeDat)){genomeDat %<>% add_column(., .after = "iRep_GC_corrected", r2_mean = NA)}
      if(!"d_prime_mean" %in% colnames(genomeDat)){genomeDat %<>% add_column(., .after = "iRep_GC_corrected", d_prime_mean = NA)}
      refDF = rbind(refDF, genomeDat)
      setdiff(colnames(refDF), colnames(genomeDat))
    }
  }
  refDF$uid = metaData3$uid[match(refDF$sample, metaData3$run_accession)]
  refDF$nucl_diversity %<>% as.numeric()
  refDF$SNV_count %<>% as.numeric()
  refDF$species = str_replace_all(refDF$genome, '\\.fasta|_|s__', ' ') %>% str_trim()
  table(refDF$species)
  refDF$type = 'NCBI'
  refDF$length %<>% as.numeric()
  refDF$SNV_count %<>% as.numeric()
  refDF$SNPs_per_bp = refDF$SNV_count/(refDF$length*refDF$breadth_minCov)
  refDF$Ethnicity = metaData3$Ethnicity[match(refDF$sample, metaData3$run_accession)]
  refDF$Race = metaData3$Race[match(refDF$sample, metaData3$run_accession)]
}
save(refDF, file = './data/RData/PRJNA288562_refDF_minCovBreath10.RData')
