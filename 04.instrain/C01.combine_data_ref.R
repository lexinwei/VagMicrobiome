source('~/workspace/library_function.R')
#### prepare the metadata for MTMG ####
if(F){
  metaData0 = read.xlsx('./public_cohort/!!!PRJNA797778_MTMG_longitudinal/!!!media-6_ethnicity_from_mgCST_paper.xlsx', 
                       startRow = 2)
  metaData1 = metaData0[metaData0$BioProject == 'PRJNA797778', ]
  table(metaData1$UID)
  sampleData0 = read.xlsx('./public_cohort/!!!PRJNA797778_MTMG_longitudinal/!!!filereport_read_run_PRJNA797778_tsv.xlsx')
  sampleData1 = sampleData0[sampleData0$library_strategy == 'WGS' , ]
  table(sampleData1$sample_alias)
  length(intersect(metaData1$UID, sampleData1$sample_alias))
  setdiff(metaData1$UID, sampleData1$sample_alias)
  setdiff(sampleData1$sample_alias, metaData1$UID)
  metaData1$seq_ID = sampleData1$run_accession[match(metaData1$UID, sampleData1$sample_alias)]
  metaData2 = metaData1[!is.na(metaData1$seq_ID), ]
  table(metaData2$Race)
  metaData2$Ethnictiy = NA
  metaData2$Ethnictiy[metaData2$Race == 'Asian'] = 'Asian'
  metaData2$Ethnictiy[metaData2$Race == 'Black or African American'] = 'Black'
  metaData2$Ethnictiy[metaData2$Race == 'Hispanic or Latino'] = 'Latina'
  metaData2$Ethnictiy[metaData2$Race == 'White or Caucasian'] = 'White'
  table(metaData2$Ethnictiy)
  metaData2$Ethnictiy[is.na(metaData2$Ethnictiy)] = 'Other'
  metaData2$Ethnictiy %<>% factor(., levels = c('White', 'Asian', 'Latina', 'Black', 'Other'))
  unique(metaData2$UID) %>% length()
  unique(metaData2$SID) %>% length()
  metaData = metaData2
  save(metaData, file = './data/RData/MTMG_metaData.RData')
}

#### combine profile of genome and gene levels ####
#================== calculate genome wide dnds / pnps ===================#
# read stb file
if(T){
  refSTB = read.delim2('./instrain_res/NCBI_ref/reference_genomes.stb', sep = '\t', header = F, col.names = c('scaffold', 'genome'))
  refDF = data.frame()
  for(ss in sampleData1$run_accession){
    print(ss)
    geneDat = read.delim2(paste0('./instrain_res/MTMG/gene_info/', ss, '.IS_gene_info.tsv'))
    dt = fread(paste0('./instrain_res/MTMG/genes_SNP_count/', ss, '_genes_SNP_count.csv.gz'), 
               select = c('gene', 'S_sites', 'N_sites'))
    dt = dt[!duplicated(dt), ]
    geneDat$S_sites = dt$S_sites[match(geneDat$gene, dt$gene)]
    geneDat$N_sites = dt$N_sites[match(geneDat$gene, dt$gene)]
    geneDat$SNS_S_count %<>% as.numeric()
    geneDat$SNS_N_count %<>% as.numeric()
    geneDat$SNV_S_count %<>% as.numeric()
    geneDat$SNV_N_count %<>% as.numeric()
    geneDat$breadth_minCov %<>% as.numeric()
    genomeDat = read.delim2(paste0('./instrain_res/MTMG/genome_info/', ss, '.IS_genome_info.tsv'))
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
      refDF = rbind(refDF, genomeDat)
    }
  }
  refDF$pt_ID.u = metaData2$UID[match(refDF$sample, metaData2$seq_ID)]
  refDF$nucl_diversity %<>% as.numeric()
  refDF$SNV_count %<>% as.numeric()
  refDF$species = str_replace_all(refDF$genome, '\\.fasta|_|s__', ' ') %>% str_trim()
  table(refDF$species)
  refDF$type = 'NCBI'
  refDF$length %<>% as.numeric()
  refDF$SNV_count %<>% as.numeric()
  refDF$SNPs_per_bp = refDF$SNV_count/(refDF$length*refDF$breadth_minCov)
  refDF$Ethnicity = metaData2$Ethnictiy[match(refDF$sample, metaData2$seq_ID)]
}
save(refDF, file = './data/RData/MTMG_refDF_minCovBreath10.RData')
##### ref plot coverage/breath/nucleotide diversity/dnds/pnps ####
load('./data/RData/MTMG_refDF_minCovBreath10.RData')
refDF$coverage_log10 = log10(refDF$coverage %>% as.numeric())

refVarDF = refDF[, c('species', 'pt_ID.u', 'coverage_log10', 'breadth', 'SNPs_per_bp', 'nucl_diversity', 'dn_ds_common', 'Ethnicity')]
refVarDF = melt(refVarDF, id.vars = c('species', 'pt_ID.u', 'Ethnicity'))
refVarDF$value %<>% as.numeric()

if(plotFlag){
  pdf('./plot/MTMG_ref_genome_info_with_dn_ds_new.pdf', width = 15, height = 8)
  labDF = data.frame(
    label = c('Coverage (log10-transformed)', 'Breadth', 'SNPs/bp', 'Nucleotide diversity', 'dN/dS'),
    variable = paste0(table(refVarDF$variable) %>% names())
  )
  
  labDF$variable %<>% factor(., labDF$variable)
  for (sp in unique(refVarDF$species)) {
    tmpDF = refVarDF[refVarDF$species == sp, ]
    tmpDF$variable %<>% factor(., labDF$variable)
    p = ggplot(tmpDF, aes(x = pt_ID.u, y = value)) +
      geom_point(shape = 21, color = 'black', stroke = 0.4, alpha = 0.5) +
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
        vjust   = -8.5, size = 4
      )
    print(p)
  }
  dev.off()
}
