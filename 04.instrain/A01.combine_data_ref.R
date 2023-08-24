source('~/workspace/library_function.R')
load('./data/RData/metaData.RData')
load('./data/RData/swabData.RData')
load('./data/RData/metaphlanDF.RData')
#### combine profile of genome and gene levels ####
#================== calculate genome wide dnds / pnps ===================#
# read stb file
if(T){
  refSTB = read.delim2('./instrain_res/NCBI_ref/reference_genomes.stb', sep = '\t', header = F, col.names = c('scaffold', 'genome'))
  refDF = data.frame()
  for(ss in swabData$fqID){
    print(ss)
    geneDat = read.delim2(paste0('./instrain_res/NCBI_ref/gene_info/', ss, '_bwa.IS_gene_info.tsv'))
    dt = fread(paste0('./instrain_res/NCBI_ref/genes_SNP_count/', ss, '_bwa_genes_SNP_count.csv.gz'), 
               select = c('gene', 'S_sites', 'N_sites'))
    dt = dt[!duplicated(dt), ]
    geneDat$S_sites = dt$S_sites[match(geneDat$gene, dt$gene)]
    geneDat$N_sites = dt$N_sites[match(geneDat$gene, dt$gene)]
    geneDat$SNS_S_count %<>% as.numeric()
    geneDat$SNS_N_count %<>% as.numeric()
    geneDat$SNV_S_count %<>% as.numeric()
    geneDat$SNV_N_count %<>% as.numeric()
    geneDat$breadth_minCov %<>% as.numeric()
    genomeDat = read.delim2(paste0('./instrain_res/NCBI_ref/genome_info/', ss, '_bwa.IS_genome_info.tsv'))
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
  refDF$pt_ID = swabData$pt_ID[match(refDF$sample, swabData$fqID)]
  refDF$pt_ID.u = swabData$pt_ID.u[match(refDF$sample, swabData$fqID)]
  refDF$nucl_diversity %<>% as.numeric()
  refDF$SNV_count %<>% as.numeric()
  refDF$species = str_replace(refDF$genome, '.fna', '')
  table(refDF$species)
  refDF$type = 'NCBI'
  refDF$trimester = swabData$trimester[match(refDF$sample, swabData$fqID)]
  refDF$week = swabData$Sample_GA[match(refDF$sample, swabData$fqID)]
  refDF$length %<>% as.numeric()
  refDF$SNV_count %<>% as.numeric()
  refDF$SNPs_per_bp = refDF$SNV_count/(refDF$length*refDF$breadth_minCov)
}
# save(refDF, file = './data/RData/refDF_minCovBreath10.RData')
# save(refDF, file = './data/RData/refDF_minCovBreath50.RData')
save(refDF, file = './data/RData/refDF_minCovBreath10.RData')
##### ref plot coverage/breath/nucleotide diversity/dnds/pnps ####
load('./data/RData/refDF_minCovBreath10.RData')
refDF$coverage_log10 = log10(refDF$coverage %>% as.numeric())
refDF$relative_abundance = metaphlanDF$relative_abundance[match(paste0(refDF$species, '|', refDF$pt_ID.u), 
                                                                paste0(metaphlanDF$species, '|', metaphlanDF$sample_id))]
refVarDF = refDF[, c('species', 'pt_ID.u', 'pt_ID', 'coverage_log10', 'breadth', 'SNPs_per_bp', 'nucl_diversity', 'dn_ds_common', 'relative_abundance', 'trimester', 'week')]
refVarDF = melt(refVarDF, id.vars = c('species', 'pt_ID.u', 'pt_ID', 'trimester', 'week'))
refVarDF$value %<>% as.numeric()
refVarDF = refVarDF[!grepl('KIT', refVarDF$pt_ID.u), ]
refVarDF$person_trimester = paste0(
  refVarDF$pt_ID, '_',
  refVarDF$trimester, '_',
  refVarDF$week
)
if(plotFlag){
  pdf('./plot/ref_genome_info_with_dn_ds_new.pdf', width = 15, height = 8)
  labDF = data.frame(
    label = c('Coverage (log10-transformed)', 'Breadth', 'SNPs/bp', 'Nucleotide diversity', 'dN/dS', 'Relative abundance'),
    variable = paste0(table(refVarDF$variable) %>% names())
  )
  
  labDF$variable %<>% factor(., labDF$variable)
  for (sp in unique(refVarDF$species)) {
    tmpDF = refVarDF[refVarDF$species == sp, ]
    tmpDF = tmpDF[order(tmpDF$pt_ID, tmpDF$week), ]
    tmpDF$person_trimester %<>% factor(., levels = unique(tmpDF$person_trimester))
    tmpDF$variable %<>% factor(., labDF$variable)
    p = ggplot(tmpDF, aes(x = person_trimester, y = value)) +
      geom_line(aes(group = pt_ID, color = pt_ID),  size = 0.3) + 
      geom_point(aes(fill = pt_ID),shape = 21, color = 'black', stroke = 0.4, alpha = 0.5) +
      scale_color_manual(values = colorParticipant) +
      scale_fill_manual(values =  colorParticipant) +
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
