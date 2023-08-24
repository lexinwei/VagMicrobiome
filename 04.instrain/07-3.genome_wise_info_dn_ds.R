source('~/workspace/library_function.R')
load('./data/metaData.RData')
load('./data/swabData.RData')
load('./data/spName.RData')
load('./data/refDF.RData')
library(randomcoloR)
n <- 35
# person_palette <- distinctColorPalette(n)
pie(rep(1, n), col=person_palette)
names(person_palette) = metaData$pt_ID
# save(person_palette, file = './data/person_palette.RData')
mycc = person_palette
#================== calculate genome wide dnds / pnps ===================#
#### ref ####
# read stb file
refSTB = read.delim2('./instrain_res/reference_genomes.stb', sep = '\t', header = F, col.names = c('scaffold', 'genome'))
selCol = c("genome", "coverage", "breadth", "nucl_diversity", "length",
           "true_scaffolds", "detected_scaffolds", "coverage_median",
           "coverage_std", "coverage_SEM", "breadth_minCov",
           "breadth_expected", "nucl_diversity_rarefied",
           "conANI_reference", "popANI_reference", "SNS_count",
           "SNV_count", "reads_mean_PID", "reads_unfiltered_reads")
refDF = data.frame()
for(ss in swabData$fqID){
  print(ss)
  geneDat = read.delim2(paste0('./instrain_res/ref_profile/gene_info/', ss, '_bwa.IS_gene_info.tsv'))

  dt = fread(paste0('./instrain_res/ref_profile/gene_info/', ss, '_bwa_genes_SNP_count.csv.gz'), select = c('gene', 'S_sites', 'N_sites'))
  dt = dt[!duplicated(dt), ]
  geneDat$S_sites = dt$S_sites[match(geneDat$gene, dt$gene)]
  geneDat$N_sites = dt$N_sites[match(geneDat$gene, dt$gene)]
  geneDat$SNS_S_count %<>% as.numeric()
  geneDat$SNS_N_count %<>% as.numeric()
  geneDat$SNV_S_count %<>% as.numeric()
  geneDat$SNV_N_count %<>% as.numeric()
  geneDat$breadth_minCov %<>% as.numeric()
  genomeDat = read.delim2(paste0('./instrain_res/ref_profile/genome_info/', ss, '_bwa.IS_genome_info.tsv'))
  genomeDat$breadth_minCov %<>% as.numeric()
  genomeDat = genomeDat[genomeDat$breadth_minCov >= 0.1, selCol]
  if(nrow(genomeDat) > 0){
    genomeDat$dnds = NA
    genomeDat$pnps = NA
    genomeDat$dn_ds_common = NA
    for (gg in genomeDat$genome) {
      scaffold_list = refSTB$scaffold[refSTB$genome == gg]
      geneDat_sub = geneDat[(geneDat$scaffold %in% scaffold_list) & geneDat$breadth_minCov >= 0.5, ]
      if(nrow(geneDat_sub) > 0){
        dn_ds = (sum(geneDat_sub$SNS_N_count, na.rm = T)/sum(geneDat_sub$N_sites, na.rm = T)) / (sum(geneDat_sub$SNS_S_count, na.rm = T)/sum(geneDat_sub$S_sites, na.rm = T))
        pn_ps = (sum(geneDat_sub$SNV_N_count, na.rm = T)/sum(geneDat_sub$N_sites, na.rm = T)) / (sum(geneDat_sub$SNV_S_count, na.rm = T)/sum(geneDat_sub$S_sites, na.rm = T))
        dn_ds_common = (sum(c(geneDat_sub$SNV_N_count, geneDat_sub$SNS_N_count), na.rm = T)/sum(geneDat_sub$N_sites, na.rm = T)) / (sum(c(geneDat_sub$SNV_S_count, geneDat_sub$SNS_S_count), na.rm = T)/sum(geneDat_sub$S_sites, na.rm = T))  
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
save(refDF, file = './data/refDF.RData')
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
  genomeDat = read.delim2(paste0('./instrain_res/mags_profile/genome_info/', ss, '_bwa.IS_genome_info.tsv'))
  genomeDat$breadth_minCov %<>% as.numeric()
  genomeDat = genomeDat[genomeDat$breadth_minCov >= 0.1, selCol]
  if(nrow(genomeDat) > 0){
    genomeDat$dnds = NA
    genomeDat$pnps = NA
    genomeDat$dn_ds_common = NA
    for (gg in genomeDat$genome) {
      scaffold_list = magSTB$scaffold[magSTB$genome == gg]
      geneDat_sub = geneDat[(geneDat$scaffold %in% scaffold_list) & geneDat$breadth_minCov >= 0.5, ]
      if(nrow(geneDat_sub) > 0){
        dn_ds = (sum(geneDat_sub$SNS_N_count, na.rm = T)/sum(geneDat_sub$N_sites, na.rm = T)) / (sum(geneDat_sub$SNS_S_count, na.rm = T)/sum(geneDat_sub$S_sites, na.rm = T))
        pn_ps = (sum(geneDat_sub$SNV_N_count, na.rm = T)/sum(geneDat_sub$N_sites, na.rm = T)) / (sum(geneDat_sub$SNV_S_count, na.rm = T)/sum(geneDat_sub$S_sites, na.rm = T))
        dn_ds_common = (sum(c(geneDat_sub$SNV_N_count, geneDat_sub$SNS_N_count), na.rm = T)/sum(geneDat_sub$N_sites, na.rm = T)) / (sum(c(geneDat_sub$SNV_S_count, geneDat_sub$SNS_S_count), na.rm = T)/sum(geneDat_sub$S_sites, na.rm = T))  
        genomeDat$dnds[genomeDat$genome == gg] = dn_ds
        genomeDat$pnps[genomeDat$genome == gg] = pn_ps
        genomeDat$dn_ds_common[genomeDat$genome == gg] = dn_ds_common
      }
    }
    genomeDat$sample = ss
    magDF = rbind(magDF, genomeDat)
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
load('./data/refDF.RData')
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
        vjust   = -8.5, size = 4
      )
    print(p)
  }
  dev.off()
}


##### mag plot coverage/breath/nucleotide diversity/dnds/pnps ####
load('./data/magDF.RData')
magDF$coverage_log10 = log10(magDF$coverage %>% as.numeric())
magDF$species = str_replace(magDF, '', '_')
magDF$relative_abundance = metaphlanDF$relative_abundance[match(paste0(magDF$species, '|', magDF$pt_ID.u), 
                                                                paste0(metaphlanDF$species, '|', metaphlanDF$sample_id))]
magVarDF = magDF[, c('species', 'pt_ID.u', 'pt_ID', 'coverage_log10', 'breadth', 'SNPs_per_bp', 'nucl_diversity', 'dn_ds_common', 'relative_abundance', 'trimester', 'week')]
magVarDF = melt(magVarDF, id.vars = c('species', 'pt_ID.u', 'pt_ID', 'trimester', 'week'))
magVarDF$value %<>% as.numeric()
magVarDF = magVarDF[!grepl('KIT', magVarDF$pt_ID.u), ]
magVarDF$person_trimester = paste0(
  magVarDF$pt_ID, '_',
  magVarDF$trimester, '_',
  magVarDF$week
)
if(plotFlag){
  pdf('./plot/mag_genome_info_with_dn_ds_new.pdf', width = 15, height = 8)
  labDF = data.frame(
    label = c('Coverage (log10-transformed)', 'Breadth', 'SNPs/bp', 'Nucleotide diversity', 'dN/dS', 'Relative abundance'),
    variable = paste0(table(magVarDF$variable) %>% names())
  )
  
  
  labDF$variable %<>% factor(., labDF$variable)
  for (sp in unique(magVarDF$species)) {
    tmpDF = magVarDF[magVarDF$species == sp, ]
    tmpDF = tmpDF[order(tmpDF$pt_ID, tmpDF$week), ]
    tmpDF$person_trimester %<>% factor(., levels = unique(tmpDF$person_trimester))
    tmpDF$variable %<>% factor(., labDF$variable)
    p = ggplot(tmpDF, aes(x = person_trimester, y = value)) +
      geom_line(aes(group = pt_ID, color = pt_ID),linetype = 'dashed', size = 0.3) + 
      geom_point(aes(fill = pt_ID),shape = 24, color = 'black', stroke = 0.4, alpha = 0.5) +
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
        vjust   = -8.5, size = 4
      )
    print(p)
  }
  dev.off()
}


##### all plot coverage/breath/nucleotide diversity/dnds/pnps ####
source('~/workspace/library_function.R')
load('./data/refDF.RData')
load('./data/magDF.RData')
load('./data/metaphlanDF.RData')
load('./data/col_vector.RData')
magDF$species %<>% str_replace_all(., ' ', '_')
impSP = intersect(magDF$species, refDF$species)
allDF = rbind(magDF[magDF$species %in% impSP, ], refDF[refDF$species %in% impSP, ])
allDF$coverage_log10 = log10(allDF$coverage %>% as.numeric())
allDF$coverage = allDF$coverage %>% as.numeric()
allDF$relative_abundance = metaphlanDF$relative_abundance[match(paste0(allDF$species, '|', allDF$pt_ID.u), 
                                                                paste0(metaphlanDF$species, '|', metaphlanDF$sample_id))]
allDF$relative_abundance[allDF$type == 'MAG'] = NA
alllVarDF = allDF[, c('species', 'pt_ID.u', 'pt_ID', 'coverage_log10', 'breadth', 'SNPs_per_bp', 'nucl_diversity', 'dn_ds_common', 'relative_abundance', 'type', 'trimester', 'week')]
alllVarDF = melt(alllVarDF, id.vars = c('species', 'pt_ID.u', 'pt_ID', 'type', 'trimester', 'week'))
alllVarDF$value %<>% as.numeric()
alllVarDF = alllVarDF[!grepl('KIT', alllVarDF$pt_ID.u), ]
alllVarDF$person_trimester = paste0(
  alllVarDF$pt_ID, '_',
  alllVarDF$trimester, '_',
  alllVarDF$week
)
alllVarDF$type = factor(alllVarDF$type, levels = c('NCBI', 'MAG'))
if(plotFlag){
  pdf('./plot/ref_and_mags_genome_info_with_dn_ds_new.pdf', width = 15, height = 8)
  labDF = data.frame(
    label = c('Coverage (log10-transformed)', 'Breadth', 'SNPs/bp', 'Nucleotide diversity', 'dN/dS', 'Relative abundance'),
    variable = paste0(table(alllVarDF$variable) %>% names())
  )
  
  library(randomcoloR)
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


###### plot association between these variables #####
refDF$coverage  %<>% as.numeric()
refDF$breadth %<>% as.numeric()
refDF$relative_abundance = metaphlanDF$relative_abundance[match(paste0(refDF$species, '|', refDF$pt_ID.u), 
                                                                paste0(metaphlanDF$species, '|', metaphlanDF$sample_id))]

methd = 'pearson'

c_cutoff = 0
b_cutoff = 0

c_cutoff = 50
b_cutoff = 0.5

c_cutoff = 10
b_cutoff = 0.5

cand_var = c('relative_abundance', 'coverage_log10', 'coverage', 'breadth', 'SNPs_per_bp', 'nucl_diversity', 'dn_ds_common')
pdf(paste0('./plot/association_between_genetic_variables_', methd, '_', c_cutoff, '_', b_cutoff, '_all.pdf'), width = 13, height = 12)

dat = refDF[refDF$coverage >= c_cutoff & refDF$breadth >= b_cutoff, cand_var]
ppp = ggpairs(dat,title = paste0(capitalize(methd), ' | ', 'Coverage >= ', c_cutoff, ' | ', 
                                 'Breadth >= ', b_cutoff),
              diag =  list(continuous = wrap("barDiag", fill="#6BAED6")),
              lower = list(continuous = my_fn1),
              upper = list(continuous = eval(parse(text = paste0('my_fn2_', methd))))) + 
  theme(panel.grid = element_blank(), strip.text = element_text(colour = 'black', size = 10), 
        axis.text = element_text(size = 10, colour = 'black'))
print(ppp)
dev.off()

#### for each species
refDF$breadth %<>% as.numeric()
methd = 'pearson'

c_cutoff = 0
b_cutoff = 0

c_cutoff = 50
b_cutoff = 0.5

c_cutoff = 10
b_cutoff = 0.5

cand_var = c('relative_abundance', 'coverage_log10', 'coverage', 'breadth', 'SNPs_per_bp', 'nucl_diversity', 'dn_ds_common')
pdf(paste0('./plot/association_between_genetic_variables_', methd, '_', c_cutoff, '_', b_cutoff, '_sp_by_sp.pdf'), width = 13, height = 12)
for(sp in unique(refDF$species)){
  dat = refDF[refDF$coverage >= c_cutoff & refDF$breadth >= b_cutoff & refDF$species == sp, cand_var]
  try({
    ppp = ggpairs(dat,title = paste0(sp, ' | ' , capitalize(methd), ' | ', 'Coverage >= ', c_cutoff, ' | ', 
                                     'Breadth >= ', b_cutoff),
                  diag =  list(continuous = wrap("barDiag", fill="#6BAED6")),
                  lower = list(continuous = my_fn1),
                  upper = list(continuous = eval(parse(text = paste0('my_fn2_', methd))))) + 
      theme(panel.grid = element_blank(), strip.text = element_text(colour = 'black', size = 10), 
            axis.text = element_text(size = 10, colour = 'black'))
    print(ppp)
  })
}
dev.off()


#### adjust pvalue? ####
CN <- combn(cand_var, 2, simplify = FALSE)
bf_adj = c()
for(z in 1:length(CN)){
  x = dat[, CN[[z]][1]]
  y = dat[, CN[[z]][2]]
  af_adj[] = c(bf_adj, cor.test(x,y, method=methd)$p.value)
}
p.adjust(bf_adj, method = 'BH')



#####================== calculate community nucleotide diversity =====================#####
source('~/workspace/library_function.R')
load('./data/refDF.RData')
load('./data/metaData.RData')
load('./data/swabData.RData')
load('./data/metaphlanDF.RData')
load('./data/col_vector.RData')
refDF$coverage_log10 = log10(refDF$coverage %>% as.numeric())
refDF$relative_abundance = metaphlanDF$relative_abundance[match(paste0(refDF$species, '|', refDF$pt_ID.u), 
                                                                paste0(metaphlanDF$species, '|', metaphlanDF$sample_id))]
sV = unique(refDF$pt_ID.u[!grepl('KIT', refDF$pt_ID.u)])
cndV = c()
cndV2 = c()
dndsV = c()
dndsV2 = c()
for(si in sV){
  siDF = refDF[refDF$pt_ID.u == si, ] %>% na.omit()
  cnd = sum(siDF$relative_abundance * siDF$nucl_diversity)/sum(siDF$relative_abundance)
  cndV = c(cndV, cnd)
  cnd2= sum( siDF$nucl_diversity)/length(siDF$nucl_diversity)
  cndV2 = c(cndV2, cnd2)
  
  dnds = sum(siDF$relative_abundance * siDF$dn_ds_common)/sum(siDF$relative_abundance)
  dnds2 = sum( siDF$dn_ds_common)/length(siDF$relative_abundance)
  dndsV = c(dndsV, dnds)
  dndsV2 = c(dndsV2, dnds2)
}
cndDF = data.frame(
  sample_id = sV,
  community_nucl_div = cndV,
  community_nucl_div2 = cndV2,
  community_dNdS = dndsV,
  community_dNdS2 = dndsV2
)
cndDF$person_id = str_split_fixed(cndDF$sample_id, '_', 2)[, 1]
cndDF = merge(cndDF, metaData, by.x = 'person_id', by.y = 'pt_ID', all.x = T)
cndDF$trimester = swabData$trimester[match(cndDF$sample_id, swabData$pt_ID.u)]
colnames(cndDF)
ggboxplot(cndDF, x = 'Ethnicity', y = 'community_nucl_div') + 
  geom_boxplot(aes(color = Ethnicity))

if(plotFlag){
  pdf('./plot/community_nucl_div_vs_metadata.pdf', width = 10, height = 7)
}
xs = c('trimester', 'Age_char', 'BMI_char', 'Ethnicity', 'Language', 'Marriage', 'FOB', 'Employed', 'Housing', 
       'Term_char', 'Gravida_char', 'Abortions_char', 'Living_children_char', 
       'Abnormal_pap', 'Depression', 'PTSD', 'PPD_pos', 'DM2', 'Asthma', 'PCN_allergy', 'Sulfa_allergy', 'NKDA',
       'PNV', 'Fe', 'Insulin', 'B6', 'Progesterone',
       'IPV_hx', 'TED_hx', 'Tobacco_hx', 'EtOH_hx', 
       'GBS', 'PCN', 'Antibiotic',
       'Induction', 'Baby_gender', 'Baby_weight_char', 
       # repeated numeric 
       'Age', 'BMI', 'Term', 'Gravida', 'Abortions', 'Living_children', 'Baby_weight')
ys = c('community_nucl_div', 'community_nucl_div2', 'community_dNdS', 'community_dNdS2')
for (x in xs) {
  for (y in ys){
  tmpDF2 = data.frame(
    'pt_ID' = cndDF$person_id,
    'pt_ID.u' = cndDF$sample_id,
    y = cndDF[, y],
    x = cndDF[, x]
  )
  if(grepl('dNdS2', y)){
    ylabb = 'dN/dS of community (no weight)'
  }
  if(grepl('dNdS$', y)){
    ylabb = 'dN/dS of community'
  }
  if(grepl('nucl_div2', y)){
    ylabb = 'Nucleotide diversity of community (no weight)'
  }
  if(grepl('nucl_div$', y)){
    ylabb = 'Nucleotide diversity of community'
  }
  #pairs I want to compare in list format for stat_compare_means
  if(plotFlag){
    if(class(tmpDF2$x) == 'factor'){
      #all pairs I want to compare
      if(is.factor(tmpDF2$x)){
        CN <- combn(na.omit(droplevels(unique(tmpDF2$x))) %>% levels, 2, simplify = FALSE)
      }else{
        CN <- combn(na.omit(unique(tmpDF2$x)), 2, simplify = FALSE)
      }
      p = ggboxplot(na.omit(tmpDF2), x = 'x', y = 'y', color = 'x', add = 'jitter',
                    add.paramsshape = 'pt_ID') +
        mytheme3 + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        stat_compare_means(comparisons=CN, method = "wilcox.test", p.adjust.method = "BH", aes(label=..p.adj..)) +
        labs(title = paste0(x), 
             color = x,
             y = ylabb, x = x) + 
        scale_color_brewer(palette = 'Dark2')
      print(p)
    }else{
      p2 = ggplot(na.omit(tmpDF2), aes(x = x, y = y)) +
        geom_smooth(method = 'lm', color = '#969696', size = 0.5, fill = "#D5D5D5") +
        geom_point(aes(color = 'Blue', fill = 'Blue'), shape = 21, size = 2, alpha = 0.5) +
        labs(x = x, y = ylabb, title = paste0(x)) +
        stat_poly_eq(formula = y ~ x,
                     aes(group=1, label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), 
                     parse = TRUE) +
        mytheme + theme(strip.text = element_text(size = 12), legend.position = 'none')
      print(p2)
    }
    }
  }
}
dev.off()








