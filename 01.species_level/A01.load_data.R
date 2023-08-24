source('~/workspace/library_function.R')
###### sequencing sample name ######
fqID = read.delim2('./data/sample.list', header = F)[, 1] %>% str_remove_all(., '_bwa')
plotFlag = F

###### metadata ######
if(T){
  metaData = openxlsx::read.xlsx('./data/20200829_metadata.xlsx')
  nrow(metaData)
  setdiff(metaData$pt_ID, swabData$pt_ID)
  setdiff(swabData$pt_ID, metaData$pt_ID)
  metaData$BMI %<>% as.numeric()
  metaData$BMI_char[metaData$BMI < 18.5] = '<18.5' # 'low'
  metaData$BMI_char[metaData$BMI >= 18.5 & metaData$BMI < 25] = '18.5-24.9' #  'normal'
  metaData$BMI_char[metaData$BMI >= 25 & metaData$BMI < 30] =  '25.0-29.9' # 'overweight'
  metaData$BMI_char[metaData$BMI >= 30 & metaData$BMI < 35] =  '30.0-34.9' # 'obese'
  metaData$BMI_char[metaData$BMI >= 35] =  '>=35' # 'extremely obese'
  metaData$BMI_char %>% table()
  metaData$BMI_char %<>% factor(., levels = c('18.5-24.9', '25.0-29.9', '30.0-34.9', '>=35'))
  
  metaData$Age %>% range()
  metaData$Age_char[metaData$Age < 35] = '<35'
  metaData$Age_char[metaData$Age >= 35] = '>=35'
  metaData$Age_char %<>% factor(., levels = c('<35', '>=35'))
  
  metaData$Marriage = NA
  metaData$Marriage[grep('Married', metaData$Marital.Status)] = 'Married'
  metaData$Marriage[grep('Unmarried|Single', metaData$Marital.Status)] = 'Unmarried'
  metaData$Marriage[grep('unknown', metaData$Marital.Status, ignore.case = T)] = NA
  metaData$Marriage %>% table()
  metaData$Marriage %<>% factor(., levels = c('Married', 'Unmarried'))
  
  metaData$FOB = NA
  metaData$FOB[grep('FOB involved', metaData$Marital.Status)] = 'Involved'
  metaData$FOB[grep('not', metaData$Marital.Status)] = 'Uninvolved'
  metaData$FOB[metaData$Marital.Status == "Single"] = 'Uninvolved'
  metaData$FOB[metaData$Marital.Status == "Single, living w/ FOB"] = 'Involved'
  metaData$FOB %<>% factor(., levels = c('Involved', 'Uninvolved'))
  table(metaData$FOB, metaData$Marriage)
  
  metaData$Pregnancy.period
  metaData$Pregnancy.period[metaData$pt_ID == 'SF1872'] = 38.1
  sum(metaData$Pregnancy.period < 37)
  metaData$term %<>% factor(., levels = c('PTB', 'ETB', 'TB'))
  metaData$term_range[metaData$term == 'PTB'] = '<37'
  metaData$term_range[metaData$term == 'ETB'] = '>=37&<39'
  metaData$term_range[metaData$term == 'TB'] = '>=39'
  metaData$term_range %<>% factor(., levels = c('<37', '>=37&<39', '>=39'))
  table(metaData$term_range)
  metaData$Term_char = metaData$term
  metaData$Term_char2 = metaData$Term_char %>% as.character()
  metaData$Term_char2[metaData$Term_char2 %in% c('PTB', 'ETB')] = 'Preterm'
  metaData$Term_char2[metaData$Term_char2 == 'TB'] = 'Full-term'
  metaData$Term_char2 %<>% factor(., levels = c('Preterm', 'Full-term'))
  table(metaData$Term_char2)
  metaData$Term = metaData$Pregnancy.period
  metaData$Employed[metaData$Employed %in% c('Partime','Employed part time')] = 'Part-time'
  metaData$Employed %>% table()
  metaData$Employed %<>% factor(., levels = c("Employed", "Part-time", "Unemployed", 'Unknown'))
  metaData$Employment = metaData$Employed 
  metaData$Housing %<>% factor(., levels = c('Housed', 'Marginally housed', 'Homeless'))
  colnames(metaData)[colnames(metaData) == 'Ethinicity'] = 'Ethnicity'
  metaData$Ethnicity %>% table()
  metaData$Ethnicity[metaData$Ethnicity == 'Caucasian'] = 'White'
  metaData$Ethnicity[metaData$Ethnicity == 'Other '] = 'Other'
  metaData$Ethnicity2 = metaData$Ethnicity
  metaData$Ethnicity %<>% factor(., levels = c('White', 'Asian', 'Latina', 'Black', 'Pacific Islander', 'Other'))
  metaData$Ethnicity2[metaData$Ethnicity2 == 'Pacific Islander'] = 'Other'
  metaData$Ethnicity2 %<>% factor(., levels = c('White', 'Asian', 'Latina', 'Black', 'Other'))
  # Parity
  metaData$Gravida = metaData$Parity %>% str_extract(., 'G\\d+') %>% str_extract(., '\\d+') %>% as.numeric()
  metaData$Gravida %>% summary()
  Gravida_char_cf = 2
  metaData$Gravida_char[metaData$Gravida > Gravida_char_cf] = paste0('>', Gravida_char_cf)
  metaData$Gravida_char[metaData$Gravida <= Gravida_char_cf] = paste0('<=', Gravida_char_cf)
  metaData$Gravida_char %<>% factor(., levels = c(paste0('<=', Gravida_char_cf), paste0('>', Gravida_char_cf)))
  metaData$Gravida_char %>% table()
  
  metaData$Full_term = metaData$Parity %>% str_extract(., 'P\\d+') %>% str_sub(., start = 2, end = 2) %>% as.numeric()
  metaData$Preterm = metaData$Parity %>% str_extract(., 'P\\d+') %>% str_sub(., start = 3, end = 3) %>% as.numeric()
  metaData$Abortions = metaData$Parity %>% str_extract(., 'P\\d+') %>% str_sub(., start = 4, end = 4) %>% as.numeric()
  metaData$Abortions %>% summary()
  
  Abortions_char_cf = 2
  metaData$Abortions_char[metaData$Abortions > Abortions_char_cf] = paste0('>', Abortions_char_cf)
  metaData$Abortions_char[metaData$Abortions <= Abortions_char_cf] = paste0('<=', Abortions_char_cf)
  metaData$Abortions_char %<>% factor(., levels = c(paste0('<=', Abortions_char_cf), paste0('>', Abortions_char_cf)))
  metaData$Abortions_char %>% table()
  
  metaData$Living_children = metaData$Parity %>% str_extract(., 'P\\d+') %>% str_sub(., start = 5, end = 5) %>% as.numeric()
  metaData$Living_children %>% summary() 
  Living_children_char_cf = 1
  metaData$Living_children_char[metaData$Living_children > Living_children_char_cf] = paste0('>', Living_children_char_cf)
  metaData$Living_children_char[metaData$Living_children <= Living_children_char_cf] = paste0('<=', Living_children_char_cf)
  metaData$Living_children_char %<>% factor(., levels = c(paste0('<=', Living_children_char_cf), paste0('>', Living_children_char_cf)))
  metaData$Living_children_char %>% table()
  metaData$Full_term + metaData$Preterm
  metaData$Insurance %>% table()
  
  # PMH
  metaData$Abnormal_pap = 'No'
  metaData$Abnormal_pap[grep('abnormal pap', metaData$PMH, ignore.case = T)] = 'Yes'
  metaData$Abnormal_pap %<>% factor(., levels = c('No', 'Yes'))
  metaData$Abnormal_pap %>% table()
  metaData$Depression = 'No'
  metaData$Depression[grep('Depression', metaData$PMH, ignore.case = T)] = 'Yes'
  metaData$Depression %<>% factor(., levels = c('No', 'Yes'))
  metaData$Depression %>% table()
  metaData$PTSD = 'No'
  metaData$PTSD[grep('PTSD', metaData$PMH, ignore.case = T)] = 'Yes'
  metaData$PTSD %<>% factor(., levels = c('No', 'Yes'))
  metaData$PTSD %>% table()
  metaData$PPD_pos = 'No'
  metaData$PPD_pos[grep('PPD', metaData$PMH)] = 'Yes'
  metaData$PPD_pos %<>% factor(., levels = c('No', 'Yes'))
  metaData$PPD_pos %>% table()
  metaData$DM2 = 'No'
  metaData$DM2[grep('DM', metaData$PMH)] = 'Yes'
  metaData$DM2 %<>% factor(., levels = c('No', 'Yes'))
  metaData$DM2 %>% table()
  metaData$Asthma = 'No'
  metaData$Asthma[grep('Asthma', metaData$PMH)] = 'Yes'
  metaData$Asthma %<>% factor(., levels = c('No', 'Yes'))
  metaData$Asthma %>% table()
  metaData$PCN_allergy = 'No'
  metaData$PCN_allergy[grep('PCN', metaData$All)] = 'Yes'
  metaData$PCN_allergy %<>% factor(., levels = c('No', 'Yes'))
  metaData$PCN_allergy %>% table()
  metaData$Sulfa_allergy = 'No'
  metaData$Sulfa_allergy[grep('Sulfa', metaData$All, ignore.case = T)] = 'Yes'
  metaData$Sulfa_allergy %<>% factor(., levels = c('No', 'Yes'))
  metaData$Sulfa_allergy %>% table()
  metaData$NKDA = 'No' # No Known Drug Allergies
  metaData$NKDA[grep('NKDA', metaData$All)] = 'Yes'
  metaData$NKDA %<>% factor(., levels = c('No', 'Yes'))
  metaData$NKDA %>% table()
  
  metaData$PNV = 'No'
  metaData$PNV[grep('PNV', metaData$Med)] = 'Yes'
  metaData$PNV %<>% factor(., levels = c('No', 'Yes'))
  metaData$PNV %>% table()
  metaData$Fe = 'No'
  metaData$Fe[grep('Fe', metaData$Med)] = 'Yes'
  metaData$Fe %<>% factor(., levels = c('No', 'Yes'))
  metaData$Fe %>% table()
  metaData$Insulin = 'No'
  metaData$Insulin[grep('Insulin', metaData$Med, ignore.case = T)] = 'Yes'
  metaData$Insulin %<>% factor(., levels = c('No', 'Yes'))
  metaData$Insulin %>% table()
  metaData$B6 = 'No'
  metaData$B6[grep('pyridoxine|B6', metaData$Med, ignore.case = T)] = 'Yes'
  metaData$B6 %<>% factor(., levels = c('No', 'Yes'))
  metaData$B6 %>% table()
  metaData$Progesterone = 'No'
  metaData$Progesterone[grep('Progesterone', metaData$Med, ignore.case = T)] = 'Yes'
  metaData$Progesterone %<>% factor(., levels = c('No', 'Yes'))
  metaData$Progesterone %>% table()
  
  # Social hx
  metaData$IPV_hx = NA
  metaData$IPV_hx[grepl('hx IPV|prior IPV|possible IPV|IPV in prior relationships', metaData$Social.hx)
                  & !grepl('no IPV|no hx IPV', metaData$Social.hx)] = 'Yes'
  metaData$IPV_hx[grepl('no IPV|no hx IPV|denies IPV', metaData$Social.hx, ignore.case = T)] = 'No'
  metaData$IPV_hx %<>% factor(., levels = c('No', 'Yes'))
  metaData$IPV_hx %>% table()
  
  metaData$TED_hx = 'Yes' # At least ever/currently take one of them
  metaData$TED_hx[metaData$TED == 'No'] = 'No' # none of T/E/D
  metaData$TED_hx %<>% factor(., levels = c('No', 'Yes'))
  metaData$TED_hx %>% table()
  
  metaData$Tobacco_hx[metaData$Tobacco == 'Never'] = 'No'
  metaData$Tobacco_hx[grep('smoker', metaData$Tobacco)] = 'Yes'
  metaData$Tobacco_hx %<>% factor(., levels = c('No', 'Yes'))
  metaData$Tobacco_hx %>% table()
  
  metaData$EtOH_hx[metaData$EtOH == 'Never'] = 'No'
  metaData$EtOH_hx[grep('drinker', metaData$EtOH)] = 'Yes'
  metaData$EtOH_hx %<>% factor(., levels = c('No', 'Yes'))
  metaData$EtOH_hx %>% table()
  
  metaData$GBS %>% table()
  metaData$GBS[metaData$GBS == 'Unk'] = NA
  metaData$GBS %<>% factor(., levels = c('Neg', 'Pos'))
  
  metaData$PCN %>% table()
  
  metaData$Antibiotic = NA
  metaData$Antibiotic[metaData$PCN %in% c('N')] = 'No'
  metaData$Antibiotic[metaData$PCN %in% c('Y', 'y')] = 'Yes'
  metaData$Antibiotic[metaData$PCN %in% c('Cefazolin (d/t PNC all)', 'N (Clindamycin)')] = 'Yes'
  metaData$Antibiotic %>% table()
  metaData$Antibiotic %<>% factor(., levels = c('No', 'Yes'))
  
  metaData$PCN[metaData$PCN %in% c('N', 'N (Clindamycin)', 'Cefazolin (d/t PNC all)')] = 'No'
  metaData$PCN[metaData$PCN %in% c('Y', 'y')] = 'Yes'
  metaData$PCN[metaData$PCN %in% c('Unk')] = NA
  metaData$PCN %<>% factor(., levels = c('No', 'Yes'))
  metaData$PCN %>% table()
  
  metaData$Induction %>% table()
  metaData$Baby_gender = metaData$Sex
  metaData$Baby_gender %<>% factor(., levels = c('M', 'F', 'M, M', 'M, F', 'F, F'))
  metaData$Baby_gender %>% table()
  
  metaData$Baby_weight = metaData$Birth.wt
  metaData$Baby_weight[grep(';', metaData$Birth.wt)] = NA
  metaData$Baby_weight %<>% as.numeric()
  metaData$Baby_weight %>% summary()
  metaData$Baby_weight_char[metaData$Baby_weight <= 3500] = '<=3.5kg' # The average birth weight for babies is around 7.5 lb (3.5 kg)
  metaData$Baby_weight_char[metaData$Baby_weight > 3500] = '>3.5kg'
  metaData$Baby_weight_char %<>% factor(., levels = c('<=3.5kg', '>3.5kg'))
  metaData$Baby_weight_char %>% table()
  
  metaData$Twin = NA
  metaData$Twin[grepl(',', metaData$Baby_gender)] = 'Yes'
  metaData$Twin[metaData$Baby_gender %in% c('M', 'F')] = 'No'
  metaData$Twin %<>% factor(., levels = c('No', 'Yes'))
  table(metaData$Twin)
  metaData$Language %<>% factor(., c('English', 'Spanish'))
  colnames(metaData)
  useful_vars = c('pt_ID', 'Age', 'Age_char', 'BMI', 'BMI_char', 'Ethnicity', 'Ethnicity2', 'Language', # Self attributes
                  'Marriage', 'FOB', 'Employment', 'Housing', # Living status
                  'Term', 'Term_char','Term_char2', 'Gravida', 'Gravida_char','Full_term', 'Preterm', 'Abortions', 'Abortions_char', 'Living_children', 'Living_children_char', # Parity
                  'Abnormal_pap', 'Depression', 'PTSD', 'PPD_pos', 'DM2', 'Asthma',  # PMH
                  'PCN_allergy', 'Sulfa_allergy', 'NKDA', # Allergy, NKDA - no known drug allergies
                  'PNV', 'Fe', 'Insulin', 'B6', 'Progesterone', # Medicine
                  'IPV_hx', 'TED', 'TED_hx','Tobacco_hx', 'EtOH_hx', # Social hx
                  'GBS', 'PCN', 'Antibiotic', 'Induction', 'Sex', 'Baby_gender', 'Birth.wt', 'Baby_weight', 'Baby_weight_char', 'Twin')
  metaData = metaData[, useful_vars]
  library(numform)
  metaData$SampleID = paste0('P', f_pad_zero(order(metaData$pt_ID)))
}

save(metaData, file = './data/RData/metaData.RData')

###### swabdata #####
swabData = read.delim2('./data/20200829_swab_ID.txt') # 30 samples don't have swabID
swabData$fqID = paste0('Ming_nova_VS_', swabData$seq_ID)
nrow(swabData)
if(plotFlag){
  sink('./data/fqID_without_swabID.txt')
  cat(paste0(setdiff(fqID, swabData$fqID), '_bwa'), sep = '\n')
  sink()
  
  sink('./data/fqID_with_swabID.txt')
  cat(paste0(intersect(fqID, swabData$fqID), '_bwa'), sep = '\n')
  sink()
  
  setdiff(swabData$fqID, fqID)
  sum(grepl('kit', swabData$pt_ID, ignore.case = T))
  
  uPID = swabData$pt_ID %>% unique()
  uPID = uPID[!grepl('kit', uPID, ignore.case = T)]
  length(uPID)
}

# add pregnancy trimester
if(T){
  swabData$Sample_GA %<>% as.numeric()
  table(swabData$Trimester)
  swabData$trimester[swabData$Sample_GA >= 0 & swabData$Sample_GA <= 14] = 'T1'
  swabData$trimester[swabData$Sample_GA > 14 & swabData$Sample_GA <= 26] = 'T2'
  swabData$trimester[swabData$Sample_GA > 26] = 'T3'
  swabData$trimester[swabData$Trimester == 'Postpartum'] = 'P'
  swabData$trimester = factor(swabData$trimester, levels = c('T1', 'T2', 'T3', 'P'))
  table(swabData$trimester)
  swabDataStage = dcast(swabData[!grepl('kit', swabData$pt_ID, ignore.case = T),], pt_ID ~ trimester)
  row.names(swabDataStage) = swabDataStage[, 1]
  swabDataStage = swabDataStage[, -1] %>% as.matrix()
  table(swabDataStage)
  
  swabData$Term = metaData$Term_char[match(swabData$pt_ID, metaData$pt_ID)]
  swabData$Term2 = metaData$Term_char2[match(swabData$pt_ID, metaData$pt_ID)]
  swabData$Ethnicity = metaData$Ethnicity[match(swabData$pt_ID, metaData$pt_ID)]
  swabData$Ethnicity2 = metaData$Ethnicity2[match(swabData$pt_ID, metaData$pt_ID)]
}

# add week for Postpartum
if(T){
  swabData = swabData[mixedorder(swabData$pt_ID.u), ]
  swabData$date = as.Date(swabData$Date_Acquired, format = "%d-%b-%y")
  for(i in 1:nrow(swabData)){
    if(swabData$Trimester[i] == 'Postpartum'){
      swabData$Sample_GA[i] = difftime(swabData$date[i], swabData$date[i-1], units="weeks") %>% as.numeric() %>% round(., 1) + swabData$Sample_GA[i-1]
    }
  }
}
if(F){
  # output filename for each individual
  for(p in unique(swabData$pt_ID[!grepl('kit', swabData$pt_ID, ignore.case = T)])){
    sink(paste0('./data/fq_name_for_each_individual/', p, '.txt'))
    cat(swabData$fqID[swabData$pt_ID == p], sep = '\n')
    sink()
  }
  # NC samples
  sink(paste0('./data/fq_name_for_each_individual/', 'NC', '.txt'))
  cat(swabData$fqID[grepl('kit', swabData$pt_ID, ignore.case = T)], sep = '\n')
  sink()
}
# keep 3_r for SF1769, remove 11_r for SF1651
swabData = swabData[!(swabData$pt_ID.u %in% c('SF1769_3', 'SF1651_11_r')), ]
swabData$pt_ID.u[swabData$pt_ID.u == 'SF1769_3_r'] = 'SF1769_3'
# add new sampleID
swabData$SampleID = metaData$SampleID[match(swabData$pt_ID, metaData$pt_ID)]
swabData$SampleID[is.na(swabData$SampleID)] = 'NC'
swabData$SampleID.u = paste0(swabData$SampleID, '-', swabData$Study_visit_number)
swabData$SampleID.u %<>% str_remove('_r$')
save(swabData, file = './data/RData/swabData.RData')


##### get taxa to NCBI taxid table #####
if(F){
  clade2taxidDF = data.frame()
  for (x in swabData$fqID){
    df0 = read.delim2(paste0('./metaphlan4_results/profiles/', x,'_bwa_profile.tsv'), sep = '\t', skip = 4)
    df1 = data.frame(clade_name = df0$X.clade_name, 
                     NCBI_tax_id = df0$NCBI_tax_id)
    clade2taxidDF = rbind(clade2taxidDF, df1)
  }
  clade2taxidDF = unique(clade2taxidDF)
  clade2taxidDF1 = str_split_fixed(clade2taxidDF$clade_name, '\\|', Inf) %>% as.data.frame()
  clade2taxidDF2 = str_split_fixed(clade2taxidDF$NCBI_tax_id, '\\|', Inf) %>% as.data.frame()
  taxa2taxidDF = data.frame()
  for(i in 1:8){
    taxa2taxidDF = rbind(taxa2taxidDF, data.frame(taxa = clade2taxidDF1[, i], taxid = clade2taxidDF2[, i]))
  }
  taxa2taxidDF = unique(taxa2taxidDF)
  taxa2taxidDF = taxa2taxidDF[!(taxa2taxidDF$taxa == '' & taxa2taxidDF$taxid == ''), ]
  save(clade2taxidDF, file = './data/RData/clade2taxidDF.RData')
  save(taxa2taxidDF, file = './data/RData/taxa2taxidDF.RData')
}

##### metaphlan4 abundance table ####
s2S = data.frame(
  S=c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Strain')
)
s2S$s = str_sub(s2S$S, 1, 1) %>% tolower()
s2S$s[8] = 't'
save(s2S, file = './data/RData/s2S.RData')
if(T){
  abuDF0 = read.table('./metaphlan4_results/merged_abundance_table.txt', 
                     comment.char = '#', sep = '\t', header = TRUE, check.names=FALSE)
  colnames(abuDF0) %<>% str_remove(., '_bwa$') %>% str_remove(., '^Ming_nova_VS_') 
  abuDF_tmp = abuDF0[abuDF0$clade_name %in% c('k__Bacteria', 'k__Eukaryota', 'UNCLASSIFIED'), ]
  colSums(abuDF_tmp[, -1]) %>% range()
  # split to get profiles for each level
  abuList = list()
  for (i in 1:8) {
    l = s2S$s[i]
    if(i != 8){
      la = s2S$s[i+1]
      abuDF1 = abuDF0[grepl(paste0(l, '__'), abuDF0$clade_name) & !grepl(paste0(la, '__'), abuDF0$clade_name), ]
    }else{
      abuDF1 = abuDF0[grepl(paste0(l, '__'), abuDF0$clade_name), ]
    }
    if(i != 1){
      abuDF1 %<>% mutate(taxa = str_split_fixed(abuDF1$clade_name, '\\|', i)[, i],
                               .after = clade_name)
    }else{
      abuDF1 %<>% mutate(taxa = abuDF1$clade_name,
                         .after = clade_name)
    }
    abuDF1 %<>% mutate(taxid = taxa2taxidDF$taxid[match(abuDF1$taxa, taxa2taxidDF$taxa)],
                       .after = taxa)
    abuDF1$taxa %<>% str_remove(., '^[kpcofgst]__') 
    abuDF2 = rbind(abuDF1, c('unclassified', 'unclassified', -1, 100-colSums(abuDF1[, 4:ncol(abuDF1)])))
    abuMat2 = apply(abuDF2[, 4:ncol(abuDF2)], 2, as.numeric)
    range(abuMat2)
    abuMat2[abuMat2 < 0] = 0
    abuDF2 = cbind(abuDF2[, 1:3], abuMat2)
    write.csv(abuDF2, file = paste0('./data/abundance_', l, '.csv'), row.names = F, quote = F)
    abuList[[l]] = abuDF2
  }
  abuS = abuList[['s']]
  abuS$taxa[abuS$taxa == 'GGB1215_SGB1581'] = 'Prevotellaceae_SGB1581'
  abuS$taxa[abuS$taxa == 'GGB12785_SGB19825'] = 'Candidatus_Nanogingivalaceae_SGB19825'
  abuS$taxa[abuS$taxa == 'GGB12796_SGB19893'] = 'Candidatus_Saccharibacteria_SGB19893'
  abuS$taxa[abuS$taxa == 'GGB1455_SGB2018'] = 'Bacteroidetes_SGB2018'
  abuS$taxa[abuS$taxa == 'GGB1456_SGB2019'] = 'Bacteroidetes_SGB2019'
  abuS$taxa[abuS$taxa == 'GGB2722_SGB3663'] = 'Lawsonellaceae_SGB3663'
  abuS$taxa[abuS$taxa == 'GGB2945_SGB3918'] = 'Firmicutes_SGB3918'
  abuS$taxa[abuS$taxa == 'GGB3012_SGB4003'] = 'Hungateiclostridiaceae_SGB4003'
  abuS$taxa[abuS$taxa == 'GGB3293_SGB4348'] = 'Clostridia_SGB4348'
  abuS$taxa[abuS$taxa == 'GGB39918_SGB47522'] = 'Actinobacteria_SGB47522'
  abuS$taxa[abuS$taxa == 'GGB4239_SGB5731'] = 'Acidaminococcaceae_SGB5731'
  abuS$taxa[abuS$taxa == 'GGB4277_SGB5832'] = 'Veillonellaceae_SGB5832'
  abuS$taxa[abuS$taxa == 'GGB753_SGB989'] = 'Actinomycetaceae_SGB989'
  abuList[['s']] = abuS
  save(abuList, file = './data/RData/abuList.RData')
}
##### output species list for strainphlan #####
if(F){
  load('./data/RData/taxa2taxidDF.RData')
  abuT = abuList[['t']]
  nrow(abuT)
  sink('./data/reference_genomes_all.list')
  refGenome = data.frame(
    species_name = str_split_fixed(abuT$clade_name[abuT$taxa != 'unclassified'], '\\|', 8)[, 7],
    SGB_name = abuT$taxa[abuT$taxa != 'unclassified']
  )
  refGenome$taxid = taxa2taxidDF$taxid[match(refGenome$species_name, taxa2taxidDF$taxa)]
  cat(refGenome$SGB_name, sep = '\n')
  sink()
  write.csv(refGenome[refGenome$taxid != '',], file = './data/SGB_species_to_taxid.csv',
            row.names = F, quote = F)
  table(duplicated(refGenome$taxid))
  table(duplicated(refGenome$species_name))
  table(duplicated(abuS$taxid))
  a = abuS[abuS$taxid %in% abuS$taxid[which(duplicated(abuS$taxid))],]
}

# scp ~/workspace/vag/data/SGB_species_to_taxid.csv weixin@10.73.29.10:~/workspace/vag/NCBI_genomes

####### output species for instrain #####
if(T){ # prepare abundance matrix
  abuS = abuList[['s']]
  abuS_mat = as.matrix(abuS[abuS$taxa != 'unclassified', 4:ncol(abuS)])
  row.names(abuS_mat) = abuS$taxa[abuS$taxa != 'unclassified']
  range(abuS_mat)
  # abuS_mat = apply(abuS_mat, 2, function(x){x/sum(x)})
  abuS_mat = abuS_mat/100
  range(abuS_mat)
  range(rowSums(abuS_mat))
  colSums(abuS_mat)
  abuS_mat_T = t(abuS_mat)
  interSample = intersect(rownames(abuS_mat_T), swabData$seq_ID[!grepl('^KIT', swabData$sample_ID)])
  abuMat = abuS_mat[, match(interSample, colnames(abuS_mat))]
  range(abuMat)
  range(rowSums(abuMat))
  abuMat = abuMat[which(rowSums(abuMat) > 0), ]
  dim(abuMat)
  criteria1 = which(rowSums(abuMat >= 0.01) >= ncol(abuMat)*0.05)
  criteria2 = which(rowSums(abuMat >= 0.001) >= ncol(abuMat)*0.15)
  setdiff(criteria2, criteria1)
  setdiff(criteria1, criteria2)
  keySpecies = row.names(abuMat)[union(criteria1, criteria2)]
}

##### define color vector ####
if(T){ # Phylum
  abuP  = abuList[['p']]
  cat(paste0(paste0('"', abuP$taxa, '" ='), collapse = '\n'))
  colorPhylum = list(
    "Firmicutes" = '#01665D',
    "Actinobacteria" = '#DA7E50',
    "Ascomycota" = '#8BB759',
    "Bacteroidetes" = '#908BC1',
    "Proteobacteria" = '#E34EA2',
    "Basidiomycota" = '#B3915B',
    "Fusobacteria" = '#E6B94D', 
    "Tenericutes" = '#6B8ABE',
    "Synergistetes" = '#878787'
  )
}
if(T){ # Kingdoom
  abuK  = abuList[['k']]
  cat(paste0(paste0('"', abuK$taxa, '" ='), collapse = '\n'))
  colorKingdom = list(
    "Bacteria" = '#B1D3E5',
    "Fungi" = '#FABF82'
  )
}
colorVagitype = c(
  "Lactobacillus crispatus" = '#90ADDF', 
  "Gardnerella vaginalis" = '#FFF180',
  "Lactobacillus iners" = '#A68CC4',
  "Lactobacillus jensenii" = '#F59A4B',
  "Lactobacillus gasseri" = '#EE4E00',
  "Fannyhessea vaginae" = '#83D55A',
  "Lactobacillus coleohominis" = '#DDB780',
  "Other" = '#2D58A3',
  "None dominant" = '#CCCFD8'
)
colorSpecies = c(
  colorVagitype[1:7],
  "Lactobacillus vaginalis" = '#E791BE',
  'Megasphaera genomosp type 1' = '#673CA4',
  'Actinomycetaceae SGB989' = '#B20017',
  'Candida albicans' = '#31A4B9',
  'Hungateiclostridiaceae SGB4003' = '#747474',
  'Coriobacteriales bacterium DNF00809' = '#9CD7B3',
  'Aerococcus christensenii' = '#85CAD8',
  'Veillonellaceae bacterium DNF00626' = '#A3756E',
  'Other' =  '#2D58A3',
  'Unclassified' = '#CCCFD8'
)
colorVagitypeShort = c(
  "L. crispatus" = '#90ADDF', 
  "G. vaginalis" = '#FFF180',
  "L. iners" = '#A68CC4',
  "L. jensenii" = '#F59A4B',
  "L. gasseri" = '#EE4E00',
  "F. vaginae" = '#83D55A',
  "L. coleohominis" = '#DDB780',
  "Other" = '#2D58A3',
  "None dominant" = '#CCCFD8'
)
colorCST = c(
  "I" = '#F57400', 
  "II" = '#E90011',
  "III" = '#80BE15',
  "IV" = '#357A60',
  "V" = '#26338F'
)

colorTerm = c(
  # 'PTB' = '#FB8072',
  # 'ETB' = '#8DA0CB',
  # 'TB' = '#C7E9B4'
  'PTB' = '#363C44',
  'ETB' = '#477C99',
  'TB' = '#B2C8C0'
)
colorTerm2 = c(
  'Preterm' = '#DF8F44',
  'Full-term' = '#374E55'
)
colorTrimester = c(
  'T1' = '#E34447', 
  'T2' = '#32A108', 
  'T3' = '#4275FF', 
  'P' = '#A046FF'
)
colorEthnicity = c(
  'White' = '#B9D392',
  'Asian' = '#D8BB7A' ,
  'Latino' = '#7DBAA4',
  'Black' = '#9A8A8B',
  'Pacific Islander' = '#F48DB7' ,
  'Other' = '#CDCDCD'
)
colorEthnicity2 = c(
  'White' = '#B9D392',
  'Asian' = '#D8BB7A' ,
  'Latino' = '#7DBAA4',
  'Black' = '#9A8A8B',
  'Other' = '#CDCDCD'
)
colorEmployment = c(
  'Employed' = '#F8766D',
  'Part-time' = '#01BFC4',
  'Unemployed' = '#7CAE00',
  'Unknown' = 'gray'
)
colorHousing = c(
  'Housed' = '#F8766D',
  'Marginally housed' = '#01BFC4',
  'Homeless' = '#7CAE00'
)
colorAge = c(
  '<35' = '#F8766D',
  '>=35' = '#01BFC4'
)
colorBMI = c(
  '18.5-24.9' = '#B9D392',
  '25.0-29.9' = '#D8BB7A' ,
  '30.0-34.9' = '#7DBAA4',
  '>=35' = '#9A8A8B'
)
colorQuality = c(
  'High' = '#F8766D',
  'Medium' =  '#00BA38'
  # 'Low' = 
)
colorMarriage = c(
  'Married' = '#F8766D',
  'Unmarried' =  '#00BA38',
  'Unknown' = 'gray'
)
colorFOB = c(
  'Involved' = '#F62B1E',
  'Uninvolved' =  '#3A94EC'
)
colorDepression = c(
  'Yes' = '#FCBA00',
  'No' =  '#2CB45C'
)
colorStrain = c(
  '1' = '#8AB4D5',
  '2' = '#F9B56B',
  '3' = '#B9E16E',
  '4' = '#F5847B',
  '5' = '#C25DD4',
  '6' = '#877C5F',
  '0' = '#D4D4D4'
)
colorPhylumMAGs = c(
  "Actinobacteriota" = '#CC683D',
  'Bacteroidota' = '#7C77B5',
  "Campylobacterota" = '#DDAB36',
  "Firmicutes" = '#01665D',
  "Firmicutes_A" = '#35978F',
  "Firmicutes_B" = '#80CDC1',
  "Firmicutes_C" = '#C7EAE5',
  "Patescibacteria" = '#5977B2',
  "Proteobacteria" = '#D43393'
)
colorFamily = c(
  'Actinomycetaceae' = '#66C2A5',
  'Aerococcaceae' = '#1F78B4',
  'Bifidobacteriaceae' = '#8DA0CB',
  'Coriobacteriales_unclassified' = '#E78AC3',
  'Debaryomycetaceae' = '#D40015',
  'Lactobacillaceae' = '#F4A34E',
  'Veillonellaceae' = '#A6D854'
)
colorLactobacillaceae= c('Lactobacillaceae' = '#F4A34E', 
                         'Non-Lactobacillaceae' = '#70A2C9')
save(colorKingdom, colorPhylum, colorVagitype, colorCST, colorVagitypeShort,
     colorTerm, colorTerm2, colorPhylumMAGs,colorFamily,
     colorTrimester, colorEthnicity, colorEthnicity2,colorFOB,colorMarriage,colorDepression,
     colorEmployment, colorHousing, colorStrain,colorLactobacillaceae,
     colorSpecies, colorParticipant,colorQuality,
     colorAge, colorBMI, file = './data/Rdata/colorList.RData')

###### count reads, base, and contigs ####
load('./data/RData/swabData.RData')
countDF0 = read.delim('./data/fastq_count_results.txt', header = F)
countDF1 = data.frame(
  seqID = str_extract(countDF0$V1, '\\d_\\d\\d'),
  base = str_extract(countDF0$V2, '\\d+,') %>% str_remove(., ',') %>% as.numeric() / 2,
  reads = str_extract(countDF0$V2, '\\d+\\}$') %>% str_remove(., '\\}') %>% as.numeric() / 2
)
countDF = countDF1[!(countDF1$seqID %in% swabData$seq_ID[grepl('KIT', swabData$pt_ID)]), ]
sum(countDF$base / 1000000000) 
range(countDF$base / 1000000000)
mean(countDF$base / 1000000000)
sd(countDF$base / 1000000000)

sum(countDF$reads / 1000000) 
range(countDF$reads / 1000000)
mean(countDF$reads / 1000000)
sd(countDF$reads / 1000000000)
