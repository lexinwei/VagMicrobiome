source('./script/library_function.R')
load('./data/RData/swabData.RData')
dat = read.delim2('./data/contig_stat_new.txt', header = T, sep = '\t')
dat$Sample %<>% str_extract(., '\\d+_\\d+')
dat$Sample.ID = swabData$SampleID[match(dat$Sample, swabData$seq_ID)]
dat1 = dat[which(dat$Sample.ID != 'NC' | is.na(dat$Sample.ID)), ]
range(dat1$Number.of.contigs)
mean(dat1$Number.of.contigs)
mean(dat1$N50)
sd(dat1$N50)

dat = read.delim2('./data/contig_stat_coReads_new.txt', header = T, sep = '\t')
dat1 = dat[which(dat$Sample != 'NC'), ]
range(dat1$Number.of.contigs)
mean(dat1$Number.of.contigs)
mean(dat1$N50)
sd(dat1$N50)
