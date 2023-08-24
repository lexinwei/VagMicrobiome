source('./script/library_function.R')
# get dRep resulting table

dat1 = read.csv('./dRep_co_drep_results/Widb.csv')]
dat2 = read.csv('./dRep_co_drep_results/Chdb.csv')

dat1 = read.csv('./dRep_indiv_drep_results/Widb.csv')
dat2 = read.csv('./dRep_indiv_drep_results/Chdb.csv', sep = ',')


dat1$number_of_contigs = dat2$X..contigs[match(dat1$genome, dat2$Bin.Id)]
range(dat1$score)
piris <- ggplot(dat1, aes(completeness, contamination)) +
  geom_point(color = '#1F78B4', aes(size = score)) +
  labs(x = 'Completeness (%)', y = 'Contamination (%)', size = 'dRep score') +
  scale_size_continuous(range = c(1, 3)) + 
  mytheme + theme(legend.position = 'bottom')
ggMarginal(piris, fill = '#F1F1F1')

pdf('./plot/completeness_vs_contamination_of_contigs.pdf', width = 4.73, height = 4.98)
ggMarginal(piris, fill = '#F1F1F1')
dev.off()

# distribution of the number of contigs
p <- ggplot(dat1, aes(x=number_of_contigs)) +
  geom_histogram(binwidth = 100, fill="#1F78B4", color="#F1F1F1", alpha=0.9) +
  scale_y_continuous(breaks = seq(0, 8, 2)) +
  labs(x = 'Number of contigs', y = 'Number of genomes') + mytheme
pdf('./plot/contigs_per_bins.pdf', width = 3.9, height = 3.60)
print(p)
dev.off()

# distribution of the genome size
dat1$size_mbp = dat1$size/1000000
p <- ggplot(dat1, aes(x=size_mbp)) +
  geom_histogram(binwidth = 0.1, fill="#1F78B4", color="#F1F1F1", alpha=0.9) +
  scale_y_continuous(breaks = seq(0, 8, 2)) +
  labs(x = 'Genome size (Mbp)', y = 'Number of genomes') + mytheme
pdf('./plot/size_per_bins.pdf', width = 3.9, height = 3.60)
print(p)
dev.off()

# distribution of the N50
dat1$n50_kbp = dat1$N50/1000
p <- ggplot(dat1, aes(x=n50_kbp)) +
  geom_histogram(binwidth = 10, fill="#1F78B4", color="#F1F1F1", alpha=0.9) +
  scale_y_continuous(breaks = seq(0, 10, 2)) +
  labs(x = 'N50 (Kbp)', y = 'Number of genomes') + mytheme
pdf('./plot/N50_per_bins.pdf', width = 3.9, height = 3.60)
print(p)
dev.off()
#### include bins filtered out ###

dat2$Group[dat2$Bin.Id %in% dat1$genome] = 'remain'
dat2$Group[!dat2$Bin.Id %in% dat1$genome] = 'filter out'

piris <- ggplot(dat2, aes(Completeness, Contamination)) +
  geom_point(aes(color = Group), size = 2) +
  labs(x = 'Completeness (%)', y = 'Contamination (%)', size = 'dRep score') +
  scale_size_continuous(range = c(1, 3)) + 
  mytheme + theme(legend.position = 'bottom')
ggMarginal(piris, fill = '#F1F1F1')


# distribution of the number of contigs
p <- ggplot(dat2, aes(x=X..contigs)) +
  geom_histogram(binwidth = 100, fill="#1F78B4", color="#F1F1F1", alpha=0.9) +
  # scale_y_continuous(breaks = seq(0, 8, 2)) +
  labs(x = 'Number of contigs', y = 'Number of genomes') + mytheme
pdf('./plot/contigs_per_bins.pdf', width = 3.9, height = 3.60)
print(p)
dev.off()

# distribution of the genome size
dat2$size_mbp = dat2$Genome.size..bp./1000000
p <- ggplot(dat2, aes(x=size_mbp)) +
  geom_histogram(binwidth = 0.1, fill="#1F78B4", color="#F1F1F1", alpha=0.9) +
  labs(x = 'Genome size (Mbp)', y = 'Number of genomes') + mytheme
pdf('./plot/size_per_bins.pdf', width = 3.9, height = 3.60)
print(p)
dev.off()






