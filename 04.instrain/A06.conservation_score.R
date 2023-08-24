BiocManager::install('bio3d')
library(bio3d)

alnFa = read.fasta('/Users/xinwei/workspace/vag/colabfold/1594_lcl_Query_32930 and 99 other sequences.aln')
conValue = conserv(x=alnFa$ali, method="similarity", sub.matrix="bio3d")

alnFa = read.fasta('/Users/xinwei/workspace/vag/colabfold/48_lcl_Query_66140 and 99 other sequences.aln')
conValue = conserv(x=alnFa$ali, method="similarity", sub.matrix="bio3d")

alnFa = read.fasta('/Users/xinwei/workspace/vag/colabfold/1_lcl_Query_127727 and 99 other sequences.aln')
conValue = conserv(x=alnFa$ali, method="similarity", sub.matrix="bio3d")

alnFa = read.fasta('/Users/xinwei/workspace/vag/colabfold/5_lcl_Query_73331 and 99 other sequences.aln')
conValue = conserv(x=alnFa$ali, method="similarity", sub.matrix="bio3d")

res = data.frame(
  AA = alnFa[["ali"]][1, ],
  Score = conValue
)
res = res[res$AA != '-', ]
res$Pos = 1:nrow(res)
range(res$Score)
plot(res$Pos, res$Score)

resOut = data.frame(
  'fixedStep  chrom=chr1  start=1  step=1' = res$Score,
  check.names = F
)
write.table(res, './colabfold/5_conserv.txt', row.names = F, quote = F)


aln <- read.fasta('./instrain_res/input_fasta/s__Lactobacillus_jensenii.fasta_NZ_CP018809.1_5_nucl.fasta')
seq = str_c(aln$ali)
length(seq)
mucList = list()
i = 1
sink('./instrain_res/input_fasta/6muc.fasta')
for (s in c(1476, 1597, 1718, 1839, 1960, 2081)) {
  startP = s*3-2
  endP = (s + 79 - 1) * 3
  mucList[[i]] = str_c(seq[startP:endP], collapse = '')
  cat('>', i, '|', startP, ':', endP, '\n', sep = '')
  cat(mucList[[i]], '\n')
  i = i + 1
}
sink()

aln <- read.fasta('./instrain_res/input_fasta/s__Lactobacillus_jensenii.fasta_NZ_CP018809.1_5_prot.fasta')
seq = str_c(aln$ali)
length(seq)
mucList = list()
i = 1
sink('./instrain_res/input_fasta/6muc_prot.fasta')
for (s in c(1476, 1597, 1718, 1839, 1960, 2081)) {
  startP = s
  endP = s + 79 - 1
  mucList[[i]] = str_c(seq[startP:endP], collapse = '')
  cat('>', i, '|', startP, ':', endP, '\n', sep = '')
  cat(mucList[[i]], '\n')
  i = i + 1
}
sink()
