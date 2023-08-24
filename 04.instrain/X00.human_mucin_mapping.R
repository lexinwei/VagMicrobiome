load('./data/RData/colorList.RData')
baseDir = '~/workspace/vag/human_mucin/vaginal_mucin/mucin'
fileList = list.files(baseDir)

resDF = read.delim(file = file.path(baseDir, fileList[1]))
colnames(resDF)
for(f in fileList[2:length(fileList)]){
  tabDF = read.delim(file = file.path(baseDir, f))
  resDF = merge(resDF, tabDF, all = T)
}
resDF = resDF[resDF$Genome != 'unmapped',]
LengthDF = resDF[, c(1, which(endsWith(colnames(resDF), '.Length'))[1])]
colnames(LengthDF) = c('Genome', 'Length')
CountDF = resDF[, c(1, which(endsWith(colnames(resDF), '.Read.Count')))]
colnames(CountDF) %<>% str_remove_all(., 'Ming_nova_VS_|_sorted.Read.Count')
Covered.BasesDF = resDF[, c(1, which(endsWith(colnames(resDF), '.Covered.Bases')))]
colnames(Covered.BasesDF) %<>% str_remove_all(., 'Ming_nova_VS_|_sorted.Covered.Bases')
Reads.per.baseDF = resDF[, c(1, which(endsWith(colnames(resDF), '.Reads.per.base')))]
colnames(Reads.per.baseDF) %<>% str_remove_all(., 'Ming_nova_VS_|_sorted.Reads.per.base')
Covered.FractionDF = resDF[, c(1, which(endsWith(colnames(resDF), '.Covered.Fraction')))]
colnames(Covered.FractionDF) %<>% str_remove_all(., 'Ming_nova_VS_|_sorted.Covered.Fraction')
cfMat = as.matrix(Covered.FractionDF[, -1])
row.names(cfMat) = Covered.FractionDF$Genome
cfMat = cfMat[, colnames(cfMat) %in% swabData$seq_ID]
colnames(cfMat) = swabData$pt_ID.u[match(colnames(cfMat), swabData$seq_ID)]
cfMatT = t(cfMat)

cfMat_kit = cfMat[, grep('^KIT', colnames(cfMat))]
annoTop = data.frame(
  Participant = swabData$pt_ID[match(colnames(cfMat), swabData$pt_ID.u)]
)
annoTop$Participant[grep('^KIT', annoTop$Participant)] = 'NC'
annoColor = list(
  Participant = c(colorParticipant, 'NC' = 'black')
)
row.names(annoTop) = colnames(cfMat)
pheatmap(cfMat, annotation_col = annoTop, 
         annotation_colors = annoColor,scale = 'none')
pdf('./plot/mucin_heatmap.pdf', width = 9.29, height = 22)
pheatmap(cfMatT, annotation_row = annoTop, fontsize_row = 5,
         annotation_colors = annoColor,scale = 'none')
dev.off()
