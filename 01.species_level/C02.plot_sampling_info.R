source('./script/library_function.R')
load('data/RData/metaData.RData')
load('data/RData/swabData.RData')
load('data/RData/colorList.RData')
plotFlag = T
#### show trimester_vs_participant ####
# heatmap
mycolors = brewer.pal(9, 'Greens')
if(plotFlag){pdf('./plot/trimester_vs_participant_heatmap_new.pdf', width = 6.5, height = 9.5)}
if(plotFlag){
  ha = rowAnnotation('Ethnicity' = metaData$Ethnicity[match(row.names(swabDataStage), metaData$pt_ID)],
                     'Term' = metaData$Term_char[match(row.names(swabDataStage), metaData$pt_ID)],
                     col = list('Ethnicity' = colorEthnicity,
                                'Term' =colorTerm))
  ph = Heatmap(swabDataStage, name = paste0("No. of sample"), cluster_rows = F, cluster_columns = F,
               column_names_side = "top", row_names_side = "left", column_title_side = 'top', 
               row_names_centered = T, column_names_centered = T,
               rect_gp = gpar(col = "lightgray", lwd = 1), na_col = 'lightgray',
               column_names_rot = T, col = mycolors, 
               width = ncol(swabDataStage)*unit(10, "mm"), 
               height = nrow(swabDataStage)*unit(5, "mm"),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if(!is.na(swabDataStage[i, j])){
                   grid.text(sprintf("%1.f", swabDataStage[i, j]), x, y, gp = gpar(fontsize = 10, col = 'black', alpha = 0.7))
                 }
               },
               row_title = "Participants", column_title = "Trimester", left_annotation = ha
  )
  draw(ph, heatmap_legend_side = "right", annotation_legend_side = "bottom", merge_legend = TRUE)
}
if(plotFlag){dev.off()}


####### stat for Table S1 #####
load('data/RData/MTMG_metaData.RData')
ptData = unique(metaData[, c('SID', 'Age_cat', 'Ethnictiy')])
tmp = metaData[, c('SID', 'UID')] %>% group_by(SID) %>% 
  summarise(sampleNum = n())
ptData$sampleNum = tmp$sampleNum[match(ptData$SID, tmp$SID)]
range(ptData$sampleNum)

table(ptData$Age_cat)
table(metaData$Age_cat)
table(metaData$Age_cat)/sum(table(metaData$Age_cat))
nrow(metaData)

table(ptData$Ethnictiy)
table(metaData$Ethnictiy)
table(metaData$Ethnictiy)/sum(table(metaData$Ethnictiy))
nrow(metaData)
