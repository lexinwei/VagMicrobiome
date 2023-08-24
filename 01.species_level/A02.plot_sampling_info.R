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


##### line point plot ####
if(plotFlag){pdf('./plot/trimester_vs_participant_line_point_4.pdf', width = 6.67, height = 7.69)}
if(plotFlag){
  pdfd = swabData[!grepl('kit', swabData$pt_ID, ignore.case = T) & swabData$trimester != 'P',]
  pdfd$pt_ID %<>% factor(., levels = metaData$pt_ID[order(metaData$Term, decreasing = T)])
  pdfd$SampleID =  metaData$SampleID[match(pdfd$pt_ID, metaData$pt_ID)]
  pdfd$SampleID %<>% factor(., levels = metaData$SampleID[order(metaData$Term, decreasing = T)])
  pdfd_P = swabData[!grepl('kit', swabData$pt_ID, ignore.case = T) & swabData$trimester == 'P',]
  pdfd_P$Sample_GA = 45
  p = ggplot() +
    geom_point(data = pdfd, mapping = aes(x = Sample_GA, y = SampleID, color = trimester), size = 2.3) +
    geom_line(data = pdfd, mapping = aes(x = Sample_GA, y = SampleID, color = trimester)) +
    geom_point(data = pdfd_P, mapping = aes(x = Sample_GA, y = SampleID, color = colorTrimester[4]), size = 2.3) +
    scale_color_manual(values = colorTrimester[1:3], name = 'Trimester') +
    new_scale("color") +
    geom_point(data = metaData, mapping = aes(x = Term, y = SampleID, color = Term_char), shape = ')', size = 3.5) +
    scale_color_manual(values = colorTerm, name = 'Term') +
    geom_vline(xintercept = c(14, 26), linetype = 'dashed', size = 0.3) +
    geom_vline(xintercept = 39, linetype = 'dashed', color = 'red', size = 0.3) +
    labs(x = 'Gestational Weeks', y = 'Participants', 
         # title = '    T1                   T2                                  T3', color = 'Trimester'
         ) +
    scale_x_continuous(breaks = c(seq(0, 40, 10), 45),
                       expand = expansion(mult = c(0, 0.035)),
                       limits = c(10, 45),
                       labels = c(seq(0, 40, 10), 'P')) +
    # facet_wrap(. ~ trimester, ncol = 4, scales = "free_x") + 
    mytheme2 + theme(axis.title = element_text(size = 16), legend.position = 'bottom', 
                     legend.box="vertical", legend.text = element_text(size = 16),
                     legend.title = element_text(size = 16),
                     legend.direction = 'horizontal', legend.margin = margin(-9, 0, 5, 0), 
                     plot.margin = unit(c(0, 0, 0, 0.1), 'cm'))
  pdfdAll = rbind(pdfd, pdfd_P)
  pTop <- ggplot(pdfdAll, aes(x = Sample_GA)) +
    geom_density(color = 'black', fill = '#FCF4D9', ) +
    scale_x_continuous(breaks = c(seq(0, 40, 10), 45),
                       expand = expansion(mult = c(0, 0.035)),
                       limits = c(10, 45),
                       labels = c(seq(0, 40, 10), 'P')) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), name = 'Density') + 
    mytheme3 +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
          axis.ticks.x = element_blank(), axis.text.y = element_text(size = 13, color = 'black'),
          axis.title.y = element_text(size = 16, color = 'black'),
          plot.margin  = unit(c(0.5, 0, 0.1, 0.09), 'cm'))
  pRight <- ggplot(pdfdAll, aes(x = pt_ID)) +
    geom_bar(stat = 'count', width = 0.75, color = 'black', fill = '#FCF4D9') + coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), name = 'Counts') + 
    mytheme3 +
    theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
          axis.ticks.y = element_blank(), axis.text.x = element_text(size = 13, color = 'black'),
          axis.title.x = element_text(size = 16, color = 'black'),
          plot.margin = unit(c(0, 0.5, 1.7, 0), 'cm'))
  pEmpty <- ggplot(pdfdAll, aes(x = pt_ID)) +
    geom_blank() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          line = element_blank(),
          panel.background = element_blank())
  grid.arrange(pTop, pEmpty, p, pRight,
               ncol = 2, nrow = 2, widths = c(5, 1), heights = c(0.9, 5))
}
if(plotFlag){dev.off()}

if(plotFlag){pdf('./plot/trimester_vs_participant_line_point_with_P.pdf', width = 15, height = 9.5)}
if(plotFlag){
  pdfd = swabData[!grepl('kit', swabData$pt_ID, ignore.case = T),]
  pdfd$pt_ID %<>% factor(., levels = metaData$pt_ID[order(metaData$Term, decreasing = T)])
  p2 = ggplot() +
    geom_point(data = pdfd, mapping = aes(x = Sample_GA, y = pt_ID, color = trimester)) +
    geom_line(data = pdfd, mapping = aes(x = Sample_GA, y = pt_ID, color = trimester)) +
    scale_color_manual(values = colorTrimester, name = 'Trimester') +
    new_scale("color") +
    geom_point(data = metaData, mapping = aes(x = Term, y = pt_ID, color = Term_char2), shape = ')', size = 4) +
    scale_color_manual(values = colorTerm2, name = 'Term') +
    geom_vline(xintercept = c(14, 26), linetype = 'dashed', size = 0.3) +
    geom_vline(xintercept = 39, linetype = 'dashed', color = 'red', size = 0.3) +
    labs(x = 'Pestational Weeks', y = 'Participants', 
         # title = '        T1                            T2                                      T3                                 P', 
         color = 'Trimester') +
    scale_x_continuous(breaks = seq(0, 200, 10), labels = seq(0, 200, 10)) +
    # facet_wrap(. ~ trimester, ncol = 4, scales = "free_x") +
    mytheme2 + theme(plot.title = element_text(size = 12))
  print(p2)
}
if(plotFlag){dev.off()}


###### number of sequenced samples for each participants #####
if(T){
  dat = data.frame(table(swabData$pt_ID[swabData$pt_ID %in% uPID]))
  mean(dat$Freq)
  sd(dat$Freq)
  # plot
  p <- ggplot(dat, aes(x=Freq)) +
    geom_histogram(binwidth = 1, fill="#FDF4D6", color="black", alpha=0.9) +
    scale_x_continuous(breaks = seq(0, max(dat$Freq), 3)) +
    scale_y_continuous() +
    labs(x = 'Number of samples', y = 'Number of participants') + mytheme2
  cairo_pdf('./plot/samples_per_participants.pdf', width = 3.07, height = 2.84)
  print(p)
  dev.off()
}
######  plot positions 96-well plate #####
if(F){
  wellInfo = data.frame(plate = paste('Plate', swabData$library_plate),
             well = swabData$library_well,
             pt_ID = swabData$pt_ID,
             pt_ID.u = swabData$pt_ID.u)
  
  wellInfo$y = swabData$library_well %>% str_extract(., '[A-Z]')
  wellInfo$y %<>% factor(., levels = c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H') %>% rev())
  table(wellInfo$y)
  wellInfo$x = swabData$library_well %>% str_extract(., '[0-9]+') %>% as.numeric()
  table(wellInfo$x)
  wellInfo$kit = 'No'
  wellInfo$kit[grepl('kit', wellInfo$pt_ID, ignore.case = T)] = 'Yes'
  wellInfo$nearkit = 'No'
  
  wellInfo$pt_ID.u %<>% str_replace(., '_plate.+', '')
  pdf('./plot/kit_96_well_plate.pdf', width = 20, height = 10)
  kip = ggplot(wellInfo, aes(x, y, label = pt_ID.u)) + 
    geom_point(shape = 21, size = 6, color = 'black', stroke = 0.5, aes(fill = kit)) +
    geom_text(size = 3, vjust = 2.5) + 
    scale_x_continuous(breaks = seq(1, 12, 1)) +
    scale_fill_manual(values = alpha(c('lightblue', '#FDB462'), 0.6)) + 
    labs(x = '', y = '') +
    facet_wrap(.~plate, nrow = 2) +
    mytheme + theme(panel.grid.minor = element_blank(),
                    strip.text = element_text(size = 12),
                    panel.grid.major = element_line(size = 0.5))
  print(kip)
  dev.off()
  display.brewer.all()
}

####### stat for Table S1 #####
load('data/RData/metaData.RData')
load('data/RData/swabData.RData')
swabData2 = swabData[!grepl('KIT', swabData$pt_ID.u), ]
load('./data/RData/sampleVagitypes.RData')
table(metaData$Age)
metaData$Age_char = NA
metaData$Age_char[metaData$Age >= 21 & metaData$Age <= 25 ] = '21-25'
metaData$Age_char[metaData$Age >= 26 & metaData$Age <= 30 ] = '26-30'
metaData$Age_char[metaData$Age >= 31 & metaData$Age <= 35 ] = '31-35'
metaData$Age_char[metaData$Age >= 36 & metaData$Age <= 40 ] = '36-40'
table(metaData$Age_char)
table(metaData$Age_char)/sum(table(metaData$Age_char))
swabData2$Age = metaData$Age_char[match(swabData2$pt_ID, metaData$pt_ID)]
table(swabData2$Age)
table(swabData2$Age)/sum(table(swabData2$Age))

table(metaData$BMI)
table(metaData$BMI_char)
swabData2$BMI = metaData$BMI_char[match(swabData2$pt_ID, metaData$pt_ID)]
table(swabData2$BMI)
table(swabData2$BMI)/sum(table(swabData2$BMI))

table(metaData$Ethnicity)
table(metaData$Ethnicity2)
swabData2$Ethnicity = metaData$Ethnicity[match(swabData2$pt_ID, metaData$pt_ID)]
table(swabData2$Ethnicity)

table(metaData$Employed)
swabData2$Employment = metaData$Employed[match(swabData2$pt_ID, metaData$pt_ID)]
table(swabData2$Employment)

table(metaData$Housing)
swabData2$Housing = metaData$Housing[match(swabData2$pt_ID, metaData$pt_ID)]
table(swabData2$Housing)

table(metaData$Term_char)
swabData2$Term_char = metaData$Term_char[match(swabData2$pt_ID, metaData$pt_ID)]
table(swabData2$Term_char)

table(metaData$Marriage)
swabData2$Marriage = metaData$Marriage[match(swabData2$pt_ID, metaData$pt_ID)]
table(swabData2$Marriage)

table(metaData$FOB)
swabData2$FOB = metaData$FOB[match(swabData2$pt_ID, metaData$pt_ID)]
table(swabData2$FOB)

table(metaData$Depression)
swabData2$Depression = metaData$Depression[match(swabData2$pt_ID, metaData$pt_ID)]
table(swabData2$Depression)

table(swabData$trimester)
