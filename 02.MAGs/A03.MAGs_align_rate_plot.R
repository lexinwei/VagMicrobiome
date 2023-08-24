source('./script/library_function.R')
load('./data/RData/swabData.RData')
# ##### Alignment rates of individually mapping, results from bowtie2 #####
# arDF = read.delim2('./MAGs/vag_mag138_align_rate.txt', sep = ' ', header = F)
# arDF$jobID = arDF$V1 %>% str_extract(., '\\d+.err') %>% str_remove_all(., '\\.err')
# arDF$alignment_rate = arDF$V1 %>% str_remove_all(., '/public/home/weixin/workspace/vag/mapping/stdout/vag_E02_\\d+_\\d+.err:\\d+:')
# arDF$alignment_rate %<>% str_remove_all(., '%') %>% as.numeric()
# swabData$jobID = order(swabData$seq_ID)
# swabData$alignment_rate = arDF$alignment_rate[match(swabData$jobID, arDF$jobID)]

##### abudance covered by MAGs ####
load('./data/RData/magAbu.RData')
map.vs.unmap = data.frame(
  SeqID = colnames(magAbu[, -1]),
  unmapped = magAbu[row.names(magAbu) == 'unmapped', -1] %>% as.numeric()
)
# map.vs.unmap = map.vs.unmap[1:10, ]
map.vs.unmap$mapped = 100 - map.vs.unmap$unmapped
map.vs.unmap$group = swabData$pt_ID[match(map.vs.unmap$SeqID, swabData$seq_ID)]
map.vs.unmap$group %<>% factor(., levels = names(colorParticipant))
map.vs.unmap = map.vs.unmap[order(map.vs.unmap$group, map.vs.unmap$mapped), ]
map.vs.unmap$sampleOrder = 1:nrow(map.vs.unmap)
# map.vs.unmap = melt(map.vs.unmap)
library(dplyr)
map.vs.unmap2 = map.vs.unmap %>%
  group_by(group) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  mutate(group = lag(group),
         sampleOrder = sampleOrder + 0.000001,
         mapped = mapped - 0.000001) %>%
  filter(row_number() %in% 2:n()) %>%
  bind_rows(map.vs.unmap)
mean(map.vs.unmap2$mapped)
sd(map.vs.unmap2$mapped)
range(map.vs.unmap2$mapped)
p = ggplot(map.vs.unmap2, aes(x = sampleOrder, y=mapped, group=group, 
                              fill = group)) +
  geom_area(aes(x = sampleOrder, y = 100), fill = '#E2E5E9', alpha = 0.5, ) +
  geom_area(alpha=1) +
  geom_line(size = 0.2, aes(color = group)) +
  # geom_col(width = 1, color = NA) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100), name = 'Percentage of mapped reads (%)') +
  scale_x_continuous(expand = c(0, 0),
                     breaks = c(0, map.vs.unmap2$sampleOrder[1:34]) +
                       (c(map.vs.unmap2$sampleOrder[1:34], 319) - c(0, map.vs.unmap2$sampleOrder[1:34]))/2,
                     labels = levels(map.vs.unmap2$group),
                     name = 'Samples of each participant') +
  scale_color_manual(values = colorParticipant,
                     aesthetics = c("color", "fill"), name = "Participant") +
  mytheme3 + theme(legend.position = 'none',aspect.ratio = 1,
                   # axis.ticks.x = element_blank(), 
                   axis.text.x =  element_blank())
cairo_pdf('./plot/Percentage of mapped reads2.pdf', width = 4, height = 4)
print(p)
dev.off()
