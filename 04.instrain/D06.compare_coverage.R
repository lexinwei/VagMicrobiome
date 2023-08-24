source('./script/library_function.R')
load('./data/RData/colorList.RData')
load('./data/RData/refDF_minCovBreath10.RData')
refDF1 = refDF
refDF1$species %<>% str_replace_all(., 's__|\\.fasta|_', ' ') %>% str_trim()
refDF1$coverage_log = log10(as.numeric(refDF1$coverage))
refDF1$breadth %<>% as.numeric()
countSP = refDF1 %>% group_by(species) %>%
  summarise(sampleNum = n())
seletedSP = countSP$species[countSP$sampleNum > 20]
refDF1$D_rev = 1 - as.numeric(refDF1$d_prime_mean)
load('./data/RData/MTMG_refDF_minCovBreath10.RData')
refDF2 = refDF
refDF2$coverage_log = log10(as.numeric(refDF2$coverage))
refDF2$D_rev = 1 - as.numeric(refDF2$d_prime_mean)
refDF2 = refDF2[refDF2$D_rev < 0.5 & refDF2$dnds < 1.5, ]
refDF2$breadth %<>% as.numeric()
overlapSP = unique(intersect(refDF1$species, refDF2$species))
length(overlapSP)
overlapSP_sel = overlapSP[overlapSP %in% seletedSP]

if(T){ # compare coverage
  yvals = c('coverage_log', "breadth", "breadth_minCov", 'nucl_diversity', 'D_rev', "pnps", "dnds")
  yLabs = list(coverage_log = 'log10(Coverage)', breadth = 'Breath', breadth_minCov = 'Breath (coverage > 5)', 
               nucl_diveristy = expression(pi), D_rev = "1 - D'", pnps = 'pN/pS', dnds = 'dN/dS')
  cairo_pdf('./plot/compare_of_two_cohorts.pdf', width = 8.79, height = 5.5, onefile = T)
  for(yval in yvals){
    comDF = rbind(
      data.frame(refDF1[refDF1$species %in% overlapSP_sel, c('species', yval)],
                 Cohort = 'Stanford'),
      data.frame(refDF2[refDF2$species %in% overlapSP_sel, c('species', yval)],
                 Cohort = 'PRJNA797778')
    )
    table(comDF$species)
    # comDF$Cohort = factor(comDF$Cohort, levels = c('Stanford', 'PRJNA797778'))
    comDF$y = comDF[, yval] %>% as.numeric()
    avgVal = comDF %>% group_by(species) %>%
      summarise(meanVal = median(y))
    comDF$species %<>% factor(., levels = avgVal$species[order(avgVal$meanVal)])
    ppG = ggboxplot(comDF, x = 'species', y = 'y', color = 'Cohort', orientation = "horizontal",
              add = 'jitter', add.params = list(size = 0.3, width = 0.1),
              ylab = yLabs[yval], xlab = 'Species') +
      scale_color_manual(values = c('Stanford' = '#F8766D', 'PRJNA797778' = '#00BFC4')) +
      mytheme + theme(legend.position = 'top', legend.margin = margin(0,0,-5,-60), aspect.ratio = 2/1) 
    print(ppG)
  }
  dev.off()
}

if(T){ # compare coverage
  yvals = c('coverage_log', "breadth", "breadth_minCov", 'nucl_diversity', 'D_rev', "pnps", "dnds")
  yvals = c('coverage_log', "breadth", "breadth_minCov", 'nucl_diversity', 'D_rev', "pnps", "dnds")
  yLabs = list(coverage_log = 'log10(Coverage)', breadth = 'Breath', breadth_minCov = 'Breath (coverage > 5)', 
               nucl_diveristy = expression(pi), D_rev = "1 - D'", pnps = 'pN/pS', dnds = 'dN/dS')
  cairo_pdf('./plot/compare_of_two_cohorts_facet.pdf', width = 20, height = 5, onefile = T)

  comDF = rbind(
    data.frame(refDF1[refDF1$species %in% overlapSP_sel, c('species', yvals)],
               Cohort = 'Stanford'),
    data.frame(refDF2[refDF2$species %in% overlapSP_sel, c('species', yvals)],
               Cohort = 'PRJNA797778')
  )
  table(comDF$species)
  comDF2 = melt(comDF)
  avgVal = comDF %>% group_by(species) %>%
    summarise(meanVal = median(coverage_log))
  comDF2$species %<>% factor(., levels = avgVal$species[order(avgVal$meanVal)])
  labeller <- function(variable,value){
    return(yLabs[value])
  }
  ppG = ggboxplot(comDF2, x = 'species', y = 'value', color = 'Cohort', orientation = "horizontal",
                  add = 'jitter', add.params = list(size = 0.3, width = 0.1),
                  ylab = 'Value', xlab = 'Species') +
    scale_color_manual(values = c('Stanford' = '#F8766D', 'PRJNA797778' = '#00BFC4')) +
    mytheme + theme(legend.position = 'top', legend.margin = margin(0,0,-5,-60), aspect.ratio = 2/1,
                    strip.text = element_text(size = 13), axis.text.x = element_text(size = 9)) +
    facet_grid(.~variable, scales = 'free', labeller = labeller)
  print(ppG)

  dev.off()
}
