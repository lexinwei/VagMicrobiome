ado1 = adonis2(abuMat_T~Age + BMI + Term + Ethnicity + CST, data = bray_curtis_pcoa_df, method="bray", by = 'margin')

anosim_result<-anosim(abuMat_T,bray_curtis_pcoa_df$Term,permutations = 999)
plot(anosim_result, col = c('gray', colorTerm))
mod <- betadisper(bray_curtis_dist, bray_curtis_pcoa_df$Term)
permutest(mod)