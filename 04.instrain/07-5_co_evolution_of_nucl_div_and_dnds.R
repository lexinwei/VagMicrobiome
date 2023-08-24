source('~/workspace/library_function.R')
load('./data/metaData.RData')
load('./data/swabData.RData')
load('./data/refDF.RData')
load('./data/magDF.RData')
load('./data/spName.RData')

###### by nucleotide diversity ref ####
# dcast df into matrix
refDF = refDF[!grepl('KIT', refDF$pt_ID), ]
refMat = reshape2::dcast(refDF, pt_ID.u ~ species, 
                         value.var = 'nucl_diversity')
refMat2 = apply(refMat[, -1], 2, as.numeric) %>% as.matrix()
row.names(refMat2) = refMat$pt_ID.u
refMat3 = refMat2[which(rowSums(!is.na(refMat2)) > 1),]
refMat4 = refMat3[, which(colSums(!is.na(refMat3)) > 1)]
rowSums(!is.na(refMat4))
colSums(!is.na(refMat4))
pheatmap(refMat4, cluster_rows = F, cluster_cols = F, scale = 'none')

# calculate paired-wise correlation
pdf('./plot/genome_wise_nucl_pairwise_corr_ref.pdf', width = 5, height = 5)
if(T){
cormethod = 'spearman'
spl1 = c()
spl2 = c()
cop = c()
corre = c()
for(ii in 1:(ncol(refMat4)-1)){
  for(jj in (ii + 1):ncol(refMat4)){
    if(length(intersect(which(!is.na(refMat4[, ii])), which(!is.na(refMat4[, jj])))) > 1){
      xycor = cor.test(refMat4[, ii], refMat4[, jj], method = 'spearman')
      cop = c(cop, xycor$p.value)
      corre = c(corre, xycor$estimate %>% as.numeric())
      spl1 = c(spl1, colnames(refMat4)[ii])
      spl2 = c(spl2, colnames(refMat4)[jj])
    }
  }
}

corrRes = data.frame(
  sp1 = spl1,
  sp2 = spl2,
  pval = cop,
  esti = corre
)
corrRes$pval.adj = p.adjust(corrRes$pval, method = "fdr")
xlsx::write.xlsx(corrRes, file = './Routput/genome_wise_nucl_dnds_pairwise_corr.xlsx',
                 row.names = F, append = F, sheetName = 'nucl_div_ref')
corrRes = corrRes[corrRes$pval.adj < 0.05, ]
corrRes_rev = data.frame(
  sp1 = spl2,
  sp2 = spl1,
  pval = cop,
  esti = corre
)
corrRes_rev$pval.adj = p.adjust(corrRes_rev$pval, method = "fdr")
corrRes_rev = corrRes_rev[corrRes_rev$pval.adj < 0.05, ]
corrMat = reshape2::dcast(rbind(corrRes, corrRes_rev), sp1 ~ sp2, value.var = 'esti')
corrMat = data.frame(corrMat, row.names = 1) %>% as.matrix()
corrMatP = reshape2::dcast(rbind(corrRes, corrRes_rev), sp1 ~ sp2, value.var = 'pval.adj')
corrMatP = data.frame(corrMatP, row.names = 1) %>% as.matrix()
corrplot(corrMat,  diag = F, method="color", type = 'upper',
         col = brewer.pal(n=10, name="PuOr"),
         outline= F,
         p.mat = corrMatP,
         addrect = 3, rect.lwd = 1.5, na.label = '.', na.label.col = 'grey', 
         sig.level = 0.05, mar = c(1, 1, 2, 1),
         insig = "blank", tl.cex = 0.4, title = cormethod,
         tl.srt= 45, tl.col="black")
for(i in 1:nrow(corrRes)){
  i1 = corrRes$sp1[i]
  i2 = corrRes$sp2[i]
  refMat4DF = data.frame(refMat4)
  refMat4DFxy = data.frame(
    x = refMat4DF[, i1],
    y = refMat4DF[, i2],
    sample = row.names(refMat4)
  )
  refMat4DFxy %<>% na.omit()
  xycor = cor.test(refMat4DFxy$x, refMat4DFxy$y, method = 'spearman')
  cop = xycor$p.value
  corre = xycor$estimate %>% as.numeric()
  corp = ggplot(refMat4DFxy, aes(x = x,  y = y)) +
    geom_point() +
    # geom_abline(linetype = 'dashed') +
    geom_smooth(method = 'lm', se = F) +
    labs(x = i1, y = i2, title = 'Nucleotide diversity') +
    annotate("text", x = max(refMat4DFxy$x) * 0.75, y = max(refMat4DFxy$y) * 0.15, 
             label = paste0('Spearman ', round(corre, 3), ', p=', p_format(cop)),) +
    mytheme3
  print(corp)
}
dev.off()
}

###### by nucleotide diversity mag ####
magDF = magDF[!grepl('KIT', magDF$pt_ID), ]
magMat = reshape2::dcast(magDF, pt_ID.u ~ genome, 
                         value.var = 'nucl_diversity')
magMat2 = apply(magMat[, -1], 2, as.numeric) %>% as.matrix()
row.names(magMat2) = magMat$pt_ID.u
colnames(magMat2) = paste0(magDF$species[match(colnames(magMat2), magDF$genome)], ' (', colnames(magMat2), ')')
magMat3 = magMat2[which(rowSums(!is.na(magMat2)) > 1),]
magMat4 = magMat3[, which(colSums(!is.na(magMat3)) > 1)]
rowSums(!is.na(magMat4))
colSums(!is.na(magMat4))
pheatmap(magMat4, cluster_rows = F, cluster_cols = F, scale = 'none')

# calculate paired-wise correlation
pdf('./plot/genome_wise_nucl_pairwise_corr_mags.pdf', width = 5, height = 5)
if(T){
  cormethod = 'spearman'
  spl1 = c()
  spl2 = c()
  cop = c()
  corre = c()
  for(ii in 1:(ncol(magMat4)-1)){
    for(jj in (ii + 1):ncol(magMat4)){
      if(length(intersect(which(!is.na(magMat4[, ii])), which(!is.na(magMat4[, jj])))) > 1){
        xycor = cor.test(magMat4[, ii], magMat4[, jj], method = 'spearman')
        cop = c(cop, xycor$p.value)
        corre = c(corre, xycor$estimate %>% as.numeric())
        spl1 = c(spl1, colnames(magMat4)[ii])
        spl2 = c(spl2, colnames(magMat4)[jj])
      }
    }
  }
  
  corrRes = data.frame(
    sp1 = spl1,
    sp2 = spl2,
    pval = cop,
    esti = corre
  )
  corrRes$pval.adj = p.adjust(corrRes$pval, method = "fdr")
  xlsx::write.xlsx(corrRes, file = './Routput/genome_wise_nucl_dnds_pairwise_corr.xlsx',
                   row.names = F, append = T, sheetName = 'nucl_div_mag')
  corrRes = corrRes[corrRes$pval.adj < 0.05, ]
  corrRes_rev = data.frame(
    sp1 = spl2,
    sp2 = spl1,
    pval = cop,
    esti = corre
  )
  corrRes_rev$pval.adj = p.adjust(corrRes_rev$pval, method = "fdr")
  corrRes_rev = corrRes_rev[corrRes_rev$pval.adj < 0.05, ]
  corrMat = reshape2::dcast(rbind(corrRes, corrRes_rev), sp1 ~ sp2, value.var = 'esti')
  corrMat = data.frame(corrMat, row.names = 1) %>% as.matrix()
  colnames(corrMat) = row.names(corrMat)
  corrMatP = reshape2::dcast(rbind(corrRes, corrRes_rev), sp1 ~ sp2, value.var = 'pval.adj')
  corrMatP = data.frame(corrMatP, row.names = 1) %>% as.matrix()
  colnames(corrMatP) = row.names(corrMatP)
  corrplot(corrMat,  diag = F, method="color", type = 'upper',
           col = brewer.pal(n=10, name="PuOr"),
           outline= F,
           p.mat = corrMatP,
           addrect = 3, rect.lwd = 1.5, na.label = '.', na.label.col = 'grey', 
           sig.level = 0.05, mar = c(1, 1, 2, 1),
           insig = "blank", tl.cex = 0.4, title = cormethod,
           tl.srt= 45, tl.col="black")
  for(i in 1:nrow(corrRes)){
    i1 = corrRes$sp1[i]
    i2 = corrRes$sp2[i]
    magMat4DF = data.frame(magMat4)
    colnames(magMat4DF) = colnames(magMat4)
    magMat4DFxy = data.frame(
      x = magMat4DF[, i1],
      y = magMat4DF[, i2],
      sample = row.names(magMat4)
    )
    magMat4DFxy %<>% na.omit()
    xycor = cor.test(magMat4DFxy$x, magMat4DFxy$y, method = 'spearman')
    cop = xycor$p.value
    corre = xycor$estimate %>% as.numeric()
    corp = ggplot(magMat4DFxy, aes(x = x,  y = y)) +
      geom_point() +
      # geom_abline(linetype = 'dashed') +
      geom_smooth(method = 'lm', se = F) +
      labs(x = i1, y = i2, title = 'Nucleotide diversity') +
      annotate("text", x = max(magMat4DFxy$x) * 0.75, y = max(magMat4DFxy$y) * 0.15, 
               label = paste0('Spearman ', round(corre, 3), ', p=', p_format(cop)),) +
      mytheme3
    print(corp)
  }
  dev.off()
}

###### by dnds ref ####
# dcast df into matrix
refDF = refDF[!grepl('KIT', refDF$pt_ID), ]
refMat = reshape2::dcast(refDF, pt_ID.u ~ species, 
                         value.var = 'dn_ds_common')
refMat2 = apply(refMat[, -1], 2, as.numeric) %>% as.matrix()
row.names(refMat2) = refMat$pt_ID.u
refMat3 = refMat2[which(rowSums(!is.na(refMat2)) > 1),]
refMat4 = refMat3[, which(colSums(!is.na(refMat3)) > 1)]
rowSums(!is.na(refMat4))
colSums(!is.na(refMat4))
pheatmap(refMat4, cluster_rows = F, cluster_cols = F, scale = 'none')

# calculate paired-wise correlation
pdf('./plot/genome_wise_dnds_pairwise_corr_ref.pdf', width = 5, height = 5)
if(T){
  cormethod = 'spearman'
  spl1 = c()
  spl2 = c()
  cop = c()
  corre = c()
  for(ii in 1:(ncol(refMat4)-1)){
    for(jj in (ii + 1):ncol(refMat4)){
      if(length(intersect(which(!is.na(refMat4[, ii])), which(!is.na(refMat4[, jj])))) > 1){
        xycor = cor.test(refMat4[, ii], refMat4[, jj], method = 'spearman')
        cop = c(cop, xycor$p.value)
        corre = c(corre, xycor$estimate %>% as.numeric())
        spl1 = c(spl1, colnames(refMat4)[ii])
        spl2 = c(spl2, colnames(refMat4)[jj])
      }
    }
  }
  
  corrRes = data.frame(
    sp1 = spl1,
    sp2 = spl2,
    pval = cop,
    esti = corre
  )
  corrRes$pval.adj = p.adjust(corrRes$pval, method = "fdr")
  xlsx::write.xlsx(corrRes, file = './Routput/genome_wise_nucl_dnds_pairwise_corr.xlsx',
                   row.names = F, append = T, sheetName = 'dnds_ref')
  corrRes = corrRes[corrRes$pval.adj < 0.05, ]
  corrRes_rev = data.frame(
    sp1 = spl2,
    sp2 = spl1,
    pval = cop,
    esti = corre
  )
  corrRes_rev$pval.adj = p.adjust(corrRes_rev$pval, method = "fdr")
  corrRes_rev = corrRes_rev[corrRes_rev$pval.adj < 0.05, ]
  corrMat = reshape2::dcast(rbind(corrRes, corrRes_rev), sp1 ~ sp2, value.var = 'esti')
  corrMat = data.frame(corrMat, row.names = 1) %>% as.matrix()
  corrMatP = reshape2::dcast(rbind(corrRes, corrRes_rev), sp1 ~ sp2, value.var = 'pval.adj')
  corrMatP = data.frame(corrMatP, row.names = 1) %>% as.matrix()
  corrplot(corrMat,  diag = F, method="color", type = 'upper',
           col = brewer.pal(n=10, name="PuOr"),
           outline= F,
           p.mat = corrMatP,
           addrect = 3, rect.lwd = 1.5, na.label = '.', na.label.col = 'grey', 
           sig.level = 0.05, mar = c(1, 1, 2, 1),
           insig = "blank", tl.cex = 0.4, title = cormethod,
           tl.srt= 45, tl.col="black")
  for(i in 1:nrow(corrRes)){
    i1 = corrRes$sp1[i]
    i2 = corrRes$sp2[i]
    refMat4DF = data.frame(refMat4)
    refMat4DFxy = data.frame(
      x = refMat4DF[, i1],
      y = refMat4DF[, i2],
      sample = row.names(refMat4)
    )
    refMat4DFxy %<>% na.omit()
    xycor = cor.test(refMat4DFxy$x, refMat4DFxy$y, method = 'spearman')
    cop = xycor$p.value
    corre = xycor$estimate %>% as.numeric()
    corp = ggplot(refMat4DFxy, aes(x = x,  y = y)) +
      geom_point() +
      # geom_abline(linetype = 'dashed') +
      geom_smooth(method = 'lm', se = F) +
      labs(x = i1, y = i2, title = 'dN/dS') +
      annotate("text", x = max(refMat4DFxy$x) * 0.75, y = max(refMat4DFxy$y) * 0.15, 
               label = paste0('Spearman ', round(corre, 3), ', p=', p_format(cop)),) +
      mytheme3
    print(corp)
  }
  dev.off()
}

###### by dnds mag ####
magDF = magDF[!grepl('KIT', magDF$pt_ID), ]
magMat = reshape2::dcast(magDF, pt_ID.u ~ genome, 
                         value.var = 'dn_ds_common')
magMat2 = apply(magMat[, -1], 2, as.numeric) %>% as.matrix()
row.names(magMat2) = magMat$pt_ID.u
colnames(magMat2) = paste0(magDF$species[match(colnames(magMat2), magDF$genome)], ' (', colnames(magMat2), ')')
magMat3 = magMat2[which(rowSums(!is.na(magMat2)) > 1),]
magMat4 = magMat3[, which(colSums(!is.na(magMat3)) > 1)]
rowSums(!is.na(magMat4))
colSums(!is.na(magMat4))
pheatmap(magMat4, cluster_rows = F, cluster_cols = F, scale = 'none')

# calculate paired-wise correlation
pdf('./plot/genome_wise_dnds_corr_mags.pdf', width = 5, height = 5)
if(T){
  cormethod = 'spearman'
  spl1 = c()
  spl2 = c()
  cop = c()
  corre = c()
  for(ii in 1:(ncol(magMat4)-1)){
    for(jj in (ii + 1):ncol(magMat4)){
      if(length(intersect(which(!is.na(magMat4[, ii])), which(!is.na(magMat4[, jj])))) > 1){
        xycor = cor.test(magMat4[, ii], magMat4[, jj], method = 'spearman')
        cop = c(cop, xycor$p.value)
        corre = c(corre, xycor$estimate %>% as.numeric())
        spl1 = c(spl1, colnames(magMat4)[ii])
        spl2 = c(spl2, colnames(magMat4)[jj])
      }
    }
  }
  
  corrRes = data.frame(
    sp1 = spl1,
    sp2 = spl2,
    pval = cop,
    esti = corre
  )
  corrRes$pval.adj = p.adjust(corrRes$pval, method = "fdr")
  xlsx::write.xlsx(corrRes, file = './Routput/genome_wise_nucl_dnds_pairwise_corr.xlsx',
                   row.names = F, append = T, sheetName = 'dnds_mag')
  corrRes = corrRes[corrRes$pval.adj < 0.05, ]
  corrRes_rev = data.frame(
    sp1 = spl2,
    sp2 = spl1,
    pval = cop,
    esti = corre
  )
  corrRes_rev$pval.adj = p.adjust(corrRes_rev$pval, method = "fdr")
  corrRes_rev = corrRes_rev[corrRes_rev$pval.adj < 0.05, ]
  corrMat = reshape2::dcast(rbind(corrRes, corrRes_rev), sp1 ~ sp2, value.var = 'esti')
  corrMat = data.frame(corrMat, row.names = 1) %>% as.matrix()
  colnames(corrMat) = row.names(corrMat)
  corrMatP = reshape2::dcast(rbind(corrRes, corrRes_rev), sp1 ~ sp2, value.var = 'pval.adj')
  corrMatP = data.frame(corrMatP, row.names = 1) %>% as.matrix()
  colnames(corrMatP) = row.names(corrMatP)
  corrplot(corrMat,  diag = F, method="color", type = 'upper',
           col = brewer.pal(n=10, name="PuOr"),
           outline= F,
           p.mat = corrMatP,
           addrect = 3, rect.lwd = 1.5, na.label = '.', na.label.col = 'grey', 
           sig.level = 0.05, mar = c(1, 1, 2, 1),
           insig = "blank", tl.cex = 0.5, title = cormethod,
           tl.srt= 45, tl.col="black")
  for(i in 1:nrow(corrRes)){
    i1 = corrRes$sp1[i]
    i2 = corrRes$sp2[i]
    magMat4DF = data.frame(magMat4)
    colnames(magMat4DF) = colnames(magMat4)
    magMat4DFxy = data.frame(
      x = magMat4DF[, i1],
      y = magMat4DF[, i2],
      sample = row.names(magMat4)
    )
    magMat4DFxy %<>% na.omit()
    xycor = cor.test(magMat4DFxy$x, magMat4DFxy$y, method = 'spearman')
    cop = xycor$p.value
    corre = xycor$estimate %>% as.numeric()
    corp = ggplot(magMat4DFxy, aes(x = x,  y = y)) +
      geom_point() +
      # geom_abline(linetype = 'dashed') +
      geom_smooth(method = 'lm', se = F) +
      labs(x = i1, y = i2, title = 'dN/dS') +
      annotate("text", x = max(magMat4DFxy$x) * 0.75, y = max(magMat4DFxy$y) * 0.15, 
               label = paste0('Spearman ', round(corre, 3), ', p=', p_format(cop)),) +
      mytheme3
    print(corp)
  }
  dev.off()
}




