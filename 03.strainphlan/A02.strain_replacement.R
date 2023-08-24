source('./script/library_function.R')
load('./data/RData/clusterDF3.RData')
load('./data/RData/metaData.RData')
load('./data/RData/swabData.RData')
load('./data/RData/colorList.RData')
###### estimate the transition of strain #####

# set window of sequential samples and filter out samples
# window = c(11, 17) # biweekly, 11, 17 about two weeks; 
# window = c(4, 10) # weekly, 4, 10, ~ about 1 weeks
window = c(3.5, 55.3) # all
# space = 'biweekly'
# space = 'weekly'
space = 'all'
if(T){
  sequentialDF = data.frame()
  for(p in metaData$pt_ID){
    subDF = swabData[swabData$pt_ID == p & swabData$trimester != 'P', ]
    subDF = subDF[naturalorder(subDF$pt_ID.u), ]
    spaceDays = (subDF$Sample_GA[2:nrow(subDF)] - subDF$Sample_GA[1:(nrow(subDF)-1)])*7
    preInd = which(spaceDays >= window[1] & spaceDays <= window[2])
    nextInd = preInd + 1
    sequentialDF = rbind(sequentialDF, data.frame(
      preSample = subDF$SampleID.u[preInd],
      nextSample = subDF$SampleID.u[nextInd]
    ))
  }
  nrow(sequentialDF)
}

cairo_pdf(paste0("./plot/strain_markov_chain_", space, "3.pdf"), 
          width=4.41, height=3.72, onefile = T)
sp = "Gardnerella vaginalis (SGB17301)"

candSP = c('Candida albicans',
  'Gardnerella vaginalis (SGB17301)', 
           'Gardnerella vaginalis (SGB17302)',
  'Gardnerella vaginalis (SGB17305)',
           "Gardnerella vaginalis (SGB17307)",
           'Gardnerella vaginalis (SGB21500)',
           'Gardnerella vaginalis (SGB7097)',
           # 'Hungateiclostridiaceae SGB4003',
           'Veillonellaceae bacterium DNF00626',
           # 'Dialister micraerophilus',
           'Megasphaera genomosp type 1',
           # 'Anaerococcus tetradius',
           # 'Peptoniphilus harei',
           # 'Peptoniphilus lacrimalis',
           'Lactobacillus iners',
  'Lactobacillus jensenii (SGB7034)',
           'Lactobacillus jensenii (SGB7035)',
           'Lactobacillus gasseri',
           'Lactobacillus crispatus',
           'Lactobacillus vaginalis',
           'Lactobacillus coleohominis',
           'Aerococcus christensenii',
           'Coriobacteriales bacterium DNF00809',
           'Actinomycetaceae SGB989',
           'Fannyhessea vaginae (SGB990)',
           'Fannyhessea vaginae (SGB991)')
candSP = unique(clusterDF3$species)
for(sp in candSP){ # transition rate matrix and Markov Chain
  print(sp)
  strainType = clusterDF3[clusterDF3$species == sp & clusterDF3$cluster != '0', ]
  using_Col = 'cluster' 
  # strainType$cluster %<>% as.factor()
  sequentialDF2 = sequentialDF
  sequentialDF2$preState = strainType[match(sequentialDF2$preSample, strainType$samples), using_Col]
  sequentialDF2$nextState = strainType[match(sequentialDF2$nextSample, strainType$samples), using_Col]
  sequentialDF2 = na.omit(sequentialDF2)
  if(nrow(sequentialDF2) == 0){
    next
  }
  nstates = length(unique(c(sequentialDF2$preState, sequentialDF2$nextState)))
  if(nstates == 1){
    next
  }
  # if(sum(sequentialDF2$nextState == sequentialDF2$preState) == nrow(sequentialDF2)){
  #   next
  # }
  CSTs = unique(c(sequentialDF2$preState, sequentialDF2$nextState))
  sequentialDF2$preState %<>% factor(., levels = CSTs)
  sequentialDF2$nextState %<>% factor(., levels = CSTs)
  ttab <- table(sequentialDF2$preState, sequentialDF2$nextState) # prevstate=row, curstate=col
  trans <- matrix(ttab, nrow=nstates)
  trans <- trans/rowSums(trans)  # Normalize row sums to 1
  CSTtrans <- trans
  colnames(CSTtrans) <- CSTs
  rownames(CSTtrans) <- CSTs
  t_persist <- -1/log(diag(CSTtrans))
  table_CSTtrans = as.data.frame(CSTtrans)
  row.names(table_CSTtrans) = row.names(CSTtrans)
  # l = list(ttab, table_CSTtrans)
  # write.xlsx(l, file = './table/CST_biweekly_trans.xlsx', row.name =T)
  # grid.table(CSTtrans %>% round(., 3)) # Paper
  t_persist # Paper
  mcPreg <- new("markovchain", states=CSTs,
                transitionMatrix = trans, name="PregCST")
  mcPreg
# }
# if(T){ # Set up igraph of the markov chain
  library(markovchain)
  library(igraph)
  netMC <- markovchain:::.getNet(mcPreg, round = TRUE)
  wts <- E(netMC)$weight/100
  
  edgel <- get.edgelist(netMC)
  elcat <- paste(edgel[,1], edgel[,2])
  elrev <- paste(edgel[,2], edgel[,1])
  edge.curved <- sapply(elcat, function(x) x %in% elrev)
  
  samdf_def <- strainType
# }
# 
# if(T){
  library(igraph)
  default.par <- par(no.readonly = TRUE)
  # Define color scale
  # Plotting function for markov chain
  plotMC <- function(object, ...) {
    netMC <- markovchain:::.getNet(object, round = TRUE)
    plot.igraph(x = netMC, ...)  
  }
  # Color bar for the markov chain visualization, gradient in strength of preterm association
  color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title=NULL) {
    scale = (length(lut)-1)/(max-min)
    
    #    dev.new(width=1.75, height=5)
    # cur.par <- par(no.sreadonly=T)
    par(mar=c(0,7,1,7)+0.1, oma=c(0,0,0,0)+0.1, mgp=c(2,0.2,0))
    # par(ps = 10, cex = 0.8)
    # par(tcl=-0.2, cex.axis=0.8, cex.lab = 0.8)
    par(tcl=-0.2)
    plot(c(min,max), c(0,10), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
    axis(1, seq(0, 1, 0.25))
    for (i in 1:(length(lut)-1)) {
      x = (i-1)/scale + min
      rect(x,0,x+1/scale,10, col=lut[i], border=NA)
    }
  }
  
  vert.clrs <- sapply(states(mcPreg), function(x) colorStrain[x])
  vert.sz <- 10 + 2*sapply(states(mcPreg), 
                          function(x) length(unique(samdf_def[samdf_def[, using_Col]==x,"pt_ID"])))
  vert.sz <- vert.sz * 0.85
  vert.font.clrs <- c("white", "white", "white", "white", "white")
  # E(netMC) to see edge list, have to define loop angles individually by the # in edge list, not vertex
  edge.loop.angle = c(0, 0, 0, 0, 3.14, 3.14, 0, 0, 0, 3.14, 3.14, 0.60, 0, 0, 0, 0.45, 0,  0, 0, 0, 0)-0.45
  
  edge.loop.angle = rep(0, length(E(netMC)))
  edgesMat = ends(netMC,es = E(netMC))
  edgesList = paste0(edgesMat[, 1], ' -> ', edgesMat[, 2])
  edge.loop.angle[edgesList == 'I -> I'] = -0.45
  edge.loop.angle[edgesList == 'II -> II'] = -0.45
  edge.loop.angle[edgesList == 'III -> III'] = 3.25
  edge.loop.angle[edgesList == 'IV -> IV'] = 0.45
  
  
  # layout <- matrix(c(0.6,0.95, 0.43,1, 0.3,0.66, 0.55,0.3, 0.75,0.65), nrow=5, ncol=2, byrow=T)
  
  # Colored by association with preterm birth
  # library(viridis)
  # pal <- magma(101)
  # layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE), heights=c(1,5))
  # color.bar(pal, min=0, max=1, nticks=6, title="Fraction preterm")
  par(mar=c(1,1,1,1)+1)
  edge.arrow.size=0.6
  edge.arrow.width=1.2
  edge.width = (15*wts + 0.1)*0.6
  edge.labels <- as.character(round(E(netMC)$weight/100, 3))
  # edge.labels[edge.labels<0.4] <- NA  # labels only for self-loops
  par(family = "serif")
  plotMC(mcPreg, 
         edge.arrow.size=edge.arrow.size, edge.arrow.width = edge.arrow.width,
         edge.label = edge.labels, edge.label.cex=1.2, edge.label.color="black",
         # FIX EDGE LABELS FOR PUBLICATION IN POST-PROCESSING
         edge.width=edge.width, edge.curved=edge.curved, edge.color = 'lightgray',
         vertex.color=vert.clrs, vertex.size=(vert.sz),
         vertex.label.font = 2, vertex.label.cex = 1, main = sp,
         vertex.label.color = 'black', vertex.frame.color = 'black', 
         edge.loop.angle = edge.loop.angle)
  #dev.off()
  par(default.par)
}
dev.off()

###### strain transition vs ethnicity #####
load('./data/RData/clusterDF3.RData')
table(clusterDF3$species) %>% length()
table(clusterDF3$species, clusterDF3$cluster) 

candSP = setdiff(unique(clusterDF3$species), c(
  'Bifidobacteriaceae bacterium NR047',
  'Candida albicans',
  'Lactobacillus gasseri'
))
candSP = c(
           'Gardnerella vaginalis (SGB17302)',
           "Gardnerella vaginalis (SGB17307)",
           'Veillonellaceae bacterium DNF00626',
           'Lactobacillus iners',
           'Lactobacillus jensenii (SGB7035)',
           'Lactobacillus crispatus',
           'Fannyhessea vaginae (SGB991)')
ifreplDF = data.frame()
for (s in candSP) {
  print(s)
  tmpDf = clusterDF3[clusterDF3$species == s,]
  tmpDf = tmpDf[tmpDf$cluster != 0, ]
  taB = table(tmpDf$pt_ID, tmpDf$cluster)
  tmpDf2 = data.frame(
    pt_ID = row.names(taB), 
    Replacement = ifelse(rowSums(taB > 1) > 1, 'Yes', 'No')
  )
  tmpDf2$Ethnicity = metaData$Ethnicity2[match(tmpDf2$pt_ID, metaData$SampleID)]
  print(table(tmpDf2$Replacement, tmpDf2$Ethnicity) )
  chiQ = chisq.test(table(tmpDf2$Replacement, tmpDf2$Ethnicity) )
  print(chiQ$p.value)
  
  tmpDf2$Term = metaData$Term_char2[match(tmpDf2$pt_ID, metaData$SampleID)]
  print(table(tmpDf2$Replacement, tmpDf2$Term) )
  chiQ2 = chisq.test(table(tmpDf2$Replacement, tmpDf2$Term) )
  print(chiQ2$p.value)
  
  ifreplDF = rbind(ifreplDF, tmpDf2)
}
print(table(ifreplDF$Replacement, ifreplDF$Ethnicity) )
chiQ1 = chisq.test(table(ifreplDF$Replacement, ifreplDF$Ethnicity),correct = F)
print(chiQ1$p.value)
print(table(ifreplDF$Replacement, ifreplDF$Term) )
chiQ2 = chisq.test(table(ifreplDF$Replacement, ifreplDF$Term),correct = F)
print(chiQ2$p.value)

