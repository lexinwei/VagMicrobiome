
current_dir <- getwd()
dir.create(file.path(current_dir, opt$out_directory), showWarnings = FALSE)

outfile <- paste(current_dir, opt$out_directory, outfile_prefix, sep="/")
taxon_separator = 's__'
tree = '/Users/xinwei/workspace/vag/reference_code/MetaPhlAn-master/metaphlan/utils/mpa_vJan21_CHOCOPhlAnSGB_202103.nwk'

### table preprocessing ###
mpa_table <- read.table('./metaphlan4_results/merged_abundance_table.txt', 
                        comment.char = '#', sep = '\t', header = TRUE, check.names=FALSE)
mpa_tableS = mpa_table[grep('s__', mpa_table$clade_name), ]


cldDF = mpa_table$clade_name %>% str_split_fixed(., '\\|', Inf) %>% as.data.frame()
unique(cldDF$V7[cldDF$V7 != ""]) %>% length()
unique(cldDF$V8) %>% length()
cldDF2 = cldDF[cldDF$V8 != "", 7:8]
unique(cldDF2$V7) %>% length()
unique(cldDF2$V7) %>% length()

# check if NCBI id is present
vec <- grepl("ncbi", colnames(mpa_table), ignore.case=TRUE)
if(any(vec)){
  # keep all the columns except the one with NCBI id
  mpa_table <- mpa_table[grep(taxon_separator, mpa_table[,1]), !vec]
} else {
  mpa_table <- mpa_table[grep(taxon_separator, mpa_table[,1]),]
}

if(taxon_separator == "t__"){
  mpa_table[,1] <- gsub(".+\\|t__SGB", "", mpa_table[,1])
  mpa_table[,1] <- gsub(".+\\|t__EUK", "", mpa_table[,1])
} else { 
  mpa_table <- mpa_table[!grepl('t__', mpa_table[,1]),]
  mpa_table[,1] <- gsub(".+\\|s__", "", mpa_table[,1])
}

mpa_table[,1] <- gsub("_group$", "", mpa_table[,1])
rownames(mpa_table) <- mpa_table[,1]
mpa_table <- mpa_table[,-1]

# remove samples with all unknowns
removed <- which(colSums(mpa_table) == 0)

if(length(removed)>0){
  
  if(length(removed)==1){
    message = "# WARNING: 1 sample with 100% unknown species was removed from the input table."
  } else {
    message = paste0("# WARNING: ", length(removed), " samples with 100% unknown species were removed from the input table.")
  }
  
  write(message, stdout())
  write(paste(names(removed), collapse='\n'), stdout())
  
  write(message, file=paste0(outfile, '_samples.log'))
  write(paste(names(removed), collapse='\n'), file=paste0(outfile, '_samples.log'), append = TRUE)
  
  # remove samples
  mpa_table <- mpa_table[, -removed]
}

### Data transformation
mpa_table <- mpa_table / 100


### Beta diversity ###

if (opt$diversity == "beta"){
  
  # Bray-Curtis
  if (opt$metric == "bray-curtis"){
    mat <- rbiom::beta.div(as.matrix(mpa_table), method="bray-curtis", weighted=TRUE)
  }
  
  # Jaccard
  if (opt$metric == "jaccard"){
    mat <- rbiom::beta.div(as.matrix(mpa_table), method="jaccard", weighted=FALSE)
  }
  
  # Unifrac
  if (grepl("unifrac", opt$metric)){
    mpa_tree <- ape::read.tree(tree)
    
    if(opt$taxon_separator == "s__"){
      mpa_tree$tip.label <- gsub(".+\\|s__", "", mpa_tree$tip.label)
    }
    
    removed <- setdiff(rownames(mpa_table), mpa_tree$tip.label)
    
    if(length(removed)){
      message = paste0("# WARNING: ", length(removed), " species not present in the tree were removed from the input profile.")
      write(message, stdout())
      write(paste(removed, collapse='\n'), stdout())
      
      write(message, file=paste0(outfile, '_species.log'))
      write(paste(removed, collapse='\n'), file=paste0(outfile, '_species.log'), append = TRUE)
    }
    filt_tree <- ape::keep.tip(mpa_tree, setdiff(rownames(mpa_table), removed))
    filt_mpa_table <- mpa_table[filt_tree$tip.label,]
    
    # check again if after species removal some samples have 0s for all the remaining species, and remove them
    removed <- which(colSums(filt_mpa_table) == 0)
    
    if(length(removed)){
      message = paste0("# WARNING: after removal of species not in the tree, ", length(removed), " samples with 0 abundance of the remaining species were removed from the input table.")
      
      write(message, stdout())
      write(paste(names(removed), collapse='\n'), stdout())
      
      if(file.exists(paste0(outfile, '_samples.log'))){
        write(message, file=paste0(outfile, '_samples.log'), append = TRUE)
      } else { 
        write(message, file=paste0(outfile, '_samples.log'))
      }
      
      write(paste(names(removed), collapse='\n'), file=paste0(outfile, '_samples.log'), append = TRUE)
      
      # remove samples
      filt_mpa_table <- filt_mpa_table[, -removed]
    }
    
    if (opt$metric == "weighted-unifrac"){
      mat <- rbiom::beta.div(as.matrix(filt_mpa_table), tree=filt_tree, method="unifrac", weighted=TRUE)
      
    } else if (opt$metric == "unweighted-unifrac"){
      mat <- rbiom::beta.div(as.matrix(filt_mpa_table), tree=filt_tree, method="unifrac", weighted=FALSE)
    } 
    
  }
  
  # CLR or Aitchison
  if (opt$metric == "clr" || opt$metric == "aitchison"){
    # Centered Log-Ratio
    ait_norm_mpa_table <- compositions::clr(mpa_table)
    mat <- as.matrix(compositions::as.data.frame.rmult(ait_norm_mpa_table))

    if (opt$metric == "aitchison"){
      # Aitchison
      mat <- rbiom::beta.div(mat, method="euclidean", weighted=TRUE)
    }
  }
  
  ### Alpha Diversity ###
  
} else if (opt$diversity == "alpha"){
  
  # Richness
  if (opt$metric == "richness"){
    mat <- microbiome::alpha(mpa_table, index = c("richness_observed"))
  }
  
  # Shannon
  if (opt$metric == "shannon"){
    mat <- microbiome::alpha(mpa_table, index = c("diversity_shannon"))
  }
  
  # Simpson
  if (opt$metric == "simpson"){
    mat <- microbiome::alpha(mpa_table, index = c("diversity_gini_simpson"))
  }
  
  # Gini
  if (opt$metric == "gini"){
    mat <- microbiome::alpha(mpa_table, index = c("dominance_gini"))
    row.names(mat) <- colnames(mpa_table)
  }
}

write.table(as.matrix(mat), paste0(outfile, '_', opt$metric, '.tsv'), sep = '\t', quote = FALSE)
