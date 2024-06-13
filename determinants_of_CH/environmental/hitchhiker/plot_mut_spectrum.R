
## 
## Mutation types:
mut_type <- function (nuc_subs) {
  muts <- nuc_subs
  types <- unlist(muts)
  types <- gsub("G>T", "C>A", types)
  types <- gsub("G>C", "C>G", types)
  types <- gsub("G>A", "C>T", types)
  types <- gsub("A>T", "T>A", types)
  types <- gsub("A>G", "T>C", types)
  types <- gsub("A>C", "T>G", types)
  return(types)
}


## Spectrum
plot_mut_spectrum <- function(mutations,
                              save=NULL,
                              tissue=NULL,
                              genomeFile = ref_genome,
                              add_mut_type="Nonsynonymous CHIP ",
                              add_to_title="", 
                              unique_var=NULL)
{
  mutations <- as.data.frame(mutations)
  
  colnames(mutations) <- c("chr","pos","ref","mut")
  
  if(!is.null(unique_var) && unique_var == TRUE){
    mutations <- unique(mutations)
  }
  
  
  mutations$pos <- as.numeric(mutations$pos)
  
  mutations$trinuc_ref <- as.character(Biostrings::getSeq(BSgenome::getBSgenome(ref_genome), 
                                                          mutations$chr, 
                                                          mutations$pos - 1, 
                                                          mutations$pos + 1))
  
  
  
  
  # as.vector(scanFa(genomeFile, 
  #                  GRanges(mutations$chr, 
  #                          IRanges(as.numeric(mutations$pos)-1, 
  #                                  as.numeric(mutations$pos)+1))))
  
  # 2. Annotating the mutation from the pyrimidine base
  ntcomp <- c(T="A",G="C",C="G",A="T")
  mutations$sub <- paste(mutations$ref,
                         mutations$mut,
                         sep=">")
  
  mutations$trinuc_ref_py <- mutations$trinuc_ref
  
  for (j in 1:nrow(mutations)) {
    if (mutations$ref[j] %in% c("A","G")) { # Purine base
      
      mutations$sub[j] <- paste(ntcomp[mutations$ref[j]],
                                ntcomp[mutations$mut[j]],
                                sep=">")
      
      mutations$trinuc_ref_py[j] <- paste(ntcomp[rev(strsplit(mutations$trinuc_ref[j],
                                                              split="")[[1]])],
                                          collapse="")
    }
  }
  
  # 3. Counting subs
  freqs <- table(paste(mutations$sub,
                       paste(substr(mutations$trinuc_ref_py, 1, 1),
                             substr(mutations$trinuc_ref_py, 3, 3),
                             sep="-"),
                       sep=","))
  
  sub_vec <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  
  ctx_vec <- paste(rep(c("A", "C", "G", "T"), each=4),
                   rep(c("A", "C", "G", "T"), times=4),
                   sep="-")
  
  full_vec <- paste(rep(sub_vec,each=16),
                    rep(ctx_vec,times=6),
                    sep=",")
  
  freqs_full <- freqs[full_vec]; freqs_full[is.na(freqs_full)] <- 0; names(freqs_full) <- full_vec
  
  xstr <- paste(substr(full_vec, 5, 5), 
                substr(full_vec, 1, 1), 
                substr(full_vec, 7, 7), 
                sep="")
  
  if(is.null(save)){dev.new(width=12, height=4)}
  colvec <- rep(c("dodgerblue", "black", "red", 
                  "grey70", "olivedrab3", "plum2"),
                each=16)
  
  y <- freqs_full; maxy <- max(y)
  
  h <- barplot(y, las=2, col=colvec, border=NA, 
               ylim=c(0, maxy*1.5), space=1, cex.names=0.6, 
               names.arg=xstr, 
               ylab="Number of mutations", 
               main=paste0("Total number of ",add_mut_type,"mutations: ",
                           nrow(mutations),
                           add_to_title))
  
  for (j in 1:length(sub_vec)) {
    
    xpos <- h[c((j-1)*16+1, j*16)]
    
    rect(xpos[1]-0.5, maxy*1.2, xpos[2]+0.5, 
         maxy*1.3, border=NA, col=colvec[j*16])
    
    text(x=mean(xpos), y=maxy*1.3, 
         pos=3, label=sub_vec[j])
  } 
  if(!is.null(save)){ 
    dev.copy(pdf,
             save,
             width=12,
             height=4); 
    dev.off()}
}

