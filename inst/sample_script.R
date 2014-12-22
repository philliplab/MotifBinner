library(MotifBinner)

x <- readFastq("~/projects/MotifBinner/data/CAP177_2040_v1merged.fastq")
x <- x@sread
x <- padAndClip(x, IRanges(start = 28, end=width(x)), 
                Lpadding.letter="+", Rpadding.letter="+")
seq_data <- x
prefix <- "CCAGCTGGTTATGCGATTCTMARGTG"
suffix <- "CTGAGCGTGTGGCAAGGCCC"
motif_length <- 9
max.mismatch <- 5
fixed = FALSE

y <- extract_motifs(seq_data, prefix, suffix, motif_length, max.mismatch, fixed)

bin_seqs <- bin_by_name(y$matched_seq)

print(y)

print(table(table(names(y$matched_seq))))

print(length(bin_seqs))

plot_bin <- function(bin){
  plot_name <- paste0(names(bin)[1], '_', length(bin), '_bin.pdf')
  names(bin) <- paste0(names(bin), 1:length(bin))
  sdz <- stringDist(bin)
  sdz[sdz>99] <- 99
  #writeXStringSet(bin, gsub(".pdf", ".FASTA", plot_name))
  #pdf(plot_name, h = 15)
  plot(bionj(sdz))
  #dev.off()
}

for (bin_name in names(which(lapply(bin_seqs, length) > 100))){
  print(bin_name)
#  bin <- bin_seqs[[bin_name]]
#  plot_bin(bin)
  sdz <- stringDist(bin_seqs[[bin_name]])
  sdzu <- stringDist(unique(bin_seqs[[bin_name]]))
  
#  par(mfrow=c(3,4))
  hist(sdz)
#  x11()
  hist(sdzu)
}

remove_contamination <- function(bin, method = 'really bad'){
  if (method == 'really bad'){
    dists <- NULL
    bin_dists <- stringDist(bin)
    dmat <- as.matrix(bin_dists)
    row.names(dmat) <- 1:nrow(dmat)
    removed_sequences <- NULL
    orig_dmat <- dmat
    pdr_psr_rat <- 1000
    threshold <- 2
    while((nrow(unique(dmat)) > 1) & (pdr_psr_rat > threshold)){
      dists <- c(dists, mean(dmat))
      dvec <- apply(dmat, 1, sum)
      max_indx <- which(dvec == max(dvec))
      new_dmat <- dmat[-max_indx, -max_indx]
      
#      print(dim(dmat))
#      print(mean(dmat))
#      print(max_indx)
      
      print('percentage sequences removed')
      psr <- length(max_indx) / nrow(dmat)
      print(psr)
      print('percentage distance reduction')
      pdr <- (mean(dmat) - mean(new_dmat)) / mean(orig_dmat)
      print(pdr)
      pdr_psr_rat <- pdr / psr
      print(pdr_psr_rat)
      print('--------------------------------------')

      if (pdr_psr_rat > threshold){
        removed_sequences <- c(removed_sequences, row.names(dmat)[max_indx])
        dmat <- new_dmat
      }
    }
    return(list(dists = dists,
                sequences = row.names(dmat),
                removed_sequences = removed_sequences))
  } else { stop('unknown method') }
}


