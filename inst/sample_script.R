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

plot_bin <- function(bin){
  plot_name <- paste0(names(bin)[1], '_', length(bin), '_bin.pdf')
  names(bin) <- paste0(names(bin), 1:length(bin))
  sdz <- stringDist(bin)
  sdz[sdz>99] <- 99
  writeXStringSet(bin, gsub(".pdf", ".FASTA", plot_name))
  pdf(plot_name, h = 15)
  plot(bionj(sdz))
  dev.off()
}

for (bin_name in names(which(lapply(bin_seqs, length) > 80))){
  print(bin_name)
  bin <- bin_seqs[[bin_name]]
  plot_bin(bin)
}

