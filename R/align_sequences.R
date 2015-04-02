# Aligns sequences using external tools
# A convenient wrapper is provided for benchmarking purposes.

#' A wrapper for a number of alignment tools
#'
#' @param seqs The sequences to align
#' @param technique Which alignment program to use
#' @param params A list of parameters to pass to the alignment program
#' @export

align_sequences <- function(seqs, technique = 'muscle', params = list()){
  if (is.null(names(seqs))){
    names(seqs) <- paste('seq', 1:length(seqs), sep='_')
  }
  seqs <- as.DNAbin(seqs)
  aligned_seqs <- muscle_par(seqs)
  seqs <- DNAStringSet(toupper(apply(as.character(aligned_seqs), 1, paste0, collapse="")))
  return(seqs)
}

#' Copy of muscle function from ape package with tweaks to allow parallel
#' execution
#'
#' @param x an object of class ‘"DNAbin"’.
#' @param pw.gapopen gap opening and gap extension penalties used by Clustal
#' during pairwise alignments. 
#' @param pw.gapext gap opening and gap extension penalties used by Clustal
#' during pairwise alignments.
#' @param gapopen idem for global alignment.
#' @param gapext idem for global alignment.
#' @param exec a character string giving the name of the program, with its path
#' if necessary. ‘clustal’ tries to guess it depending on the operating system
#' (see details).
#' @param MoreArgs a character string giving additional options.
#' @param quiet a logical: the default is to not print on R's console the
#' messages from the external program.
#' @param original.ordering a logical specifying whether to return the aligned
#' sequences in the same order than in ‘x’.
#' @export

muscle_par <- function(x, exec = "muscle", MoreArgs = "", quiet = TRUE, 
                       original.ordering = TRUE){
  if (missing(x)){
    system(exec)
    return(invisible(NULL))
  }
  d <- tempdir()
  random_fname <- paste(sample(c(LETTERS, letters, 0:9), 20), collapse="")
  inf_name <- paste("input_muscle_", random_fname, ".fas", sep = "")
  outf_name <- paste("output_muscle_", random_fname, ".fas", sep = "")
  inf <- paste(d, inf_name, sep = "/")
  outf <- paste(d, outf_name, sep = "/")
  write.dna(x, inf, "fasta")
  opts <- paste("-in", inf, "-out", outf)
  if (quiet){
    opts <- paste(opts, "-quiet")
  }
  opts <- paste(opts, MoreArgs)
  system(paste(exec, opts))
  res <- read.dna(outf, "fasta")
  if (original.ordering){
    res <- res[labels(x), ]
  }
  return(res)
}

