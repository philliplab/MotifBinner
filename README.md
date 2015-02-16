MotifBinner
===========

Bins reads from NGS sequencing machines by a motif. Motif are identified by
special affixes. Binned sequences are aligned and consensus sequences
extracted.

The code below is greatly  simplified, but should serve as a high level map to
allow one to get an idea for the flow of the code.

```r
extract_motif <- function(the_seq, prefix, suffix,
                          extra_transformation_options,
                          extra_location_options){
  motif <- subseq(vmatchPattern(...)) # quite complex to actually implement 
  names(the_seq) <- motif 
  named_seq = list(matched_seq = the_seq[!is.na(motif)],
                   unmatched_seq = the_seq[is.na(motif)])
  return(named_seq)
}

bin_by_name <- function(seq_dat){
  bins <- list()
  for (motif in unique(names(seq_dat))){
    bins[[motif]] <- seq_dat[names(seq_dat) == motif]
  }
  return(bins)
}

process_bin <- function(bin){
  x <- classify_bin()
  x <- align_sequences(x$src)
  x <- construct_consensus(x)
  return(x)
}

classify_bin <- function(bin, technique, params){
  # The bin will be split into a list with 2 element 'src' and 'out'.
  # The elements of 'out' are to be discarded as mislabelled sequences
  # Several techniques are implemented
}

align_sequences <- function(seqs, technique, params){
  # just a wrapper for various different multiple alignment programs
  # Returns a DNAStringSet of aligned sequences
}

construct_consensus <- function(seqs, technique , params){
  # just a wrapper for various consensus string functions
  # Given a MSA, return the consensus sequence
  # Would be great to use qualities in this process
}

file_to_consensus <- function(file_name, ...){
  binned <- bin_file(file_name, ...)
  binned <- lapply(binned, process_bin, ...)
  return(binned)
}

bin_file <- function(file_name, ...){
  seq_dat <- read_file()
  seq_dat <- extract_motifs(seq_dat, ...)
  binned <- bin_by_name(seq_dat$matched_seq, ...)
  return(binned)
}
```

## Benchmarking System

In order to benchmark the system, simulate reads from a sequence. Possibly
contaminate the reads. Score how accurately the input sequence was
reconstructed.

Generating reads are achieved through two functions:
```r
mut_seq <- function(the_seq, n_reads, error_profile){
  new_seq <- rep(the_seq, n_reads)
  for (i in 1:nchar(the_seq)){
    for (j in 1:n_reas){
      x <- runif()
      if (x > error_profile[i, the_seq[i]]){
        mutate_at_position(new_seq, i, j)
      }
    }
  }
  return(new_seq) 
}

gen_error_profile <- function(...){
  # Generates an error profile indicating how likely it is to mutate from base
  # x to base y at position i based on a variety of inputs
}
```

Next various scenarios and setups must be constructed and ran so that
performance can be assessed. See the consensus_investigation.Rmd document in
the inst folder.
