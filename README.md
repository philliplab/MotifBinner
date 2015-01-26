MotifBinner
===========

Bins reads from NGS sequencing machines by a motif. Motif are identified by
special affixes. Binned sequences are aligned and consensus sequences
extracted.

```r
extract_motif <- function(the_seq, prefix, suffix,
                          extra_transformation_options,
                          extra_location_options){
  the_seq <- substr(the_seq) # TODO can more information be know about the
    # location of the motif? should this be used in this function?
    # Issue number 7
  motif <- subseq(vmatchPattern(...)) # More complicated version of this
    # Need to think carefully about the data structure. Should we look at
    # IRanges::List?
  names(the_seq) <- motif # Currently passing motif info via the names
    # attribute. No original data from the FASTA headers are required
    # to be included in the final result.
  return(motif)
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
```
