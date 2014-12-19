MotifBinner
===========

Bins reads from NGS sequencing machines by a motif. Motif are identified by
special affixes. Binned sequences can be aligned and consensus sequences
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

# Now check for mislabeled seqences. This is done by checking for outliers in
# the distance matrix
# See Issue number 9

align_bins <- function(bins){
  aligned_bins <- list()
  for (motif in names(bins)){
    alignment <- fancy_aligner(bins[[motif]])
    aligned_bins[[motif]] <- alignment
  }
  return(aligned_bins)
}

# A function that computes distances based on a MSA and report on poorly
#   aligned bins

# A function that returns only the consensus sequences of each bin
```
