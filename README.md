MotifBinner
===========

Development to start in Dec 2014

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
  return(motif)
}

bin_by_motif <- function(seq_data, ...){
  bins <- list() # TODO more efficient data structure however, much of this we
# will want to run in parallel, so the data structure must have things like
# lapply defined for it in packages like snowfall.
  for (the_seq in seq_dat){
    motif <- extract_motif(the_seq, ...)
    bins[[motif]] <- c(bins[[motif]], the_seq)
  }
}

# should contaminated sequences be removed before alignment?
# This allows a simple distance criteria to be used instead
#   of performing a complicated alignment due to poorly matched sequences?

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
