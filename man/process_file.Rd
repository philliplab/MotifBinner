% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/batch_utilities.R
\name{process_file}
\alias{process_file}
\title{Processes a file into consensus bins}
\usage{
process_file(file_name, output, prefix = "CCAGCTGGTTATGCGATTCTMARGTG",
  suffix = "CTGAGCGTGTGGCAAGGCCC", motif_length = 9,
  max.mismatch_start = 0, max.mismatch = 5, threshold = 8/600,
  start_threshold = 8/600, max_sequences = 400, remove_gaps = TRUE,
  strip_uids = TRUE, n_bins_to_process = 0, verbose = TRUE,
  prefix_for_names = "")
}
\arguments{
\item{file_name}{The file name}

\item{prefix}{The prefix that is used to identify the motif}

\item{suffix}{The suffix that is used to identify the motif}

\item{motif_length}{The length of the motif that forms the pid.}

\item{max.mismatch}{The maximum number of mismatches to allow when searching
for the pid}

\item{threshold}{Outlier sequences are removed from the bin until the
maximum distance between any two sequences drops below this threshold.}

\item{start_threshold}{Only start the classification if the maximum
distance between and two sequences in the bin is greater than this.}

\item{max_sequences}{The maximum number of sequences to use for the
computation of the distance matrix. If more sequences than this is present,
then randomly select this many sequences and run the classification
algorithm on them. This is only to improve the computation speed.}

\item{remove_gaps}{If set to TRUE (the default, then gaps will be removed
from the consensus sequences)}

\item{strip_uids}{Remove the unique identifiers from the sequence. It is not
intelligent. The names will be split on '_' and the first and last pieces
will be kept.}

\item{n_bins_to_process}{The number of bins to process through the outlier
detection, alignment and consensus generation. If smaller than or equal to
0, all bins will be processed.}

\item{verbose}{Progress information will be provided if set to TRUE}

\item{prefix_for_names}{Add this bit of text to the front of each
sequence in the resulting consensus sequences.}

\item{output_dir}{The directory where the output must be stored}
}
\description{
Processes a file into consensus bins
}

