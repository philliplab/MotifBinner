% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/align_sequences.R
\name{muscle_par}
\alias{muscle_par}
\title{Copy of muscle function from ape package with tweaks to allow parallel
execution}
\usage{
muscle_par(x, exec = "muscle", MoreArgs = "", quiet = TRUE,
  original.ordering = TRUE, keep_alignment_files = FALSE)
}
\arguments{
\item{x}{an object of class ‘"DNAbin"’.}

\item{exec}{a character string giving the name of the program, with its path
if necessary. ‘clustal’ tries to guess it depending on the operating system
(see details).}

\item{MoreArgs}{a character string giving additional options.}

\item{quiet}{a logical: the default is to not print on R's console the
messages from the external program.}

\item{original.ordering}{a logical specifying whether to return the aligned
sequences in the same order than in ‘x’.}

\item{keep_alignment_files}{Default is FALSE. If set to TRUE, then the files
produced by muscle will not be deleted.}

\item{pw.gapopen}{gap opening and gap extension penalties used by Clustal
during pairwise alignments.}

\item{pw.gapext}{gap opening and gap extension penalties used by Clustal
during pairwise alignments.}

\item{gapopen}{idem for global alignment.}

\item{gapext}{idem for global alignment.}
}
\description{
Copy of muscle function from ape package with tweaks to allow parallel
execution
}

