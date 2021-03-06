% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/classification_checker.R
\name{check_classification}
\alias{check_classification}
\title{Checks that the output from a classifier satisfies basic constraints}
\usage{
check_classification(bin, classified)
}
\arguments{
\item{bin}{The input bin that was classified}

\item{classied}{The result from a classifier}
}
\description{
Makes sure that the following holds:
\itemize{
 \item{Output from the classifier is a list with 'in' and 'out' elements}
 \item{Each of the in and out elements are DNAStringSets}
 \item{None of the sequences in $in is in $out}
 \item{None of the sequences in $out is in $in}
 \item{Check that all the input sequences are present in the output data at
 the right frequencies}
 \item{Check that no new sequences were introduced during the classification
 step}
}
}

