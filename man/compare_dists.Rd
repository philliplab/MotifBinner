% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/distances.R
\name{compare_dists}
\alias{compare_dists}
\title{Compare two distance matrices.}
\usage{
compare_dists(d1, d2)
}
\arguments{
\item{d1}{The first distance matrix}

\item{d2}{The second distance matrix}
}
\description{
Given two distance matrices (or dist objects) for the same elements in the
same ordering, return a data.frame with 2 columns, and with each row giving the
2 distances given to the current pairwise distance
}

