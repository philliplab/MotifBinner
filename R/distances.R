#' Computes a distance matrix between sequences by using relative distances to
#' a few carefully selected sequences
#'
#' The first sequence in the input data is chosen and the pairwise distances
#' between all sequences and this sequence is computed. The sequence most distant
#' from the first sequence is then chosen as the next anchor. The process
#' continues by selecting the sequence furthest from all previous anchors as the
#' next anchor.
#'
#' Once all the distances to the chosen number of anchors have been computed, a
#' distance matrix is constructed using the distances between the sequences and
#' the anchors. For very large numbers of sequences (>500) this realizes a big
#' speedup.
#' @param x The input sequences
#' @param anchors The number of sequences to use as anchors
#' @export
fast_stringDist <- function(x, anchors = 3){
  stopifnot(anchors < length(unique(x))-5)
  distances <- matrix(0, ncol = anchors, nrow=length(x))
  row.names(distances) <- names(x)
  all_anchors <- NULL
  for (i in 1:anchors){
    if (i == 1){
      anchor <- x[1]
      all_anchors <- DNAStringSet(x[1])
    }
    for (j in 1:length(x)){
      distances[j,i] <- stringDist(DNAStringSet(c(anchor, x[j])))
    }
    if (i < anchors){
      total_dists <- apply(distances, 1, sum)
      anchor <- x[which(total_dists == max(total_dists))[1]]
      all_anchors <- DNAStringSet(c(all_anchors,
                                  anchor))
    }
  }
  return(dist(distances))
}

#' Compare two distance matrices.
#'
#' Given two distance matrices (or dist objects) for the same elements in the
#' same ordering, return a data.frame with 2 columns, and with each row giving the
#' 2 distances given to the current pairwise distance
#'
#' @param d1 The first distance matrix
#' @param d2 The second distance matrix
#' @export
compare_dists <- function(d1, d2){
  stopifnot(attr(d1, 'names') == attr(d2, 'names'))
  d1 <- as.matrix(d1)
  d2 <- as.matrix(d2)
  comp_dist <- data.frame(d1=numeric(0),
                          d2=numeric(0))
  for (i in 1:(nrow(d1)-1)){
    new_rows <- data.frame(d1=rep(0, nrow(d1)-i+1),
                           d2=rep(0, nrow(d1)-i+1))
    for (j in i:nrow(d1)){
      new_rows[j-i+1, 1] <- d1[i, j]
      new_rows[j-i+1, 2] <- d2[i, j]
    }
    comp_dist <- rbind(comp_dist, 
                       new_rows)
  }
  return(comp_dist)
}
