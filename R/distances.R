i <- 1
i <- 2

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

compare_dists <- function(d1, d2){
  d1 <- as.matrix(d1)
  d2 <- as.matrix(d2)
  comp_dist <- data.frame(d1=numeric(0),
                          d2=numeric(0))
  for (i in 1:(nrow(d1)-1)){
    new_rows <- data.frame(d1=numeric(0),
                           d2=numeric(0))
    for (j in i:nrow(d1)){
      new_rows <- rbind(new_rows, data.frame(d1 = d1[i,j],
                                             d2 = d2[i,j]))
    }
    comp_dist <- rbind(comp_dist, 
                       new_rows)
  }
  return(comp_dist)
}

easy_compare_dists <- function(x, anchors = 3){
  fast <- system.time(yo <- fast_stringDist(x, anchors))
  normal <- system.time(yp <- stringDist(x))
  dists <- compare_dists(yo, yp)
  rsq <- summary(lm(dists[,1] ~ dists[,2]))$adj.r.squared
  plot(dists, main = paste("n seq = ", length(x), 
                           "; fast = ", round(fast[1], 1), 
                           "\nnormal = ", round(normal[1], 1), 
                           "; Rsq = ", round(rsq, 4), 
                           sep = ""))
}
