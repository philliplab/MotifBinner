


### Arguments


```r
x <- pb_dat
x$pb_out <- NULL
x$seqs <- NULL
print(x)
```

```
## $classification_technique
## [1] "absolute"
## 
## $classification_params
## $classification_params$threshold
## [1] 0.01333333
## 
## $classification_params$start_threshold
## [1] 0.02333333
## 
## $classification_params$max_sequences
## [1] 100
## 
## 
## $alignment_technique
## [1] "muscle"
## 
## $alignment_params
## list()
## 
## $consensus_technique
## [1] "mostConsensusString"
## 
## $consensus_params
## list()
## 
## $remove_gaps
## [1] TRUE
```

### Compute metrics


```r
metrics <- NULL
i <- 1
for (i in seq_along(pb_out)){
  print(i)
  input_size <- length(pb_out[[i]]$src) + length(pb_out[[i]]$out)
  output_size <- length(pb_out[[i]]$src)
  out <- pb_out[[i]]$out
  src <- pb_out[[i]]$src
  dmat <- pb_out[[i]]$dmat
  min_dist <- min(as.dist(dmat))
  in_max_dist <- max(pb_out[[i]]$dmat)
  out_dmat <- dmat[match(names(src), colnames(dmat)), 
                   match(names(src), colnames(dmat))]
  out_max_dist <- max(out_dmat)
  min_out_dist <- min(as.dist(out_dmat))
  if (is.null(pb_out[[i]]$alignment)){
    gaps <- NA
    pos_no_mismatch <- NA
    pos_mismatch <- NA
    total_mismatches <- NA
    non_acgt <- NA
  } else {
    conMat <- consensusMatrix(pb_out[[i]]$alignment)
    gaps <- sum(conMat[-(1:4),])
    mismatches <- apply(conMat, 2, sum) - apply(conMat, 2, max)
    pos_no_mismatch <- sum(mismatches==0)
    pos_mismatch <- sum(mismatches!=0)
    total_mismatches <- sum(mismatches)
    non_acgt <- sum(consensusMatrix(DNAStringSet(pb_out[[i]]$consensus[[1]]))[-(1:4),])
  }
  metrics <- rbind(metrics,
    data.frame(input_size = input_size,
               output_size = output_size,
               min_dist = min_dist,
               in_max_dist = in_max_dist,
               out_max_dist = out_max_dist,
               min_out_dist = min_out_dist,
               gaps = gaps,
               pos_no_mismatch = pos_no_mismatch,
               pos_mismatch = pos_mismatch,
               total_mismatches = total_mismatches,
               non_acgt))
}
```

```
## [1] 1
## [1] 2
## [1] 3
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 4
## [1] 5
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 6
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 7
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 8
## [1] 9
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 10
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 11
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 12
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 13
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 14
## [1] 15
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 16
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 17
## [1] 18
## [1] 19
## [1] 20
## [1] 21
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 22
## [1] 23
## [1] 24
## [1] 25
## [1] 26
## [1] 27
## [1] 28
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 29
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 30
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 31
## [1] 32
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 33
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 34
## [1] 35
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 36
## [1] 37
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 38
```

```
## Warning in max(out_dmat): no non-missing arguments to max; returning -Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 39
## [1] 40
## [1] 41
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 42
## [1] 43
## [1] 44
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 45
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 46
## [1] 47
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 48
## [1] 49
## [1] 50
## [1] 51
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 52
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 53
## [1] 54
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 55
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 56
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 57
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 58
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 59
## [1] 60
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 61
## [1] 62
## [1] 63
## [1] 64
## [1] 65
## [1] 66
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 67
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 68
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 69
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 70
## [1] 71
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 72
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 73
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 74
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 75
## [1] 76
## [1] 77
## [1] 78
## [1] 79
## [1] 80
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 81
## [1] 82
## [1] 83
## [1] 84
## [1] 85
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 86
## [1] 87
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 88
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 89
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 90
## [1] 91
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 92
## [1] 93
## [1] 94
## [1] 95
## [1] 96
## [1] 97
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 98
## [1] 99
## [1] 100
## [1] 101
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 102
## [1] 103
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 104
## [1] 105
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 106
## [1] 107
## [1] 108
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 109
## [1] 110
## [1] 111
## [1] 112
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 113
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 114
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 115
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 116
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 117
## [1] 118
## [1] 119
## [1] 120
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 121
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 122
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 123
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 124
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 125
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 126
## [1] 127
## [1] 128
## [1] 129
## [1] 130
## [1] 131
## [1] 132
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 133
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 134
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 135
## [1] 136
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 137
## [1] 138
## [1] 139
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 140
## [1] 141
## [1] 142
## [1] 143
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 144
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 145
## [1] 146
## [1] 147
## [1] 148
## [1] 149
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 150
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 151
## [1] 152
## [1] 153
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 154
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 155
## [1] 156
## [1] 157
## [1] 158
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 159
## [1] 160
## [1] 161
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 162
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 163
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 164
## [1] 165
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 166
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 167
## [1] 168
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 169
## [1] 170
## [1] 171
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 172
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 173
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 174
## [1] 175
## [1] 176
## [1] 177
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 178
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 179
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 180
## [1] 181
## [1] 182
## [1] 183
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 184
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 185
## [1] 186
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 187
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 188
## [1] 189
## [1] 190
## [1] 191
## [1] 192
## [1] 193
## [1] 194
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 195
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 196
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 197
## [1] 198
## [1] 199
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 200
## [1] 201
## [1] 202
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 203
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 204
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 205
## [1] 206
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 207
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 208
## [1] 209
## [1] 210
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 211
## [1] 212
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 213
## [1] 214
## [1] 215
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 216
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 217
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 218
## [1] 219
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 220
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 221
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 222
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 223
## [1] 224
## [1] 225
## [1] 226
## [1] 227
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 228
## [1] 229
## [1] 230
## [1] 231
## [1] 232
## [1] 233
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 234
## [1] 235
## [1] 236
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 237
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 238
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 239
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 240
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 241
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 242
## [1] 243
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 244
## [1] 245
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 246
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 247
## [1] 248
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 249
## [1] 250
## [1] 251
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 252
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 253
## [1] 254
## [1] 255
## [1] 256
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 257
## [1] 258
## [1] 259
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 260
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 261
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 262
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 263
## [1] 264
## [1] 265
## [1] 266
## [1] 267
## [1] 268
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 269
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 270
## [1] 271
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 272
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 273
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 274
## [1] 275
## [1] 276
## [1] 277
## [1] 278
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 279
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 280
## [1] 281
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 282
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 283
## [1] 284
## [1] 285
## [1] 286
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 287
## [1] 288
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 289
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 290
## [1] 291
## [1] 292
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 293
## [1] 294
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 295
## [1] 296
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 297
## [1] 298
## [1] 299
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 300
## [1] 301
## [1] 302
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 303
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 304
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 305
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 306
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 307
## [1] 308
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 309
```

```
## Warning in max(out_dmat): no non-missing arguments to max; returning -Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 310
## [1] 311
## [1] 312
## [1] 313
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 314
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 315
## [1] 316
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 317
## [1] 318
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 319
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 320
## [1] 321
## [1] 322
## [1] 323
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 324
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 325
## [1] 326
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 327
## [1] 328
## [1] 329
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 330
## [1] 331
## [1] 332
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 333
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 334
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 335
## [1] 336
## [1] 337
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 338
## [1] 339
## [1] 340
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 341
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 342
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 343
## [1] 344
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 345
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 346
## [1] 347
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 348
## [1] 349
## [1] 350
## [1] 351
## [1] 352
## [1] 353
## [1] 354
## [1] 355
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 356
## [1] 357
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 358
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 359
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 360
## [1] 361
## [1] 362
## [1] 363
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 364
## [1] 365
## [1] 366
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 367
## [1] 368
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 369
## [1] 370
## [1] 371
## [1] 372
## [1] 373
## [1] 374
## [1] 375
## [1] 376
## [1] 377
## [1] 378
## [1] 379
## [1] 380
## [1] 381
## [1] 382
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 383
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 384
## [1] 385
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 386
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 387
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 388
## [1] 389
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 390
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 391
## [1] 392
## [1] 393
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 394
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 395
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 396
## [1] 397
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 398
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 399
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 400
## [1] 401
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 402
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 403
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 404
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 405
## [1] 406
## [1] 407
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 408
## [1] 409
## [1] 410
## [1] 411
## [1] 412
## [1] 413
## [1] 414
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 415
## [1] 416
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 417
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 418
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 419
## [1] 420
## [1] 421
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 422
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 423
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 424
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 425
## [1] 426
## [1] 427
## [1] 428
## [1] 429
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 430
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 431
## [1] 432
## [1] 433
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 434
## [1] 435
## [1] 436
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 437
## [1] 438
## [1] 439
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 440
## [1] 441
## [1] 442
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 443
## [1] 444
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 445
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 446
## [1] 447
## [1] 448
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 449
## [1] 450
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 451
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 452
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 453
## [1] 454
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 455
## [1] 456
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 457
## [1] 458
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 459
## [1] 460
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 461
## [1] 462
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 463
## [1] 464
## [1] 465
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 466
## [1] 467
## [1] 468
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 469
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 470
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 471
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 472
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 473
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 474
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 475
## [1] 476
## [1] 477
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 478
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 479
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 480
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 481
## [1] 482
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 483
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 484
## [1] 485
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 486
## [1] 487
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 488
## [1] 489
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 490
## [1] 491
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 492
## [1] 493
## [1] 494
## [1] 495
## [1] 496
## [1] 497
## [1] 498
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 499
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 500
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 501
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 502
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 503
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 504
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 505
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 506
```

```
## Warning in max(out_dmat): no non-missing arguments to max; returning -Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 507
## [1] 508
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 509
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 510
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 511
## [1] 512
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 513
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 514
## [1] 515
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 516
## [1] 517
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 518
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 519
## [1] 520
## [1] 521
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 522
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 523
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 524
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 525
## [1] 526
## [1] 527
## [1] 528
## [1] 529
## [1] 530
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 531
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 532
## [1] 533
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 534
## [1] 535
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 536
## [1] 537
## [1] 538
## [1] 539
## [1] 540
## [1] 541
## [1] 542
## [1] 543
## [1] 544
## [1] 545
## [1] 546
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 547
## [1] 548
## [1] 549
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 550
## [1] 551
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 552
## [1] 553
## [1] 554
## [1] 555
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 556
## [1] 557
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 558
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 559
## [1] 560
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 561
## [1] 562
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 563
## [1] 564
## [1] 565
## [1] 566
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 567
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 568
## [1] 569
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 570
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 571
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 572
## [1] 573
## [1] 574
## [1] 575
## [1] 576
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 577
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 578
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 579
## [1] 580
## [1] 581
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 582
## [1] 583
## [1] 584
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 585
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 586
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 587
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 588
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 589
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 590
## [1] 591
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 592
## [1] 593
## [1] 594
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 595
## [1] 596
## [1] 597
## [1] 598
## [1] 599
## [1] 600
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 601
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 602
## [1] 603
## [1] 604
## [1] 605
## [1] 606
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 607
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 608
## [1] 609
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 610
## [1] 611
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 612
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 613
## [1] 614
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 615
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 616
## [1] 617
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 618
## [1] 619
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 620
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 621
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 622
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 623
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 624
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 625
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 626
## [1] 627
## [1] 628
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 629
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 630
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 631
## [1] 632
## [1] 633
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 634
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 635
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 636
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 637
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 638
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 639
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 640
## [1] 641
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 642
## [1] 643
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 644
## [1] 645
## [1] 646
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 647
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 648
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 649
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 650
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 651
## [1] 652
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 653
## [1] 654
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 655
## [1] 656
## [1] 657
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 658
## [1] 659
## [1] 660
## [1] 661
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 662
## [1] 663
## [1] 664
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 665
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 666
## [1] 667
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 668
## [1] 669
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 670
## [1] 671
## [1] 672
## [1] 673
## [1] 674
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 675
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 676
## [1] 677
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 678
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 679
## [1] 680
## [1] 681
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 682
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 683
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 684
## [1] 685
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 686
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 687
## [1] 688
## [1] 689
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 690
## [1] 691
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 692
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 693
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 694
```

```
## Warning in max(out_dmat): no non-missing arguments to max; returning -Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 695
## [1] 696
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 697
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 698
## [1] 699
## [1] 700
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 701
## [1] 702
## [1] 703
## [1] 704
## [1] 705
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 706
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 707
## [1] 708
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 709
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 710
## [1] 711
## [1] 712
## [1] 713
## [1] 714
```

```
## Warning in max(out_dmat): no non-missing arguments to max; returning -Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 715
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 716
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 717
## [1] 718
## [1] 719
## [1] 720
## [1] 721
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 722
## [1] 723
## [1] 724
## [1] 725
## [1] 726
## [1] 727
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 728
## [1] 729
## [1] 730
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 731
## [1] 732
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 733
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 734
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 735
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 736
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 737
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 738
## [1] 739
## [1] 740
## [1] 741
## [1] 742
## [1] 743
## [1] 744
## [1] 745
## [1] 746
## [1] 747
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 748
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 749
## [1] 750
## [1] 751
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 752
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 753
## [1] 754
## [1] 755
## [1] 756
## [1] 757
## [1] 758
## [1] 759
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 760
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 761
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 762
## [1] 763
## [1] 764
## [1] 765
## [1] 766
## [1] 767
## [1] 768
## [1] 769
## [1] 770
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 771
## [1] 772
## [1] 773
## [1] 774
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 775
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 776
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 777
## [1] 778
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 779
## [1] 780
## [1] 781
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 782
## [1] 783
## [1] 784
## [1] 785
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 786
## [1] 787
## [1] 788
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 789
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 790
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 791
## [1] 792
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 793
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 794
## [1] 795
## [1] 796
## [1] 797
## [1] 798
## [1] 799
## [1] 800
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 801
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 802
## [1] 803
## [1] 804
## [1] 805
## [1] 806
## [1] 807
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 808
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 809
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 810
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 811
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 812
## [1] 813
## [1] 814
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 815
## [1] 816
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 817
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 818
## [1] 819
## [1] 820
## [1] 821
## [1] 822
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 823
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 824
## [1] 825
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 826
## [1] 827
## [1] 828
## [1] 829
## [1] 830
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 831
## [1] 832
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 833
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 834
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 835
```

```
## Warning in max(out_dmat): no non-missing arguments to max; returning -Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 836
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 837
## [1] 838
## [1] 839
## [1] 840
## [1] 841
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 842
```

```
## Warning in max(out_dmat): no non-missing arguments to max; returning -Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 843
## [1] 844
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 845
## [1] 846
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 847
## [1] 848
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 849
## [1] 850
## [1] 851
## [1] 852
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 853
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 854
## [1] 855
## [1] 856
## [1] 857
## [1] 858
## [1] 859
## [1] 860
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 861
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 862
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 863
## [1] 864
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 865
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 866
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 867
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 868
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 869
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 870
## [1] 871
## [1] 872
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 873
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 874
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 875
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 876
## [1] 877
## [1] 878
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 879
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 880
## [1] 881
## [1] 882
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 883
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 884
## [1] 885
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 886
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 887
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 888
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 889
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 890
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 891
## [1] 892
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 893
## [1] 894
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 895
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 896
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 897
## [1] 898
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 899
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 900
## [1] 901
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 902
## [1] 903
## [1] 904
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 905
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 906
## [1] 907
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 908
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 909
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 910
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 911
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 912
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 913
## [1] 914
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 915
## [1] 916
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 917
## [1] 918
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 919
## [1] 920
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 921
## [1] 922
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 923
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 924
## [1] 925
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 926
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 927
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 928
## [1] 929
## [1] 930
## [1] 931
## [1] 932
## [1] 933
## [1] 934
## [1] 935
## [1] 936
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 937
## [1] 938
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 939
## [1] 940
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 941
## [1] 942
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 943
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 944
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 945
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 946
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 947
## [1] 948
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 949
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 950
```

```
## Warning in max(out_dmat): no non-missing arguments to max; returning -Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 951
## [1] 952
## [1] 953
## [1] 954
## [1] 955
## [1] 956
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 957
## [1] 958
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 959
## [1] 960
## [1] 961
## [1] 962
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 963
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 964
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 965
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 966
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 967
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 968
## [1] 969
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 970
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 971
## [1] 972
## [1] 973
## [1] 974
## [1] 975
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 976
## [1] 977
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 978
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 979
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 980
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 981
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 982
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 983
## [1] 984
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 985
## [1] 986
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 987
## [1] 988
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 989
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 990
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 991
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 992
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 993
```

```
## Warning in min(as.dist(dmat)): no non-missing arguments to min; returning
## Inf
```

```
## Warning in min(as.dist(out_dmat)): no non-missing arguments to min;
## returning Inf
```

```
## [1] 994
## [1] 995
## [1] 996
## [1] 997
## [1] 998
```

### Effects of outlier removal

#### Before and after bin sizes


```r
hist(metrics$input_size)
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png) 

```r
hist(metrics$output_size)
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-2.png) 

#### Before and after bin distances


```r
hist(metrics$in_max_dist)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png) 

```r
hist(metrics$out_max_dist)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-2.png) 

```r
hist(metrics$min_out_dist)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-3.png) 

```r
hist(metrics$in_max_dist - metrics$out_max_dist)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-4.png) 

#### Where are the small bins coming from?


```r
plot(metrics$output_size ~ metrics$input_size)
abline(a=0, b=1)
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png) 


```r
small_bins <- subset(metrics, output_size < 6)
plot(jitter(small_bins$output_size) ~ jitter(small_bins$input_size))
abline(a=0, b=1)
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png) 

### Alignment step

#### Number of gaps inserted


```r
hist(metrics$gaps)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png) 

```r
kable(data.frame(num_gaps = as.numeric(names(table(metrics$gaps))),
                 num_alignments = as.numeric(table(metrics$gaps))))
```



| num_gaps| num_alignments|
|--------:|--------------:|
|        0|            359|
|        1|             76|
|        2|             26|
|        3|              8|
|        4|              5|
|        5|              2|
|        6|              2|
|        7|              2|
|        8|              1|
|        9|              2|
|       15|              2|
|       17|              1|
|       18|              1|
|       19|              1|
|       20|              1|
|       21|              1|
|       24|              1|
|       25|              1|
|       26|              1|
|       27|              1|
|       28|              1|
|       32|              1|
|       35|              2|
|       38|              1|
|       42|              1|
|       44|              1|
|       45|              2|
|       47|              1|
|       51|              2|
|       54|              1|
|      118|              1|

#### Number positions with mismatches in the alignment


```r
hist(metrics$pos_no_mismatch)
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png) 

```r
hist(metrics$pos_mismatch)
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-2.png) 

#### Total number of mismatches


```r
hist(metrics$total_mismatches)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-1.png) 

#### NON-ACGT characters in consensus


```r
hist(metrics$non_acgt)
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-1.png) 

```r
kable(data.frame(num_non_acgt = as.numeric(names(table(metrics$non_acgt))),
                 num_con_seq = as.numeric(table(metrics$non_acgt))))
```



| num_non_acgt| num_con_seq|
|------------:|-----------:|
|            0|         326|
|            1|          74|
|            2|          42|
|            3|          19|
|            4|          17|
|            5|          14|
|            6|           9|
|            8|           2|
|            9|           1|
|           10|           1|
|           11|           1|
|           13|           2|
