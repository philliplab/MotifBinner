

## Motif Searches


```r
Sys.time()
```

```
## [1] "2015-03-19 13:09:27 SAST"
```

### Input parameters for the searches

```r
print(prefix)
```

```
## [1] "CCAGCTGGTTATGCGATTCTMARGTG"
```

```r
print(suffix)
```

```
## [1] "CTGAGCGTGTGGCAAGGCCC"
```

```r
print(motif_length)
```

```
## [1] 9
```

### Matched vs Unmatched sequences

```r
n_matched <- length(motif_dat$matched_seq)
n_unmatched <- length(motif_dat$unmatched_seq)
n_total <- n_matched + n_unmatched
kable(data.frame(seq_matched_count = n_matched,
                 seq_matched_perc = 100*round(n_matched / n_total, 4),
                 seq_unmatched_count = n_unmatched,
                 seq_unmatched_perc = 100*round(n_unmatched / n_total, 4)))
```



| seq_matched_count| seq_matched_perc| seq_unmatched_count| seq_unmatched_perc|
|-----------------:|----------------:|-------------------:|------------------:|
|             11561|            79.32|                3014|              20.68|

### Matched Sequences

#### Summary of sequences


```r
print(motif_dat$matched_seq)
```

```
##   A DNAStringSet instance of length 11561
##         width seq                                      names               
##     [1]   299 TATGGGATCAAAGTCTAAA...CCTATACATTATTGTGCT CTTACCTAA
##     [2]   299 TATGGGACGAAAGCCTCAA...CCTATACATTATTGTGCT TTAGATGTT
##     [3]   308 TATGGGACGAAAGTCTAAA...CCTATACATTATTGTGCT TATTGGGGC
##     [4]   299 TATGGGACGAAAGCCTCAA...CCTATACATTATTGTGCT ATTGTGCCA
##     [5]   299 TATGGGATGAAAGCCTCAA...CCTATACATTATTGTGCT GTAGACTGT
##     ...   ... ...
## [11557]   299 TATGGGACGAAAGTCTAAA...CCTATACATTATTGTGCT GTACACAAG
## [11558]   299 TATGGGACGAAAGCCTCAA...CCTATACATTATTGTGCT CATTTAAAA
## [11559]   299 TATGGGACGAAAGTCTCAA...CCTATACATTATTGTGCT GCTCAAGTA
## [11560]   308 TATGGGATCAAAGTCTAAA...CCTATACATTATTGTGCT AATCCTCTA
## [11561]   299 TATGGGACGAAAGTCTAAA...CCTATACATTATGGTGCT TAATTGCGA
```

#### The letter frequencies

Counts of the letters in the matched sequences.

```r
x <- apply(consensusMatrix(motif_dat$matched_seq), 1, sum)
x <- x[x != 0]
kable(data.frame(letter = names(x),
                 count = x,
                 count_per_seq = round(x/length(motif_dat$matched_seq),2),
                 row.names = 1:length(x)))
```



|letter |   count| count_per_seq|
|:------|-------:|-------------:|
|A      | 1296178|        112.12|
|C      |  675110|         58.40|
|G      |  536220|         46.38|
|T      |  978310|         84.62|
|+      |       2|          0.00|

Normalized counts of the letters in the matched sequences.

```r
x_freq <- round(x/sum(x), 4)
kable(data.frame(letter = names(x_freq),
                 count = x_freq,
                 row.names = 1:length(x_freq)))
```



|letter |  count|
|:------|------:|
|A      | 0.3718|
|C      | 0.1937|
|G      | 0.1538|
|T      | 0.2807|
|+      | 0.0000|


```r
barplot(x/sum(x),
        main = 'Normalized counts of the letters\nin the matched sequences')
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png) 

#### The sequence lengths


```r
hist(width(motif_dat$matched_seq),
     main = "Histogram of the lengths of the matched sequences",
     xlab = "Sequence Length",
     ylab = "Number of Sequences")
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png) 

### Unmatched Sequences
#### Summary of sequences


```r
print(motif_dat$unmatched_seq)
```

```
##   A DNAStringSet instance of length 3014
##        width seq                                       names               
##    [1]   419 AATGTCAGCACAGTACAAT...CCTGAGCGTGTGGCAAGGC seq_4
##    [2]   351 TATGGGACCAAAGTCTAAA...ACTGAGCGTGTGGCAAGGC seq_5
##    [3]   353 TATGGGACGAAAGCCTCAA...CCTGAGCGTGTGGCAAGGC seq_10
##    [4]   351 TATGGGACGAAAGCCTCAA...CCTGAGCGTGTGGCAAGGC seq_15
##    [5]   343 TATGGGATGAAAGCCTAAA...CCTGAGCGTGTGGCAAGGC seq_17
##    ...   ... ...
## [3010]   419 AATGTCAGTACAGTACAAT...ACTGAGCGTGTGGCAAGGC seq_14561
## [3011]   419 AATGTCAGTACAGTACAAT...ACTGAGCGTGTGGCAAGGC seq_14562
## [3012]   418 AATGTCAGTACAGTACAAT...CCTGAGCGTGTGGCAAGGC seq_14563
## [3013]   419 AATGTAAGCACAGTACAAT...CCTGAGCGTGTGGCAAGGC seq_14564
## [3014]   419 AATGTCAGTACAGTACAAT...GCTGAGCTTGTGGCAAGGC seq_14566
```

#### The letter frequencies

Counts of the letters in the unmatched sequences.

```r
x <- apply(consensusMatrix(motif_dat$unmatched_seq), 1, sum)
x <- x[x != 0]
kable(data.frame(letter = names(x),
                 count = x,
                 count_per_seq = round(x/length(motif_dat$unmatched_seq),2),
                 row.names = 1:length(x)))
```



|letter |  count| count_per_seq|
|:------|------:|-------------:|
|A      | 466241|        154.69|
|C      | 213997|         71.00|
|G      | 241376|         80.08|
|T      | 266341|         88.37|
|+      |      8|          0.00|

Normalized counts of the letters in the unmatched sequences.

```r
x_freq <- round(x/sum(x), 4)
kable(data.frame(letter = names(x_freq),
                 count = x_freq,
                 row.names = 1:length(x_freq)))
```



|letter |  count|
|:------|------:|
|A      | 0.3925|
|C      | 0.1801|
|G      | 0.2032|
|T      | 0.2242|
|+      | 0.0000|


```r
barplot(x/sum(x),
        main = 'Normalized counts of the letters\nin the unmatched sequences')
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13-1.png) 

#### The sequence lengths


```r
hist(width(motif_dat$unmatched_seq),
     main = "Histogram of the lengths of the unmatched sequences",
     xlab = "Sequence Length",
     ylab = "Number of Sequences")
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-1.png) 
