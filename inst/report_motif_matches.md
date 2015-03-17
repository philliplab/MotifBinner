

## A Basic report on search for PIDs in the sequences


```r
Sys.time()
```

```
## [1] "2015-03-17 16:28:59 SAST"
```

### Input parameters

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


```r
x <- apply(consensusMatrix(motif_dat$matched_seq), 1, sum)
x <- x[x != 0]
print(x)
```

```
##       A       C       G       T       + 
## 1296178  675110  536220  978310       2
```

```r
print(round(x/sum(x), 4))
```

```
##      A      C      G      T      + 
## 0.3718 0.1937 0.1538 0.2807 0.0000
```

```r
barplot(x/sum(x))
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png) 

#### The sequence lengths


```r
hist(width(motif_dat$matched_seq))
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png) 

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


```r
x <- apply(consensusMatrix(motif_dat$unmatched_seq), 1, sum)
x <- x[x != 0]
print(x)
```

```
##      A      C      G      T      + 
## 466241 213997 241376 266341      8
```

```r
print(round(x/sum(x), 4))
```

```
##      A      C      G      T      + 
## 0.3925 0.1801 0.2032 0.2242 0.0000
```

```r
barplot(x/sum(x))
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png) 

#### The sequence lengths


```r
hist(width(motif_dat$unmatched_seq))
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png) 
