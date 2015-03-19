

## Input File


```r
Sys.time()
```

```
## [1] "2015-03-19 12:51:11 SAST"
```

### A summary of the input sequences


```r
print(seq_dat)
```

```
##   A DNAStringSet instance of length 14575
##         width seq
##     [1]   352 TATGGGATCAAAGTCTAAAACCATGTGTA...TGCTTACCTAACTGAGCGTGTGGCAAGGC
##     [2]   352 TATGGGACGAAAGCCTCAAGCCATGTGTA...TGTTAGATGTTCTGAGCGTGTGGCAAGGC
##     [3]   361 TATGGGACGAAAGTCTAAAGCCATGTGTA...TGTATTGGGGCCTGAGCGTGTGGCAAGGC
##     [4]   419 AATGTCAGCACAGTACAATGTACACATGG...ACTCATGCAACCTGAGCGTGTGGCAAGGC
##     [5]   351 TATGGGACCAAAGTCTAAAGCCATGTGTA...GTGCTATCGAACTGAGCGTGTGGCAAGGC
##     ...   ... ...
## [14571]   352 TATGGGACGAAAGTCTAAAGCCATGTGTA...TGGTACACAAGCTGAGCGTGTGGCAAGGC
## [14572]   352 TATGGGACGAAAGCCTCAAGCTATGTGTA...TGCATTTAAAACTGAGCGTGTGGCAAGGC
## [14573]   352 TATGGGACGAAAGTCTCAAACCATGTGTA...TGGCTCAAGTACTGAGCGTGTGGCAAGGC
## [14574]   361 TATGGGATCAAAGTCTAAAACCATGTGTA...TGAATCCTCTACTGAGCGTGTGGCAAGGC
## [14575]   352 TATGGGACGAAAGTCTAAAGCCATGTGTA...TGTAATTGCGACCGAGCGTGTGGCAAGGC
```

### The letter frequencies

Counts of the letters in the input sequences.

```r
x <- apply(consensusMatrix(seq_dat), 1, sum)
x <- x[x != 0]
kable(data.frame(letter = names(x),
                 count = x,
                 row.names = 1:length(x)))
```



|letter |   count|
|:------|-------:|
|A      | 1889682|
|C      | 1028655|
|G      |  973999|
|T      | 1394116|
|N      |      12|

Normalized counts of the letters in the input sequences.

```r
x_freq <- round(x/sum(x), 4)
kable(data.frame(letter = names(x_freq),
                 count = x_freq,
                 row.names = 1:length(x_freq)))
```



|letter |  count|
|:------|------:|
|A      | 0.3575|
|C      | 0.1946|
|G      | 0.1842|
|T      | 0.2637|
|N      | 0.0000|


```r
barplot(x/sum(x),
        main = 'Normalized counts of the letters\nin the input sequences')
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png) 

### The sequence lengths


```r
hist(width(seq_dat),
     main = "Histogram of the lengths of the input sequences",
     xlab = "Sequence Length",
     ylab = "Number of Sequences")
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png) 
