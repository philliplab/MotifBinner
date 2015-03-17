

## A Basic report on the input file that was provided


```r
Sys.time()
```

```
## [1] "2015-03-17 16:12:40 SAST"
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


```r
x <- apply(consensusMatrix(seq_dat), 1, sum)
x <- x[x != 0]
print(x)
```

```
##       A       C       G       T       N 
## 1889682 1028655  973999 1394116      12
```

```r
print(round(x/sum(x), 4))
```

```
##      A      C      G      T      N 
## 0.3575 0.1946 0.1842 0.2637 0.0000
```

```r
barplot(x/sum(x))
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png) 

### The sequence lengths


```r
hist(width(seq_dat))
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png) 
