

## The unprocessed bins produced


```r
Sys.time()
```

```
## [1] "2015-03-19 13:25:20 SAST"
```

### The number of bins


```r
print(length(bin_seqs))
```

```
## [1] 998
```

### The sizes of the bins


```r
bin_sizes <- unlist(lapply(bin_seqs, length))
hist(bin_sizes)
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png) 

```r
kable(data.frame(bin_size = as.numeric(names(table(bin_sizes))),
                 num_bins = as.numeric(table(bin_sizes))))
```



| bin_size| num_bins|
|--------:|--------:|
|        1|      482|
|        2|       76|
|        3|       35|
|        4|       32|
|        5|       24|
|        6|       20|
|        7|        9|
|        8|       13|
|        9|        9|
|       10|       17|
|       11|       15|
|       12|       10|
|       13|        8|
|       14|        6|
|       15|        6|
|       16|       10|
|       17|        9|
|       18|        8|
|       19|       11|
|       20|        9|
|       21|        5|
|       22|        8|
|       23|        5|
|       24|        5|
|       25|        4|
|       26|        6|
|       27|        6|
|       28|        4|
|       29|        8|
|       30|        8|
|       31|        5|
|       32|        4|
|       33|        4|
|       34|        4|
|       35|        5|
|       36|        3|
|       37|        7|
|       38|        3|
|       39|        5|
|       40|        4|
|       41|        3|
|       43|        3|
|       44|        5|
|       45|        2|
|       46|        5|
|       47|        1|
|       49|        4|
|       51|        1|
|       52|        1|
|       53|        1|
|       54|        1|
|       55|        6|
|       56|        1|
|       57|        3|
|       58|        1|
|       59|        1|
|       60|        3|
|       61|        1|
|       62|        2|
|       63|        1|
|       65|        5|
|       66|        3|
|       67|        3|
|       68|        3|
|       70|        1|
|       73|        4|
|       75|        1|
|       76|        1|
|       77|        1|
|       80|        1|
|       81|        2|
|       82|        1|
|       85|        1|
|       86|        2|
|       87|        1|
|       94|        2|
|       96|        1|
|       97|        2|
|      100|        1|
|      101|        1|
|      103|        1|
|      113|        1|
|      119|        1|

### scatter and log log plot of the bin sizes


```r
loglog <- data.frame(bin_size = as.numeric(names(table(bin_sizes))),
                     num_bins = as.numeric(table(bin_sizes)))
plot(num_bins ~ bin_size, data = loglog)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png) 

```r
plot(log10(num_bins) ~ log10(bin_size), data = loglog)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-2.png) 

### Categories and count frequencies


```r
cats <- hist(bin_sizes, breaks = c(0,1,2,3,4,10,20,400,1000000), plot=FALSE)
counts <- cats$counts
names(counts) <- paste('(', cats$breaks[1:(length(cats$breaks)-1)], 
                  '-', cats$breaks[2:(length(cats$breaks))], ']',
                  sep='')
names(counts)[length(counts)] <- paste(cats$breaks[length(cats$breaks)-1], '+',
                                       sep = '')
barplot(counts, main = 'Histogram of Bin Sizes')
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png) 

```r
kable(data.frame(category = names(counts),
                 count = counts,
                 row.names = 1:length(counts)))
```



|category | count|
|:--------|-----:|
|(0-1]    |   482|
|(1-2]    |    76|
|(2-3]    |    35|
|(3-4]    |    32|
|(4-10]   |    92|
|(10-20]  |    92|
|(20-400] |   189|
|400+     |     0|
