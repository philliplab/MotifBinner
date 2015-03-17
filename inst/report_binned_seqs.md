

## A Basic report on the bins produced


```r
Sys.time()
```

```
## [1] "2015-03-17 17:00:28 SAST"
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
table(bin_sizes)
```

```
## bin_sizes
##   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18 
## 482  76  35  32  24  20   9  13   9  17  15  10   8   6   6  10   9   8 
##  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35  36 
##  11   9   5   8   5   5   4   6   6   4   8   8   5   4   4   4   5   3 
##  37  38  39  40  41  43  44  45  46  47  49  51  52  53  54  55  56  57 
##   7   3   5   4   3   3   5   2   5   1   4   1   1   1   1   6   1   3 
##  58  59  60  61  62  63  65  66  67  68  70  73  75  76  77  80  81  82 
##   1   1   3   1   2   1   5   3   3   3   1   4   1   1   1   1   2   1 
##  85  86  87  94  96  97 100 101 103 113 119 
##   1   2   1   2   1   2   1   1   1   1   1
```

```r
hist(bin_sizes)
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png) 

### log log plot of the bin sizes


```r
loglog <- data.frame(bin_size = as.numeric(names(table(bin_sizes))),
                     num_bins = as.numeric(table(bin_sizes)))
plot(num_bins ~ bin_size, data = loglog)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png) 

```r
plot(log(num_bins) ~ bin_size, data = loglog)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-2.png) 

```r
plot(log(num_bins) ~ log(bin_size), data = loglog)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-3.png) 

### The sequence lengths


```r
hist(width(seq_dat))
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png) 
