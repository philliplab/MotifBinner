


## A Basic report on the unprocessed bins produced


```r
Sys.time()
```

```
## [1] "2015-03-19 13:43:46 SAST"
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
plot(log10(num_bins) ~ log10(bin_size), data = loglog)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-2.png) 

### Detailed investigations of the distances in the bins



#### Sequence lengths


```r
hist(comp_bins$min_seq_len)
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png) 

```r
hist(comp_bins$mean_seq_len)
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-2.png) 

```r
hist(comp_bins$max_seq_len)
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-3.png) 

#### Bins of size 2


```r
c_bins <- subset(comp_bins, bin_size == 2)
hist(c_bins$mean_dist)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png) 

```r
kable(data.frame(bin_dist = as.numeric(names(table(c_bins$mean_dist))),
                 num_bins = as.numeric(table(c_bins$mean_dist))))
```



| bin_dist| num_bins|
|--------:|--------:|
|        0|        3|
|        1|        6|
|        2|        7|
|        3|       15|
|        4|       12|
|        5|       13|
|        6|        6|
|        8|        2|
|        9|        1|
|       10|        1|
|       11|        1|
|       13|        2|
|       16|        1|
|       17|        1|
|       18|        3|
|       19|        1|
|       21|        1|

#### Bins of size 3


```r
c_bins <- subset(comp_bins, bin_size == 3)
hist(c_bins$min_dist)
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png) 

```r
kable(data.frame(min_dist = as.numeric(names(table(c_bins$min_dist))),
                 num_bins = as.numeric(table(c_bins$min_dist))))
```



| min_dist| num_bins|
|--------:|--------:|
|        0|        4|
|        1|        5|
|        2|        5|
|        3|        9|
|        4|        6|
|        5|        2|
|        6|        4|

```r
hist(c_bins$max_dist)
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-2.png) 

```r
kable(data.frame(max_dist = as.numeric(names(table(c_bins$max_dist))),
                 num_bins = as.numeric(table(c_bins$max_dist))))
```



| max_dist| num_bins|
|--------:|--------:|
|        1|        1|
|        2|        2|
|        3|        2|
|        4|        6|
|        5|        4|
|        6|        5|
|        7|        2|
|        8|        1|
|        9|        2|
|       12|        1|
|       14|        2|
|       15|        1|
|       16|        2|
|       18|        3|
|       19|        1|

```r
plot(c_bins$min_dist ~ c_bins$max_dist)
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-3.png) 

#### Bins of size 4 and above


```r
c_bins <- subset(comp_bins, bin_size > 3)
hist(c_bins$bin_size)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-1.png) 

```r
hist(c_bins$min_dist)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-2.png) 

```r
kable(data.frame(min_dist = as.numeric(names(table(c_bins$min_dist))),
                 num_bins = as.numeric(table(c_bins$min_dist))))
```



| min_dist| num_bins|
|--------:|--------:|
|        0|      276|
|        1|       69|
|        2|       37|
|        3|       15|
|        4|        5|
|        5|        2|
|        7|        1|

```r
hist(c_bins$mean_dist)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-3.png) 

```r
kable(data.frame(mean_dist = as.numeric(names(table(round(c_bins$mean_dist,0)))),
                 num_bins = as.numeric(table(round(c_bins$mean_dist,0)))))
```



| mean_dist| num_bins|
|---------:|--------:|
|         2|        2|
|         3|       27|
|         4|       56|
|         5|       59|
|         6|      107|
|         7|       71|
|         8|       37|
|         9|       16|
|        10|       19|
|        11|        5|
|        12|        1|
|        13|        3|
|        20|        1|
|        42|        1|

```r
hist(c_bins$max_dist)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-4.png) 

```r
kable(data.frame(max_dist = as.numeric(names(table(c_bins$max_dist))),
                 num_bins = as.numeric(table(c_bins$max_dist))))
```



| max_dist| num_bins|
|--------:|--------:|
|        3|        1|
|        4|        7|
|        5|       12|
|        6|       16|
|        7|       19|
|        8|       10|
|        9|       15|
|       10|        8|
|       11|        7|
|       12|        6|
|       13|        5|
|       14|        9|
|       15|       10|
|       16|        6|
|       17|       13|
|       18|       16|
|       19|       13|
|       20|       18|
|       21|       18|
|       22|       17|
|       23|       18|
|       24|       15|
|       25|       13|
|       26|       13|
|       27|       14|
|       28|        8|
|       29|       14|
|       30|        7|
|       31|        8|
|       32|        7|
|       33|        5|
|       34|        4|
|       35|        9|
|       36|        5|
|       37|        1|
|       38|        5|
|       39|        3|
|       40|        3|
|       41|        2|
|       42|        4|
|       43|        4|
|       44|        3|
|       48|        2|
|       49|        4|
|       51|        1|
|       53|        1|
|       54|        1|
|       55|        1|
|       59|        1|
|       62|        1|
|       80|        1|
|      154|        1|

```r
plot(c_bins$min_dist ~ c_bins$max_dist)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-5.png) 

```r
hist(c_bins$max_occ)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-6.png) 

```r
kable(data.frame(max_occ = as.numeric(names(table(c_bins$max_occ))),
                 num_bins = as.numeric(table(c_bins$max_occ))))
```



| max_occ| num_bins|
|-------:|--------:|
|       1|      129|
|       2|      124|
|       3|       85|
|       4|       42|
|       5|       14|
|       6|        6|
|       7|        1|
|       8|        1|
|       9|        2|
|      10|        1|
