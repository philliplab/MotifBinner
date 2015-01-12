
```r
library(MotifBinner)
```

```
## Loading required package: ShortRead
## Loading required package: BiocGenerics
## Loading required package: parallel
## 
## Attaching package: 'BiocGenerics'
## 
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
## 
## The following object is masked from 'package:stats':
## 
##     xtabs
## 
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, as.vector, cbind,
##     colnames, do.call, duplicated, eval, evalq, Filter, Find, get,
##     intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rep.int, rownames, sapply, setdiff, sort,
##     table, tapply, union, unique, unlist
## 
## Loading required package: BiocParallel
## Loading required package: Biostrings
## Loading required package: IRanges
## Loading required package: XVector
## Loading required package: Rsamtools
## Loading required package: GenomicRanges
## Loading required package: GenomeInfoDb
## Loading required package: GenomicAlignments
## Loading required package: BSgenome
## Loading required package: ape
## 
## Attaching package: 'ape'
## 
## The following object is masked from 'package:ShortRead':
## 
##     zoom
## 
## Loading required package: shiny
## Loading required package: testthat
```

# Investigation of various approach to treating mislabel detection

## Overview

Given a bin of sequences that are theoretically from only a single virus
molecule, remove all of those that are actually from a different molecule.

This is accomplished by looking at the distances between each of the sequences
and removing the outliers.

The difficult parts are deciding how to compute the distances efficiently and
setting the thresholds that start and stop the process.

This process is further complicated by these concerns:
1) bin sizes varies from 1 to several hundred
2) We have to deal with indels
3) This process must be completely automated
4) This process must be very computationally efficient

## Testing strategy

We need to test the system using bins with known answers. 

### Test Data

As always, the design data structure for the test data is important. Use a
structure like this:


```r
list('test1' = list('in' = DNAStringSet(...),
                    'out' = DNAStringSet(...)),
     'test2' = ...)
```

The data is available in the package:

```r
test_dat <- get_mislabel_test_data()
```

Using this data for a 'bin' by putting the 'src' and 'out' data together and check
that only the 'out' data is removed and all the 'src' data is kept.

### Basic Code for running tests

Basic code to run the tests is available from the package:


```r
score_classification()
```

### Metrics of interest

A number of metrics must be considered when looking at the accuracy of
classification.

Keep it basic. 

Consider these standard classification metrics:

#### Sensitivity:
Number of true 'in' classifications / (Total size of true 'in' population)

#### Specificity:
Number of true 'out' classifications / (Total size of the 'out' population)

Maximize both simultaneously

Now, in addition to these two, there are other metrics of interest that we can
derive based on our knowledge of the system.

#### Maximum distance in the final dataset: 
We know that the only source of errors
should be the sequencing process. We have access to data about the error rates
of the sequencing process. We can use this to make a statement like:

The sequencing process is accurate to such a degree that no two reads of the
same molecule should differ by more than one base per 100 bases. Build a metric
around this information.

#### Speed

The time it took to classify the reads in the bin.

## The different strategies

Design and implement a who set of strategies and then benchmark them to find
the best ones.

### Remove None

This strategy keeps all the data

Use the random strategy but set the parameter 'n' to 0. So that 0 percent of
the data will be randomly removed.


```r
kable(score_all_classifications(test_dat, 'random', params = list(n=0)))
```



|name   | sn| sp| max_dist| time_taken|
|:------|--:|--:|--------:|----------:|
|test1  |  1|  0|       29|      0.031|
|test2  |  1|  0|        8|      0.015|
|test3  |  1|  0|       13|      0.014|
|test4  |  1|  0|       23|      0.012|
|test5  |  1|  1|        3|      0.013|
|test6  |  1|  0|       13|      0.014|
|test7  |  1|  0|       26|      0.013|
|test8  |  1|  0|       28|      0.017|
|test9  |  1|  0|       17|      0.012|
|test10 |  1|  0|       37|      0.018|

### Remove Random

This strategy removes x% of the sample at random

### Information Variance Balance

This strategy will keep on removing the most outlying sequence as long as is
leads to a percentage reduction in variance that is x times larger than the
percentage of information that was discarded

### The silver bullet

To be devised

