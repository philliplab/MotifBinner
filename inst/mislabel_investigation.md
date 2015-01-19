

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

#### Last words about the metrics

The euclidean distance from the sensitivity and specificity from (1, 1) will be
reports as 'combo' for each test. If the average of this column for all test
datasets is 0, then we have a perfect classifier. This is the metric of
interest to be on the lookout for.

## The different strategies

Design and implement a who set of strategies and then benchmark them to find
the best ones.

### Remove None

This strategy keeps all the data

Use the random strategy but set the parameter 'n' to 0. So that 0 percent of
the data will be randomly removed.


```r
kable(score_all_classifications(test_dat, 'random', params = list(n=0)), digits = 2)
```



|name    | sn|   sp| combined| max_dist| time_taken| input_n| output_n| true_n|
|:-------|--:|----:|--------:|--------:|----------:|-------:|--------:|------:|
|test1   |  1| 0.00|     1.00|   100.00|       0.03|       3|        3|      2|
|test2   |  1| 0.00|     1.00|     2.94|       0.01|       3|        3|      2|
|test3   |  1| 0.00|     1.00|     4.78|       0.01|       3|        3|      0|
|test4   |  1| 0.00|     1.00|     8.19|       0.01|       4|        4|      3|
|test5   |  1| 1.00|     0.00|     1.10|       0.01|       4|        4|      4|
|test6   |  1| 0.00|     1.00|     4.78|       0.01|       8|        8|      7|
|test7   |  1| 0.00|     1.00|     9.25|       0.01|       8|        8|      7|
|test8   |  1| 0.00|     1.00|    10.29|       0.01|      20|       20|     19|
|test9   |  1| 0.00|     1.00|     6.25|       0.01|      20|       20|     18|
|test10  |  1| 0.00|     1.00|    11.78|       0.01|      30|       30|     25|
|test11  |  1| 0.00|     1.00|    59.88|       0.01|     113|      113|    107|
|test12  |  1| 0.00|     1.00|    12.46|       0.01|     119|      119|    107|
|summary |  1| 0.08|     0.92|       NA|       0.01|      NA|       NA|     NA|

### Remove Random

This strategy removes 40% of the sample at random


```r
kable(score_all_classifications(test_dat, 'random', params = list(n=0.4)), digits = 2)
```



|name    |   sn|   sp| combined| max_dist| time_taken| input_n| output_n| true_n|
|:-------|----:|----:|--------:|--------:|----------:|-------:|--------:|------:|
|test1   | 1.00| 0.00|     1.00|   100.00|       0.01|       3|        3|      2|
|test2   | 1.00| 0.00|     1.00|     2.94|       0.01|       3|        3|      2|
|test3   | 1.00| 0.33|     0.67|     4.78|       0.02|       3|        2|      0|
|test4   | 0.67| 0.00|     1.05|     8.19|       0.02|       4|        3|      3|
|test5   | 0.75| 1.00|     0.25|     0.37|       0.02|       4|        3|      4|
|test6   | 0.14| 0.00|     1.32|     4.78|       0.02|       8|        2|      7|
|test7   | 0.71| 0.00|     1.04|     8.90|       0.02|       8|        6|      7|
|test8   | 0.84| 0.00|     1.01|     9.93|       0.02|      20|       17|     19|
|test9   | 0.39| 0.00|     1.17|     6.25|       0.02|      20|        9|     18|
|test10  | 0.32| 0.40|     0.91|     8.92|       0.02|      30|       11|     25|
|test11  | 0.85| 0.17|     0.85|    59.88|       0.02|     113|       96|    107|
|test12  | 0.32| 0.58|     0.80|    11.74|       0.02|     119|       39|    107|
|summary | 0.67| 0.21|     0.92|       NA|       0.02|      NA|       NA|     NA|
### Information Variance Balance

This strategy will keep on removing the most outlying sequence as long as is
leads to a percentage reduction in variance that is x times larger than the
percentage of information that was discarded

#### Threshold of 1

```r
kable(score_all_classifications(test_dat, 'infovar_balance', params = list(threshold = 1)), digits = 2)
```



|name    |   sn|   sp| combined| max_dist| time_taken| input_n| output_n| true_n|
|:-------|----:|----:|--------:|--------:|----------:|-------:|--------:|------:|
|test1   | 1.00| 1.00|     0.00|     0.00|       0.01|       3|        2|      2|
|test2   | 1.00| 1.00|     0.00|     0.00|       0.01|       3|        2|      2|
|test3   | 1.00| 0.33|     0.67|     0.74|       0.01|       3|        2|      0|
|test4   | 1.00| 1.00|     0.00|     1.78|       0.01|       4|        3|      3|
|test5   | 0.50| 1.00|     0.50|     0.00|       0.01|       4|        2|      4|
|test6   | 1.00| 1.00|     0.00|     0.37|       0.02|       8|        7|      7|
|test7   | 1.00| 1.00|     0.00|     1.07|       0.02|       8|        7|      7|
|test8   | 0.95| 1.00|     0.05|     1.47|       0.10|      20|       18|     19|
|test9   | 0.72| 1.00|     0.28|     1.10|       0.10|      20|       13|     18|
|test10  | 0.96| 1.00|     0.04|     2.87|       0.29|      30|       24|     25|
|test11  | 0.98| 1.00|     0.02|     1.77|       4.83|     113|      105|    107|
|test12  | 0.96| 1.00|     0.04|     1.42|       3.94|     119|      103|    107|
|summary | 0.92| 0.94|     0.13|       NA|       0.78|      NA|       NA|     NA|

#### Threshold of 2

```r
kable(score_all_classifications(test_dat, 'infovar_balance', params = list(threshold = 2)), digits = 2)
```



|name    |   sn|   sp| combined| max_dist| time_taken| input_n| output_n| true_n|
|:-------|----:|----:|--------:|--------:|----------:|-------:|--------:|------:|
|test1   | 1.00| 1.00|     0.00|     0.00|       0.01|       3|        2|      2|
|test2   | 1.00| 1.00|     0.00|     0.00|       0.01|       3|        2|      2|
|test3   | 1.00| 0.33|     0.67|     0.74|       0.01|       3|        2|      0|
|test4   | 1.00| 1.00|     0.00|     1.78|       0.01|       4|        3|      3|
|test5   | 0.75| 1.00|     0.25|     0.37|       0.01|       4|        3|      4|
|test6   | 1.00| 1.00|     0.00|     0.37|       0.02|       8|        7|      7|
|test7   | 1.00| 1.00|     0.00|     1.07|       0.02|       8|        7|      7|
|test8   | 1.00| 1.00|     0.00|     1.84|       0.10|      20|       19|     19|
|test9   | 1.00| 0.50|     0.50|     2.94|       0.10|      20|       19|     18|
|test10  | 1.00| 0.80|     0.20|     4.78|       0.29|      30|       26|     25|
|test11  | 1.00| 0.50|     0.50|     3.54|       4.77|     113|      110|    107|
|test12  | 1.00| 1.00|     0.00|     2.49|       3.75|     119|      107|    107|
|summary | 0.98| 0.84|     0.18|       NA|       0.76|      NA|       NA|     NA|

#### Threshold of 3

```r
kable(score_all_classifications(test_dat, 'infovar_balance', params = list(threshold = 3)), digits = 2)
```



|name    | sn|   sp| combined| max_dist| time_taken| input_n| output_n| true_n|
|:-------|--:|----:|--------:|--------:|----------:|-------:|--------:|------:|
|test1   |  1| 0.00|     1.00|   100.00|       0.01|       3|        3|      2|
|test2   |  1| 0.00|     1.00|     2.94|       0.01|       3|        3|      2|
|test3   |  1| 0.00|     1.00|     4.78|       0.01|       3|        3|      0|
|test4   |  1| 1.00|     0.00|     1.78|       0.02|       4|        3|      3|
|test5   |  1| 1.00|     0.00|     1.10|       0.01|       4|        4|      4|
|test6   |  1| 1.00|     0.00|     0.37|       0.02|       8|        7|      7|
|test7   |  1| 1.00|     0.00|     1.07|       0.02|       8|        7|      7|
|test8   |  1| 1.00|     0.00|     1.84|       0.10|      20|       19|     19|
|test9   |  1| 0.50|     0.50|     2.94|       0.10|      20|       19|     18|
|test10  |  1| 0.80|     0.20|     4.78|       0.30|      30|       26|     25|
|test11  |  1| 0.33|     0.67|     4.13|       4.76|     113|      111|    107|
|test12  |  1| 0.67|     0.33|     3.91|       3.72|     119|      111|    107|
|summary |  1| 0.61|     0.39|       NA|       0.76|      NA|       NA|     NA|

### Most Frequent Sequence

```r
kable(score_all_classifications(test_dat, 'most_frequent', params = list()), digits = 2)
```



|name    |   sn|   sp| combined| max_dist| time_taken| input_n| output_n| true_n|
|:-------|----:|----:|--------:|--------:|----------:|-------:|--------:|------:|
|test1   | 1.00| 1.00|     0.00|        0|       0.03|       3|        2|      2|
|test2   | 1.00| 1.00|     0.00|        0|       0.03|       3|        2|      2|
|test3   | 1.00| 0.67|     0.33|        0|       0.02|       3|        1|      0|
|test4   | 0.33| 1.00|     0.67|        0|       0.03|       4|        1|      3|
|test5   | 0.50| 1.00|     0.50|        0|       0.02|       4|        2|      4|
|test6   | 0.86| 1.00|     0.14|        0|       0.03|       8|        6|      7|
|test7   | 0.57| 1.00|     0.43|        0|       0.02|       8|        4|      7|
|test8   | 0.63| 1.00|     0.37|        0|       0.02|      20|       12|     19|
|test9   | 0.39| 1.00|     0.61|        0|       0.02|      20|        7|     18|
|test10  | 0.52| 1.00|     0.48|        0|       0.02|      30|       13|     25|
|test11  | 0.62| 1.00|     0.38|        0|       0.03|     113|       66|    107|
|test12  | 0.55| 1.00|     0.45|        0|       0.02|     119|       59|    107|
|summary | 0.66| 0.97|     0.36|       NA|       0.03|      NA|       NA|     NA|

### The silver bullet

To be devised

## Meeting Agenda

- Simulation of test data for MotifBinner
  - grinder
  - GemSim
  - Evolveagene
  - poorly understood and supremely complex MiSeq error profiles
  - Why am I trying to duplicate the entire process, I should consider a more focussed simulation approach (Available software forces me to, but I should rather adapt polyester?)
- Testing strategies for MotifBinner
  - Benchmarking system demonstration
  - Why am I benchmarking the mislabel detection? Should I just benchmark the entire process?
- Mislabel detection strategies for MotifBinner
  - Why kinky diagrams failed
  - Distance computation inefficiencies caused by duplicates and possible hacks
  - Describe the information variance balance strategy
  - How relevant are the metrics Im using for mislabel detection? Since
    technically we only care about the consensus sequence. Is all that
    matters that none of the 'out' sequences makes it into the 'src' bin?
  - The most "frequent" strategy and some brainstorming about what can go wrong
    - Two+ variant versions
- Gaps and consensus sequences and homopolymer errors
