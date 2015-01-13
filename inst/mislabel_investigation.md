

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



|name    | sn| sp| combined|   max_dist| time_taken| input_n| output_n| true_n|
|:-------|--:|--:|--------:|----------:|----------:|-------:|--------:|------:|
|test1   |  1|  0|      1.0| 100.000000|      0.026|       3|        3|      2|
|test2   |  1|  0|      1.0|   2.941177|      0.015|       3|        3|      2|
|test3   |  1|  0|      1.0|   4.779412|      0.011|       3|        3|      0|
|test4   |  1|  0|      1.0|   8.185053|      0.014|       4|        4|      3|
|test5   |  1|  1|      0.0|   1.102941|      0.012|       4|        4|      4|
|test6   |  1|  0|      1.0|   4.779412|      0.012|       8|        8|      7|
|test7   |  1|  0|      1.0|   9.252669|      0.012|       8|        8|      7|
|test8   |  1|  0|      1.0|  10.294118|      0.013|      20|       20|     19|
|test9   |  1|  0|      1.0|   6.250000|      0.012|      20|       20|     18|
|test10  |  1|  0|      1.0|  11.783440|      0.013|      30|       30|     25|
|summary | NA| NA|      0.9|         NA|      0.014|      NA|       NA|     NA|

### Remove Random

This strategy removes 40% of the sample at random


```r
kable(score_all_classifications(test_dat, 'random', params = list(n=0.4)))
```



|name    |        sn|        sp|  combined|    max_dist| time_taken| input_n| output_n| true_n|
|:-------|---------:|---------:|---------:|-----------:|----------:|-------:|--------:|------:|
|test1   | 1.0000000| 0.0000000| 1.0000000| 100.0000000|     0.0130|       3|        3|      2|
|test2   | 1.0000000| 0.0000000| 1.0000000|   2.9411765|     0.0130|       3|        3|      2|
|test3   | 1.0000000| 0.3333333| 0.6666667|   0.7352941|     0.0200|       3|        2|      0|
|test4   | 0.6666667| 0.0000000| 1.0540926|   8.1850534|     0.0200|       4|        3|      3|
|test5   | 0.5000000| 1.0000000| 0.5000000|   1.1029412|     0.0180|       4|        2|      4|
|test6   | 0.8571429| 0.0000000| 1.0101525|   4.4117647|     0.0200|       8|        7|      7|
|test7   | 0.8571429| 1.0000000| 0.1428571|   0.7117438|     0.0180|       8|        6|      7|
|test8   | 0.8947368| 1.0000000| 0.1052632|   1.8382353|     0.0180|      20|       17|     19|
|test9   | 0.7222222| 0.0000000| 1.0378634|   6.2500000|     0.0190|      20|       15|     18|
|test10  | 0.8400000| 0.6000000| 0.4308132|   8.2802548|     0.0180|      30|       23|     25|
|summary |        NA|        NA| 0.6947709|          NA|     0.0177|      NA|       NA|     NA|
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
|test9   | 0.72| 1.00|     0.28|     1.10|       0.11|      20|       13|     18|
|test10  | 0.96| 1.00|     0.04|     2.87|       0.30|      30|       24|     25|
|summary |   NA|   NA|     0.15|       NA|       0.06|      NA|       NA|     NA|

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
|test7   | 1.00| 1.00|     0.00|     1.07|       0.03|       8|        7|      7|
|test8   | 1.00| 1.00|     0.00|     1.84|       0.10|      20|       19|     19|
|test9   | 1.00| 0.50|     0.50|     2.94|       0.10|      20|       19|     18|
|test10  | 1.00| 0.80|     0.20|     4.78|       0.30|      30|       26|     25|
|summary |   NA|   NA|     0.16|       NA|       0.06|      NA|       NA|     NA|

#### Threshold of 3

```r
kable(score_all_classifications(test_dat, 'infovar_balance', params = list(threshold = 3)), digits = 2)
```



|name    | sn|  sp| combined| max_dist| time_taken| input_n| output_n| true_n|
|:-------|--:|---:|--------:|--------:|----------:|-------:|--------:|------:|
|test1   |  1| 0.0|     1.00|   100.00|       0.01|       3|        3|      2|
|test2   |  1| 0.0|     1.00|     2.94|       0.01|       3|        3|      2|
|test3   |  1| 0.0|     1.00|     4.78|       0.01|       3|        3|      0|
|test4   |  1| 1.0|     0.00|     1.78|       0.01|       4|        3|      3|
|test5   |  1| 1.0|     0.00|     1.10|       0.01|       4|        4|      4|
|test6   |  1| 1.0|     0.00|     0.37|       0.02|       8|        7|      7|
|test7   |  1| 1.0|     0.00|     1.07|       0.02|       8|        7|      7|
|test8   |  1| 1.0|     0.00|     1.84|       0.10|      20|       19|     19|
|test9   |  1| 0.5|     0.50|     2.94|       0.10|      20|       19|     18|
|test10  |  1| 0.8|     0.20|     4.78|       0.30|      30|       26|     25|
|summary | NA|  NA|     0.37|       NA|       0.06|      NA|       NA|     NA|

### The silver bullet

To be devised

