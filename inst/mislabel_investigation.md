

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



|name   | sn| sp|   max_dist| time_taken|
|:------|--:|--:|----------:|----------:|
|test1  |  1|  0| 100.000000|      0.027|
|test2  |  1|  0|   2.941177|      0.016|
|test3  |  1|  0|   4.779412|      0.013|
|test4  |  1|  0|   8.185053|      0.014|
|test5  |  1|  1|   1.102941|      0.013|
|test6  |  1|  0|   4.779412|      0.012|
|test7  |  1|  0|   9.252669|      0.012|
|test8  |  1|  0|  10.294118|      0.013|
|test9  |  1|  0|   6.250000|      0.012|
|test10 |  1|  0|  11.783440|      0.012|

### Remove Random

This strategy removes 40% of the sample at random


```r
kable(score_all_classifications(test_dat, 'random', params = list(n=0.4)))
```



|name   |        sn|        sp|    max_dist| time_taken|
|:------|---------:|---------:|-----------:|----------:|
|test1  | 1.0000000| 0.0000000| 100.0000000|      0.012|
|test2  | 1.0000000| 0.0000000|   2.9411765|      0.013|
|test3  | 1.0000000| 0.3333333|   4.7794118|      0.021|
|test4  | 0.6666667| 0.0000000|   8.1850534|      0.019|
|test5  | 0.5000000| 1.0000000|   1.1029412|      0.019|
|test6  | 0.1428571| 0.0000000|   4.7794118|      0.020|
|test7  | 0.8571429| 1.0000000|   0.7117438|      0.018|
|test8  | 0.8947368| 1.0000000|   1.4705882|      0.020|
|test9  | 0.7777778| 0.5000000|   6.2500000|      0.018|
|test10 | 0.8400000| 0.6000000|   8.9171975|      0.020|
### Information Variance Balance

This strategy will keep on removing the most outlying sequence as long as is
leads to a percentage reduction in variance that is x times larger than the
percentage of information that was discarded

#### Threshold of 1

```r
kable(score_all_classifications(test_dat, 'infovar_balance', params = list(threshold = 1)), digits = 2)
```



|name   |   sn|   sp| max_dist| time_taken|
|:------|----:|----:|--------:|----------:|
|test1  | 1.00| 1.00|     0.00|       0.01|
|test2  | 1.00| 1.00|     0.00|       0.01|
|test3  | 1.00| 0.33|     0.74|       0.01|
|test4  | 1.00| 1.00|     1.78|       0.01|
|test5  | 0.50| 1.00|     0.00|       0.01|
|test6  | 1.00| 1.00|     0.37|       0.02|
|test7  | 1.00| 1.00|     1.07|       0.02|
|test8  | 0.95| 1.00|     1.47|       0.10|
|test9  | 0.72| 1.00|     1.10|       0.11|
|test10 | 0.96| 1.00|     2.87|       0.30|

#### Threshold of 2

```r
kable(score_all_classifications(test_dat, 'infovar_balance', params = list(threshold = 2)), digits = 2)
```



|name   |   sn|   sp| max_dist| time_taken|
|:------|----:|----:|--------:|----------:|
|test1  | 1.00| 1.00|     0.00|       0.01|
|test2  | 1.00| 1.00|     0.00|       0.01|
|test3  | 1.00| 0.33|     0.74|       0.01|
|test4  | 1.00| 1.00|     1.78|       0.01|
|test5  | 0.75| 1.00|     0.37|       0.01|
|test6  | 1.00| 1.00|     0.37|       0.02|
|test7  | 1.00| 1.00|     1.07|       0.02|
|test8  | 1.00| 1.00|     1.84|       0.11|
|test9  | 1.00| 0.50|     2.94|       0.10|
|test10 | 1.00| 0.80|     4.78|       0.31|

#### Threshold of 3

```r
kable(score_all_classifications(test_dat, 'infovar_balance', params = list(threshold = 3)), digits = 2)
```



|name   | sn|  sp| max_dist| time_taken|
|:------|--:|---:|--------:|----------:|
|test1  |  1| 0.0|   100.00|       0.01|
|test2  |  1| 0.0|     2.94|       0.01|
|test3  |  1| 0.0|     4.78|       0.01|
|test4  |  1| 1.0|     1.78|       0.01|
|test5  |  1| 1.0|     1.10|       0.01|
|test6  |  1| 1.0|     0.37|       0.02|
|test7  |  1| 1.0|     1.07|       0.03|
|test8  |  1| 1.0|     1.84|       0.10|
|test9  |  1| 0.5|     2.94|       0.10|
|test10 |  1| 0.8|     4.78|       0.30|

### The silver bullet

To be devised

