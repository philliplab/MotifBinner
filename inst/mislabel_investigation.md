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

We need to test the system using bins with known answers. As always, the design
data structure for the test data is important. Use a structure like this:


```r
list('test1' = list('in' = DNAStringSet(...),
                    'out' = DNAStringSet(...)),
     'test2' = ...)
```

Using this data for a 'bin' by putting the in and out data together and check
that only the 'out' data is removed and all the 'in' data is kept.

Basic code to run the tests:


```r
test_dat <- get_mislabel_test_data()
results <- list()
for (test_num in names(test_dat)){
  all_dat <- c(test_dat[[test_num]]$in,
               test_dat[[test_num]]$out)
  pure_dat <- remove_contamination(all_dat)
  results(test_num) <- c(sum(test_dat[[test_num]]$in %in% pure_dat)/length(test_dat[[test_num]]$in),
                         sum(test_dat[[test_num]]$out %in% pure_dat)/min(1,length(test_dat[[test_num]]$out)))

}
```

## The different strategies

Design and implement a who set of strategies and then benchmark them to find
the best ones.

### Remove None

This strategy keeps all the data

### Remove Random

This strategy removes x% of the sample at random

### Information Variance Balance

This strategy will keep on removing the most outlying sequence as long as is
leads to a percentage reduction in variance that is x times larger than the
percentage of information that was discarded

### The silver bullet

To be devised


