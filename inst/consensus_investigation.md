

# Investigation of various approach to constructing a consensus string


```r
Sys.time()
```

```
## [1] "2015-03-11 13:27:40 SAST"
```

## Overview

Three important steps needs to be optimized:
* mislabel classification
* alignment
* consensus string construction

This document investigates all three steps simultaneously by a simulation
process:
* Given an input string and a read error profile
* Potentatially given a contamination string
* Simulate reads from both strings and form a bin
* Use mislabel classification to throw out outliers
* Align the remaining sequences
* Construct a final sequence
* Compare the reconstructed sequence to the input sequence with edit distance

We need very clear terminology to keep track of the different aspects of the
simulations:
* All the parameters that goes into the three steps are called a 'setup'
* A input sequence, number of reads, contamination sequence, number of
  contamination reads and the error profile is called a 'scenario' or 'test
  scenario'.
* A test scenario, setup and a seed is called a 'case' or 'test case'.

## The setups


```r
setups <- list()
setups[['base']] <- list(classification_technique = 'infovar_balance', 
                         classification_params = list(threshold = 1), 
                         alignment_technique = 'muscle', 
                         alignment_params = list(), 
                         consensus_technique = 'Biostrings::consensusString', 
                         consensus_params = list(ambiguityMap = 'N'))

setups[['base1']] <- list(classification_technique = 'infovar_balance', 
                         classification_params = list(threshold = 1,
                                          start_threshold = 0.02), 
                         alignment_technique = 'muscle', 
                         alignment_params = list(), 
                         consensus_technique = 'Biostrings::consensusString', 
                         consensus_params = list(ambiguityMap = 'N'))


setups[['most_con']] <- list(classification_technique = 'infovar_balance', 
                         classification_params = list(threshold = 1,
                                          start_threshold = 0.02), 
                         alignment_technique = 'muscle', 
                         alignment_params = list(), 
                         consensus_technique = 'mostConsensusString', 
                         consensus_params = list())
```

## The error profiles

### Three very basic error profiles


```r
err_profiles <- list()

x <- list(read_len = 1000,
          A_err_rates = list('A' = 0.997, 'C' = 0.001, 'G' = 0.001, 'T' = 0.001),
          C_err_rates = list('C' = 0.997, 'A' = 0.001, 'G' = 0.001, 'T' = 0.001),
          G_err_rates = list('G' = 0.997, 'C' = 0.001, 'A' = 0.001, 'T' = 0.001),
          T_err_rates = list('T' = 0.997, 'C' = 0.001, 'G' = 0.001, 'A' = 0.001))
params <- list()
params[['params']] <- x
params[['technique']] <- 'uniform'
err_profiles[['unif_1_1000']] <- do.call(gen_error_profile, params)


x <- list(read_len = 1000,
          A_err_rates = list('A' = 0.97, 'C' = 0.01, 'G' = 0.01, 'T' = 0.01),
          C_err_rates = list('C' = 0.97, 'A' = 0.01, 'G' = 0.01, 'T' = 0.01),
          G_err_rates = list('G' = 0.97, 'C' = 0.01, 'A' = 0.01, 'T' = 0.01),
          T_err_rates = list('T' = 0.97, 'C' = 0.01, 'G' = 0.01, 'A' = 0.01))
params <- list()
params[['params']] <- x
params[['technique']] <- 'uniform'
err_profiles[['unif_1_100']] <- do.call(gen_error_profile, params)

x <- list(read_len = 1000,
          A_err_rates = list('A' = 0.7, 'C' = 0.1, 'G' = 0.1, 'T' = 0.1),
          C_err_rates = list('C' = 0.7, 'A' = 0.1, 'G' = 0.1, 'T' = 0.1),
          G_err_rates = list('G' = 0.7, 'C' = 0.1, 'A' = 0.1, 'T' = 0.1),
          T_err_rates = list('T' = 0.7, 'C' = 0.1, 'G' = 0.1, 'A' = 0.1))
params <- list()
params[['params']] <- x
params[['technique']] <- 'uniform'
err_profiles[['unif_1_10']] <- do.call(gen_error_profile, params)
```

## The test scenarios


```r
scenarios <- list()
scenarios[['unif_read_1']] <- list(name = 'unif_read_1',
                                   ref_seq = paste(rep('A', 500), collapse = ""), 
                                   n_reads = 10, 
                                   error_rates = err_profiles[['unif_1_1000']], 
                                   contam_seq = NULL, 
                                   n_contam = 0)
scenarios[['unif_read_2']] <- list(name = 'unif_read_2',
                                   ref_seq = paste(rep('A', 500), collapse = ""), 
                                   n_reads = 10, 
                                   error_rates = err_profiles[['unif_1_100']], 
                                   contam_seq = NULL, 
                                   n_contam = 0)
scenarios[['unif_read_3']] <- list(name = 'unif_read_3',
                                   ref_seq = paste(rep('A', 500), collapse = ""), 
                                   n_reads = 10, 
                                   error_rates = err_profiles[['unif_1_10']], 
                                   contam_seq = NULL, 
                                   n_contam = 0)

scenarios[['unif_contam_1']] <- list(name = 'unif_contam_1',
                                     ref_seq = paste(rep('A', 500), collapse = ""), 
                                     n_reads = 10, 
                                     error_rates = err_profiles[['unif_1_100']], 
                                     contam_seq = paste(c(rep('A', 200), rep('C', 100), rep('A', 200)), collapse = ""), 
                                     n_contam = 1)
scenarios[['unif_contam_2']] <- list(name = 'unif_contam_2',
                                     ref_seq = paste(rep('A', 500), collapse = ""), 
                                     n_reads = 10, 
                                     error_rates = err_profiles[['unif_1_100']], 
                                     contam_seq = paste(c(rep('A', 200), rep('C', 100), rep('A', 200)), collapse = ""), 
                                     n_contam = 5)
scenarios[['unif_contam_3']] <- list(name = 'unif_contam_3',
                                     ref_seq = paste(rep('A', 500), collapse = ""), 
                                     n_reads = 10, 
                                     error_rates = err_profiles[['unif_1_100']], 
                                     contam_seq = paste(c(rep('A', 200), rep('C', 100), rep('A', 200)), collapse = ""), 
                                     n_contam = 9)
```

## The test cases


```r
cases <- list()
seeds <- 1:3

# Note this code is gonna become complex and is going to interact with the case
# specification is weird ways - so leave it here
for (setup in names(setups)){
  for (scenario in names(scenarios)){
    for (seed in seeds){
      case_name <- paste(setup, scenario, seed, sep = '-')
      cases[[case_name]] <- list(scenario = scenarios[[scenario]],
                                 seed = seed,
                                 setup = setups[[setup]])
    }
  }
}
```

## The test runner


```r
unique_scenarios <- list_unique_scenarios(cases)

cache_file <- '~/projects/MotifBinner/code/MotifBinner/inst/scenario_cache.rdata'

scenario_cache <- load_or_initialize_cache(cache_file)

scenario_cache <- create_scenario_data(unique_scenarios, scenario_cache)

save(scenario_cache, file = cache_file)
```

## Running the test cases.


```r
results <- run_all_tests(cases, scenario_cache)
```

## The results

```r
kable(results)
```



|setup    |scenario      |seed | mismatch| input_len| output_len| mismatch_rate|
|:--------|:-------------|:----|--------:|---------:|----------:|-------------:|
|base     |unif_read_1   |1    |        1|       500|        500|         0.002|
|base     |unif_read_1   |2    |        0|       500|          0|         0.000|
|base     |unif_read_1   |3    |        0|       500|          0|         0.000|
|base     |unif_read_2   |1    |        0|       500|          0|         0.000|
|base     |unif_read_2   |2    |        0|       500|          0|         0.000|
|base     |unif_read_2   |3    |        0|       500|          0|         0.000|
|base     |unif_read_3   |1    |      171|       500|        507|         0.342|
|base     |unif_read_3   |2    |      173|       500|        500|         0.346|
|base     |unif_read_3   |3    |      160|       500|        519|         0.320|
|base     |unif_contam_1 |1    |        0|       500|          0|         0.000|
|base     |unif_contam_1 |2    |        0|       500|          0|         0.000|
|base     |unif_contam_1 |3    |        0|       500|          0|         0.000|
|base     |unif_contam_2 |1    |      197|       500|        597|         0.394|
|base     |unif_contam_2 |2    |      199|       500|        599|         0.398|
|base     |unif_contam_2 |3    |      200|       500|        600|         0.400|
|base     |unif_contam_3 |1    |      200|       500|        600|         0.400|
|base     |unif_contam_3 |2    |      200|       500|        600|         0.400|
|base     |unif_contam_3 |3    |      200|       500|        600|         0.400|
|base1    |unif_read_1   |1    |        0|       500|          0|         0.000|
|base1    |unif_read_1   |2    |        0|       500|          0|         0.000|
|base1    |unif_read_1   |3    |        0|       500|          0|         0.000|
|base1    |unif_read_2   |1    |        0|       500|          0|         0.000|
|base1    |unif_read_2   |2    |        0|       500|          0|         0.000|
|base1    |unif_read_2   |3    |        0|       500|          0|         0.000|
|base1    |unif_read_3   |1    |      171|       500|        507|         0.342|
|base1    |unif_read_3   |2    |      173|       500|        500|         0.346|
|base1    |unif_read_3   |3    |      160|       500|        519|         0.320|
|base1    |unif_contam_1 |1    |        0|       500|          0|         0.000|
|base1    |unif_contam_1 |2    |        0|       500|          0|         0.000|
|base1    |unif_contam_1 |3    |        0|       500|          0|         0.000|
|base1    |unif_contam_2 |1    |      197|       500|        597|         0.394|
|base1    |unif_contam_2 |2    |      199|       500|        599|         0.398|
|base1    |unif_contam_2 |3    |      200|       500|        600|         0.400|
|base1    |unif_contam_3 |1    |      200|       500|        600|         0.400|
|base1    |unif_contam_3 |2    |      200|       500|        600|         0.400|
|base1    |unif_contam_3 |3    |      200|       500|        600|         0.400|
|most_con |unif_read_1   |1    |        0|       500|          0|         0.000|
|most_con |unif_read_1   |2    |        0|       500|          0|         0.000|
|most_con |unif_read_1   |3    |        0|       500|          0|         0.000|
|most_con |unif_read_2   |1    |        0|       500|          0|         0.000|
|most_con |unif_read_2   |2    |        0|       500|          0|         0.000|
|most_con |unif_read_2   |3    |        0|       500|          0|         0.000|
|most_con |unif_read_3   |1    |       79|       500|        500|         0.158|
|most_con |unif_read_3   |2    |       87|       500|        500|         0.174|
|most_con |unif_read_3   |3    |       76|       500|        500|         0.152|
|most_con |unif_contam_1 |1    |        0|       500|          0|         0.000|
|most_con |unif_contam_1 |2    |        0|       500|          0|         0.000|
|most_con |unif_contam_1 |3    |        0|       500|          0|         0.000|
|most_con |unif_contam_2 |1    |        0|       500|          0|         0.000|
|most_con |unif_contam_2 |2    |        0|       500|          0|         0.000|
|most_con |unif_contam_2 |3    |        0|       500|          0|         0.000|
|most_con |unif_contam_3 |1    |       27|       500|        500|         0.054|
|most_con |unif_contam_3 |2    |       38|       500|        500|         0.076|
|most_con |unif_contam_3 |3    |       36|       500|        500|         0.072|

## Fixes resulting from benchmarking

### ur1_b_1

A single mismatch was found on this super simple and easy case. The primary
cause is that the infovar_balance classification technique tries to remove
contamination from a bin with no contamination. Set a starting criteria for the
infovar_balance technique so that it will only start if the maximum distance
in the distance matrix is above some threshold.


```r
# Generate test data
test_bin <- do.call(gen_and_contaminate_reads, c(scenarios[['unif_read_1']], list(seed=1)))

# See how easy the problem is
consensusString(test_bin$src) == test_bin$true_consensus
```

```
## [1] TRUE
```

```r
# See how the basic technique fails

params <- setups[['base']]
params$test_bin <- test_bin
result <- do.call(score_consensus, params)
result$edit_dist
```

```
## [1] 1
```

```r
# Fix it with the starting threshold:
params[['classification_params']] <- list(threshold = 1,
                                          start_threshold = 0.02)
result <- do.call(score_consensus, params)
result$edit_dist
```

```
## [1] 0
```

## Utility Functions


```r
list_to_env <- function(x){
  for (i in names(x)){
    p <- list(x = i, value = x[[i]], envir = .GlobalEnv)
    do.call(assign, p)
  }
}
```

## How does the consensusString Function in Biostrings actually work?

When and how does ambigueity characters get added to the data. Should I write
my own version of the function to deal with it?

### It behaves as expected for matching sequences

```r
consensusString(DNAStringSet(c('AAA', 'AAA')))
```

```
## [1] "AAA"
```

```r
consensusString(DNAStringSet(c('ACGT', 'ACGT')))
```

```
## [1] "ACGT"
```

### For mismatches - a very strict threshold must be met

```r
consensusString(DNAStringSet(c('AAC', 'AAA')))
```

```
## [1] "AAM"
```

```r
consensusString(DNAStringSet(c('AAC', 'AAA', 'AAA')))
```

```
## [1] "AAM"
```

```r
consensusString(DNAStringSet(c('AAC', 'AAA', 'AAA', 'AAA')))
```

```
## [1] "AAM"
```

```r
consensusString(DNAStringSet(c('AAC', 'AAA', 'AAA', 'AAA', 'AAA')))
```

```
## [1] "AAA"
```

```r
consensusString(DNAStringSet(c('AAC', 'AAT', 'AAA', 'AAA', 'AAA', 'AAA')))
```

```
## [1] "AAA"
```

### Can this threshold be lowered?

```r
consensusString(DNAStringSet(c('AAC', 'AAA')), threshold = 0.5)
```

```
## Error in .local(x, ...): 'threshold' must be a numeric in (0, 1/sum(nchar(ambiguityMap) == 1)]
```

No because they are using some strange restrictions.  Specifically the
threshold is that no other character may occur more than 25% of the time. This
is stricter than the current behaviour. Is this a problem or can the stricter
behaviour be used?

If you want to use the 50% criteria again, then I would have to either write my
own consensusString generator or modify theirs - it might take a bit of time.

I think it would probably be better to reduce that restriction because what
this means is that any bin with 2, 3 or 4 sequences that are not identical will
have degeneracy (unless the sequences with mismatches were removed by the
outlier detector)

Added easyConsensusString to handle this. Which later morphed into
mostConsensusString.

## Benchmarks for contamination

The purpose here is to create a benchmark that will test the effectiveness of
outlier removal and explore what happens in some edge cases. Some of these
benchmarks should make it into the unit tests eventually.

Testing the interaction between contamination and read errors is not so
important, so only use 1 in 100 read error rate.

What is important to explore is what happens when the level of contamination
changes. so try with 10%, 25% and 45% and 50% contamination and check what
happens. Note that the mislabel detector will really struggle on 45% and 50% in
its current form so this is probably where I will hit the crossroads and have
to implement the mislabel detector that is based on absolute thresholds.

The sequence to use for the contamination should not be too rediculous
otherwise it will make the aligner go nuts if the mislabel detector fails to
remove it. The current test string is 500 A's So for contamination try using
200 As followed by 100 Cs then another 200 As. This is 1 mutation / read error
for every 5 bases which is very obviously a contaminant. Also try a lower level
of errors. Try 245 As, 10 Cs, 245 As. This is still a read error / mutation
rate of 1 in 50 which is higher than the simulated read error rate, so it
should be detectable.

The problem is still that the aligner chooses to insert gaps with these thus
greatly inflating the mismatch scores. Well, maybe this is not that bad since
the benchmark shows that the mislabel detector failed and that is kind the only
thing that matters.

## Read Error Rate Computations

Assuming that the read error rate is described as 1 errorneous base for each
100 bases sequenced, how many mismatches can you expect between two reads in a
bin? What kind of a distribution is this? This is a binomial distribution.

So, what is the probability that there will be 0, 1, 2, 3, ... read errors in a
read of length 500 if the read error rate is 1 in 100?


```r
plot(dbinom(0:20, 500, 1/100))
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-1.png) 

How do we use this to decide when a bin contains outliers? If the distance
between reads are such that it is unlikely that the only source of errors is
the read errors. So the question now becomes:

What is the probability that the distance between two reads is 0 given the only
source of discrepancies is the read errors. This is the probability that
neither read contains any errors plus the chance that the two reads has exactly
the same read errors. The chance of exactly the same read errors occurring are
complicating the matter and the chance of that happening is small, so lets
ignore that for now. 

Hence the question becomes: "What is the chance of there being 0 read errors in
a read of length 1000 when the chance of a read error is 1/100.


```r
plot(dbinom(0:20, 1000, 1/100))
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-15-1.png) 

Now, how many read errors can be allowed in a bin before we consider it to be
suspicious? If the chance is less than 5% of there being that many read errors,
then we should suspect alternative error sources. So lets assume a bin of size
two. We just have to find the 95th percentile of the distribution and we are
done.


```r
sample_20 <- pbinom(0:20, 1000, 1/100)
plot(x = 0:20, y=sample_20)
abline(h=0.95)
dist_95 <- abs(sample_20-0.95)
closest_95 <- (0:20)[which(dist_95 == min(dist_95))]
abline(v=closest_95)
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16-1.png) 

Now, consider the problematic case where there are more than two sequences in
the bins. Now the question becomes, What is the probability that the largest
and second largest order statistic of a random sample from a binomial
distribution is distance x from each other. The distributions of the order
statistics are significantly more complex than the binomial distribution.

How can you prove this point? Consider the cases where the only source of
distance between sequences are the read errors. Now think what will happen as
the bin size approaches infinity? For any given number x, there will exist a
number y so that if the bin size exceeds y, the probability of there being a
read with more than x errors will exceed say 50%. So this demonstrates that
even if the only source of distance between sequences is read errors, there
will be a big bin size at which any absolute threshold is exceeded.

Now the question is how significant is this effect?

First we need to compute some order statistics. See
http://en.wikipedia.org/wiki/Order_statistic


```r
p1f <- function(x, n, p){
  pbinom(x, n, p) - dbinom(x, n, p)
}

p2f <- function(x, n, p){
  dbinom(x, n, p)
}

p3f <- function(x, n, p){
  1 - pbinom(x, n, p)
}

ord_F <- function(x, k, o_n, n, p){
  p1 <- p1f(x, n, p)
  p2 <- p2f(x, n, p)
  p3 <- p3f(x, n, p)
  summ <- 0
  for (j in 0:(o_n-k)){
    summ <- summ + choose(o_n, j) * p3^j * (p1+p2)^(o_n-j)
  }
  summ
}
```

Lets just start of simple. What happens to the probability that a single read
contains 5 or more read errors given that the sequence is of length 500 and the
read error rate is 1 in 100?

So lets just look at the maximum order statistic and systematically increase
the number of reads in the bin

```r
probs <- rep(0, 50)
for (i in 1:50){
  probs[i] <- ord_F(10, i, i, 500, 1/100)
}
plot(probs)
```

![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-18-1.png) 




















