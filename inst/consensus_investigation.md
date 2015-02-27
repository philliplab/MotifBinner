

# Investigation of various approach to constructing a consensus string

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
scenarios[['unif_read_1']] <- list(ref_seq = paste(rep('A', 500), collapse = ""), 
                                   n_reads = 10, 
                                   error_rates = err_profiles[['unif_1_1000']], 
                                   contam_seq = NULL, 
                                   n_contam = 0)
scenarios[['unif_read_2']] <- list(ref_seq = paste(rep('A', 500), collapse = ""), 
                                   n_reads = 10, 
                                   error_rates = err_profiles[['unif_1_100']], 
                                   contam_seq = NULL, 
                                   n_contam = 0)
scenarios[['unif_read_3']] <- list(ref_seq = paste(rep('A', 500), collapse = ""), 
                                   n_reads = 10, 
                                   error_rates = err_profiles[['unif_1_10']], 
                                   contam_seq = NULL, 
                                   n_contam = 0)

scenarios[['unif_contam_1']] <- list(ref_seq = paste(rep('A', 500), collapse = ""), 
                                   n_reads = 10, 
                                   error_rates = err_profiles[['unif_1_100']], 
                                   contam_seq = paste(rep('C', 500), collapse = ""), 
                                   n_contam = 1)
scenarios[['unif_contam_2']] <- list(ref_seq = paste(rep('A', 500), collapse = ""), 
                                   n_reads = 10, 
                                   error_rates = err_profiles[['unif_1_100']], 
                                   contam_seq = paste(rep('C', 500), collapse = ""), 
                                   n_contam = 5)
scenarios[['unif_contam_3']] <- list(ref_seq = paste(rep('A', 500), collapse = ""), 
                                   n_reads = 10, 
                                   error_rates = err_profiles[['unif_1_100']], 
                                   contam_seq = paste(rep('C', 500), collapse = ""), 
                                   n_contam = 9)
```

## The test cases


```r
cases <- list()
seeds <- 1:3

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
gacr_m <- memoise(gen_and_contaminate_reads)

run_test <- function(scenario, seed, setup){
  params <- scenario
  params[['seed']] <- seed

  test_bin <- do.call( gacr_m, params)

  params <- setup
  params$test_bin <- test_bin
  result <- do.call(score_consensus, params)
  return(result$edit_dist)
}
```

## Running the test cases.


```r
results <- data.frame(setup = character(0),
                      scenario = character(0),
                      seed = character(0),
                      mismatch = numeric(0),
                      input_len = numeric(0),
                      mismatch_rate = numeric(0))

for (tc in names(cases)){
  print(tc)
  split_names <- strsplit(tc, split='-')[[1]]
  setup <- split_names[1]
  scenario <- split_names[2]
  seed <- split_names[3]
  mismatch <- do.call(run_test, cases[[tc]])
  input_len <- nchar(cases[[tc]][['scenario']][['ref_seq']])
  results <- rbind(results,
                   data.frame(setup = setup,
                              scenario = scenario,
                              seed = seed,
                              mismatch = mismatch,
                              input_len = input_len,
                              mismatch_rate = mismatch / input_len))
}
```

## The results

```r
kable(results)
```



|setup    |scenario      |seed | mismatch| input_len| mismatch_rate|
|:--------|:-------------|:----|--------:|---------:|-------------:|
|base     |unif_read_1   |1    |        1|       500|         0.002|
|base     |unif_read_1   |2    |        0|       500|         0.000|
|base     |unif_read_1   |3    |        0|       500|         0.000|
|base     |unif_read_2   |1    |        0|       500|         0.000|
|base     |unif_read_2   |2    |        0|       500|         0.000|
|base     |unif_read_2   |3    |        0|       500|         0.000|
|base     |unif_read_3   |1    |      220|       500|         0.440|
|base     |unif_read_3   |2    |      225|       500|         0.450|
|base     |unif_read_3   |3    |      209|       500|         0.418|
|base     |unif_contam_1 |1    |        0|       500|         0.000|
|base     |unif_contam_1 |2    |        0|       500|         0.000|
|base     |unif_contam_1 |3    |        0|       500|         0.000|
|base     |unif_contam_2 |1    |        0|       500|         0.000|
|base     |unif_contam_2 |2    |        0|       500|         0.000|
|base     |unif_contam_2 |3    |        0|       500|         0.000|
|base     |unif_contam_3 |1    |      999|       500|         1.998|
|base     |unif_contam_3 |2    |      999|       500|         1.998|
|base     |unif_contam_3 |3    |      999|       500|         1.998|
|base1    |unif_read_1   |1    |        0|       500|         0.000|
|base1    |unif_read_1   |2    |        0|       500|         0.000|
|base1    |unif_read_1   |3    |        0|       500|         0.000|
|base1    |unif_read_2   |1    |        0|       500|         0.000|
|base1    |unif_read_2   |2    |        0|       500|         0.000|
|base1    |unif_read_2   |3    |        0|       500|         0.000|
|base1    |unif_read_3   |1    |      220|       500|         0.440|
|base1    |unif_read_3   |2    |      225|       500|         0.450|
|base1    |unif_read_3   |3    |      209|       500|         0.418|
|base1    |unif_contam_1 |1    |        0|       500|         0.000|
|base1    |unif_contam_1 |2    |        0|       500|         0.000|
|base1    |unif_contam_1 |3    |        0|       500|         0.000|
|base1    |unif_contam_2 |1    |        0|       500|         0.000|
|base1    |unif_contam_2 |2    |        0|       500|         0.000|
|base1    |unif_contam_2 |3    |        0|       500|         0.000|
|base1    |unif_contam_3 |1    |      999|       500|         1.998|
|base1    |unif_contam_3 |2    |      999|       500|         1.998|
|base1    |unif_contam_3 |3    |      999|       500|         1.998|
|most_con |unif_read_1   |1    |        0|       500|         0.000|
|most_con |unif_read_1   |2    |        0|       500|         0.000|
|most_con |unif_read_1   |3    |        0|       500|         0.000|
|most_con |unif_read_2   |1    |        0|       500|         0.000|
|most_con |unif_read_2   |2    |        0|       500|         0.000|
|most_con |unif_read_2   |3    |        0|       500|         0.000|
|most_con |unif_read_3   |1    |      138|       500|         0.276|
|most_con |unif_read_3   |2    |      168|       500|         0.336|
|most_con |unif_read_3   |3    |      147|       500|         0.294|
|most_con |unif_contam_1 |1    |        0|       500|         0.000|
|most_con |unif_contam_1 |2    |        0|       500|         0.000|
|most_con |unif_contam_1 |3    |        0|       500|         0.000|
|most_con |unif_contam_2 |1    |        0|       500|         0.000|
|most_con |unif_contam_2 |2    |        0|       500|         0.000|
|most_con |unif_contam_2 |3    |        0|       500|         0.000|
|most_con |unif_contam_3 |1    |      138|       500|         0.276|
|most_con |unif_contam_3 |2    |      149|       500|         0.298|
|most_con |unif_contam_3 |3    |      141|       500|         0.282|

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




























