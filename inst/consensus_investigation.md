

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
                         consensus_params = list())

setups[['base1']] <- list(classification_technique = 'infovar_balance', 
                         classification_params = list(threshold = 1,
                                          start_threshold = 0.02), 
                         alignment_technique = 'muscle', 
                         alignment_params = list(), 
                         consensus_technique = 'Biostrings::consensusString', 
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
```

## The test cases


```r
cases <- list()
cases[['ur1_b_1']] <- list(scenario = scenarios[['unif_read_1']],
                           seed = 1,
                           setup = setups[['base']])
cases[['ur2_b_1']] <- list(scenario = scenarios[['unif_read_2']],
                           seed = 1,
                           setup = setups[['base']])
cases[['ur3_b_1']] <- list(scenario = scenarios[['unif_read_3']],
                           seed = 1,
                           setup = setups[['base']])
cases[['ur1_b1_1']] <- list(scenario = scenarios[['unif_read_1']],
                            seed = 1,
                            setup = setups[['base1']])
cases[['ur2_b1_1']] <- list(scenario = scenarios[['unif_read_2']],
                            seed = 1,
                            setup = setups[['base1']])
cases[['ur3_b1_1']] <- list(scenario = scenarios[['unif_read_3']],
                            seed = 1,
                            setup = setups[['base1']])
```

## The test runner


```r
run_test <- function(scenario, seed, setup){
  params <- scenario
  params[['seed']] <- seed

  test_bin <- do.call( gen_and_contaminate_reads, params)

  params <- setup
  params$test_bin <- test_bin
  result <- do.call(score_consensus, params)
  return(result$edit_dist)
}
```

## Running the test cases.


```r
results <- data.frame(case = character(0),
                      mismatch = numeric(0),
                      input_len = numeric(0),
                      mismatch_rate = numeric(0))

for (tc in names(cases)){
  mismatch <- do.call(run_test, cases[[tc]])
  input_len <- nchar(cases[[tc]][['scenario']][['ref_seq']])
  results <- rbind(results,
                   data.frame(case = tc,
                              mismatch = mismatch,
                              input_len = input_len,
                              mismatch_rate = mismatch / input_len))
}
```

## The results

```r
kable(results)
```



|case     | mismatch| input_len| mismatch_rate|
|:--------|--------:|---------:|-------------:|
|ur1_b_1  |        1|       500|         0.002|
|ur2_b_1  |        0|       500|         0.000|
|ur3_b_1  |      220|       500|         0.440|
|ur1_b1_1 |        0|       500|         0.000|
|ur2_b1_1 |        0|       500|         0.000|
|ur3_b1_1 |      220|       500|         0.440|

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