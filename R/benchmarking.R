#' Given all the test cases list the unique scenarios
#' @param cases list of lists of the cases
#' @export

list_unique_scenarios <- function(cases){

  unique_scenarios <- list()
  
  for (tc in names(cases)){
    params <- cases[[tc]]$scenario
    params[['seed']] <- cases[[tc]]$seed
    hash <- digest(params)
    unique_scenarios[[hash]] <- params
  }
  return(unique_scenarios)
}

#' Given the location of a cache file, intialize the cache
#' @param cache_file The absolute path to the cache file
#' @export

load_or_initialize_cache <- function(cache_file){
  if (file.exists(cache_file)){
    load(cache_file)
  } else {
    scenario_cache <- list()
  }
  return(scenario_cache)
}

#' Compute the scenarios that are not already in the scenario cache
#' @param unique_scenarios A list of the unique scenarios as produced by
#' list_unique_scenarios
#' @param scenario_cache The scenario cache list
#' @export

create_scenario_data <- function(unique_scenarios, scenario_cache){
  x <- foreach (hash = names(unique_scenarios)) %dopar% {
    if (hash %in% names(scenario_cache)){
      result <- list(res = scenario_cache[[hash]], hash = hash)
    } else {
      params <- unique_scenarios[[hash]]
      params$name <- NULL
      res <- do.call( gen_and_contaminate_reads, params)
      result <- list(res = res, hash = hash)
    }
    result
  }

  for (i in seq_along(x)){
    hash <- x[[i]]$hash
    res <- x[[i]]$res
    scenario_cache[[hash]] <- res
  }
  return(scenario_cache)
}

#' Run a test case
#' @param scenario The computed scenario to run the test on
#' @param seed The seed that was used to generate the scenario
#' @param setup A list containing parameters that will control how the bin is
#' processed
#' @param scenario_cache The cache with the pre computed scenarios
#' @export

run_test <- function(scenario, seed, setup, scenario_cache){
  params <- scenario
  params[['seed']] <- seed

  test_bin <- scenario_cache[[digest(params)]]

  params <- setup
  params$test_bin <- test_bin
  result <- do.call(score_consensus, params)
  return(list(mismatch = result$edit_dist,
              output_len = nchar(result$alignment)))
}

#' Runs all the tests
#' @param cases All the test cases to run
#' @param scenario_cache The cache with the pre computed scenarios
#' @export

run_all_tests <- function(cases, scenario_cache){
  results <- data.frame(setup = character(0),
                        scenario = character(0),
                        seed = character(0),
                        mismatch = numeric(0),
                        input_len = numeric(0),
                        output_len = numeric(0),
                        mismatch_rate = numeric(0))
  
  for (tc in names(cases)){
    print(tc)
    split_names <- strsplit(tc, split='-')[[1]]
    setup <- split_names[1]
    scenario <- split_names[2]
    seed <- split_names[3]
    params <- cases[[tc]]
    params[['scenario_cache']] <- scenario_cache
    test_result <- do.call(run_test, params)
    output_len <- test_result$output_len
    mismatch <- test_result$mismatch
    input_len <- nchar(cases[[tc]][['scenario']][['ref_seq']])
    results <- rbind(results,
                     data.frame(setup = setup,
                                scenario = scenario,
                                seed = seed,
                                mismatch = mismatch,
                                input_len = input_len,
                                output_len = output_len,
                                mismatch_rate = mismatch / input_len))
  }
  return(results)
}
