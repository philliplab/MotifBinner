# Given some input sequences the functions in this script will produce
# simulated reads from those sequences

#' Generates a read error profile
#' @param technique What technique should be used to generate the error
#' profile?
#' @param params A list of parameters to pass to the specific profile generator
#' @export

gen_error_profile <- function(technique = 'uniform', 
                              params = list(read_len = 1000,
                                            A_mut_rates = list('A' = 0.997, 'C' = 0.001,
                                                               'G' = 0.001, 'T' = 0.001),
                                            C_mut_rates = list('C' = 0.997, 'A' = 0.001,
                                                               'G' = 0.001, 'T' = 0.001),
                                            G_mut_rates = list('G' = 0.997, 'C' = 0.001,
                                                               'A' = 0.001, 'T' = 0.001),
                                            T_mut_rates = list('T' = 0.997, 'C' = 0.001,
                                                               'G' = 0.001, 'T' = 0.001))){
  if (technique == 'uniform'){
    err_prof <- do.call(gen_error_profile_uniform, params)
  }
  return(err_prof)
}

