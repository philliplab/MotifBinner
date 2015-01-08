#' Returns that dataset used for testing the mislabel detection
#' @export

get_mislabel_test_data <- function(){
  x <- list('test1' = list('src' = DNAStringSet(c('AAAAAAAAAAAAAAAAAAAAAAAAAAAAA',
                                               'AAAAAAAAAAAAAAAAAAAAAAAAAAAAA')),
                           'out' = DNAStringSet('CCCCCCCCCCCCCCCCCCCCCCCCCCCCC')
                           )
  )
  return(x)
}
