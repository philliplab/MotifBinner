#' Given a test bin of sequences, benchmark how well the input consensus
#' sequence can be reconstructed
#' @param test_bin The test bin to benchmark
#' @param classification_technique The technique to use to search for
#' mislabelled sequences
#' @param classification_params The parameters for the classification technique
#' @param alignment_technique The alignment technique to use
#' @param alignment_params The alignment parameters to use
#' @param consensus_technique The technique to use when generating the
#' consensus sequence
#' @param consensus_params The parameters for the consensus generator
#' @export

score_consensus <- function(test_bin, 
                            classification_technique = 'infovar_balance',
                            classification_params = list(threshold = 1),
                            alignment_technique = 'muscle',
                            alignment_params = list(),
                            consensus_technique = 'Biostrings::consensusString',
                            consensus_params = list()){
  input_seq <- test_bin$true_consensus
  result <- process_bin(list(src = test_bin$src, out = test_bin$out),
                        classification_technique,
                        classification_params,
                        alignment_technique,
                        alignment_params,
                        consensus_technique,
                        consensus_params)
  if (length(result) == 0){
    result <- list(a = paste(rep('-', nchar(input_seq)), sep="", collapse=""))
  }
  aligned <- pairwiseAlignment(input_seq, result[[1]])
  return(list(result = result,
              input_seq = input_seq, 
              score = score(aligned),
              edit_dist = nedit(aligned),
              alignment = aligned))
}
