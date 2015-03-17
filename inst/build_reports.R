library(knitr)
args <- commandArgs(TRUE)
if (length(args) < 2) stop("Bad args, usage report_script_name input_data_file_name")

report_script_name <- args[1]
input_data_file_name <- args[2]
load(input_data_file_name)

knit2html(report_script_name)
