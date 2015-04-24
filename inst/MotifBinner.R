#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("MotifBinner"))

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to
# make_option(c("-h", "--help"), action="store_true", default=FALSE,
# help="Show this help message and exit")

option_list <- list(
make_option("--file_name", 
            help = "The file name"),
make_option("--output", 
            help = "The directory where the output must be stored"),
make_option("--prefix", 
            help = "The prefix that is used to identify the motif"),
make_option("--suffix", 
            help = "The suffix that is used to identify the motif"),
make_option("--motif_length", 
            type = 'numeric',
            help = "The length of the motif that forms the pid"),
make_option("--max_mismatch_start", 
            dest = 'max.mismatch_start', type = 'integer', default = 0,
            help = "Maximum number of mismatches are iteratively increased when search - start at this number"),
make_option("--max_mismatch", 
            dest = 'max.mismatch', type = 'integer',
            help = "The maximum number of mismatches to allow when searching for the pid - iteratively incease the number of mismatches allowed until this number is reached."),
make_option("--threshold", 
            type = 'numeric',
            help = "Outlier sequences are removed from the bin until the maximum distance between any two sequences drops below this threshold."),
make_option("--start_threshold", 
            type = 'numeric',
            help = "Only start the classification if the maximum distance between and two sequences in the bin is greater than this."),
make_option("--max_sequences", 
            type = 'numeric', default = 400,
            help = "Warning: values larger than 500 might incur the wrath of the sysadmin. 400 is a good middle ground. The maximum number of sequences to use for the computation of the distance matrix. If more sequences than this is present, then randomly select this many sequences and run the classification algorithm on them. This is only to improve the computation speed."),
make_option("--dont_remove_gaps", 
            dest = 'remove_gaps', action = 'store_false', default = TRUE,
            help = "use this option to prevent gaps from being removed from the consensus sequences"),
make_option("--dont_strip_uids", 
            dest = 'strip_uids', action = 'store_false', default = TRUE,
            help = "Use this option to prevent the removal of the unique identifiers from the final consensus sequences. It is not intelligent. The names will be split on '_' and the first and last pieces will be kept."),
make_option("--n_bins_to_process", 
            default = 0, type = 'integer',
            help = "The number of bins to process through the outlier detection, alignment and consensus generation. If smaller than or equal to 0, all bins will be processed."),
make_option("--verbose", 
            default = FALSE, action = 'store_true',
            help = "Adds fancy progress reporting via printouts to the screen and creation of files in the R'sessions tempory folder.")
                      )

opt <- parse_args(OptionParser(option_list=option_list))

opt$help <- NULL
print(opt)
dput(opt)
x <- do.call(process_file, opt)
rm(x)

# ./MotifBinner.R --file_name=~/projects/MotifBinner/data/CAP177_2040_v1merged.fastq --output=/tmp/MotifBinner --prefix=CCAGCTGGTTATGCGATTCTMARGTG --suffix=CTGAGCGTGTGGCAAGGCCC --motif_length=9 --max_mismatch=5 --threshold=0.01333 --start_threshold=0.01333 --max_sequences=100

