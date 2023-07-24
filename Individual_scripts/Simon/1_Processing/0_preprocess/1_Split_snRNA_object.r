
library(tidyverse)
library(Seurat)
library(Signac)
library(qs)
library(optparse)

# Define the options
option_list <- list(
    make_option(c("-i", "--input"), type="character", default=NULL,
                            help="Path to the input file"),
    make_option(c("-o", "--output"), type="character", default=NULL,
                            help="Path to the output directory"),
    make_option(c("-s", "--sample_column"), type="character", default="Sample",
                            help="Name of the column containing the sample information"),
    make_option(c("-n", "--num_cpu"), type="integer", default=50,
                            help="Number of CPUs to use")
)

# Parse the command line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Set variable
sample_column <- opt$sample_column
out_folder <- opt$output

# Load the data
# Read the file based on its extension
input <- opt$input
if (grepl("\\.rds$", input)) {
    sn <- readRDS(input)
} else if (grepl("\\.qs$", input)) {
    sn <- qs::qread(input)
} else {
    stop("File extension not recognized")
}

# # save as qs object
# Split object and save by Sample
dir.create(paste0(out_folder, '/individual/'), recursive = TRUE)
for(sample_use in unique(sn@meta.data[[sample_column]])){
        message("Processing and Saving:", sample_use)
        cells_keep = sn@meta.data %>% filter(.data[[sample_column]] == sample_use) %>% rownames
        sn_sample = subset(sn, cells = cells_keep)
        qs::qsave(sn_sample, file = paste0(out_folder, '/individual/', sample_use, '.qs'), nthreads = opt$num_cpu)
}

message('Step1 Split snRNA objects Done!')