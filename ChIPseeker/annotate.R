#!/usr/bin/env Rscript

# ~~~~~ PACKAGES ~~~~~ # 
library("optparse")



# ~~~~~ FUNCTIONS ~~~~~ # 
msprintf <- function(fmt, ...) {
    message(sprintf(fmt, ...))
}

make_annotation_filename <- function (input_file, suffix = '_annotations.tsv') {
    # Convert '/path/to/file.bed' to '/path/to/file_annotations.tsv'
    return(file.path(dirname(input_file), gsub(pattern = '.bed', replacement = suffix, x = basename(input_file))))
}

check_numlines <- function(input_file, min_value = 0) {
    # make sure a file has >0 lines
    has_enough_lines <- FALSE
    if (length(readLines(input_file)) > min_value) has_enough_lines <- TRUE
    return(has_enough_lines)
}

validate_file <- function(input_file) {
    # make sure that all files are .bed, and that they have >0 lines
    # validation passes if all files are .bed
    all_exist <- all(file.exists(input_file))
    if ( ! isTRUE(all_exist)) {
        msprintf("WARNING: Input file do not exist:\n%s\nFile will not be processed\n\n", input_file)
        return(FALSE)
    }
    all_bed <- all(grepl(pattern = '*.bed$', x = basename(input_file)))
    if ( ! isTRUE(all_bed)) {
        msprintf("WARNING: Input file is not .bed:\n%s\nFile will not be processed\n\n", input_file)
        return(FALSE)
    }
    all_min_linenum <- all(sapply(input_file, check_numlines))
    if ( ! isTRUE(all_min_linenum)) {
        msprintf("WARNING: Input file does not have enough lines:\n%s\nFile will not be processed\n\n", input_file)
        return(FALSE)
    }
    return(TRUE)
}

annotate_beds <- function(bed_files, suffix, tss_dist = 3000, annoDb = "org.Hs.eg.db") {
    # annotate all .bed files with ChIPSeeker 
    
    
    # ~~~~~ VALIDATION ~~~~~ # 
    # check to make sure at least one files has >0 lines before we try to load data, because it takes a while to load
    any_min_linenum <- any(sapply(names(bed_files), check_numlines))
    if ( ! isTRUE(any_min_linenum)) {
        msprintf("ERROR: No input files have enough lines to be processed\nExiting...\n\n")
        quit()
    }
    
    # ~~~~~ LOAD DATA ~~~~~ # 
    message("\nLoading packages and data...\n")
    # source("http://bioconductor.org/biocLite.R")
    # biocLite("ChIPseeker")
    suppressPackageStartupMessages(library("ChIPseeker"))
    suppressPackageStartupMessages(library("clusterProfiler"))
    suppressPackageStartupMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene"))
    txdb <- get("TxDb.Hsapiens.UCSC.hg19.knownGene")
    
    # ~~~~~ RUN ~~~~~ # 
    # iterate over bed files
    msprintf('\n------------------------------\n')
    msprintf('\n------------------------------\n')
    for(i in seq_along(bed_files)){
        bed_file <- names(bed_files[i])
        process_file <- bed_files[i] # TRUE or FALSE
        output_file <- make_annotation_filename(input_file = bed_file, suffix = suffix)
        
        msprintf("Input File:\n%s\n\n\nFile will be processed:\n%s\n\n", bed_file, process_file)
        if(isTRUE(as.logical(process_file))){
            
            msprintf("Reading peaks file...\n\n")
            peak <- readPeakFile(bed_file)
            
            msprintf("Getting peak annotations...\n\n")
            peakAnno <- annotatePeak(peak, tssRegion = c(-tss_dist, tss_dist), 
                                     TxDb = txdb, 
                                     annoDb = annoDb)
            
            msprintf("Saving tables to file:\n%s\n\n", output_file)
            write.table(peakAnno, quote=FALSE, sep="\t", row.names =FALSE, file=output_file)
        }
        msprintf('\n------------------------------\n')
    }
}


find_all_beds <- function (input_dirs) {
    # find all .bed files in the supplied dirs
    return(dir(input_dirs, pattern = '.bed', full.names = TRUE, recursive = TRUE))
}



# ~~~~~ SCRIPT ARGS ~~~~~ # 
option_list <- list(
    make_option(c("-d", "--dir"), action="store_true", default=FALSE,
                dest="dir_mode", help="Directories with peaks to annotate"),
    # make_option(c("--genome"), type="character", default="hg19",
    #             dest = "genome", help="Genome version to use for annotations [default %default]",
    #             metavar="genome"),
    make_option(c("--suffix"), type="character", default='_annotations.tsv',
                dest = "suffix", help="File suffix to append [default %default]",
                metavar="suffix"),
    make_option(c("--tss-dist"), type="numeric", default=3000,
                dest = "tss_dist", help="TSS distance to use [default %default]",
                metavar="tss-dist")
)
opt <- parse_args(OptionParser(option_list=option_list), positional_arguments = TRUE)

# genome <- opt$options$genome
dir_mode <- opt$options$dir_mode 
suffix <- opt$options$suffix
tss_dist <- opt$options$tss_dist

input_items <- opt$args

# get script dir
args <- commandArgs(trailingOnly = F)  
scriptPath <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))
save.image(file.path(scriptPath, "loaded_args.Rdata"))



# ~~~~~ RUN ~~~~~ # 
if (isTRUE(dir_mode)) input_items <- find_all_beds(input_items)

validated_items <- sapply(input_items, validate_file)

annotate_beds(bed_files = validated_items, suffix = suffix, tss_dist = tss_dist)
