#!/usr/bin/env Rscript

# ~~~~~ PACKAGES ~~~~~ # 
library("optparse")



# ~~~~~ FUNCTIONS ~~~~~ # 
msprintf <- function(fmt, ...) {
    message(sprintf(fmt, ...))
}

make_annotation_filename <- function (input_file) {
    # Convert '/path/to/file.bed' to '/path/to/file_annotations.tsv'
    return(file.path(dirname(input_file), gsub(pattern = '.bed', replacement = '_annotations.tsv', x = basename(input_file))))
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

annotate_beds <- function(bed_files, genome, scriptPath, biodat = FALSE) {
    # annotate all .bed files
    
    # ~~~~~ GENOME ACCESS DATA ~~~~~ # 
    # items needed to access the genome data from the biomart servers
    genome_key <- list("hg19" = list(
        host = "grch37.ensembl.org", 
        biomart = "ENSEMBL_MART_ENSEMBL", 
        dataset = "hsapiens_gene_ensembl"
    )
    # , # add more genomes here !! 
    # "mm10" = list(
    #     
    #     biomart = "ENSEMBL_MART_ENSEMBL", 
    #     dataset = "mmusculus_gene_ensembl"
    # )
    )
    
    # ~~~~~ VALIDATION ~~~~~ # 
    # check to make sure at least one files has >0 lines before we try to load data, because it takes a while to load
    if ( ! genome %in% names(genome_key)) {
        msprintf("ERROR: Supplied genome '%s' is not recognized by this script\n", genome)
        msprintf("Genome options are:\n")
        msprintf(names(genome_key))
        quit()
    }
    any_min_linenum <- any(sapply(names(bed_files), check_numlines))
    if ( ! isTRUE(any_min_linenum)) {
        msprintf("ERROR: No input files have enough lines to be processed\nExiting...\n\n")
        quit()
    }
    
    # ~~~~~ LOAD DATA ~~~~~ # 
    message("\nLoading packages for annotation...\n")
    # source("https://bioconductor.org/biocLite.R")
    # biocLite("ChIPpeakAnno")
    suppressPackageStartupMessages(library("ChIPpeakAnno"))
    suppressPackageStartupMessages(library("biomaRt"))
    
    if(biodat != FALSE){
        biomart_data_file <- as.character(biodat)
    } else {
        biomart_data_file <- file.path(scriptPath, "data", genome, "biomart_data.RData") 
    }

    msprintf("Looking for previously saved biomaRt data for %s in location:\n%s\n\n", genome, biomart_data_file)
    if(file.exists(biomart_data_file)){
        msprintf("Found biomaRt data file:\n%s\nLoading data from file...\n\n", biomart_data_file)
        load(biomart_data_file)
    } else {
        msprintf("Saved biomaRt data file not found!")
        msprintf("Retreiving reference information for %s from biomaRt, this might take a few minutes...", genome)
        martEns <- useMart(host = genome_key[[genome]][["host"]], 
                           biomart = genome_key[[genome]][["biomart"]], 
                           dataset = genome_key[[genome]][["dataset"]], verbose=F)
        martEnsTSS <- getAnnotation(mart=martEns, featureType="TSS")
        martEnsDF <- getBM(attributes=c("ensembl_gene_id", "external_gene_name", "gene_biotype"), mart=martEns)
        msprintf("Successfully retrieved genome data from biomaRt")
        msprintf("Saving a copy of the biomaRt data to file:\n%s\n", biomart_data_file)
        save(martEns, martEnsTSS, martEnsDF, file = biomart_data_file)
    }
    
    # ~~~~~ RUN ~~~~~ # 
    # iterate over bed files
    msprintf('\n------------------------------\n')
    msprintf('\n------------------------------\n')
    for(i in seq_along(bed_files)){
        bed_file <- names(bed_files[i])
        process_file <- bed_files[i] # TRUE or FALSE
        output_file <- make_annotation_filename(bed_file)
        msprintf("Input File:\n%s\n\n\nFile will be processed:\n%s\n\n", bed_file, process_file)
        if(isTRUE(as.logical(process_file))){
            msprintf("Reading in the BED file...\n\n")
            peaks_granges <- toGRanges(bed_file, format="BED", header=FALSE)
            
            # get the annotations
            msprintf("\nGetting annotations...\n")
            peaks_granges <- annotatePeakInBatch(peaks_granges, AnnotationData = martEnsTSS, PeakLocForDistance = "middle", FeatureLocForDistance = "TSS", output = "shortestDistance", multiple = TRUE)
            
            # merge the annotations with the peaks
            msprintf("\nMerging annotations...\n")
            peaks_granges_df <- merge(as.data.frame(peaks_granges) , martEnsDF , by.x=c("feature"), by.y=c("ensembl_gene_id") , all.x=TRUE)
            
            # save the output
            msprintf("\nSaving the output to file:\n%s\n\n", output_file)
            write.table(peaks_granges_df, row.names = FALSE, sep = '\t', quote = FALSE,
                        file = output_file)
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
    make_option(c("--genome"), type="character", default="hg19",
                dest = "genome", help="Genome version to use for annotations [default %default]",
                metavar="genome"),
    make_option(c("-b", "--biodat"), type="character", default=FALSE,
                dest = "biodat", help="Path to the biomaRt data file to use for annotations [default %default]",
                metavar="biodat")
)
opt <- parse_args(OptionParser(option_list=option_list), positional_arguments = TRUE)

genome <- opt$options$genome
dir_mode <- opt$options$dir_mode 
biodat <- opt$options$biodat
input_items <- opt$args

# get script dir
args <- commandArgs(trailingOnly = F)  
scriptPath <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))
save.image(file.path(scriptPath, "loaded_args.Rdata"))



# ~~~~~ RUN ~~~~~ # 
if (isTRUE(dir_mode)) input_items <- find_all_beds(input_items)

validated_items <- sapply(input_items, validate_file)

annotate_beds(bed_files = validated_items, genome =  genome, scriptPath = scriptPath, biodat = biodat)
