#!/usr/bin/env Rscript

# ~~~~~ PACKAGES ~~~~~ # 
library("optparse")
library("tools")



# ~~~~~ FUNCTIONS ~~~~~ # 
msprintf <- function(fmt, ...) {
    message(sprintf(fmt, ...))
}

remove_ext <- function(input_file){
    old_ext <- file_ext(input_file)
    filename_base <- gsub(pattern = sprintf('.%s$', old_ext), replacement = '', x = basename(input_file))
    return(filename_base)
}

make_filename <- function (input_file, new_ext) {
    # Convert '/path/to/file.bed' to '/path/to/file_annotations.tsv'
    old_ext <- file_ext(input_file)
    filename_base <- gsub(pattern = sprintf('.%s$', old_ext), replacement = '', x = basename(input_file))
    filename_new <- sprintf('%s.%s', filename_base, new_ext)
    return(file.path(dirname(input_file), filename_new))
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



find_all_beds <- function (input_dirs) {
    # find all .bed files in the supplied dirs
    return(dir(input_dirs, pattern = '.bed', full.names = TRUE, recursive = TRUE))
}



sort_bed_df <- function(df){
    # sort the bed and remove duplicate entries
    df <- df[1:3]
    df <- df[ order(df[,1], df[,2]), ]
    df <- df[! duplicated(df), ]
    return(df)
}

read_bed <- function(bedfile){
    # read a bedfile into a dataframe
    df <- read.delim(file = bedfile, header = FALSE, sep = '\t', stringsAsFactors = FALSE)
    return(df)
}

write_bed <- function(df, output_file){
    write.table(x = df, file = output_file, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
}


default_bin_dir <- function(){
    user_home <- system(command = "echo $HOME", intern = TRUE)
    bin_dir <- file.path(user_home, "annovar")
    return(bin_dir)
}

default_db_dir <- function(){
    db_dir <- file.path(default_bin_dir(), "db")
    return(db_dir)
}

annotate_beds <- function(bed_files, genome, bin_dir, db_dir) {
    # ~~~~~ VALIDATION ~~~~~ # 
    # check to make sure at least one files has >0 lines before we try to load data, because it takes a while to load
    any_min_linenum <- any(sapply(names(bed_files), check_numlines))
    if ( ! isTRUE(any_min_linenum)) {
        msprintf("ERROR: No input files have enough lines to be processed\nExiting...\n\n")
        quit()
    }
    
    # ~~~~~ RUN ~~~~~ # 
    # iterate over bed files
    msprintf('\n------------------------------\n')
    msprintf('\n------------------------------\n')
    for(i in seq_along(bed_files)){
        bed_file <- names(bed_files[i])
        process_file <- bed_files[i] # TRUE or FALSE
        
        msprintf("Input File:\n%s\n\n\nFile will be processed:\n%s\n\n", bed_file, process_file)
        if(isTRUE(as.logical(process_file))){
            
            # read in the .bed file
            bed_df <- read_bed(bedfile = bed_file)
            
            # sort & uniq the bed
            bed_df <- sort_bed_df(df = bed_df)
            
            # add blank columns for ANNOVAR input format
            bed_df[[4]] <- 0
            bed_df[[5]] <- 0
            
            # write the sorted bed ANNOVAR input file
            avinput_file <- make_filename(input_file = bed_file, new_ext = "avinput")
            write_bed(df = bed_df, output_file = avinput_file)
            
            # annotate with ANNOVAR
            multianno_output <- table_annovar(avinput = avinput_file, bin_dir = bin_dir, db_dir = db_dir, buildver = genome)
            
            # annotations_df <- read_multianno(input_file = multianno_output)
            # return(annotations_df)
        }
        msprintf('\n------------------------------\n')
    }
}


table_annovar <- function(avinput, bin_dir = default_bin_dir(), db_dir = default_db_dir(), buildver = "hg19"){
    protocol <- "cytoBand,refGene"
    operation <- "r,g"
    
    output_prefix <- remove_ext(input_file = avinput)
    output_file <- file.path(dirname(avinput), output_prefix) # example-data/Sample1
    
    multianno_output_suffix <- sprintf("%s_multianno.txt", buildver)
    multianno_output <- make_filename(input_file = avinput, new_ext = multianno_output_suffix)
    
    
    genome_db_dir = file.path(db_dir, buildver)
    
    annovar_command = sprintf('
perl "%s/annovar/table_annovar.pl" "%s" "%s" --outfile "%s" --buildver "%s" --protocol "%s" --operation "%s" --nastring "%s" --remove
', bin_dir, avinput, genome_db_dir, output_file, buildver, protocol, operation, na_string
    )
    msprintf("Input file:\n\n%s\n\n", avinput) # example-data/Sample1.bed
    msprintf("Output file:\n\n%s\n\n", multianno_output) # example-data/Sample1.hg19_multianno.txt
    msprintf("ANNOVAR command is:\n\n%s\n\n", annovar_command)
    system(command = annovar_command)
    return(multianno_output)
}


read_multianno <- function(input_file){
    df <- read.delim(file = input_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE, na.strings = na_string)
    return(df)
}

# ~~~~~ GLOBALS ~~~~~ # 
na_string <- "."

# ~~~~~ SCRIPT ARGS ~~~~~ # 
option_list <- list(
    make_option(c("-d", "--dir"), action="store_true", default=FALSE,
                dest="dir_mode", help="Directories with peaks to annotate"),
    make_option(c("--bin-dir"), type="character", default=default_bin_dir(),
                dest = "bin_dir", help="Default directory for ANNOVAR installed binaries [default %default]",
                metavar="bin-dir"),
    make_option(c("--db-dir"), type="character", default=default_db_dir(),
                dest = "db_dir", help="Default directory for ANNOVAR databases [default %default]",
                metavar="db-dir"),
    make_option(c("--genome"), type="character", default="hg19",
                dest = "genome", help="Genome version to use for annotations [default %default]",
                metavar="genome")
)
opt <- parse_args(OptionParser(option_list=option_list), positional_arguments = TRUE)

genome <- opt$options$genome
dir_mode <- opt$options$dir_mode 
bin_dir <- opt$options$bin_dir
db_dir <- opt$options$db_dir

input_items <- opt$args

# get script dir
args <- commandArgs(trailingOnly = F)  
scriptPath <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))
save.image(file.path(scriptPath, "loaded_args.Rdata"))




# ~~~~~ RUN ~~~~~ # 
if (isTRUE(dir_mode)) input_items <- find_all_beds(input_items)

validated_items <- sapply(input_items, validate_file)

annotate_beds(bed_files = validated_items, bin_dir = bin_dir, db_dir = db_dir, genome =  genome)
# save.image()