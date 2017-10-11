#!/bin/bash

# This script will install ANNOVAR with the required databases
# http://annovar.openbioinformatics.org/en/latest/user-guide/download/

bin_dir="${1:-$HOME/annovar}" # default to home dir
db_dir="${2:-$bin_dir/db}"

printf "Bin dir:\n%s\n\n" "$bin_dir"
printf "Databse dir:\n%s\n\n" "$db_dir"


install () {
    local bin_dir="$bin_dir"
    local db_dir="$db_dir"
    
    mkdir -p "$bin_dir" && mkdir -p "$db_dir" && (
    
    cd "$bin_dir"
    
    # download ANNOVAR and extract 
    if [ ! -f annovar.latest.tar.gz ]; then
        wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz
        tar -zxvf annovar.latest.tar.gz
    fi

    # install ANNOVAR db's
    # hg19
    hg19_db_dir="${db_dir}/hg19/" 
    mkdir -p "$hg19_db_dir" 

    "${bin_dir}/annovar/annotate_variation.pl" -downdb -buildver hg19 -webfrom annovar refGene "$hg19_db_dir"
    "${bin_dir}/annovar/annotate_variation.pl" -buildver hg19 -downdb cytoBand "$hg19_db_dir"

    # mm10
    mm10_db_dir="${db_dir}/mm10/" 
    mkdir -p "$mm10_db_dir" 

    "${bin_dir}/annovar/annotate_variation.pl" -downdb -buildver mm10 -webfrom annovar refGene "$mm10_db_dir"
    "${bin_dir}/annovar/annotate_variation.pl" -buildver mm10 -downdb cytoBand "$mm10_db_dir"
    
    ) || printf "ERROR: Could not make directory for installation:\n%s\n\nExiting..." "$bin_dir" "$db_dir" && exit 1

}

install "$bin_dir" "$db_dir"