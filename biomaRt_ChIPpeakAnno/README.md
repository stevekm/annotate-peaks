Annotate peaks using the R packages `biomaRt` and `ChIPpeakAnno`. 

# Usage
Annotate all supplied `.bed` files:
```
./annotate.R test1.bed test2.bed
```

Annotate all `.bed` files in a directory:
```
./annotate.R bed_dir_1 bed_dir_2 -d
```

Specify a different genome to use:
```
./annotate.R test1.bed test2.bed --genome hg19
```
- __NOTE:__ Only `hg19` is currently supported.
- __NOTE:__ May output duplicate annotations for different gene isoforms

## Example

```
$ ./annotate.R test_beds/test1.bed test_beds/test3.bed
WARNING: Input file does not have enough lines:
test_beds/test3.bed
File will not be processed



Loading packages for annotation...

Looking for previously saved biomaRt data for hg19 in location:
/Users/kellys04/projects/annotate-peaks/data/hg19/biomart_data.RData


Found biomaRt data file:
/Users/kellys04/projects/annotate-peaks/data/hg19/biomart_data.RData
Loading data from file...



------------------------------


------------------------------

Input File:
test_beds/test1.bed


File will be processed:
TRUE


Reading in the BED file...


duplicated or NA names found. Rename all the names by numbers.

Getting annotations...


Merging annotations...


Saving the output to file:
test_beds/test1_annotations.tsv



------------------------------

Input File:
test_beds/test3.bed


File will be processed:
FALSE



------------------------------
```

# Software Requirements
Tested with:
- R version 3.3.0
- [biomaRt_2.28.0](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)
- [ChIPpeakAnno_3.6.4](https://bioconductor.org/packages/release/bioc/html/ChIPpeakAnno.html)
- optparse_1.3.2
