Annotate peaks using the R package `ChIPSeeker`

R script to run the ChIPSeeker peak annotation and summary plotting pipeline on .bed files. 

For a full `ChIPSeeker` pipeline with plot and peak type output, see [peak-type-summary](https://github.com/stevekm/peak-type-summary)

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

## Options

- `-d`, `--dir`: Dir mode; treat input items as directories in which to search for .bed files
- `--tss-dist`: TSS region distance
- `--suffix`: Suffix to use for the annotation file output

# Examples

Annotate & plot .bed files

```bash
$ ./annotate.R example-data/Sample1.bed

Loading packages and data...


------------------------------


------------------------------

Input File:
example-data/Sample1.bed


File will be processed:
TRUE


Reading peaks file...


Getting peak annotations...


>> preparing features information...		 2017-07-17 14:13:34
>> identifying nearest features...		 2017-07-17 14:13:35
>> calculating distance from peak to TSS...	 2017-07-17 14:13:36
>> assigning genomic annotation...		 2017-07-17 14:13:36
>> adding gene annotation...			 2017-07-17 14:14:17
Loading required package: org.Hs.eg.db

'select()' returned many:many mapping between keys and columns
>> assigning chromosome lengths			 2017-07-17 14:14:19
>> done...					 2017-07-17 14:14:19
Saving tables to file:
example-data/Sample1_annotations.tsv



------------------------------

```

# Software
- Tested with R version 3.2.3 and 3.3.0, with the following packages:
  - `ChIPseeker_1.6.7`
  - `clusterProfiler_2.4.3`
  - `TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2`
  - `optparse_1.3.2 `
