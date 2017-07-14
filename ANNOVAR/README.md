Annotate peaks with [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/) using `refGene` database.

# Install ANNOVAR

Install ANNOVAR using the `install.sh` script provided:

```bash
./install.sh
```

By default, it will install binaries to `~/annovar` and the databases to `~/annovar/db`. To change this, add arguments like this:

```bash
install.sh /path/to/bin_dir /path/to/db_dir
```

# Usage
Annotate all supplied `.bed` files:
```
./annotate.R test1.bed test2.bed
```

# Software Requirements
Tested with:
- R version 3.3.0
- ANNOVAR version: 2017-06-01 23:08:16 -0400 (Thu,  1 Jun 2017)


# Example 

```bash
$ ./annotate.R example-data/Sample1.bed

------------------------------


------------------------------

Input File:
example-data/Sample1.bed


File will be processed:
TRUE


Input file:

example-data/Sample1.avinput


Output file:

example-data/Sample1.hg19_multianno.txt


ANNOVAR command is:


perl "/ifs/home/kellys04/annovar/annovar/table_annovar.pl" "example-data/Sample1.avinput" "/ifs/home/kellys04/annovar/db/hg19" --outfile "example-data/Sample1" --buildver "hg19" --protocol "cytoBand,refGene" --operation "r,g" --nastring "." --remove



-----------------------------------------------------------------
NOTICE: Processing operation=r protocol=cytoBand

NOTICE: Running with system command <annotate_variation.pl -regionanno -dbtype cytoBand -buildver hg19 -outfile example-data/Sample1 example-data/Sample1.avinput /ifs/home/kellys04/annovar/db/hg19>
NOTICE: Output file is written to example-data/Sample1.hg19_cytoBand
NOTICE: Reading annotation database /ifs/home/kellys04/annovar/db/hg19/hg19_cytoBand.txt ... Done with 862 regions
NOTICE: Finished region-based annotation on 5000 genetic variants
-----------------------------------------------------------------
NOTICE: Processing operation=g protocol=refGene

NOTICE: Running with system command <annotate_variation.pl -geneanno -buildver hg19 -dbtype refGene -outfile example-data/Sample1.refGene -exonsort example-data/Sample1.avinput /ifs/home/kellys04/annovar/db/hg19>
NOTICE: Output files were written to example-data/Sample1.refGene.variant_function, example-data/Sample1.refGene.exonic_variant_function
NOTICE: Reading gene annotation from /ifs/home/kellys04/annovar/db/hg19/hg19_refGene.txt ... Done with 63481 transcripts (including 15216 without coding sequence annotation) for 27720 unique genes
NOTICE: Processing next batch with 5000 unique variants in 5000 input lines
-----------------------------------------------------------------
NOTICE: Multianno output file is written to example-data/Sample1.hg19_multianno.txt

------------------------------
```
