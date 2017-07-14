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
