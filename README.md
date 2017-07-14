# annotate-peaks

Scripts to simplify the annotation of .bed formatted files containing peaks or genomic regions. 

See each directory's `README` file for descriptions of each script. 

These scripts are designed for easy usage with GNU `find`, when you have a lot of .bed files that you want to annotate at once. For example:

```bash
find /path/to/peaks_dir -name "*.bed" ! -path "*/dir_I_dont_like/*" | xargs ./annotate.R
```
