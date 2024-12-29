# NucFlag Paper
Code for NucFlag paper. WIP.


### Usage

```bash
snakemake -np -j 40 -s misasim/Snakefile \
--use-conda \
--rerun-triggers mtime \
--rerun-incomplete \
--conda-frontend conda \
--executor lsf \
--default-resources lsf_extra="'-M 100GB'" lsf_queue="epistasis_long" lsf_project="default"
```