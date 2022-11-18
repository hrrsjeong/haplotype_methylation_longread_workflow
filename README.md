# haplotype-specific methylation analysis using long-read sequencing data
## Description

This is a modular Snakemake workflow for processing long-read sequencing data (pacbio HiFi or ONT) to generate haplotype assigned methylation and to perform differential methylation analysis. To run this pipeline, it requires phased genome assemblies (`hap1` and `hap2`; doesn't have to be trio phased) and `modbam` files including methylation with `MM` and `ML` tags.

## Loading in **Snakemake**

try to run the pipeline something like:

```
mkdir -p log; snakemake --ri --jobname "{rulename}.{jobid}" \
--drmaa " -l centos=7 -V -cwd -j y -o ./log -e ./log -l h_rt={resources.hrs}:00:00 \
 -l mfree={resources.mem}G -pe serial {threads} -w n -S /bin/bash" -w 60 \
 -k -j 180 -p --resources load=100 --resources pbcpg=300 --reason
```
