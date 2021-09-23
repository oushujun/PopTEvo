## Files explained
family_sum_combined_counts_allNAM_15Jun21.txt.gz, raw count table for all TE families. Read counts for all genes is summarized in the entry `Gene`.

## RNAseq mapping and counts

### Reads 

The RNAseq data was softlinked, creating respecitve NAM line folders using [`linker-new.sh`](scripts/linker-new.sh) script.

```bash
./linker-new.sh
```

### Genomes

A file listing absolute paths for the NAM genomes was created using the `find` comamnd ([`genomes.fofn`](scripts/genomes.fofn))

```bash
find $(pwd -P) -name "*pseudomolecules.fasta" > genomes.fofn
```

### Index

Indexing was performed using the [`hisat2-index.sh`](scripts/hisat2-index.sh), which was submitted as a slurm job [`index.sub`](scripts/index.sub).

```bash
sbatch --array=0-26 index.sub
```

### Mapping

Mapping was performed using the [`runHISAT2.sh`](scripts/runHISAT2.sh), which was submitted as a slurm job using [`hisat2.sub`](scripts/hisat2.sub) as template and [`job-submit.sh`](scripts/job-submit.sh) as a submitter script.

```bash
while read line; do ./job-submit.sh $line; done<genomes.fofn
```

### From BAM to raw count table
Follow scripts and protocal here: https://github.com/SNAnderson/maizeTEexpression/tree/master/NAM_2021
