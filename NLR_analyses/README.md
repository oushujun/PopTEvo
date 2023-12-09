# Analysis of LTR neighborhoods of NLR genes 

This analysis is composed of multiple steps, each documented in a separate subdirectory with `data` and `scripts` directories.

1) `nlrannotator`: Annotating NLR genes in the NAM genomes
2) `ltr`: LTR annotation for each NAM line used for the analysis
3) `clustering`: Clustering of NLR genes and background genes based on physical distance
4) `nlr_neighborhood_age`: Analysis of young and old LTRs in the neighborhoods of foreground and background gene sets
5) `nlr_neighborhood_amplified`: Analysis of tropical-amplifying and not tropical-amplifying LTRs in the neighborhoods of foreground and background gene sets in tropical and temperate NAM lines
6) `plots`: Pre-computed results tables and plotting script for visualising the main results

Note that all large files in these subdirectories are compressed and need to be uncompressed using `gzip` for reproducing most analyses.

## Software dependencies

The scripts used for this analysis rely on bash, python and R as well as `bedtools`.

### Python libraries

```
scipy
```

### R libraries 

```
ggplot2
reshape2
dplyr
cowplot
```
