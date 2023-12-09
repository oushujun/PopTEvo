## Splitting tandem and singleton genes

Tandem and singleton genes are determined using different gene annotation data in `data/gene_info` including tandem gene annotations from the [Hufford et al. (2021)](https://www.science.org/doi/10.1126/science.abg5289) NAM paper [pan_gene_matrix_v3_cyverse.csv](https://datacommons.cyverse.org/browse/iplant/home/shared/NAM/NAM_genome_and_annotation_Jan2021_release/SUPPLEMENTAL_DATA/pangene-files/pan_gene_matrix_v3_cyverse.csv). The set of NLR genes `data/gene_info/nlrannotator_nlrome_unique.bed` is the non-redundant combination of the NLR-annotator set and the NLRome set described in `1_nlrannotator`. For tandem genes, a physical clustering threshold in bp is used to detect whether tandem duplicates are clustered. A python script ouputs bed files with the various gene types.

```
scripts/find_tandem_dup.py singleton nonnlr 100000 > singleton_nonnlr.bed
scripts/find_tandem_dup.py duplicate nonnlr 100000 > duplicate_nonnlr.bed
scripts/find_tandem_dup.py singleton nlr 100000 > singleton_nlr.bed
scripts/find_tandem_dup.py duplicate nlr 100000 > duplicate_nlr.bed
```

Files can then be split by NAM line (column 5) using awk.

```
awk '{print > $5}' singleton_nonnlr.bed
```

The resulting files per NAM were renamed and are provided in the relevant directories:

```
nam_singleton_nlr/
nam_singleton_nonnlr/
nam_tandem_nlr/
nam_tandem_nonnlr/
```

