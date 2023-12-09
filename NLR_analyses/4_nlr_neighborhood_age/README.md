# Analysing the LTR neighborhood of target gene clusters

An exploratory analysis of LTR age distribution in neighborhoods with k LTRs was conducted, where k is a number from 1 to 50 and indicates the number of closest LTRs to the target gene cluster selected (including those that overlap the region). All NAM lines as listed in `data/groups.txt` were analysed.

For example, we can analyse the k=10 LTR neighborhood of the NAM line CML69, with a background set of tandem non-NLR genes and a foreground set of tandem NLR genes as follows.

```
scripts/get_nlr_stats_k.sh CML69 10 ../3_clustering/data/nam_tandem_nonnlr/CML69_nonnlr_100kb_merged.bed ../3_clustering/data/nam_tandem_nlr/CML69_nlr_100kb_merged.bed tNLR_tnonNLR > summary_stats_age.txt
```

Summary statistics will be printed to `summary_stats_age.txt` with the following columns.

```
Col1: Group name
Col2: Number of LTRs in neighborhood (k)
Col3: NAM line
Col4: Background young divided by background old
Col5: Foreground young divided by foreground old
Col6: Count of background young LTRs
Col7: Count of background old LTRs
Col8: Count of foreground young LTRs
Col9: Count of foreground old LTRs
```

This will also generate an output script `OUT_tNLR_tnonNLR_10` containing the following tab-separated rows from which the summary stats where calculated:

```
CML69   7152    Moderate        BACKGROUND
CML69   7548    Old     BACKGROUND
CML69   8257    Young   BACKGROUND
CML69   1619    unknown BACKGROUND
CML69   51      Moderate        NLR
CML69   52      Old     NLR
CML69   83      Young   NLR
CML69   20      unknown NLR
```

This analysis can be conducted for all NAM lines and for all relevant foreground and background combinations using a bash wrapper that takes only the value of k as input.

```
./scripts/run_tests.sh 5 > all_NAM_summary_stats_age.txt
```

The results for this analysis have been precomputed and compiled here: `6_plots/data/ltr_neighborhood_age.tsv`.

Next we can explore the maximum distance between the target cluster and the LTRs in the neighborhood. This will allow us to determine the relationship between k and the physical size of the analysed LTR neighborhood.

```
scripts/get_nlr_stats_k_dist_overlap.sh CML69 10 ../3_clustering/data/nam_tandem_nonnlr/CML69_nonnlr_100kb_merged.bed ../3_clustering/data/nam_tandem_nlr/CML69_nlr_100kb_merged.bed tNLR_tnonNLR
```

The output file `OUT_tNLR_tnonNLR_10_OVERLAP_DIST` contains the maximum physical distance in basepairs for the user-specified k, in this case k=10 and the physical distance from the cluster is ~160kb with the background cluster having a slightly greater distance.

```
CML69   172363  BACKGROUND      10
CML69   153280  NLR     10
```

We can also extract all LTR families within the k=10 LTR neighborhood of target genes or gene clusters and then test for enrichment of specific LTR families against the NAM-wide global LTR family background (`ltr_family_NAM_counts.tsv`). Test results are written to `test_result_tNLR_tnonNLR_10.txt`.

```
scripts/get_nlr_stats_k_family_all_overlap.sh CML69 10 ../3_clustering/data/nam_tandem_nonnlr/CML69_nonnlr_100kb_merged.bed ../3_clustering/data/nam_tandem_nlr/CML69_nlr_100kb_merged.bed tNLR_tnonNLR
scripts/te_family_enrichment_test.py OUT_ALL_tNLR_tnonNLR_10_OVERLAPFAM > test_result_tNLR_tnonNLR_10.txt
```

By modifying the python script to use `data/ltr_family_tropical_counts.tsv` and supplying LTR families from the neighborhood of tropical NAM lines only, the same test can be conducted specifically for tropical lines.
