## Testing for enrichment of tropical-amplifying LTRs in neighborhoods of NLR clusters

A test was conducted for the tropical and temperate NAM lines to determine whether tropical-amplifying LTRs were enriched in the neighborhoods of NLR clusters. Statistics can be generated with a bash wrapper as shown below.

```
./scripts/run_tests.sh 5 tropical >> tropical_amplifying_test_statistics.txt
./scripts/run_tests.sh 5 temperate >> tropical_amplifying_test_statistics.txt
```

Statistics are printed to `tropical_amplifying_test_statistics.txt` and the columns in the statistics contain the following information.

```
Col1: Group name
Col2: NAM line
Col3: Number of LTRs in neighborhood (k)
Col4 : Background tropical-amplifying divided by background not tropical-amplifying
Col5 : Foreground tropical-amplifying divided by background not tropical-amplifying
Col6 : Count of background not tropical-amplifying
Col7 : Count of background tropical-amplifying
Col8 : Count of foreground not tropical-amplifying
Col9 : Count of foreground tropical-amplifying
```

To provide additional raw outputs, the wrapper als writes individual output files for each test and group. For instance, for the tropical lines the following output files are generated.

* `OUT_AMP_tropical_tNLR100kb_sNLR_5`: NLR clusters (foreground) versus singleton NLRs (background)
* `OUT_AMP_tropical_tNLR100kb_tnonNLR100kb_5`: NLR clusters (foreground) versus non-NLR clusters (background)
* `OUT_AMP_tropical_tnonNLR100kb_snonNLR_5`: non-NLR clusters (foreground) versus singleton non-NLRs (background)

Pre-computed and processed results are provided in `6_plots/data/ltr_neighborhood_temperate_tropical-amplifying.tsv` and `6_plots/data/ltr_neighborhood_tropical_tropical-amplifying.tsv`.

