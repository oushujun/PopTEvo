#!/usr/bin/env bash

k=$1
# tandem nlr versus tandem non-nlr
cat data/groups.txt | cut -f1| while read l; do scripts/get_nlr_stats_k.sh $l ${k} ../3_clustering/data/nam_tandem_nonnlr/${l}_nonnlr_100kb_merged.bed ../3_clustering/data/nam_tandem_nlr/${l}_nlr_100kb_merged.bed tNLR100kb_tnonNLR100kb;done

# tandem nlr versus singleton nlr
cat data/groups.txt | cut -f1| while read l; do scripts/get_nlr_stats_k.sh $l ${k} ../3_clustering/data/nam_singleton_nlr/${l}_nlr_singleton.bed ../3_clustering/data/nam_tandem_nlr/${l}_nlr_100kb_merged.bed tNLR100kb_sNLR;done

# tandem nlr versus singleton non-nlr
cat data/groups.txt | cut -f1| while read l; do scripts/get_nlr_stats_k.sh $l ${k} ../3_clustering/data/nam_singleton_nonnlr/${l}_nonnlr_singleton_sampled_5k.bed ../3_clustering/data/nam_tandem_nlr/${l}_nlr_100kb_merged.bed tNLR100kb_snonNLR;done

# tandem non-nlr versus singleton non-nlr
cat data/groups.txt | cut -f1| while read l; do scripts/get_nlr_stats_k.sh $l ${k} ../3_clustering/data/nam_singleton_nonnlr/${l}_nonnlr_singleton_sampled_5k.bed ../3_clustering/data/nam_tandem_nonnlr/${l}_nonnlr_100kb_merged.bed  tnonNLR100kb_snonNLR;done

# singleton non-nlr versus singleton nlr
cat data/groups.txt | cut -f1| while read l; do scripts/get_nlr_stats_k.sh $l ${k} ../3_clustering/data/nam_singleton_nonnlr/${l}_nonnlr_singleton_sampled_5k.bed ../3_clustering/data/nam_singleton_nlr/${l}_nlr_singleton.bed sNLR_snonNLR;done

cat data/groups.txt | cut -f1| while read l; do scripts/get_nlr_stats_k.sh $l ${k} ../3_clustering/data/nam_tandem_nonnlr/${l}_nonnlr_100kb_merged.bed ../3_clustering/data/nam_tandem_nlr/${l}_nlr_100kb_merged.bed tNLR100kb_tnonNLR100kb;done

# tandem nlr versus singleton nlr
cat data/groups.txt | cut -f1| while read l; do scripts/get_nlr_stats_k.sh $l ${k} ../3_clustering/data/nam_singleton_nlr/${l}_nlr_singleton.bed ../3_clustering/data/nam_tandem_nlr/${l}_nlr_100kb_merged.bed tNLR100kb_sNLR;done

# tandem nlr versus singleton non-nlr
cat data/groups.txt | cut -f1| while read l; do scripts/get_nlr_stats_k.sh $l ${k} ../3_clustering/data/nam_singleton_nonnlr/${l}_nonnlr_singleton_sampled_5k.bed ../3_clustering/data/nam_tandem_nlr/${l}_nlr_100kb_merged.bed tNLR100kb_snonNLR;done

# tandem non nlur verus singleton non nlr 100kb
cat data/groups.txt | cut -f1| while read l; do scripts/get_nlr_stats_k.sh $l ${k} ../3_clustering/data/nam_singleton_nonnlr/${l}_nonnlr_singleton_sampled_5k.bed ../3_clustering/data/nam_tandem_nonnlr/${l}_nonnlr_100kb_merged.bed tnonNLR100kb_snonNLR;done


