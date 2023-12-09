#!/usr/bin/env bash

k=$1
group=$2
# tandem nlr versus tandem non-nlr
grep ${group} ../4_nlr_neighborhood_age/data/groups.txt|  cut -f1| while read l; do scripts/get_nlr_stats_k_amplified.sh $l ${k} ../3_clustering/data/nam_tandem_nonnlr/${l}_nonnlr_100kb_merged.bed ../3_clustering/data/nam_tandem_nlr/${l}_nlr_100kb_merged.bed ${group}_tNLR100kb_tnonNLR100kb;done

# tandem nlr versus singleton nlr
grep ${group} ../4_nlr_neighborhood_age/data/groups.txt|  cut -f1| while read l; do scripts/get_nlr_stats_k_amplified.sh $l ${k} ../3_clustering/data/nam_singleton_nlr/${l}_nlr_singleton.bed ../3_clustering/data/nam_tandem_nlr/${l}_nlr_100kb_merged.bed ${group}_tNLR100kb_sNLR;done

# tandem non-nlr versus singleton non-nlr
grep ${group} ../4_nlr_neighborhood_age/data/groups.txt|  cut -f1| while read l; do scripts/get_nlr_stats_k_amplified.sh $l ${k} ../3_clustering/data/nam_singleton_nonnlr/${l}_nonnlr_singleton_sampled_5k.bed ../3_clustering/data/nam_tandem_nonnlr/${l}_nonnlr_100kb_merged.bed  ${group}_tnonNLR100kb_snonNLR;done

