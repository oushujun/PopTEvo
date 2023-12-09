#!/usr/bin/env python3

import sys
import scipy.stats as stats

backdict = {}
foredict = {}

tropicals = ["CML103","CML228","CML247","CML277","CML322","CML333","CML52","CML69","Ki11","Ki3","NC350","NC358","Tzi8"]

total_count = 0
with open("data/ltr_family_NAM_counts.tsv") as background:
    for l in background:
        l = l.strip().split(' ')
        count = int(l[0])
        family = l[1]
        if count >= 100:
            backdict[family] = {'this_family':count,'not_this_family':0}
            total_count += count

for f in backdict:
    famcount = backdict[f]['this_family']
    backdict[f] = {'this_family':famcount,'not_this_family':total_count - famcount}

kdict = {}
with open(sys.argv[1]) as foreground:
    # for each k, count the NLR associated TE families by summing over (tropical) NAM lines
    for l in foreground:
        l = l.strip().split('\t')

        group = l[0].split('_')[3]
        overlap = l[0].split('_')[5]
        k = l[0].split('_')[4]
        # BACKGROUND or NLR
        clust_group = l[4]
        count = int(l[2])
        family = l[3]
        nam_line = l[1]
        if k not in kdict:
            kdict[k] = {}
        if nam_line in tropicals and clust_group == "NLR" and group == "100kb" and overlap == "OVERLAPFAM":
            if family not in kdict[k]:
                kdict[k][family] = 0
            kdict[k][family] += count

for k,v in kdict.items():
    total_count = 0
    for fam,count in v.items():
        total_count += count

    for fam,count in v.items():
        if fam in backdict:
            this_family = count
            not_this_family = total_count-count
            back_this_family = backdict[fam]['this_family']
            back_not_this_family = backdict[fam]['not_this_family']
            contingency_tab = [[this_family,not_this_family],[back_this_family,back_not_this_family]]
            odd_ratio, p_value = stats.fisher_exact(contingency_tab)
            outline = [k,fam,round(odd_ratio,4),p_value*len(kdict[k]),p_value,this_family,not_this_family,back_this_family,back_not_this_family]
            print(*outline,sep='\t')

