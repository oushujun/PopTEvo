#!/usr/bin/env python3

import sys

if sys.argv[1] == "singleton":
    dup = False
elif sys.argv[1] == "duplicate":
    dup = True

if sys.argv[2] == "nonnlr":
    nlrgene = False
elif sys.argv[2] == "nlr":
    nlrgene = True

# max_dist = 50000
max_dist = int(sys.argv[3])

def interval_dist(r1, r2):
     # sort the two ranges such that the range with smaller first element
     # is assigned to x and the bigger one is assigned to y
     x, y = sorted((r1, r2))
     #now if x[1] lies between x[0] and y[0](x[1] != y[0] but can be equal to x[0])
     #then the ranges are not overlapping and return the differnce of y[0] and x[1]
     #otherwise return 0 
     if x[0] <= x[1] < y[0] and all( y[0] <= y[1] for y in (r1,r2)):
        return y[0] - x[1]
     return 0


nlr = {}
with open('../data/gene_info/nlrannotator_nlrome_unique.bed','r') as nlrlist:
    for l in nlrlist:
        l = l.strip().split('\t')
        namline = l[6]
        gene = l[3]
        if namline not in nlr:
            nlr[namline] = []
        nlr[namline].append(gene)

# store all genes per line
gpos = {}
with open('../data/gene_info/nam_gene_locations.bed','r') as pos:
    for l in pos:
        l = l.strip().split('\t')
        namline = l[0].split('-')[1]
        gene = l[4]
        chrom = l[1]
        start = int(l[2])
        end = int(l[3])
        if namline not in gpos:
            gpos[namline] = {}
        #if gene not in gpos[namline]:
        #    gpos[namline][gene] = []
        gpos[namline][gene] = [chrom,start,end]

#Zm-Ki11-REFERENCE-NAM-1.0_Zm00030ab.1.gff3	chr1	54156	59780	Zm00030ab000020

# encode pos as [chr,start,end]

with open('../data/gene_info/pangene-files_pan_gene_matrix_v3_cyverse.csv','r') as pan:
    header = next(pan)
    header = header.strip().split(',')
    for line in pan:
        line = line.strip().split(',')
        for i,namline in enumerate(header[3:29]):
            idx = i + 3
            genes = line[idx].split(';')
            good_genes = []
            good_neighbors = set()
            # singleton
            if not dup and len(genes) == 1 and "gmap" not in genes[0] and "NA" not in genes[0]:
                gname = genes[0].split('_')[0]
                gpos_single = gpos[namline][gname]
                outline = gpos_single + [gname,namline]
                if nlrgene and gname in nlr[namline]:
                    good_genes.append(gname)
                    print(*outline,sep='\t')
                elif not nlrgene and gname not in nlr[namline]:
                    good_genes.append(gname)
                    print(*outline,sep='\t')
            elif dup and len(genes) > 1:
                for g in genes:
                    if "gmap" not in g and "NA" not in g:
                        gname = g.split('_')[0]
                        if nlrgene and gname in nlr[namline]:
                            good_genes.append(gname)
                        elif not nlrgene and gname not in nlr[namline]:
                            good_genes.append(gname)
                # keep only neighboring genes
                for gg in good_genes:
                    for m in good_genes:
                        if gg != m:
                            gpos1 = gpos[namline][gg]
                            gpos2 = gpos[namline][m]
                            # same chrom
                            if gpos1[0] == gpos2[0]:
                                dist = interval_dist([gpos1[1],gpos1[2]],[gpos2[1],gpos2[2]])
                                if dist <= max_dist:
                                    good_neighbors.add(gg)
                                    good_neighbors.add(m)
                # print good neighbors
                for gn in good_neighbors:
                    gn_gpos = gpos[namline][gn]
                    out_list = ','.join(list(good_neighbors))
                    outline = gn_gpos + [gn,namline,out_list]
                    print(*outline,sep='\t')
