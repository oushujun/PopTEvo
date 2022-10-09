# PopTEvo
Study of TE evolution in a population of genomes

## Contents
 ### TE_annotation
 panEDTA annotation of NAM genomes
  * bin: scripts to summarize panEDTA annotations
  * data: 
    + B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.homo.chr.gt.pres.25k.h: 25k NAM SNPs using B73v5 as the reference
    + NAM.EDTA1.9.0.MTEC02052020.TE.v1.anno.intact.LTR.sort.gz: Intact LTR superfamily classification based on TEsorter
    + NAM.EDTA2.0.0.MTEC02052020.TElib.fa: panEDTA library generated for NAM genomes
    + NAM.intact.LTR.genedist.gz: physical distance of intact LTRs to the nearest genes.
    + pan_TE_bootstrap1000.summary26.txt: the occurrance of 1000 bootstrap resampling of panTE families in NAM genomes
    + soloLTR.txt: All solo LTRs found in NAM
    + TE_Fam_stats.txt: Statistics and descriptions of pangenome TE families
    + individual_genomes: panEDTA annotation results generated for each NAM genome
 ### panTE
 Population genomics of intact LTR elements
  * 1.pairwise: identification of syntenic LTRs between any two NAM genomes
    + bin: scripts to identify pairwise syntenic LTRs
    + data: syntenic LTR information for each genome pair
  * 2.combine: combine pairwise information to a pangenome table
    + bin: scripts to combine pairwise syntenic LTRs
    + data/final_27_genome_TE_resolved_B73_added_v2.txt.gz: syntenic LTR information in 26 NAM genomes
  * 3.rooting: identify ancestral states of LTR insertions using a teosinte genome
    + bin: scripts to root syntenic LTRs
    + data: syntenic LTR information in a teosinte genome
  * 4.site_frequency_spectrum: SFS studies for both SNP and syntenic LTRs
    + SNP_calling: scripts to call SNPs
    + data: VCF and SFS files for both SNPs nad syntenic LTRs
 ### methylation
 Methylation analysis of intact LTR elements
  * data: UMR information of each NAM genomes using B73v5 as the reference
 ### expression
 Expression analysis of TE families
  * scripts: scripts to map short reads to each genome
  * TEexpression: scripts to generate per-family read counts
  * data: raw count data for each TE family in each RNA library
 ### TE_fam_clean2.R
 R commands to compile all data (except expression) and generate most figures
 ### TE_fam_exp_clean.R
 R commands for expression-related analyses


