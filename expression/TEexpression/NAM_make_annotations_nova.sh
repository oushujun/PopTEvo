#!/bin/bash

# Copy/paste this job script into a text file and submit with the command: 
#sbatch --array=2-26 --constraint=AVX2 NAM_make_annotations_nova.sh

#SBATCH --time=1:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --mail-type=BEGIN 
#SBATCH --mail-type=END 
#SBATCH --mail-type=FAIL


#Load modules
module load samtools/1.10-py3-xuj7ylj
module load bedtools2/2.27.1-s2mtpsu
module load gcc/7.3.0-xegsmw4

#Run Code

readarray -t FILES < NAM_genomes.txt

SAMPLE=${FILES[$SLURM_ARRAY_TASK_ID - 1]}
IFS=' '
read -ra INFO <<< "$SAMPLE"

FULL="$(cut -f 1 <<< $INFO)"
GENO="$(cut -f 2 <<< $INFO)"
ID="$(cut -f 3 <<< $INFO)"
VER="$(cut -f 4 <<< $INFO)"

grep exon ${FULL}_${ID}.1.gff3 > ${FULL}_${ID}.1.exon.gff3

grep gene ${FULL}_${ID}.1.gff3 > ${FULL}_${ID}.1.gene.gff3

sed -i 's/?/./g' ${GENO}.pseudomolecules-${VER}.fasta.mod.EDTA.TEanno.split.gff3

#Go in and manually run chr replacement commands
#sed -i 's/${GENO}_//g' ${GENO}.pseudomolecules-${VER}.fasta.mod.EDTA.TEanno.split.gff3

bedtools subtract -a ${GENO}.pseudomolecules-${VER}.fasta.mod.EDTA.TEanno.split.gff3 -b ${FULL}_${ID}.1.exon.gff3 > ${GENO}.allTE.subtractexon.gff3

sed -i 's/ID=/Name=gene:/g' ${FULL}_${ID}.1.gene.gff3

cat ${GENO}.allTE.subtractexon.gff3 ${FULL}_${ID}.1.gene.gff3 > ${GENO}.TE.subtractexon.plusgenes.gff3 

bedtools sort -i ${GENO}.TE.subtractexon.plusgenes.gff3 | awk '{OFS = "\t"} {if ($7 == "*") {$7 = "."; print;} else {print}}' - | awk '{OFS = "\t"} {$3 = "all"; print}' - > ${GENO}.TE.subtractexon.plusgenes.sort.gff3


sed 's/ID=/ID=${GENO}_/g' ${GENO}.TE.subtractexon.plusgenes.sort.gff3 > ${GENO}.TE.subtractexon.plusgenes.sort.geno.gff3 


