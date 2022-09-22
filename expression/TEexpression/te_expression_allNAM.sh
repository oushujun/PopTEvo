#!/bin/bash

# Copy/paste this job script into a text file and submit with the command: 
#sbatch --array=1-480 --constraint=AVX2 te_expression_allNAM.sh

#SBATCH --time=18:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --job-name="NAM" 
#SBATCH --mail-type=BEGIN 
#SBATCH --mail-type=END 
#SBATCH --mail-type=FAIL


#Load modules
module load samtools/1.10-py3-xuj7ylj
module load bedtools2/2.27.1-s2mtpsu
module load gcc/7.3.0-xegsmw4

#Run Code

readarray -t FILES < samples_allNAM.txt

SAMPLE=${FILES[$SLURM_ARRAY_TASK_ID - 1]}
IFS=' '
read -ra INFO <<< "$SAMPLE"

BAM="$(cut -f 1 <<< $INFO)"
NAME="$(cut -f 2 <<< $INFO)"
GENO="$(cut -f 3 <<< $INFO)"

samtools sort -n ${BAM} > ${NAME}.hisat.sortedByName.out.bam

module load py-pip/20.2-py3-gibulwf

python -m HTSeq.scripts.count -s no -t all -i Name -m union -a 0 --nonunique all -o ${NAME}.hisat.out.sam ${NAME}.hisat.sortedByName.out.bam ${GENO}.TE.subtractexon.plusgenes.sort.gff3  > counts_intermediate_${NAME}.hisat.sortedByName.out.txt

perl  te_family_mapping_ver8.2_NAM.pl ${NAME}.hisat.out.sam ${NAME}
