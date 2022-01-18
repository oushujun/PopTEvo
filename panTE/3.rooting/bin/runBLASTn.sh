#!/bin/bash
ref="/home/oushujun/jfw/TE/MaizeNAM/divergence/teosinte/Zx-PI566673-REFERENCE-YAN-1.0.fa";
#bam=$2
query=$1
cpus=$SLURM_JOB_CPUS_PER_NODE
TMPDIR="/mnt/bgfs/$USER/$SLURM_JOBID"
out=${query}.teosinte.out
donefile=${query}_blast.done

ml purge
#source /work/LAS/mhufford-lab/shared_dir/minconda/20181213/etc/profile.d/conda.sh
conda activate mapping

if [ ! -f ${donefile} ]; then

blastn -query $query -subject $ref -out $out \
   -perc_identity 0.8 -qcov_hsp_perc 0.8 -outfmt 6 -num_alignments 1000 -num_threads 36 || {

echo >&2 ERROR: blast failed for $query
exit 1
}
fi
touch $donefile
#rsync -av ${TMPDIR}/${out} /ptmp/LAS/oushujun/
