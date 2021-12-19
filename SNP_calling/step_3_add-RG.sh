#!/bin/bash
module purge
module load jdk
source /work/LAS/mhufford-lab/shared_dir/minconda/20181213/etc/profile.d/conda.sh
conda activate star
ulimit -c unlimited
cpus=$SLURM_JOB_CPUS_PER_NODE
TMPDIR="/scratch/$USER/$SLURM_JOBID"
threads=$(echo "$SLURM_JOB_CPUS_PER_NODE - 2" | bc)
R1=$3
R2=$4
OUT=$(basename ${3} |sed 's/_R1_/_/g' | sed 's/.fastq.gz//g')
PLT=ILLUMINA
RG=$2
LIB=$1
TDATE=$(date '+%Y-%m-%d %H:%M:%S' |sed 's/ /T/g')
REF=/ptmp/LAS/arnstrm/shujun/mexicana/zmexicana.fa
picard="java -Xmx10g -Djava.io.tmpdir=$TMPDIR -jar /home/$USER/bin/picard.jar"
echo $R1
echo $R2
echo $OUT
echo $PLT
echo $TDATE
if [ ! -f "${OUT}.rg.done" ]; then
$picard AddOrReplaceReadGroups I=${OUT}_final.bam O=${OUT}_final2.bam CREATE_INDEX=TRUE RGID=${RG} RGSM=${RG} RGPU=${OUT} RGLB=${LIB} RGPL=${PLT} || {
echo >&2 ERROR: RG failed failed for $OUT
exit 1
}
fi
touch ${OUT}.rg.done
