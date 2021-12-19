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
picard="java -Xmx50g -Djava.io.tmpdir=$TMPDIR -jar /home/$USER/bin/picard.jar"
echo $R1
echo $R2
echo $OUT
echo $PLT
echo $TDATE
if [ ! -f "allsteps.done" ]; then
# convert fastq to sam and add readgroups
if [ ! -f "${OUT}.fq2sam.done" ]; then
$picard FastqToSam FASTQ=${R1} FASTQ2=${R2} OUTPUT=${OUT}_fastqtosam.bam READ_GROUP_NAME=${RG} SAMPLE_NAME=${OUT} LIBRARY_NAME=${LIB} PLATFORM_UNIT=${PLT} PLATFORM=illumina SEQUENCING_CENTER=ISU RUN_DATE=${TDATE}  || {
echo >&2 ERROR: FastqToSam failed for $OUT
exit 1
}
fi
touch ${OUT}.fq2sam.done
# marking adapters
if [ ! -f "${OUT}.markadapters.done" ]; then
$picard MarkIlluminaAdapters I=${OUT}_fastqtosam.bam O=${OUT}_markilluminaadapters.bam M=${OUT}_markilluminaadapters_metrics.txt TMP_DIR=${TMPDIR}  || {
echo >&2 ERROR: MarkIlluminaAdapters failed for $OUT
exit 1
}
fi
touch ${OUT}.markadapters.done
# convert bam back to fastq for mapping
if [ ! -f "${OUT}.sam2fq.done" ]; then
$picard SamToFastq I=${OUT}_markilluminaadapters.bam FASTQ=${OUT}_samtofastq_interleaved.fq CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true TMP_DIR=${TMPDIR} || {
echo >&2 ERROR: SamToFastq failed for $OUT
exit 1
}
fi
touch ${OUT}.sam2fq.done
# mapping reads to indexed genome
if [ ! -f "${OUT}.bwamem.done" ]; then
bwa-mem2 mem -M -t $threads -p ${REF%.*} ${OUT}_samtofastq_interleaved.fq | samtools view -buS - > ${OUT}_bwa_mem.bam || {
echo >&2 ERROR: BWA failed for $OUT
exit 1
}
fi
touch ${OUT}.bwamem.done
# merging alignments
if [ ! -f "${OUT}.mergealignments.done" ]; then
$picard MergeBamAlignment R=$REF UNMAPPED_BAM=${OUT}_fastqtosam.bam ALIGNED_BAM=${OUT}_bwa_mem.bam O=${OUT}_snippet_mergebamalignment.bam CREATE_INDEX=true ADD_MATE_CIGAR=true CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS TMP_DIR="${TMPDIR}" || {
echo >&2 ERROR: MergeBamAlignment failed for $OUT
exit 1
}
fi
touch ${OUT}.mergealignments.done
# mark duplicates
if [ ! -f "${OUT}.markdups.done" ]; then
$picard MarkDuplicates INPUT=${OUT}_snippet_mergebamalignment.bam OUTPUT=${OUT}_final.bam METRICS_FILE=${OUT}_mergebamalignment_markduplicates_metrics.txt OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 CREATE_INDEX=true TMP_DIR=$TMPDIR || {
echo >&2 ERROR: MarkDuplicates failed for $OUT
exit 1
}
fi
touch "${OUT}.markdups.done" "${OUT}.allsteps.done"
rm ${OUT}_fastqtosam.bam ${OUT}_markilluminaadapters.bam ${OUT}_samtofastq_interleaved.fq ${OUT}_bwa_mem.bam ${OUT}_snippet_mergebamalignment.bam ${OUT}_snippet_mergebamalignment.bai
fi
