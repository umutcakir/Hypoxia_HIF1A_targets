#!/bin/bash
#SBATCH -p soeding
#SBATCH -N 1
#SBATCH -n 60
#SBATCH -t 2-00:00:00
#SBATCH --mem=500G
#SBATCH --job-name=STAR_bulkRNA
#SBATCH --output=STAR_bulkRNA_%j.out
#SBATCH --error=STAR_bulkRNA_%j.err

set -euo pipefail

############################
# Activate environment
############################
source activate hg38_STAR_Salmon

############################
# Variables
############################
THREADS=60

STAR_INDEX=/cbscratch/umut.cakir/HIF_targets/hg38_genome_index/STAR_hg38
GTF=/cbscratch/umut.cakir/HIF_targets/hg38_genome_index/genes_plus_TE.sorted_hg38.gtf

OUTDIR=STAR_out
COUNTDIR=counts

mkdir -p ${OUTDIR}
mkdir -p ${COUNTDIR}

############################
# STEP 1: STAR alignment (near-default)
############################
echo "[INFO] Starting STAR alignment"

for fq in *.fastq.gz; do
    sample=$(basename ${fq} .fastq.gz)

    echo "[INFO] Aligning ${sample}"

    STAR \
        --runThreadN ${THREADS} \
        --genomeDir ${STAR_INDEX} \
        --readFilesIn ${fq} \
        --readFilesCommand zcat \
        --outFileNamePrefix ${OUTDIR}/${sample}. \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS nM
done

echo "[INFO] STAR alignment completed"

############################
# STEP 2: featureCounts (your rule enforced)
############################
echo "[INFO] Running featureCounts"

featureCounts \
    -T ${THREADS} \
    -s 0 \
    -t exon \
    -g gene_id \
    -a ${GTF} \
    -o ${COUNTDIR}/gene_counts.txt \
    ${OUTDIR}/*.bam


echo "[INFO] featureCounts completed"

