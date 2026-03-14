#!/bin/bash
#SBATCH -p soeding
#SBATCH -N 1
#SBATCH -n 96
#SBATCH -t 2-00:00:00
#SBATCH --mem=900G

# ====================
# Environment
# ====================
source ~/.bashrc
conda activate sra_env

# ====================
# Paths
# ====================
WORKDIR=/cbscratch/umut.cakir/HIF_targets/RNA_Seq/paired_end_libraries
TMPDIR=/cbscratch/umut.cakir/HIF_targets/temp/sra_tmp
CACHE=/cbscratch/umut.cakir/HIF_targets/temp/sra_cache

export TMPDIR

mkdir -p ${WORKDIR}
cd ${WORKDIR}

# ====================
# Configure SRA Toolkit
# ====================
vdb-config --set /repository/user/cache=${CACHE}
vdb-config --set /repository/user/main/public/root=${CACHE}

# ====================
# Download + compress
# ====================
while read SRR; do
    echo "[$(date)] Downloading paired-end ${SRR}"

    fasterq-dump ${SRR} \
        --split-files \
        --threads 96 \
        --temp ${TMPDIR} \
        --outdir .

    if [[ $? -ne 0 ]]; then
        echo "ERROR downloading ${SRR}" >&2
        exit 1
    fi

    echo "[$(date)] Compressing ${SRR}_1.fastq and ${SRR}_2.fastq"
    pigz -p 96 ${SRR}_1.fastq ${SRR}_2.fastq

    echo "[$(date)] Finished ${SRR}"
done < srr_list.txt

echo "All paired-end SRR downloads completed successfully."

