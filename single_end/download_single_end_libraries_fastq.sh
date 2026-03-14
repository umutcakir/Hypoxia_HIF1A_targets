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
WORKDIR=/cbscratch/umut.cakir/HIF_targets/RNA_Seq/single_end_libraries
TMPDIR=/cbscratch/umut.cakir/HIF_targets/temp/sra_tmp
CACHE=/cbscratch/umut.cakir/HIF_targets/temp/sra_cache

export TMPDIR

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
    echo "[$(date)] Downloading ${SRR}"

    fasterq-dump ${SRR} \
        --threads 96 \
        --temp ${TMPDIR} \
        --outdir .

    if [[ $? -ne 0 ]]; then
        echo "ERROR downloading ${SRR}" >&2
        exit 1
    fi

    echo "[$(date)] Compressing ${SRR}.fastq"
    pigz -p 96 ${SRR}.fastq

    echo "[$(date)] Finished ${SRR}"
done < srr_list.txt

echo "All SRR downloads and compression completed successfully."

