#!/bin/bash
#SBATCH -p soeding
#SBATCH -N 1
#SBATCH -n 96
#SBATCH -t 2-00:00:00
#SBATCH --mem=500G
#SBATCH --exclusive

source activate bowtie2_env

for i in *.fastq.gz; \
      do \
      echo "$i"
      bowtie2 \
                --phred33 -p 90 \
                --very-sensitive-local \
                -x human_index/hg38 \
                -U "$i" \
                -S ${i%.fastq.gz}.sam; \
      done
 
echo "Completed"
 
