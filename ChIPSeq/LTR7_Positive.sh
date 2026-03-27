############################################################
# INPUT FILES
############################################################
# HIF1A_all_combined.bed
# hg38_RepeatMasker.bed
# ../../../hg38_genome_index/GRCh38.primary_assembly.genome.fa

############################################################
# STEP 0 — Remove duplicate summits
############################################################

sort -k1,1 -k2,2n -k3,3n HIF1A_all_combined.bed | \
awk '!seen[$1,$2,$3]++' > HIF1A_unique.bed

############################################################
# STEP 1 — Keep only standard chromosomes
############################################################

awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/' HIF1A_unique.bed \
> HIF1A_standard.bed

awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/' hg38_RepeatMasker.bed \
> TE_standard.bed

############################################################
# STEP 2 — Convert summits to 1 bp and extend ±50 bp
############################################################

awk '{print $1"\t"$2"\t"$2+1}' HIF1A_standard.bed \
> HIF1A_1bp.bed

awk 'BEGIN{OFS="\t"}{
  start = $2 - 50;
  if (start < 0) start = 0;
  end = $3 + 50;
  print $1, start, end
}' HIF1A_1bp.bed > HIF1A_50bp.bed

############################################################
# STEP 3 — Extract only exact LTR7
############################################################

awk '$4 == "LTR7"' TE_standard.bed > LTR7_only.bed

############################################################
# STEP 4 — Sort
############################################################

sort -k1,1 -k2,2n HIF1A_50bp.bed > HIF1A_50bp_sorted.bed
sort -k1,1 -k2,2n LTR7_only.bed > LTR7_only_sorted.bed

############################################################
# STEP 5 — POSITIVE: summit windows overlapping LTR7
############################################################

bedtools intersect \
  -a HIF1A_50bp_sorted.bed \
  -b LTR7_only_sorted.bed \
  -u > LTR7_positive_50bp.bed

############################################################
# STEP 6 — NEGATIVE: LTR7 loci NOT overlapped by any summit
############################################################

bedtools intersect \
  -a LTR7_only_sorted.bed \
  -b HIF1A_50bp_sorted.bed \
  -v > LTR7_negative.bed

############################################################
# STEP 7 — Extract FASTA
############################################################

bedtools getfasta \
  -fi ../../../hg38_genome_index/GRCh38.primary_assembly.genome.fa \
  -bed LTR7_positive_50bp.bed \
  -fo LTR7_positive_50bp.fa

bedtools getfasta \
  -fi ../../../hg38_genome_index/GRCh38.primary_assembly.genome.fa \
  -bed LTR7_negative.bed \
  -fo LTR7_negative.fa \
  -name -s

############################################################
# STEP 8 — QC
############################################################

wc -l HIF1A_unique.bed
wc -l HIF1A_50bp_sorted.bed
wc -l LTR7_positive_50bp.bed
wc -l LTR7_negative.bed



