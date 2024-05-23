#Align Hybrid to Steller's Jay reference genome

```
###Align hyb to steller's reference genome (conda activate alignment)

#index reference (c_stelleri_au.fasta)
bowtie2-build raw_sequences/c_stelleri_au.fasta raw_sequences/cs_ref

#align and sort
bowtie2 --threads 24 -q -x raw_sequences/cs_ref -U prep_hybrid/raw_hybrid_trimmed.fq.gz --no-unal | samtools view -bS > prep_hybrid/hyb.bam
samtools sort prep_hybrid/hyb.bam -o prep_hybrid/sort.hyb.bam
samtools index prep_hybrid/sort.hyb.bam

#extract mitochondrial scaffold
samtools view -b prep_hybrid/sort.hyb.bam JANXIP010000352.1 -o prep_hybrid/hyb_mt.bam
samtools sort prep_hybrid/hyb_mt.bam -o prep_hybrid/sort.hyb_mt.bam
samtools index prep_hybrid/sort.hyb_mt.bam

#create mask text regions file to cover unmapped/low quality positions
samtools depth -a -Q 10 prep_hybrid/sort.hyb_mt.bam > prep_hybrid/mask.txt
awk '$3==0 { print }' prep_hybrid/mask.txt > prep_hybrid/mask2.txt
cut -d$'\t' -f1,2 prep_hybrid/mask2.txt > prep_hybrid/mask3.txt
#this is the locations to cover with "N"

samtools faidx raw_sequences/c_stelleri_au.fasta
bcftools mpileup -Ou -q 20 -r JANXIP010000352.1 -f raw_sequences/c_stelleri_au.fasta prep_hybrid/sort.hyb_mt.bam | bcftools call -Ou -c | bcftools norm -f raw_sequences/c_stelleri_au.fasta -Oz -o prep_hybrid/hyb_mt.vcf.gz
bcftools index prep_hybrid/hyb_mt.vcf.gz



###MAKE consensus of the hybrid mitochondria

samtools faidx raw_sequences/c_stelleri_au.fasta JANXIP010000352.1 | bcftools consensus -m prep_hybrid/mask3.txt prep_hybrid/hyb_mt.vcf.gz > hyb_mt_FINAL.fasta

```
