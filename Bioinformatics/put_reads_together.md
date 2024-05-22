#Align Hybrid to Steller's Jay reference genome

```
conda activate alignment

#Align trimmed hybrid reads to Steller's Jay reference genome
bowtie2 --threads 24 -f -x /scratch/08209/brian97/hybrid/mt_raw/grja_mt_prep/stja_ref/GCA_026167965.1_bCyaSte1.0.p_genomic -U /scratch/08209/brian97/hybrid/hyb_filt.fasta --no-unal | samtools view -bS > /scratch/08209/brian97/hybrid/mt_new/hyb.bam
samtools sort hyb.bam -o sort.hyb.bam

samtools view -b sort.hyb.bam JANXIP010000352.1 -o hyb_mt.bam
samtools sort hyb_mt.bam -o sort.hyb_mt.bam
samtools index sort.hyb_mt.bam

#create mask text regions file to cover unmapped/low quality positions
samtools depth -a -Q 10 sort.hyb_mt.bam > mask.txt
awk '$3==0 { print }' mask.txt > mask2.txt
cut -d$'\t' -f1,2 mask2.txt > mask3.txt
#this is the locations to cover with "N"



bcftools mpileup -Ou -q 20 -r JANXIP010000352.1 -f /scratch/08209/brian97/hybrid/mt_raw/grja_mt_prep/stja_ref/GCA_026167965.1_bCyaSte1.0.p_genomic.fna sort.hyb_mt.bam | bcftools call -Ou -mv | bcftools norm -f /scratch/08209/brian97/hybrid/mt_raw/grja_mt_prep/stja_ref/GCA_026167965.1_bCyaSte1.0.p_genomic.fna -Oz -o hyb_mt.vcf.gz
bcftools index hyb_mt.vcf.gz



###MAKE conesnsus of the hybrid mitochondria

samtools faidx /scratch/08209/brian97/hybrid/mt_raw/grja_mt_prep/stja_ref/GCA_026167965.1_bCyaSte1.0.p_genomic.fna JANXIP010000352.1 | bcftools consensus -m mask3.txt hyb_mt.vcf.gz > hyb_mt.fasta





```
