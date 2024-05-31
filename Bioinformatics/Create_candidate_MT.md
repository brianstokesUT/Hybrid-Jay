# P. pica
Leave it as 
```
raw_sequences/p_pica_mt.fasta
```

# C. stelleri
Leave it as 
```
raw_sequences/c_stelleri_mt.fasta
```

# A. woodhousei
Leave it as 
```
raw_sequences/a_californica_mt.fasta
```

# C. Cristata
Need to work on the raw sequences


```
#index mitochondrial reference (c_stelleri_mt.fasta) - no need to waste computing aligning against full reference because we know where these sequences come from
samtools faidx raw_sequences/c_stelleri_mt.fasta
bowtie2-build -f raw_sequences/c_stelleri_mt.fasta raw_sequences/cs_mt_ref

#Create Alignment using bowtie2 --local preset
bowtie2 --local -f -x raw_sequences/cs_mt_ref -f -U raw_sequences/c_cristata_mt.fasta --no-unal | samtools view -bS > prep_mt/c_cristata_mt.bam

samtools sort prep_mt/c_cristata_mt.bam -o prep_mt/sort.c_cristata_mt.bam
samtools index prep_mt/sort.c_cristata_mt.bam


(i think this creates a single fasta and fills in gaps with the steller's jay)
#Then we create a single C. cristata mitochondrial fasta file which all gaps not represented by the raw sequence files being filled in by the Steller's Jay refernce
bcftools mpileup -Ou -f raw_sequences/c_stelleri_mt.fasta prep_mt/sort.c_cristata_mt.bam | bcftools call -Ou -mv | bcftools norm -f raw_sequences/c_cristata_mt.fasta -Oz -o prep_mt/c_cristata_mt.vcf.gz


```

# Origional verion for C. Cristata

```
###########ORIGIONAL VERSION
######BLUE JAY BLAST CREATION######
conda activate alignment

#make stellars jay mt index for bowtie
cd /scratch/08209/brian97/hybrid/mt_raw/cya_ste$
bowtie2-build -f JANXIP010000352.1.fasta JANXIP010000352.1
#make faidx index for bcftools mpileup
samtools faidx JANXIP010000352.1.fasta

#align bljay to stellar's ref
cd /scratch/08209/brian97/hybrid/mt_raw/cya_cri
(USING SENSITIVE BC I KNOW THESE ALIGN SOMEWHERE)
###bowtie2 --very-sensitive -x /scratch/08209/brian97/hybrid/mt_raw/cya_ste/JANXIP010000352.1 -U cc_ncbi.fq | samtools view -bS > cya_cri.bam###
#found out can't use the fq bc it brings in the super weird low phred scores!!!
bowtie2 --very-sensitive-local -x /scratch/08209/brian97/hybrid/mt_raw/cya_ste/JANXIP010000352.1 -f -U cc_ncbi.fasta | samtools view -bS > cya_cri.bam




samtools sort cya_cri.bam -o sort_cya_cri.bam
samtools index sort_cya_cri.bam 

samtools depth -a sort_cya_cri.bam > depth.txt

###Produce quick quality report of the alignments so we have the stats for paper
conda deactivate
conda activate qualimap
qualimap bamqc -bam sort_cya_cri.bam -outdir qualimap_results
###75.04% of the ref mt genome is covered by at least one mapped basepair


#lets make variant calls, don't need any quality filter because these are NCBI and should all be high quality?
conda deactivate
conda activate bcftools

bcftools mpileup -Ou -f /scratch/08209/brian97/hybrid/mt_raw/cya_ste/JANXIP010000352.1.fasta sort_cya_cri.bam | bcftools call -Ou -mv | bcftools norm -f /scratch/08209/brian97/hybrid/mt_raw/cya_ste/JANXIP010000352.1.fasta -Oz -o cyacri_MT.vcf.gz
bcftools index cyacri_MT.vcf.gz
bcftools consensus -f /scratch/08209/brian97/hybrid/mt_raw/cya_ste/JANXIP010000352.1.fasta cyacri_MT.vcf.gz > cyacri_mt_con1.fasta













bcftools index cya_cri.vcf.gz

bcftools consensus -f /scratch/08209/brian97/hybrid/mt_raw/pica_ref/pica_mt.fasta cya_cri.vcf.gz > cya_cri.fa


```

















# C. yncas
Need to do lots of work on the raw paired end read files.

```
for R1_FILE in raw_sequences/c_yncas/*_R1_001.fastq.gz 
do
	#Extract R1 File Name
	BASE_NAME=$(basename $R1_FILE _R1_001.fastq.gz)
	PREFIX=${BASE_NAME:0:2}

	#Reconstruct R2 File Name
	R2_FILE="${BASE_NAME}_R2_001.fastq.gz"
	
	#BAM output file name
	OUTPUT_BAM="cy${PREFIX}.bam"

	bowtie2 -q -x raw_sequences/cs_ref -1 raw_sequences/c_yncas/$R1_FILE -2 raw_sequences/c_yncas/$R2_FILE --no-unal | samtools view -bS > raw_sequences/c_yncas/$OUTPUT_BAM


done
```

