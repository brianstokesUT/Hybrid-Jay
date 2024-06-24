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
bowtie2-build -f raw_sequences/c_stelleri_mt.fasta raw_sequences/c_stelleri_mt
samtools faidx raw_sequences/c_stelleri_mt.fasta

bowtie2 --very-sensitive-local -x raw_sequences/c_stelleri_mt -f -U raw_sequences/c_cristata_mt.fasta | samtools view -bS > prep_mt/c_cristata_mt.bam


#####################################3


samtools sort prep_mt/c_cristata_mt.bam -o prep_mt/sort.c_cristata_mt.bam
samtools index prep_mt/sort.c_cristata_mt.bam



########################

bcftools mpileup -Ou -f raw_sequences/c_stelleri_mt.fasta prep_mt/sort.c_cristata_mt.bam | bcftools call -Ou -mv | bcftools norm -f raw_sequences/c_stelleri_mt.fasta -Oz -o prep_mt/c_cristata_mt.vcf.gz
bcftools index prep_mt/c_cristata_mt.vcf.gz

######################################


conda deactivate
conda activate gatk4

gatk CreateSequenceDictionary -R raw_sequences/c_stelleri_mt.fasta

gatk IndexFeatureFile -I prep_mt/c_cristata_mt.vcf.gz

gatk FastaAlternateReferenceMaker -R raw_sequences/c_stelleri_mt.fasta -O c_cristata_mt_FINAL.fasta -V prep_mt/c_cristata_mt.vcf.gz

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


###########(Make this a loop or wildcard)

samtools sort raw_sequences/c_yncas/cy01.bam -o prep_mt/sort.cy01.bam
samtools index prep_mt/sort.cy01.bam

samtools sort raw_sequences/c_yncas/cy02.bam -o prep_mt/sort.cy02.bam
samtools index prep_mt/sort.cy02.bam

samtools sort raw_sequences/c_yncas/cy03.bam -o prep_mt/sort.cy03.bam
samtools index prep_mt/sort.cy03.bam

samtools sort raw_sequences/c_yncas/cy04.bam -o prep_mt/sort.cy04.bam
samtools index prep_mt/sort.cy04.bam


############


bcftools mpileup -Ou -f raw_sequences/c_stelleri_au.fasta prep_mt/sort.cy01.bam prep_mt/sort.cy02.bam prep_mt/sort.cy03.bam prep_mt/sort.cy04.bam | bcftools call -Ou -mv | bcftools norm -f raw_sequences/c_stelleri_au.fasta -Oz -o prep_mt/c_yncas.vcf.gz
bcftools index prep_mt/c_yncas.vcf.gz

#################

bcftools view prep_mt/c_yncas.vcf.gz --regions JANXIP010000352.1 -o prep_mt/c_yncas_mt.vcf.gz -Oz
bcftools index prep_mt/c_yncas_mt.vcf.gz

##################

conda deactivate
conda activate gatk4

gatk CreateSequenceDictionary -R raw_sequences/c_stelleri_au.fasta

gatk IndexFeatureFile -I prep_mt/c_yncas_mt.vcf.gz

gatk FastaAlternateReferenceMaker -R raw_sequences/c_stelleri_au.fasta -O c_yncas_mt_FINAL.fasta -V prep_mt/c_yncas_mt.vcf.gz -L JANXIP010000352.1

```


# Create MT Blast database

```
#########Tidy up our previously made *mt.fasta sequences
mkdir mt_blast
cp *mt_FINAL* ./mt_blast
cp raw_sequences/*mt.fasta ./mt_blast

###edit names of c_yncas & c_cristata because we used the c_stelleri refernce to make the fastas
sed -i '1s/.*/>c_yncas_mt/' mt_blast/c_yncas_mt_FINAL.fasta
sed -i '1s/.*/>c_cristata_mt/' mt_blast/c_cristata_mt_FINAL.fasta


#create our database by combining all fastas into a single fasta
cat mt_blast/c_yncas_mt_FINAL.fasta mt_blast/c_cristata_mt_FINAL.fasta mt_blast/a_californica_mt.fasta mt_blast/c_stelleri_mt.fasta mt_blast/p_pica_mt.fasta > mt_blast/jayz_mt.fasta

#make database
/work/08209/brian97/ls6/tools/ncbi-blast-2.14.0+/bin/makeblastdb -in mt_blast/jayz_mt.fasta -out mt_blast/jayz_mt_db -dbtype nucl -title jayz_mt_db

```

# Query database and check results
```
/work/08209/brian97/ls6/tools/ncbi-blast-2.14.0+/bin/blastn -query hyb_mt_FINAL.fasta -db mt_blast/jayz_mt_db -out jayz_mt_blast.out

#VIEW BY REMOVING "#" IN LINE BELOW
#nano jayz_mt_blast.out
```
