

First we'll pull a selected autosomal sequence of interest from the *C. cristata* dataset
In this example we'll focus on AY443280.1 which is one of the longest sequences available and includes partial coverage of the RAG1 gene region. Other regions will follow similar steps but will need some edited scripts to replace which region is used.

```
awk '/^>AY443280.1/{flag=1;print;next}/^>/{flag=0}flag' raw_sequences/c_cristata_au.fasta > raw_sequences/c_cristata_AY443280.1.fasta
```

Now we want to find the homologous region within the Steller's Jay genome - this should also be verified by checking through NCBI.
```
#make blast db
~PATH/tools/ncbi-blast-2.14.0+/bin/makeblastdb -in raw_sequences/c_stelleri_au.fasta -out rag1/c_stelleri_au_db -dbtype nucl -title rag1/c_stelleri_au_db

#search the the RAG1 sequence of C. cristata against this db to find the homologous region
~PATH/tools/ncbi-blast-2.14.0+/bin/blastn -query raw_sequences/c_cristata_AY443280.1.fasta -db rag1/c_stelleri_au_db -out rag1/AY443280.1_rag1_blast.out

#nano rag1/AY443280.1_rag1_blast.out
#JANXIP010000005.1:37834122-37836993

#now create a fasta of just this region from C. stelleri 

mamba activate data_prep
samtools faidx raw_sequences/c_stelleri_au.fasta JANXIP010000005.1:37834122-37836993 > rag1/c_stelleri_rag1.fasta
```
Then we do the same query and split process for all other samples
```
#########################A. californica
#make blast db
~PATH/tools/ncbi-blast-2.14.0+/bin/makeblastdb -in raw_sequences/a_californica_au.fasta -out rag1/a_californica_au_db -dbtype nucl -title rag1/a_californica_au_db

#search the the RAG1 sequence of A. californica against this db to find the homologous region
~PATH/tools/ncbi-blast-2.14.0+/bin/blastn -query raw_sequences/c_cristata_AY443280.1.fasta -db rag1/a_californica_au_db -out rag1/a_californica_rag1_blast.out

#nano rag1/a_californica_rag1_blast.out
#JAQMYR010000007.1:21718278-21721149

#now create a fasta of just this region from A. californica 

mamba activate data_prep
samtools faidx raw_sequences/a_californica_au.fasta JAQMYR010000007.1:21718278-21721149 > rag1/a_californica_rag1.fasta



#########################P. pica
#make blast db
~PATH/tools/ncbi-blast-2.14.0+/bin/makeblastdb -in raw_sequences/p_pica_au.fasta -out rag1/p_pica_au_db -dbtype nucl -title rag1/p_pica_au_db

#search the the RAG1 sequence of A. californica against this db to find the homologous region
~PATH/tools/ncbi-blast-2.14.0+/bin/blastn -query raw_sequences/c_cristata_AY443280.1.fasta -db rag1/p_pica_au_db -out rag1/p_pica_rag1_blast.out

#nano rag1/p_pica_rag1_blast.out
#JAOYNA010000039.1:44057536-44060407

#now create a fasta of just this region from A. californica 

mamba activate data_prep
samtools faidx raw_sequences/p_pica_au.fasta JAOYNA010000039.1:44057536-44060407 > rag1/p_pica_rag1.fasta
```


*C. yncas* alignment will have the same homologous positionality as the Steller's Jay reference.
```
#call high confidencevariants at the RAG1 gene from our previously sorted bam files
bcftools mpileup -Ou -q 20 -r JANXIP010000005.1:37834112-37836993 -f raw_sequences/c_stelleri_au.fasta prep_mt/sort.cy01.bam prep_mt/sort.cy02.bam prep_mt/sort.cy03.bam prep_mt/sort.cy04.bam | bcftools call -Ou -mv | bcftools norm -f raw_sequences/c_stelleri_au.fasta -Oz -o rag1/c_yncas_rag1.vcf.gz 
bcftools index rag1/c_yncas_rag1.vcf.gz

#create consensus sequence at the same homologous region of C. stelleri because that served as the C. yncas reference
mamba deactivate
mamba activte gatk4

gatk IndexFeatureFile -I rag1/c_yncas_rag1.vcf.gz

gatk FastaAlternateReferenceMaker -R raw_sequences/c_stelleri_au.fasta -O rag1/c_yncas_rag1.fasta -V rag1/c_yncas_rag1.vcf.gz -L JANXIP010000005.1:37834122-37836993

```

Next homologous regions of the hybrid reads along with 100bp flanking on each end.
```
###########Call aligned sequences from hybrid which were mapped to specific genes (IMPORTANT NOT TO CALL CONSENSUS)

mamba deactivate
mamba activate data_prep
#check depth at the RAG1 gene 
#samtools depth -a -r JANXIP010000005.1:37834012-37837003 prep_hybrid/sort.hyb.bam | more
#depth is really solid here

#subset highquality mapped reads from this region (min MAPQ score of 10)
samtools view -q 20 -b prep_hybrid/sort.hyb.bam JANXIP010000005.1:37834112-37836993 > rag1/hyb_rag1.bam

samtools sort rag1/hyb_rag1.bam -o rag1/sort.hyb_rag1.bam
samtools index rag1/sort.hyb_rag1.bam

samtools fasta rag1/sort.hyb_rag1.bam -0 rag1/hyb_rag1.fasta

```

Prep for BLAST database

```
#quickly edit RAG1 fasta files headers to have easier to use names
sed -i '1s/.*/>c_yncas_rag1/' rag1/c_yncas_rag1.fasta
sed -i '1s/.*/>c_stelleri_rag1/' rag1/c_stelleri_rag1.fasta
sed -i '1s/.*/>a_californica_rag1/' rag1/a_californica_rag1.fasta
sed -i '1s/.*/>p_pica_rag1/' rag1/p_pica_rag1.fasta



####Put all RAG1's into a single fasta
cat raw_sequences/c_cristata_AY443280.1.fasta rag1/c_yncas_rag1.fasta rag1/c_stelleri_rag1.fasta rag1/a_californica_rag1.fasta rag1/p_pica_rag1.fasta > rag1/candidate_rag1.fasta

#use BLAST makeblastdb tool
~PATH/tools/ncbi-blast-2.14.0+/bin/makeblastdb -in rag1/candidate_rag1.fasta -out rag1/candidate_rag1_db -dbtype nucl -title rag1/candidate_rag1_db


#Run against our hybird just to see how it looks - more accurate results will come from phasing as described in "phasing.md"
~PATH/tools/ncbi-blast-2.14.0+/bin/blastn -query ~PATH/autosomes/homolog/hyb_rag1.fasta -db rag1/candidate_rag1_db -out rag1/rag1_blast.out

#view results
nano rag1/rag1_blast.out

```






