

First we'll pull a selected autosomal sequence of interest from the *C. cristata* dataset

```
awk '/^>AY443280.1/{flag=1;print;next}/^>/{flag=0}flag' raw_sequences/c_cristata_au.fasta > raw_sequences/c_cristata_AY443280.1.fasta
```

```
#make blast db
/work/08209/brian97/ls6/tools/ncbi-blast-2.14.0+/bin/makeblastdb -in raw_sequences/c_stelleri_au.fasta -out rag1/c_stelleri_au_db -dbtype nucl -title rag1/c_stelleri_au_db

#search the the RAG1 sequence of C. cristata against this db to find the homologous region
/work/08209/brian97/ls6/tools/ncbi-blast-2.14.0+/bin/blastn -query raw_sequences/c_cristata_AY443280.1.fasta -db rag1/c_stelleri_au_db -out rag1/AY443280.1_rag1_blast.out

#nano rag1/AY443280.1_rag1_blast.out
#JANXIP010000005.1:37834122-37836993

#now create a fasta of just this region from C. stelleri 

conda activate samtools
samtools faidx raw_sequences/c_stelleri_au.fasta JANXIP010000005.1:37834122-37836993 > rag1/c_stelleri_rag1.fasta
```
Then we do the same query and split process for all other samples
```
#########################A. californica
#make blast db
/work/08209/brian97/ls6/tools/ncbi-blast-2.14.0+/bin/makeblastdb -in raw_sequences/a_californica_au.fasta -out rag1/a_californica_au_db -dbtype nucl -title rag1/a_californica_au_db

#search the the RAG1 sequence of A. californica against this db to find the homologous region
/work/08209/brian97/ls6/tools/ncbi-blast-2.14.0+/bin/blastn -query raw_sequences/c_cristata_AY443280.1.fasta -db rag1/a_californica_au_db -out rag1/a_californica_rag1_blast.out

#nano rag1/a_californica_rag1_blast.out
#JAQMYR010000007.1:21718278-21721149

#now create a fasta of just this region from A. californica 

conda activate samtools
samtools faidx raw_sequences/a_californica_au.fasta JAQMYR010000007.1:21718278-21721149 > rag1/a_californica_rag1.fasta



#########################P. pica
#make blast db
/work/08209/brian97/ls6/tools/ncbi-blast-2.14.0+/bin/makeblastdb -in raw_sequences/p_pica_au.fasta -out rag1/p_pica_au_db -dbtype nucl -title rag1/p_pica_au_db

#search the the RAG1 sequence of A. californica against this db to find the homologous region
/work/08209/brian97/ls6/tools/ncbi-blast-2.14.0+/bin/blastn -query raw_sequences/c_cristata_AY443280.1.fasta -db rag1/p_pica_au_db -out rag1/p_pica_rag1_blast.out

#nano rag1/p_pica_rag1_blast.out
#JAOYNA010000039.1:44057536-44060407

#now create a fasta of just this region from A. californica 

conda activate samtools
samtools faidx raw_sequences/p_pica_au.fasta JAOYNA010000039.1:44057536-44060407 > rag1/p_pica_rag1.fasta
```


Make *C. yncas* slightly differently
```
#call high confidencevariants at the RAG1 gene from our previously sorted bam files
bcftools mpileup -Ou -q 20 -r JANXIP010000005.1:37834112-37836993 -f raw_sequences/c_stelleri_au.fasta prep_mt/sort.cy01.bam prep_mt/sort.cy02.bam prep_mt/sort.cy03.bam prep_mt/sort.cy04.bam | bcftools call -Ou -mv | bcftools norm -f raw_sequences/c_stelleri_au.fasta -Oz -o rag1/c_yncas_rag1.vcf.gz 
bcftools index rag1/c_yncas_rag1.vcf.gz

#create consensus sequence at the same homologous region of C. stelleri because that served as the C. yncas reference
conda deactivate
conda activte gatk4

gatk IndexFeatureFile -I rag1/c_yncas_rag1.vcf.gz

gatk FastaAlternateReferenceMaker -R raw_sequences/c_stelleri_au.fasta -O rag1/c_yncas_rag1.fasta -V rag1/c_yncas_rag1.vcf.gz -L JANXIP010000005.1:37834122-37836993

```

Blastoff

```

####Put all RAG1's into a single fasta

cat raw_sequences/c_cristata_AY443280.1.fasta rag1/c_yncas_rag1.fasta rag1/c_stelleri_rag1.fasta rag1/a_californica_rag1.fasta rag1/p_pica_rag1.fasta > rag1/rag1.fasta

/work/08209/brian97/ls6/tools/ncbi-blast-2.14.0+/bin/makeblastdb -in rag1/rag1.fasta -out rag1/rag1_db -dbtype nucl -title rag1/rag1_db


/work/08209/brian97/ls6/tools/ncbi-blast-2.14.0+/bin/blastn -query /scratch/08209/brian97/hybrid/autosomes/homolog/hyb_rag1.fasta -db rag1/rag1_db -out rag1/rag1_blast.out



```






