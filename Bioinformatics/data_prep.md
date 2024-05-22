# Hybrid 
```
###HYBRID
#whatever data download steps from NCBI go here


conda create --name trim-galore
conda activate trim-galore
conda install bioconda::trim-galore

#combine our two fastq files
cat raw_hybrid/*.fastq.gz > raw_hybrid/raw_hybrid.fastq.gz

#trim adapters and low quality reads with default setting (trim_galore)
trim_galore raw_hybrid/raw_hybrid.fastq.gz -o prep_hybrid/
```




# data_prep.bash
```
###P. pica
#P. pica mitochondrial
esearch -db nucleotide -query "MT792356 [ACCN]" | efetch -format fasta > p_pica_mt.fasta
#P. pica autosomal 
datasets download genome accession GCA_025802055.1 --include gff3,rna,cds,protein,genome,seq-report --filename p_pica.zip
unzip p_pica.zip -d p_pica
#move P. pica sequences
cp p_pica/ncbi_dataset/data/GCA_025802055.1/GCA_025802055.1_ASM2580205v1_genomic.fna p_pica_au.fasta


###C. stelleri 
datasets download genome accession GCA_026167965.1 --include gff3,rna,cds,protein,genome,seq-report --filename c_stelleri.zip
unzip c_stelleri.zip -d c_stelleri
#extract C. stelleri mitochondrial sequence
samtools faidx c_stelleri/ncbi_dataset/data/GCA_026167965.1/GCA_026167965.1_bCyaSte1.0.p_genomic.fna JANXIP010000352 > c_selleri_mt.fasta
#move C. stelleri sequences
cp p_pica/ncbi_dataset/data/GCA_025802055.1/GCA_025802055.1_ASM2580205v1_genomic.fna p_pica_au.fasta

###A. woodhousei 
datasets download genome accession GCA_028536675.1 --include gff3,rna,cds,protein,genome,seq-report --filename a_californica.zip
unzip a_californica.zip -d a_californica
#extract A. californica  mitochondrial sequence
samtools faidx a_californica/ncbi_dataset/data/GCA_028536675.1/GCA_028536675.1_bAphCal1.0.hap1_genomic.fna JAQMYR010001393.1 > a_californica_mt.fasta
#move A. californica sequences
cp a_californica/ncbi_dataset/data/GCA_028536675.1/GCA_028536675.1_bAphCal1.0.hap1_genomic.fna a_californica_au.fasta


###C. cristata 
#C. cristata mitochondrial
efilter -db nucleotide -query "txid28727[organism:exp] AND mitochondrion[filter]" | efetch -format fasta > c_cristata_mt.fasta
#C. cristata autosomal
efilter -db nucleotide -query "txid28727[organism:exp] NOT mitochondrion[filter]" | efetch -format fasta > c_cristata_au.fasta
```


```
###C. yncas
#C. yncas mitochondrial
#C. yncas autosomal
```

