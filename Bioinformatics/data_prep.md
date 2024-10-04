# Tools and mamba environments

Download the BLAST+ tool based on the [directions found on NCBI](https://www.ncbi.nlm.nih.gov/books/NBK569861/#intro_Installation.RedHat_Linux) to a folder titled "~PATH/tools/"


# Next ensure that your machine can run mamba and create the following environemnts

```
# Make necessary conda channels available
mamba config --add channels bioconda
mamba config --add channels conda-forge


# Create environment "data_prep" and install associated packages
mamba create --name data_prep trim-galore=0.6.10 entrez-direct=22.1 samtools=1.20-0 ncbi-datasets-cli=16.22.1 bowtie2=2.5.4 bcftools=1.20 shapeit5=5.1.1 picard=3.1.1

# Create environment "gatk4" and install associated packages
# this version of gatk4 seemed to create conflict when installed to data_prep environment
mamba create --name gatk4 gatk4=4.5.0.0
```

# Hybrid 
Retrieval of SRA data by BioProject accession is not supported by NCBI Datasets at this time. Refer to the README.md to download relevant file. These files should be placed into ~PATH/raw_hybird/
```

mamba create --name trim-galore
mamba activate trim-galore
mamba install bioconda::trim-galore

#combine our two fastq files
cat raw_hybrid/*.fastq.gz > raw_hybrid/raw_hybrid.fastq.gz

#trim adapters and low quality reads with default setting (trim_galore)
trim_galore raw_hybrid/raw_hybrid.fastq.gz -o prep_hybrid/
mamba deactivate
```

# Download and basic prep of candidate data
```
###P. pica
#P. pica mitochondrial
esearch -db nucleotide -query "MT792356 [ACCN]" | efetch -format fasta > raw_sequences/p_pica_mt.fasta
#P. pica autosomal 
datasets download genome accession GCA_025802055.1 --include gff3,rna,cds,protein,genome,seq-report --filename raw_sequences/p_pica.zip
unzip raw_sequences/p_pica.zip -d raw_sequences/p_pica
#move P. pica sequences
cp raw_sequences/p_pica/ncbi_dataset/data/GCA_025802055.1/GCA_025802055.1_ASM2580205v1_genomic.fna raw_sequences/p_pica_au.fasta

###C. stelleri 
datasets download genome accession GCA_026167965.1 --include gff3,rna,cds,protein,genome,seq-report --filename raw_sequences/c_stelleri.zip
unzip raw_sequences/c_stelleri.zip -d raw_sequences/c_stelleri
#extract C. stelleri mitochondrial sequence
samtools faidx raw_sequences/c_stelleri/ncbi_dataset/data/GCA_026167965.1/GCA_026167965.1_bCyaSte1.0.p_genomic.fna JANXIP010000352.1 > raw_sequences/c_stelleri_mt.fasta
#move C. stelleri sequences
cp raw_sequences/c_stelleri/ncbi_dataset/data/GCA_026167965.1/GCA_026167965.1_bCyaSte1.0.p_genomic.fna raw_sequences/c_stelleri_au.fasta

###A. woodhousei 
datasets download genome accession GCA_028536675.1 --include gff3,rna,cds,protein,genome,seq-report --filename raw_sequences/a_californica.zip
unzip raw_sequences/a_californica.zip -d raw_sequences/a_californica
#extract A. californica  mitochondrial sequence
samtools faidx raw_sequences/a_californica/ncbi_dataset/data/GCA_028536675.1/GCA_028536675.1_bAphCal1.0.hap1_genomic.fna JAQMYR010001393.1 > raw_sequences/a_californica_mt.fasta
#move A. californica sequences
cp raw_sequences/a_californica/ncbi_dataset/data/GCA_028536675.1/GCA_028536675.1_bAphCal1.0.hap1_genomic.fna raw_sequences/a_californica_au.fasta

###C. cristata 
#C. cristata mitochondrial
efilter -db nucleotide -query "txid28727[organism:exp] AND mitochondrion[filter]" | efetch -format fasta > raw_sequences/c_cristata_mt.fasta
#C. cristata autosomal
efilter -db nucleotide -query "txid28727[organism:exp] NOT mitochondrion[filter]" | efetch -format fasta > raw_sequences/c_cristata_au.fasta
```
# Prep C. stelleri to serve as reference genome
```
bowtie2-build -f ~PATH/raw_sequences/c_stelleri/ncbi_dataset/data/GCA_026167965.1/GCA_026167965.1_bCyaSte1.0.p_genomic.fna raw_sequences/cs_ref
```

# Green Jay
Retrieval of SRA data by BioProject accession is not supported by NCBI Datasets at this time. Refer to the README.md to download relevant files including BioSample accessions SAMN44062996, SAMN44062997, SAMN44062998, SAMN44062999. Files should be placed into ~PATH/raw_sequences/c_yncas/

```
mamba activate data_prep

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

# Sort and index
samtools sort raw_sequences/c_yncas/cy01.bam -o prep_mt/sort.cy01.bam
samtools index prep_mt/sort.cy01.bam

samtools sort raw_sequences/c_yncas/cy02.bam -o prep_mt/sort.cy02.bam
samtools index prep_mt/sort.cy02.bam

samtools sort raw_sequences/c_yncas/cy03.bam -o prep_mt/sort.cy03.bam
samtools index prep_mt/sort.cy03.bam

samtools sort raw_sequences/c_yncas/cy04.bam -o prep_mt/sort.cy04.bam
samtools index prep_mt/sort.cy04.bam

# Create VCF
bcftools mpileup -Ou -f raw_sequences/c_stelleri_au.fasta prep_mt/sort.cy01.bam prep_mt/sort.cy02.bam prep_mt/sort.cy03.bam prep_mt/sort.cy04.bam | bcftools call -Ou -mv | bcftools norm -f raw_sequences/c_stelleri_au.fasta -Oz -o prep_mt/c_yncas.vcf.gz
bcftools index prep_mt/c_yncas.vcf.gz


```
