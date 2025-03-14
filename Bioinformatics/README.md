# Genomic Methods
## Summary
We used WGS data from the hybrid individual to determine ancestry using BLAST+ methodology. We assumed the majority of mitochondrial sequences were passed down by the hybrid's maternal species while autosomal sequences were passed by both maternal and paternal species.

All Bioinformatic processing was performed on [Texas Advanced Computing Center's (TACC) Lonestar6 cluster](https://tacc.utexas.edu/systems/lonestar6/) which is a Linux-64 based system. If you desire to run included code on your own machine or cluster you should adjust commands as necessary for your machine's opperating capacities.

## Conceptual Overview & Candidate Parents
Based on location and plumage morphology we assumed the individual was the offspring of a Green Jay and Blue Jay, but we considered all possible sources of parental ancestry within jay species found in the state of Texas to maintain unbiased analysis: 
+ Eurasian Magpie (*Pica Pica*) - **OUTGROUP**
+ Steller's Jay (*Cyanocitta stelleri*)
+ Woodhouse's Scrub Jay (*Aphelocoma woodhousei*)
+ Blue Jay (*Cyanocitta cristata*)
+ Green Jay (*Cyanocorax yncas*)

We aimed to compare reads from the hybrid individual against avaavailable representative sequences of the candidate species using BLAST+ to determine paternal ancestry. Because we used non-targeted Whole Genome Sequencing and some of the candidate parental species had limited sequencing data available, we needed to develop a method to compare homologous genes with eachother while minimizing sampling bias. 

Data collection and prep is described breifly within appropriate sections and may be split between mitochondria and autosomal sources.

# Hybrid Data Generation
Raw fastq files along with library prep/sequencing details of the putative hybrid invividual can be found in [NIH BioProject#1114044](http://www.ncbi.nlm.nih.gov/bioproject/1114044)


The raw reads are split into two fastq files named ```061323A_S1_L001_R1_001.fastq.gz``` and ```061323A_S1_L002_R1_001.fastq.gz``` 

   *note:* sequencing data for this individual is single-ended

## Hybrid Data Prep
We wanted to ensure homologous regions of each  genome were being compared during BLAST+ analysis so we first ensured the hybrid and candidate parental species sequence data were aligned to the same reference genome.

After merging and trimming the two raw hybrid fastq files, we aligned these to a Steller's Jay (*Cyanocitta stelleri*) refernece genome. We then created a representative genotype for both the hybrid mitochondrial and autosomal genomes by calling variants and then creating a individual-level masked consensus sequence where "*N*" is placed at basepair with insufficient high-confidence coverage of hybrid reads. This process results in a consensus fasta file where both REF and ALT (variants) are represented are presenent, while basepairs with low coverage are masked.



# Candidate Species Data Generation
## Eurasian Magpie (*Pica Pica*) Data Generation
### *P. pica* Mitochondrial Data
We used a *Pica pica melanotos* mitochondrion assembly: [Genbank assembly accession #MT792356.1](https://www.ncbi.nlm.nih.gov/nuccore/1899896744) ([Kryukov et al., 2020](https://doi.org/10.1080%2F23802359.2020.1838354)) to represent Eurasian Magpie mitochondrial genomes.

### *P. pica* Autosomal Data
We used [Genome assembly ASM2580205v1 (GenBank: GCA_025802055.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_025802055.1/) from the [WGS project JAOYNA01](https://www.ncbi.nlm.nih.gov/nuccore/JAOYNA000000000.1) to represent Eurasian Magpie autosomal genomes.


## Steller's Jay (*Cyanocitta stelleri*) Data Generation
We used a *Cyanocitta stelleri* mitochondrion assembly: [Genbank assembly accession #bCyaSte1.0.p](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_026167965.1/) ([Benham et al., 2023](https://doi.org/10.1093/jhered/esad042)) from [WGS project JANXIP01](https://www.ncbi.nlm.nih.gov/nuccore/JANXIP000000000.1) to represent *C. stelleri* autosomal and mitochondrial genomes.


## Woodhouse's Scrub Jay (*Aphelocoma woodhousei*) Data Generation
At the time of analysis no *A. woodhousei* sequencing data was available on any public repository. We used a Western Scrub Jay (*Aphelocoma californica*) assembly: [Genome assembly bAphCal1.0.hap1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_028536675.1/) ([DeRaad et al., 2023](https://doi.org/10.1093%2Fjhered%2Fesad047)) from [BioProject PRJNA904314](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA904314/) as representative for both *A. woodhousei* autosomal and mitochondrial genomes.

## Blue Jay (*Cyanocitta cristata*) Data Generation
We used a *Cyanocitta cristata* mitochondrion assembly: [Genbank assembly accession bCyaCrs1.hap1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_046129655.1/) ([Rhie et al., 2021](https://doi.org/10.1038/s41586-021-03451-0)) from [BioProject PRJNA1181931](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1181931/) to represent *C. cristata* autosomal and mitochondrial genomes.


## Green Jay (*Cyanocorax yncas*)Data Generation
We used WGS data from 4 Green Jay samples captured throughout Texas. Samples were sequenced using PE, 150bp, 10x coverage. We aligned these to a Steller's Jay Refernce Genome for comparison with hybrid sequencing data.

These four samples are stored in Bioproject [PRJNA1168985: WGS of Green Jay (Cyanocorax yncas) in Texas](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1168985)



