# Hybrid-Jay (*REPLACE WITH PUBLICATION NAME*)

*CITATION HERE*

## Summary
<p align="center">
  <img src="https://github.com/brianstokesUT/Hybrid-Jay/assets/91159511/e082186b-21ec-4f59-b68e-81bb3b770686" width="250">
</p>

We found a putative hybrid Green Jay x Blue Jay individual in central Texas during the Summer of 2023. Using bioinformatic methods we desscribe how we determined parental ancestry. 

This repository includes all code necessary to produce results of *INSERT CITATION HERE* . The reposity first describes how we determined genetic ancestry and generation of *figure number*. Second we describe the generation of species range maps and associated plots shown in *FIGURE 2*. 

1. Bioinformatics
   - README.md
     - Describes concetual overview of bioinformatic methodology and data sources.
   - setup.md
   - dataprep.md
     - Describes how to download and prepare raw sequence data via command line for mitochondrial and autosomal analyses
   - Create_Hybrid_Consensus_MT.md
   - Create_Hybrid_Consensus_AU.md
   - Create_candidate_MT.md
   - Create_candidate_AU.md
   - phasing.md
2. Range Models
   - README.md
     - described datasources for Range Models
   - bioclim_dataprep.R
     - Preperation of BioClim data for our climate model, some steps may not be necessary if you have a sufficiently powerful computer
   - ebird_and_maxent.R
     - Preperation of eBird data and production of maxent models for *FIGURE 2*
   - Figures.R
     - Code to produce each portion of *FIGURE 2*


# Genomic Methods and Data Sources
We used WGS data from the hybrid individual to determine paternal ancestry based on BLAST methodology. We assumed the majority of mitochondrial sequences were passed down by the hybrid's maternal species while autosomal sequences were passed by both maternal and paternal speices.


## Data Generation
Raw fastq files along with library prep/sequencing details of the putative hybrid invividual can be found [in NIH BioProject#1114044](http://www.ncbi.nlm.nih.gov/bioproject/1114044)

### Possible Parents
During analysis we considered all possible sources of ancestry within potential jay species found in the state of Texas: 
+ Green Jay (*Cyanocorax yncas*)
+ Blue Jay (*Cyanocitta cristata*)
+ Steller's Jay (*Cyanocitta stelleri*)
+ Woodhouse's Scrub Jay (*Aphelocoma woodhousei*)
+ Eurasian Magpie (*Pica Pica*) - **OUTGROUP**

Each candidate species uses data available on NCBI - data collection and prep is described within the Bioinformatics directory.
