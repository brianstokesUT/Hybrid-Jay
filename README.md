# Hybrid-Jay (*REPLACE WITH PUBLICATION NAME*)

*CITATION HERE*

## Summary
<p align="center">
  <img src="https://github.com/brianstokesUT/Hybrid-Jay/assets/91159511/e082186b-21ec-4f59-b68e-81bb3b770686" width="250">
</p>

We found a putative hybrid Green Jay x Blue Jay individual in central Texas during the Summer of 2023. Using bioinformatic methods we desscribe how we determined ancestry. 

This repository includes all code necessary to produce results of *INSERT CITATION HERE* . The reposity first describes how we determined genetic ancestry. Second we describe the generation of species range maps shown in *FIGURE NUMBER*. 

1. Genomic Methods
   - First nested list item
     - Second nested list item
2. Range Models
   - First nested list item
     - Second nested list item



# Genomic Methods
##Summary
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

Each candidate species uses data available on NCBI - data collection is described breifly within appropriate sections and may be split between mitochondria and chromosomal sources.



## Blue Jay (*Cyanocitta cristata*)
All databases are limited by what is avaialable on NCBI for Blue Jay, potential sequences are listed in the subsections below.

### Blue Jay Mitochondria
#### Data Sources

| **GenBank**    | **GI**         | **Length (bp)** | **Gene**                 |
|------------|------------|--------|----------------------|
| AF197832.1 | 6469720    | 1020   | MT-CO1               |
| AF218921.1 | 6979773    | 1293   | GB072 Control Region |
| AY509629.1 | 45655454   | 358    | MT-CYB               |
| AY666449.1 | 51101994   | 658    | MT-CO1               |
| AY666459.1 | 51102014   | 657    | MT-CO1               |
| AY666552.1 | 51102200   | 692    | MT-CO1               |
| AY666576.1 | 116832127  | 651    | MT-CO1               |
| DQ432877.1 | 116832127  | 651    | MT-CO1               |
| DQ434557.1 | 117372137  | 697    | MT-CO1               |
| DQ434558.1 | 117372139  | 697    | MT-CO1               |
| DQ434559.1 | 117372141  | 667    | MT-CO1               |
| DQ912604.1 | 114214838  | 1041   | MT-ND2               |
| KY492618.1 | 1214941616 | 354    | MT-ND3               |
| KY495244.1 | 1214941702 | 1013   | MT-ND2               |
| OM817599.1 | 2321932993 | 1041   | MT-ND2               |
| OM817600.1 | 2321932995 | 1041   | MT-ND2               |
| OQ002164.1 | 2412883929 | 1560   | MT-CO1               |
| OQ035627.1 | 2445559238 | 519    | MT-ND6               |
| OQ036935.1 | 2445561854 | 684    | MT-ATP6              |
| OQ038053.1 | 2445564090 | 168    | MT-ATP8              |
| OQ038973.1 | 2445565930 | 684    | MT-CO2               |
| OQ039913.1 | 2445567810 | 1143   | MT-CYB               |
| OQ041035.1 | 2445570054 | 978    | MT-ND1               |
| OQ042157.1 | 2445572298 | 351    | MT-ND3               |
| OQ043368.1 | 2445574720 | 297    | MT-ND4L              |
| OQ044490.1 | 2445576964 | 1818   | MT-ND5               |
| OQ045108.1 | 2445578200 | 785    | MT-CO3               |
| OQ046820.1 | 2445581624 | 1040   | MT-ND2               |
| OQ048151.1 | 2445584286 | 1378   | MT-ND4               |
| X74258.1   | 396719     | 1143   | MT-CYB               |

Because mitochondria have fairly low levels of intraspecific variation (at least compared to interspecific variation) we will create a consensus Blue Jay mitochondrial sequence based on the available Steller's Jay refernece genome which has a high quality mitochondrial scaffold.
1. Download all the sequences listed in the table above from NCBI using prefered method.
2. Align Blue Jay MT sequences to the Steller's Jay refernce MT scaffold.

``` 
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
```
3. Take a quick look at the quality and mapping depth of the resultant alignment.
```
samtools sort cya_cri.bam -o sort_cya_cri.bam
samtools index sort_cya_cri.bam 

samtools depth -a sort_cya_cri.bam > depth.txt

###Produce quick quality report of the alignments so we have the stats for paper
conda deactivate
conda activate qualimap
qualimap bamqc -bam sort_cya_cri.bam -outdir qualimap_results
###75.04% of the ref mt genome is covered by at least one mapped basepair
```
4. Call variants and create a consensus sequence
```
#lets make variant calls, don't need any quality filter because these are NCBI and should all be high quality?
conda deactivate
conda activate bcftools

bcftools mpileup -Ou -f /scratch/08209/brian97/hybrid/mt_raw/cya_ste/JANXIP010000352.1.fasta sort_cya_cri.bam | bcftools call -Ou -mv | bcftools norm -f /scratch/08209/brian97/hybrid/mt_raw/cya_ste/JANXIP010000352.1.fasta -Oz -o cyacri_MT.vcf.gz
bcftools index cyacri_MT.vcf.gz
bcftools consensus -f /scratch/08209/brian97/hybrid/mt_raw/cya_ste/JANXIP010000352.1.fasta cyacri_MT.vcf.gz > cyacri_mt_con1.fasta
```

### Blue Jay Autosomes

| GenBank    | GI         | Length (bp) | Gene                                    |
|------------|------------|-------------|-----------------------------------------|
| AY082422.1 | 28916038   | 852         | beta-fibrinogen, intron 7               |
| AY395603.1 | 39979941   | 891         | fibrinogen, intron 7                    |
| AY443137.1 | 38324817   | 1152        | RAG2                                    |
| AY443280.1 | 38324377   | 2872        | RAG1                                    |
| DQ320585.1 | 85822649   | 559         | beta-fibrinogen, intron 5               |
| DQ912621.1 | 114214869  | 507         | AK1, intron 5                           |
| DQ912641.1 | 114214889  | 853         | beta-fibrinogen, intron 7               |
| FJ598306.1 | 224593463  | 592         | TGFB2, intron 5                         |
| HM623931.1 | 339773338  | 74          | short-wavelength-sensitive 1 opsin TRM2 |
| HQ391560.1 | 319739720  | 620         | CHD1Z                                   |
| KY492640.1 | 1214941656 | 670         | myoglobin, intron 2                     |
| KY492659.1 | 1214941675 | 505         | TGFB2, intron 5                         |
| KY495265.1 | 1214941736 | 732         | fibrinogen, intron 7                    |










## Green Jay (*Cyanocorax yncas*)
We used WGS data from 7 Green Jay samples captured throughout Texas. Samples were sequenced using PE, 150bp, 10x coverage. We aligned these to a Steller's Jay Refernce Genome for comparison with hybrid sequencing data.

### Green Jay Mitochondria

'''
  this is code text?
'''


### Green Jay Autosomes

##########Create GRJA RAG1 gene fasta

```
#call variants at the RAG1 gene
bcftools mpileup -Ou -q 20 -r JANXIP010000005.1:37834112-37836993 -f /scratch/08209/brian97/hybrid/mt_raw/grja_mt_prep/stja_ref/GCA_026167965.1_bCyaSte1.0.p_genomic.fna sort.2000.bam sort.2001.bam sort.2002.bam sort.2003.bam sort.2004.bam sort.2005.bam | bcftools call -Ou -mv | bcftools norm -f /scratch/08209/brian97/hybrid/mt_raw/grja_mt_prep/stja_ref/GCA_026167965.1_bCyaSte1.0.p_genomic.fna -Oz -o /scratch/08209/brian97/hybrid/autosomes/homolog/grja_s5_rag1.vcf.gz
bcftools index /scratch/08209/brian97/hybrid/autosomes/homolog/grja_s5_rag1.vcf.gz

#create conesnsus, for just the well covered RAG1 gene region
samtools faidx /scratch/08209/brian97/hybrid/mt_raw/grja_mt_prep/stja_ref/GCA_026167965.1_bCyaSte1.0.p_genomic.fna JANXIP010000005.1:37834112-37836993 | bcftools consensus /scratch/08209/brian97/hybrid/autosomes/homolog/grja_s5_rag1.vcf.gz > /scratch/08209/brian97/hybrid/autosomes/homolog/grja_s5_rag1.fasta
```






## Steller's Jay (*Cyanocitta stelleri*)
We used a Steller's Jay genome assembly (GCA_026167965.1)

### Steller's Jay Mitochondria

### Steller's Jay Autosomes
RAG1:
```




```




