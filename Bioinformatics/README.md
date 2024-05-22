# Genomic Methods
## Summary
We used WGS data from the hybrid individual to determine ancestry using BLAST+ methodology. We assumed the majority of mitochondrial sequences were passed down by the hybrid's maternal species while autosomal sequences were passed by both maternal and paternal species.


## Hybrid Data Generation
Raw fastq files along with library prep/sequencing details of the putative hybrid invividual can be found in [NIH BioProject#1114044](http://www.ncbi.nlm.nih.gov/bioproject/1114044)


The raw reads are split into two fastq files named ```061323A_S1_L001_R1_001.fastq.gz``` and ```061323A_S1_L002_R1_001.fastq.gz``` 


### Candidate Parents
Based on location and plumage morphology we assumed the indivudal was resultant of the breeding between a Green Jay and Blue Jay, but we considered all possible sources of ancestry within potential jay species found in the state of Texas to maintain unbiased analysis: 
+ Eurasian Magpie (*Pica Pica*) - **OUTGROUP**
+ Steller's Jay (*Cyanocitta stelleri*)
+ Woodhouse's Scrub Jay (*Aphelocoma woodhousei*)
+ Blue Jay (*Cyanocitta cristata*)
+ Green Jay (*Cyanocorax yncas*)


Reads from the hybrid individual were compared against available sequences of the candidate species in NCBI Databases using BLAST+. Data collection and prep is described breifly within appropriate sections and may be split between mitochondria and autosomal sources.

## Eurasian Magpie (*Pica Pica*) Data Generation
### *P. pica* Mitochondrial Data
We used a *Pica pica melanotos* mitochondrion assembly: [Genbank assembly accession #MT792356.1](https://www.ncbi.nlm.nih.gov/nuccore/1899896744) ([Kryukov et al., 2020](https://doi.org/10.1080%2F23802359.2020.1838354)) to represent Eurasian Magpie mitochondrial genomes.

### *P. pica* Autosomal Data
We used [Genome assembly ASM2580205v1 (GenBank: GCA_025802055.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_025802055.1/) from the [WGS project JAOYNA01](https://www.ncbi.nlm.nih.gov/nuccore/JAOYNA000000000.1) to represent Eurasian Magpie autosomal genomes.


## Steller's Jay (*Cyanocitta stelleri*) Data Generation
We used a *Cyanocitta stelleri* mitochondrion assembly: [Genbank assembly accession #bCyaSte1.0.p](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_026167965.1/) ([Benham et al., 2023](https://doi.org/10.1093/jhered/esad042)) from [WGS project JANXIP01](https://www.ncbi.nlm.nih.gov/nuccore/JANXIP000000000.1) to represent both autosomal and mitochondrial genomes.


## Woodhouse's Scrub Jay (*Aphelocoma woodhousei*) Data Generation
At the time of analysis no *A. woodhousei* sequencing data was available on any public repository. We used a Western Scrub Jay (*Aphelocoma californica*) assembly: [Genome assembly bAphCal1.0.hap1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_028536675.1/) ([DeRaad et al., 2023](https://doi.org/10.1093%2Fjhered%2Fesad047)) from [BioProject PRJNA904314](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA904314/) as representative for both *A. woodhousei* autosomal and mitochondrial genomes.

## Blue Jay (*Cyanocitta cristata*) Data Generation
At the time of analysis no *C. cristata* assembly was available on any public repository.
### *C. cristata* Mitochondrial Data
#### Data Sources

The table below describe the mitochondrial sequences of *C. cristata* avaiable within NCBI datasets at the time of analysis.

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

### *C. cristata* Autosomal Data
#### Data Sources

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


## Green Jay (*Cyanocorax yncas*)Data Generation
At the time of analysis no *C. yncas* assembly was available on any public repository.
### *C. yncas* Mitochondrial Data

### *C. yncas* Autosomal Data



