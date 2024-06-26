#FOR BRIAN (Need to have setup for each folder eventually used, all environments, and other tools used (BLAST, etc)
## Export and Create by requirements.txt

**Export requirements.txt**

```
$ conda list -e > requirements.txt
```

**Create Environment by requirements.txt**

```
$ conda create [-n <environment-name>|-p /path/to/your/environment] --file requirements.txt
```


CONDA ENVS

conda config --add channels bioconda
conda config --add channels conda-forge


mamba create --name data_prep 
mamba activate data_prep
mamba install trim-galore=0.6.10 entrez-direct=22.1 samtools=1.20-0 ncbi-datasets-cli=16.22.1 bowtie2=2.5.4-1




# SETUP


##Before data_prep.bash
```
conda create --name data_prep
conda activate data_prep
conda install bioconda::entrez-direct
conda install -c conda-forge ncbi-datasets-cli
conda install bioconda::samtools

mkdir raw_sequences
```


# ALIGNMENT

```
conda creat --name alignment
conda activate alignment
conda install bioconda::bowtie2

```
