#FOR BRIAN
## Export and Create by requirements.txt

**Export requirements.txt**

```
$ conda list -e > requirements.txt
```

**Create Environment by requirements.txt**

```
$ conda create [-n <environment-name>|-p /path/to/your/environment] --file requirements.txt
```




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
