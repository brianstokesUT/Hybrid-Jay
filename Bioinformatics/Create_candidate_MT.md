# P. pica
Leave it as 
```
raw_sequences/p_pica_mt.fasta
```

# C. stelleri
Leave it as 
```
raw_sequences/c_stelleri_mt.fasta
```

# A. woodhousei
Leave it as 
```
raw_sequences/a_californica_mt.fasta
```

# C. Cristata
Leave it as 
```
raw_sequences/c_cristata_mt.fasta
```

# C. yncas

```
mamba activate data_prep
bcftools view prep_mt/c_yncas.vcf.gz --regions JANXIP010000352.1 -o prep_mt/c_yncas_mt.vcf.gz -Oz
bcftools index prep_mt/c_yncas_mt.vcf.gz

##################

mamba deactivate
mamba activate gatk4

gatk CreateSequenceDictionary -R raw_sequences/c_stelleri_au.fasta

gatk IndexFeatureFile -I prep_mt/c_yncas_mt.vcf.gz

gatk FastaAlternateReferenceMaker -R raw_sequences/c_stelleri_au.fasta -O c_yncas_mt_FINAL.fasta -V prep_mt/c_yncas_mt.vcf.gz -L JANXIP010000352.1

mamba deactivate
```


# Create MT Blast database

```
# Tidy up our previously made *mt.fasta sequences
mkdir mt_blast
cp *mt_FINAL* ./mt_blast
cp raw_sequences/*mt.fasta ./mt_blast

###edit names of c_yncas & c_cristata because we used the c_stelleri refernce to make the fastas
sed -i '1s/.*/>c_yncas_mt/' mt_blast/c_yncas_mt_FINAL.fasta
sed -i '1s/.*/>c_cristata_mt/' mt_blast/c_cristata_mt_FINAL.fasta


#create our database by combining all fastas into a single fasta
cat mt_blast/c_yncas_mt_FINAL.fasta mt_blast/c_cristata_mt_FINAL.fasta mt_blast/a_californica_mt.fasta mt_blast/c_stelleri_mt.fasta mt_blast/p_pica_mt.fasta > mt_blast/jayz_mt.fasta

#make database
~PATH/tools/ncbi-blast-2.14.0+/bin/makeblastdb -in mt_blast/jayz_mt.fasta -out mt_blast/jayz_mt_db -dbtype nucl -title jayz_mt_db

```

# Query database and check results
```
~PATHtools/ncbi-blast-2.14.0+/bin/blastn -query hyb_mt_FINAL.fasta -db mt_blast/jayz_mt_db -out jay_mt_blast.out

#VIEW RESULTS BY REMOVING "#" IN LINE BELOW
#nano jay_mt_blast.out
```
