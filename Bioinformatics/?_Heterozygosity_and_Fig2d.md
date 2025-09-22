# Download New Data for Heterozygosity Analysis

*note* we use C. monedula* becasue the *P. pica* assembly lacked an alternate haplotype.
```
mkdir het

###Download C_cristata
datasets download genome accession GCA_046129655.1 --include gff3,rna,cds,protein,genome,seq-report --filename het/c_cristata_hap1.zip
unzip het/c_cristata_hap1.zip -d het/c_cristata_hap1

datasets download genome accession GCA_046129645.1 --include gff3,rna,cds,protein,genome,seq-report --filename het/c_cristata_hap2.zip
unzip het/c_cristata_hap2.zip -d het/c_cristata_hap2

###Download A_californica
datasets download genome accession GCA_028536675.1 --include gff3,rna,cds,protein,genome,seq-report --filename het/a_califonica_hap1.zip
unzip het/a_califonica_hap1.zip -d het/a_califonica_hap1

datasets download genome accession GCA_028536645.1 --include gff3,rna,cds,protein,genome,seq-report --filename het/a_califonica_hap2.zip
unzip het/a_califonica_hap2.zip -d het/a_califonica_hap2

###Download C_stelleri
datasets download genome accession GCA_026167965.1 --include gff3,rna,cds,protein,genome,seq-report --filename het/c_stelleri_hap1.zip
unzip het/c_stelleri_hap1.zip -d het/c_stelleri_hap1

datasets download genome accession GCA_026168045.1 --include gff3,rna,cds,protein,genome,seq-report --filename het/c_stelleri_hap2.zip
unzip het/c_stelleri_hap2.zip -d het/c_stelleri_hap2

###Download C_monedula
datasets download genome accession GCA_965178545.1 --include gff3,rna,cds,protein,genome,seq-report --filename het/c_monedula_hap1.zip
unzip het/c_monedula_hap1.zip -d het/c_monedula_hap1

datasets download genome accession GCA_965178535.1 --include gff3,rna,cds,protein,genome,seq-report --filename het/c_monedula_hap2.zip
unzip het/c_monedula_hap2.zip -d het/c_monedula_hap2
```


