# Download New Data for Heterozygosity Analysis (ACTB GENE ONLY)

*note* we use *C. monedula* becasue the *P. pica* assembly lacked an alternate haplotype.
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






# Find the ACTB Gene for each reference genome based on a *A_coerulescens* refernce genome (becasue its the only one with full annotation)
We pulled this region directly from NCBI and the fasta file is found within the Bioinformatics directory (a_coerulescens_ACTB.fasta) but can be downloaded from NCBI direclty if you wish - just pull this region: NC_091028.1:11281498-11286329

```
###Download C_coerulescens and pull the ACTB Region


###make BLAST database
~PATH/tools/ncbi-blast-2.14.0+/bin/makeblastdb -in ACTB/a_coerulescens_ACTB.fasta -out ACTB/c_coerulescens_ACTB_db -dbtype nucl -title ACTB/c_coerulescens_ACTB_db


###c_cristata
~PATH/tools/ncbi-blast-2.14.0+/bin/blastn -query het/c_cristata_hap1/ncbi_dataset/data/GCA_046129655.1/GCA_046129655.1_bCyaCrs1.hap1_genomic.fna -db ACTB/a_coerulescens_ACTB_db -out ACTB/c_cristata_hap1_ACTB.out

~PATH/tools/ncbi-blast-2.14.0+/bin/blastn -query het/c_cristata_hap2/ncbi_dataset/data/GCA_046129645.1/GCA_046129645.1_bCyaCrs1.hap2_genomic.fna -db ACTB/a_coerulescens_ACTB_db -out ACTB/c_cristata_hap2_ACTB.out

###a_californica
~PATH/tools/ncbi-blast-2.14.0+/bin/blastn -query het/a_califonica_hap1/ncbi_dataset/data/GCA_028536675.1/GCA_028536675.1_bAphCal1.0.hap1_genomic.fna -db ACTB/a_coerulescens_ACTB_db -out ACTB/a_califonica_hap1_ACTB.out

~PATH/tools/ncbi-blast-2.14.0+/bin/blastn -query het/a_califonica_hap2/ncbi_dataset/data/GCA_028536645.1/GCA_028536645.1_bAphCal1.0.hap2_genomic.fna -db ACTB/a_coerulescens_ACTB_db -out ACTB/a_califonica_hap2_ACTB.out

###c_stelleri
~PATH/tools/ncbi-blast-2.14.0+/bin/blastn -query het/c_stelleri_hap1/ncbi_dataset/data/GCA_026167965.1/GCA_026167965.1_bCyaSte1.0.p_genomic.fna -db ACTB/a_coerulescens_ACTB_db -out ACTB/c_stelleri_hap1_ACTB.out

~PATH/tools/ncbi-blast-2.14.0+/bin/blastn -query het/c_stelleri_hap2/ncbi_dataset/data/GCA_026168045.1/GCA_026168045.1_bCyaSte1.0.a_genomic.fna -db ACTB/a_coerulescens_ACTB_db -out ACTB/c_stelleri_hap2_ACTB.out

###c_monedula
~PATH/tools/ncbi-blast-2.14.0+/bin/blastn -query het/c_monedula_hap1/ncbi_dataset/data/GCA_965178545.1/GCA_965178545.1_bColMon1.hap1.1_genomic.fna -db ACTB/a_coerulescens_ACTB_db -out ACTB/c_monedula_hap1_ACTB.out

~PATH/tools/ncbi-blast-2.14.0+/bin/blastn -query het/c_monedula_hap2/ncbi_dataset/data/GCA_965178535.1/GCA_965178535.1_bColMon1.hap2.1_genomic.fna -db ACTB/a_coerulescens_ACTB_db -out ACTB/c_monedula_hap2_ACTB.out
```

# Extract the proper region from each reference haplotype
```
###Index each fasta file
samtools faidx het/c_cristata_hap1/ncbi_dataset/data/GCA_046129655.1/GCA_046129655.1_bCyaCrs1.hap1_genomic.fna
samtools faidx het/c_cristata_hap2/ncbi_dataset/data/GCA_046129645.1/GCA_046129645.1_bCyaCrs1.hap2_genomic.fna

samtools faidx het/a_califonica_hap1/ncbi_dataset/data/GCA_028536675.1/GCA_028536675.1_bAphCal1.0.hap1_genomic.fna
samtools faidx het/a_califonica_hap2/ncbi_dataset/data/GCA_028536645.1/GCA_028536645.1_bAphCal1.0.hap2_genomic.fna

samtools faidx het/c_stelleri_hap1/ncbi_dataset/data/GCA_026167965.1/GCA_026167965.1_bCyaSte1.0.p_genomic.fna
samtools faidx het/c_stelleri_hap2/ncbi_dataset/data/GCA_026168045.1/GCA_026168045.1_bCyaSte1.0.a_genomic.fna

samtools faidx het/c_monedula_hap1/ncbi_dataset/data/GCA_965178545.1/GCA_965178545.1_bColMon1.hap1.1_genomic.fna
samtools faidx het/c_monedula_hap2/ncbi_dataset/data/GCA_965178535.1/GCA_965178535.1_bColMon1.hap2.1_genomic.fna

###extract sections based on blast results and store as new fastas
samtools faidx het/c_cristata_hap1/ncbi_dataset/data/GCA_046129655.1/GCA_046129655.1_bCyaCrs1.hap1_genomic.fna  CM100545.1:13429427-13436266 > ACTB/c_cristata_hap1_ACTB.fasta

samtools faidx het/c_cristata_hap2/ncbi_dataset/data/GCA_046129645.1/GCA_046129645.1_bCyaCrs1.hap2_genomic.fna   CM100589.1:6002322-6009156 > ACTB/c_cristata_hap2_ACTB.fasta

samtools faidx het/a_califonica_hap1/ncbi_dataset/data/GCA_028536675.1/GCA_028536675.1_bAphCal1.0.hap1_genomic.fna   JAQMYR010000017.1:13706762-13713597 > ACTB/a_californica_hap1_ACTB.fasta

samtools faidx het/a_califonica_hap2/ncbi_dataset/data/GCA_028536645.1/GCA_028536645.1_bAphCal1.0.hap2_genomic.fna JAQMYS010000017.1:10816510-10823346 > ACTB/a_californica_hap2_ACTB.fasta

samtools faidx het/c_stelleri_hap1/ncbi_dataset/data/GCA_026167965.1/GCA_026167965.1_bCyaSte1.0.p_genomic.fna JANXIP010000018.1:6039845-6046687 > ACTB/c_stelleri_hap1_ACTB.fasta

samtools faidx het/c_stelleri_hap2/ncbi_dataset/data/GCA_026168045.1/GCA_026168045.1_bCyaSte1.0.a_genomic.fna JANXIQ010000017.1:12376919-12383759 > ACTB/c_stelleri_hap2_ACTB.fasta

samtools faidx het/c_monedula_hap1/ncbi_dataset/data/GCA_965178545.1/GCA_965178545.1_bColMon1.hap1.1_genomic.fna OZ238480.1:6237410-6242272 > ACTB/c_monedula_hap1_ACTB.fasta

samtools faidx het/c_monedula_hap2/ncbi_dataset/data/GCA_965178535.1/GCA_965178535.1_bColMon1.hap2.1_genomic.fna OZ238440.1:11081779-11086643 > ACTB/c_monedula_hap2_ACTB.fasta
```

# We have to handle *C_yncas* uniquely becasue its not a reference genome
```
# Reference genome and region
REF="raw_sequences/c_stelleri_au.fasta"
REGION="JANXIP010000018.1:6040845-6045687"

# Output directory
OUTDIR="ACTB"
mkdir -p $OUTDIR

# Samples
for SAMPLE in cy01 cy02 cy03 cy04; do
    echo ">>> Processing $SAMPLE"

    # Add read groups
    picard AddOrReplaceReadGroups \
        I=prep_mt/sort.${SAMPLE}.bam \
        O=prep_mt/sort.${SAMPLE}.withRG.bam \
        RGID=${SAMPLE} \
        RGLB=lib1 \
        RGPL=illumina \
        RGPU=unit1 \
        RGSM=${SAMPLE}

    # Index BAM
    samtools index prep_mt/sort.${SAMPLE}.withRG.bam

    # Call variants in region
    bcftools mpileup -Ou -f $REF \
        --regions $REGION \
        prep_mt/sort.${SAMPLE}.withRG.bam | \
        bcftools call -Ou -m | \
        bcftools view -Oz -o $OUTDIR/${SAMPLE}_ACTB_region.vcf.gz

    # Normalize variants
    bcftools norm -m +any -f $REF \
        $OUTDIR/${SAMPLE}_ACTB_region.vcf.gz -Oz -o $OUTDIR/${SAMPLE}_ACTB_region_normalized.vcf.gz

    # Index normalized VCF (useful for bcftools query)
    bcftools index -f $OUTDIR/${SAMPLE}_ACTB_region_normalized.vcf.gz

    # Count variants (sanity check)
    bcftools +counts $OUTDIR/${SAMPLE}_ACTB_region_normalized.vcf.gz

    # Calculate heterozygosity
    OUTPUT="${OUTDIR}/${SAMPLE}_heterozygosity_summary.txt"
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' $OUTDIR/${SAMPLE}_ACTB_region_normalized.vcf.gz | \
    awk 'BEGIN { heterozygous=0; total=0 } 
    {
        total++; 
        if ($5 == "1/0" || $5 == "0/1") {
            heterozygous++
        }
    } 
    END { 
        if (total > 0) {
            heterozygosity = heterozygous / total
            print "Sample:", "'$SAMPLE'"
            print "Total sites:", total
            print "Heterozygous sites:", heterozygous
            print "Heterozygosity:", heterozygosity
        } else {
            print "Sample:", "'$SAMPLE'"
            print "No variants found in the VCF file"
        }
    }' > $OUTPUT

    cat $OUTPUT
    echo ">>> Finished $SAMPLE"
    echo
done

```
