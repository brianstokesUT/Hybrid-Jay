mkdir phasing

bcftools mpileup -Ou -q 20 -r JANXIP010000005.1 -f raw_sequences/c_stelleri_au.fasta prep_hybrid/sort.hyb.bam | bcftools call -Ou -c | bcftools norm -f raw_sequences/c_stelleri_au.fasta -Oz -o phasing/hyb_JANXIP010000005.1.vcf.gz
bcftools index phasing/hyb_JANXIP010000005.1.vcf.gz

