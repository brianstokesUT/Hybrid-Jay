mkdir phasing

bcftools mpileup -Ou -q 20 -r JANXIP010000005.1 -f raw_sequences/c_stelleri_au.fasta prep_hybrid/sort.hyb.bam | bcftools call -Ou -c | bcftools norm -f raw_sequences/c_stelleri_au.fasta -Oz -o phasing/hyb_JANXIP010000005.1.vcf.gz
bcftools index phasing/hyb_JANXIP010000005.1.vcf.gz


apptainer build shapeit5.sif docker://abelean/shapeit5:5.1.1


wget https://github.com/odelaneau/shapeit5/releases/download/v5.1.1/phase_common_static
chmod +x phase_common_static

phase_common --seed 01000010 --I phasing/hyb_JANXIP010000005.1.vcf.gz --output phasing/target.phased.bcf
