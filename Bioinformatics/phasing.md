mkdir phasing

mamba create --name gatk4 gatk4=4.5.0.0
mamba activate gatk4

gatk HaplotypeCaller -R raw_sequences/c_stelleri_au.fasta -I prep_hybrid/sort.hyb1.bam -O phasing/phase_hyb_JANXIP010000005.1.g.vcf.gz -ERC GVCF --intervals JANXIP010000005.1

mamba deactivate
mamba activate data_prep

bcftools consensus -H 1 -f raw_sequences/c_stelleri_au.fasta phasing/phase_hyb_JANXIP010000005.1.g.vcf.gz > hyb_haplotype1.fasta
bcftools consensus -H 2 -f raw_sequences/c_stelleri_au.fasta phasing/phase_hyb_JANXIP010000005.1.g.vcf.gz > hyb_haplotype2.fasta


samtools faidx raw_sequences/c_stelleri_au.fasta JANXIP010000005.1 | bcftools consensus -H 1 phasing/phase_hyb_JANXIP010000005.1.g.vcf.gz -o hyb_JANXIP010000005.1_haplotype1.fasta
samtools faidx raw_sequences/c_stelleri_au.fasta JANXIP010000005.1 | bcftools consensus -H 2 phasing/phase_hyb_JANXIP010000005.1.g.vcf.gz -o hyb_JANXIP010000005.1_haplotype2.fasta
