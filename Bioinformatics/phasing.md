mkdir phasing

mamba create --name gatk4 gatk4=4.5.0.0
mamba activate gatk4

gatk HaplotypeCaller -R raw_sequences/c_stelleri_au.fasta -I prep_hybrid/sort.hyb1.bam -O phasing/phase_hyb_JANXIP010000005.1.g.vcf.gz -ERC GVCF --intervals JANXIP010000005.1
