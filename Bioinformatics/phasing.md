mkdir phasing

#Filter the .bam file

mamba activate gatk4

gatk HaplotypeCaller -R raw_sequences/c_stelleri_au.fasta -I prep_hybrid/sort.hyb1.bam -O phasing/phase_hyb_JANXIP010000005.1.g.vcf.gz -ERC GVCF --intervals JANXIP010000005.1

mamba deactivate
mamba activate data_prep

bcftools consensus -H 1 -f raw_sequences/c_stelleri_au.fasta phasing/phase_hyb_JANXIP010000005.1.g.vcf.gz > hyb_haplotype1.fasta
bcftools consensus -H 2 -f raw_sequences/c_stelleri_au.fasta phasing/phase_hyb_JANXIP010000005.1.g.vcf.gz > hyb_haplotype2.fasta


samtools faidx raw_sequences/c_stelleri_au.fasta JANXIP010000005.1:37833122-37837993| bcftools consensus -H 1 phasing/phase_hyb_JANXIP010000005.1.g.vcf.gz -o hyb_rag1_haplotype1.fasta
samtools faidx raw_sequences/c_stelleri_au.fasta JANXIP010000005.1 | bcftools consensus -H 2 phasing/phase_hyb_JANXIP010000005.1.g.vcf.gz -o hyb_JANXIP010000005.1_haplotype2.fasta




#candidate_rag1_db
~PATH/tools/ncbi-blast-2.14.0+/bin/blastn -query hyb_JANXIP010000005.1_haplotype1.fasta -db rag1/candidate_rag1_db -out rag1/hap1_rag1_blast.out

~PATH/tools/ncbi-blast-2.14.0+/bin/blastn -query hyb_JANXIP010000005.1_haplotype2.fasta -db rag1/candidate_rag1_db -out rag1/hap2_rag1_blast.out
