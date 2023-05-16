# brew install samtools bcftools bedtools bedops sratoolkit

cutadapt -q 15 -e 0.12 -a TGGAATTCTCGGGTGCCAAGG -m 16 --discard-untrimmed SRR18488115.fastq.gz --output trimmed_SRR18488115.fastq.gz

# download GCF_000001405.40 genome version (golang implementation)

# download mirna annotations for genome (golang implementation)

# download mirna mature seqs (golang implementation)

# download mirna hairpins seqs (golang implementation)


# reheader genome file by https://genome.ucsc.edu/cgi-bin/hgTracks?chromInfoPage= (golang implementation)

gff2bed < hsa.gff3 > hsa.bed


bowtie-build genomes/GCF_000001405.40.reheadered.fasta bowtie_indexes/GCF_000001405.40

bowtie -S -v 1 -a --best --strata -x bowtie_indexes/GCF_000001405.40 results/SRR18488115/trimmed_SRR18488115.fastq.gz results/SRR18488115/result.sam

samtools view -bo results/SRR18488115/result.bam results/SRR18488115/result.sam

samtools sort -o results/SRR18488115/result.sorted.bam results/SRR18488115/result.bam

samtools index results/SRR18488115/result.sorted.bam

bcftools mpileup --seed 42 -Ou --max-depth 100000 -f genomes/GCF_000001405.40.reheadered.fasta -R mirna/hsa.bed results/SRR18488115/result.sorted.bam \
    | bcftools call -Ou -mv \
    | bcftools filter -g3 -G10 -e 'QUAL<15 || (DP4[0]+DP4[1] == 0) || DP < 10' > results/SRR18488115/result.sorted.filtered.vcf

bedtools intersect -wa -wb -a mirna/hsa.gff3 -b results/SRR18488115/result.sorted.filtered.vcf > result.sorted.filtered.annotated.tsv

# Generate table with ref and alt nucleotides for particular position (golang implementation)

python3 scripts/plot.py results/result-SRR18488115.tsv results/result-SRR18488115.png
