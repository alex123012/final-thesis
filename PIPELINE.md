# Installation step
```bash
brew install samtools bcftools bedtools bedops sratoolkit
```

# Preparation step
* Download sequences
```bash
# golang implementation
```

* Run adapter trimming
```bash
cutadapt -q 15 -e 0.12 -a TGGAATTCTCGGGTGCCAAGG -m 16 --discard-untrimmed <Accession>.fastq.gz --output trimmed_<Accession>.fastq.gz
```

* Download genome by version (default `GCF_000001405.40`)
```bash
# golang implementation
```

* download mirna annotations for genome
```bash
# golang implementation
```

* download mirna mature seqs
```bash
# golang implementation
```

* download mirna hairpins seqs
```bash
# golang implementation
```

* reheader genome file by https://genome.ucsc.edu/cgi-bin/hgTracks?chromInfoPage=
```bash
# golang implementation
```

* Generate BED annotations from GFF3
```bash
gff2bed < hsa.gff3 > hsa.bed
```

# Run step
* Build bowtie index for genome
```bash
bowtie-build genomes/<Genome version>.reheadered.fasta bowtie_indexes/<Genome version>
```

* Map reads with bowtie
```bash
bowtie -S -v 1 -a --best --strata -x bowtie_indexes/<Genome version> results/<Accession>/trimmed_<Accession>.fastq.gz results/<Accession>/result.sam
```

* Convert SAM bowtie output to BAM
```bash
samtools view -bo results/<Accession>/result.bam results/<Accession>/result.sam
```

* Sort Bam file
```bash
samtools sort -o results/<Accession>/result.sorted.bam results/<Accession>/result.bam
```

* Index Bam file
```bash
samtools index results/<Accession>/result.sorted.bam
```

* Run variant calling and filtering
```bash
bcftools mpileup --seed 42 -Ou --max-depth 100000 -f genomes/<Genome version>.reheadered.fasta -R mirna/hsa.bed results/<Accession>/result.sorted.bam \
    | bcftools call -Ou -mv \
    | bcftools filter -g3 -G10 -e 'QUAL<15 || (DP4[0]+DP4[1] == 0) || DP < 10' > results/<Accession>/result.sorted.filtered.vcf
```

* Intersect annotations with generated variants
```bash
bedtools intersect -wa -wb -a mirna/hsa.gff3 -b results/<Accession>/result.sorted.filtered.vcf > results/<Accession>/result.sorted.filtered.annotated.tsv
```

* Generate table with ref and alt nucleotides for particular position
```bash
# golang implementation
```

# Plot figures
* Plot heatmap for editing events
```bash
python3 scripts/plot.py results/result-<Accession>.tsv results/result-heatmap-<Accession>.png
```
