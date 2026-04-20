# mixed-viral-lineages-genome-assembly

This Snakemake pipeline analyses sequencing data to assess whether viral samples contain evidence of multiple viral lineages and assembles genomes, writes a report file about variable positions and calculates entropy of each genome segment. It is designed for viruses with segmented genomes.

## Input data

The pipeline expects these inputs: 

- Genome reference FASTA file  located in the `references/` directory

- A comma separated file with SRA accessions including information on whether each sample is paired-end or single-end

- A comma separated file with sample names and direct .fastq file downloading links

## Processing and filtering criteria

Pipeline parameters are defined in `config.yaml` file. The following filtering thresholds are applied by default, but can be modified:

- Reads shorter than 50 nucleotides are discarded

- Reads with a Phred quality score below 20 are excluded

- A minimum sequencing depth of 10 reads is required for allele frequency calculations

- Alternative alleles are considered only if they represent at least 90% of reads at a given position

These thresholds are intended to reduce noise from low-quality reads and low-confidence variants. All thresholds can be adjusted.

## Output files

The pipeline produces:

- Tab-separated reports summarizing allele frequencies that differ from reference sequence 

- Assembled consensus FASTA sequences in which nucleotides differing from the reference genome are incorporated

- Tab-separated report summarizing mean entropy and mean read depth for each genome segment 