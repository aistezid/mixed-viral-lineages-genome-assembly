#!/usr/bin/env python
# coding: utf-8

import os
import gzip
import csv
import math
from Bio import SeqIO

sample_name = snakemake.wildcards.sample

vcf_file = snakemake.input.vcf
reference_file = snakemake.input.ref
depth_file = snakemake.input.depth
fasta_output = snakemake.output.fasta
report_file = snakemake.output.report
entropy_report_file = snakemake.output.entropy_report

min_depth = snakemake.params.min_depth
min_alt_freq = snakemake.params.min_alt_freq
max_freq_for_N = snakemake.params.max_freq_for_N      ## if max nucleotide freq < this -> N

reference_seqs = {}
ref_index = {}
for record in SeqIO.parse(reference_file, "fasta"): ## read reference file
    segment = record.id.split("|")[1]
    reference_seqs[segment] = str(record.seq)
    ref_index[segment] = {}
    for i, nt in enumerate(record.seq, start=1):
        ref_index[segment][i] = nt


depth_data = {}
with open(depth_file) as f:
    for line in f:
        if line.startswith("#"): ## skip header line
            continue
        parts = line.rstrip().split("\t")
        chrom = parts[0]
        segment = chrom.split("|")[1] ## split depth file data to segment position and depth
        pos = int(parts[1])
        dp = int(parts[2])
        if segment not in depth_data:
            depth_data[segment] = {}
        depth_data[segment][pos] = dp

vcf_data = {}

with gzip.open(vcf_file, "rt") as f:
    for line in f:
        if line.startswith("#"): ## skip header line
            continue
        fields = line.rstrip().split("\t") ## split vcf line to columns, extract info
        chrom = fields[0]
        segment = chrom.split("|")[1]
        pos = int(fields[1])
        ref_nt = fields[3]
        alt_nt = fields[4].split(",")
        format_keys = fields[8].split(":")
        sample_values = fields[9].split(":")
        fmt = dict(zip(format_keys, sample_values)) ## parse the FORMAT column and save it to a dictionary
        ad = list(map(int, fmt["AD"].split(",")))

        alleles = [ref_nt] + alt_nt ## get list of all alleles in the site
        total_depth = sum(ad)

        low_freq_alleles = {}

        allele_counts = dict(zip(alleles, ad)) ## calcylate frequencies
        allele_freqs = {}
        for allele, count in allele_counts.items():
            if count > 0:
                allele_freqs[allele] = count / total_depth

        if segment not in vcf_data:
            vcf_data[segment] = {}

        entropy = 0
        for freq in allele_freqs.values(): ## entropy calculation
            entropy -= (freq * math.log2(freq)) / math.log2(4)

        vcf_data[segment][pos] = { ## store all calculated info of vcfs in one dict
            "ref": ref_nt,
            "depth": total_depth,
            "allele_counts": allele_counts,
            "allele_freqs": allele_freqs,
            "entropy": entropy
        }

## Entropy/depth accumulation — done over every reference position, not just VCF-called ones,
## so segments with no VCF sites at all still get a row in the entropy report
seg_entropy_sum = {seg: 0 for seg in reference_seqs}
seg_depth_sum = {seg: 0 for seg in reference_seqs}
seg_pos_count = {seg: 0 for seg in reference_seqs}

for segment, positions in ref_index.items():
    for position in positions:
        if segment in vcf_data and position in vcf_data[segment]:
            entropy = vcf_data[segment][position]["entropy"]
            depth = vcf_data[segment][position]["depth"]
        else:
            depth = depth_data.get(segment, {}).get(position, 0)
            entropy = 0  # no VCF call = treated as invariant at this position
        seg_entropy_sum[segment] += entropy
        seg_depth_sum[segment] += depth
        seg_pos_count[segment] += 1

## Sequence assembly 
new_seq = {}
ref_diff_positions = []
for segment, positions in ref_index.items():
    seq = []
    for position in sorted(positions):
        ref_nt = ref_index[segment][position]

        if segment not in vcf_data or position not in vcf_data[segment]: ## if nucleotide position is not in vcf
            actual_depth = depth_data.get(segment, {}).get(position, 0) ## get depth from file
            called_nt = "N" if actual_depth < min_depth else ref_nt
            seq.append(called_nt) ## take reference nt if depth is too low
            if actual_depth < min_depth: ## always report low-depth positions, with whatever read counts they have
                freqs = {ref_nt: 1.0} if actual_depth > 0 else {} ## save freq as 1 if there are reads, since this saves only entirely reference positions
                counts = {ref_nt: actual_depth} if actual_depth > 0 else {}
                ref_diff_positions.append((segment, position, ref_nt, called_nt, actual_depth, freqs, counts))  
            continue

        freqs = vcf_data[segment][position]["allele_freqs"]
        counts = vcf_data[segment][position]["allele_counts"]
        depth = vcf_data[segment][position]["depth"]

        if depth < min_depth:  ## get low depth positions - call N but report frequencies
            seq.append("N")
            ref_diff_positions.append((segment, position, ref_nt, "N", depth, freqs, counts))
            continue
        
        best_allele = max(freqs, key=freqs.get)

        if freqs[best_allele] < max_freq_for_N:
            ref_diff_positions.append((segment, position, ref_nt, "N", depth, freqs, counts))
            seq.append("N")
            continue

        alt_alleles = {a: f for a, f in freqs.items() if a != ref_nt and f >= min_alt_freq}  ## get alt alleles above min_alt_freq
        if alt_alleles:
            best_alt = max(alt_alleles, key=alt_alleles.get)  ## take the most frequent alt allele if it passes threshold
            seq.append(best_alt)
            ref_diff_positions.append((segment, position, ref_nt, best_alt, depth, freqs, counts))
        else:
            sub_threshold_alts = {a: f for a, f in freqs.items() if a != ref_nt and f > 0}
            if sub_threshold_alts:  ## variable but no alt dominates
                seq.append("N")
                counts = vcf_data[segment][position]["allele_counts"]
                ref_diff_positions.append((segment, position, ref_nt, "N", depth, freqs, counts))
            else:  ## purely reference
                seq.append(ref_nt)

    new_seq[segment] = "".join(seq)

## save files

with open(fasta_output, "w") as f:
    for segment, seq in new_seq.items():
        header = f"{sample_name}|{segment}"
        f.write(f">{header}\n{seq}\n")

with open(report_file, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["sample", "segment", "position", "depth", "reference", "called_nt",
                     "A_freq", "A_depth", "T_freq", "T_depth", "G_freq", "G_depth", "C_freq", "C_depth"])
    for segment, pos, ref, called, depth, alleles, counts in ref_diff_positions:
        writer.writerow([
            sample_name, segment, pos, depth, ref, called,
            round(alleles.get("A", 0), 4), counts.get("A", 0),
            round(alleles.get("T", 0), 4), counts.get("T", 0),
            round(alleles.get("G", 0), 4), counts.get("G", 0),
            round(alleles.get("C", 0), 4), counts.get("C", 0),
        ])

with open(entropy_report_file, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["sample", "segment", "mean_entropy", "mean_depth"])
    for segment in seg_entropy_sum:
        mean_entropy = seg_entropy_sum[segment] / seg_pos_count[segment]
        mean_depth = seg_depth_sum[segment] / seg_pos_count[segment]
        writer.writerow([sample_name, segment, mean_entropy, mean_depth])

# command line: jupyter nbconvert --to script Combined.ipynb

