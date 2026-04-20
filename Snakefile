import os
import pandas as pd

configfile: "config.yaml" ## somehow needed this part

REF = "references/wmv6-cms001_028_alco-ref.fasta"
VCF_DIR = "vcf"
sra_samples = pd.read_csv(
	config["sra_list"],
	sep=",",
	header=None,
	names=["SAMPLES", "type"],
	dtype=str,
	keep_default_na=False)

curl_samples = pd.read_csv(
	config["curl_list"],
	sep=",",
	header=None,
	names=["SAMPLES", "type", "fastq1", "fastq2"],
	dtype=str,
	keep_default_na=False)

sra_paired = sra_samples[sra_samples["type"] == "paired"]["SAMPLES"].tolist()
sra_single = sra_samples[sra_samples["type"] == "single"]["SAMPLES"].tolist()
curl_paired = curl_samples[curl_samples["type"] == "paired"]["SAMPLES"].tolist()

paired_samples = sra_paired + curl_paired
single_samples = sra_single
SAMPLES = paired_samples + single_samples

rule all:
	input:
		expand("assembled_seq/{sample}_report.tsv", sample=SAMPLES),
		expand("assembled_seq/{sample}.fasta", sample=SAMPLES),
		expand("assembled_seq/{sample}_entropy_report.tsv", sample=SAMPLES)

rule download_fastq_pe:
	output:
		fq1="samples/paired/curl/{sample}/{sample}_1.fastq",
		fq2="samples/paired/curl/{sample}/{sample}_2.fastq"
	params:
		url1=lambda wildcards: curl_samples.loc[curl_samples["SAMPLES"] == wildcards.sample, "fastq1"].values[0],
		url2=lambda wildcards: curl_samples.loc[curl_samples["SAMPLES"] == wildcards.sample, "fastq2"].values[0]
	conda:
		"envs/download_unzip.yaml"
	shell:
		"""
		curl -L --retry 5 --retry-delay 10 --max-time 600 {params.url1} -o {output.fq1}.gz
		zcat {output.fq1}.gz > {output.fq1}
		rm {output.fq1}.gz

		curl -L --retry 5 --retry-delay 10 --max-time 600 {params.url2} -o {output.fq2}.gz
		zcat {output.fq2}.gz > {output.fq2}
		rm {output.fq2}.gz
		"""

# rule download_sra_pe:
# 	output:
# 		"samples/paired/{sample}/{sample}.sra"
# 	conda:
# 		"envs/sra_tools.yaml"
# 	shell:
# 		"prefetch {wildcards.sample} -O samples/paired"

# rule download_sra_se:
# 	output:
# 		"samples/single/{sample}/{sample}.sra"
# 	conda:
# 		"envs/sra_tools.yaml"
# 	shell:
# 		"prefetch {wildcards.sample} -O samples/single"

rule sra_fastq_pe:
	output:
		fq1="samples/paired/sra/{sample}/{sample}_1.fastq",
		fq2="samples/paired/sra/{sample}/{sample}_2.fastq"
	conda:
		"envs/sra_tools.yaml"
	shell:
		"""
		fasterq-dump {wildcards.sample} \
			--split-files \
			--outdir samples/paired/sra/{wildcards.sample} \
			--temp samples/paired/sra/{wildcards.sample}
		"""

rule sra_fastq_se:
	output:
		"samples/single/sra/{sample}/{sample}.fastq"
	conda:
		"envs/sra_tools.yaml"
	shell:
		"""
		fasterq-dump {wildcards.sample} \
			--outdir samples/single/sra/{wildcards.sample} \
			--temp samples/single/sra/{wildcards.sample}
		"""

def get_paired_fastq(wildcards):
	if wildcards.sample in curl_paired:
		return {
			"fq1": f"samples/paired/curl/{wildcards.sample}/{wildcards.sample}_1.fastq",
			"fq2": f"samples/paired/curl/{wildcards.sample}/{wildcards.sample}_2.fastq"
		}
	else:
		return {
			"fq1": f"samples/paired/sra/{wildcards.sample}/{wildcards.sample}_1.fastq",
			"fq2": f"samples/paired/sra/{wildcards.sample}/{wildcards.sample}_2.fastq"
		}

# rule fastqdump_pe:
# 	input:
# 		"samples/paired/sra/{sample}/{sample}.sra"
# 	output:
# 		fq1="samples/paired/sra/{sample}/{sample}_1.fastq",
# 		fq2="samples/paired/sra/{sample}/{sample}_2.fastq"
# 	conda:
# 		"envs/sra_tools.yaml"
# 	shell:
# 		"fastq-dump --split-files {input} --outdir samples/paired/{wildcards.sample}"

# rule fastqdump_se:
# 	input:
# 		"samples/single/{sample}/{sample}.sra"
# 	output:
# 		"samples/single/{sample}/{sample}.fastq"
# 	conda:
# 		"envs/sra_tools.yaml"
# 	shell:
# 		"fastq-dump {input} --outdir samples/single/{wildcards.sample}"

rule trim_reads_pe:
	input:
		# r1="samples/paired/{sample}/{sample}_1.fastq",
		# r2="samples/paired/{sample}/{sample}_2.fastq",
		unpack(get_paired_fastq)
	output:
		trimmed_r1="trimmed/paired/{sample}/{sample}_1_trimmed.fastq",
		trimmed_r2="trimmed/paired/{sample}/{sample}_2_trimmed.fastq",
		report_html="trimmed/paired/{sample}/fastq_report.html",
		report_json="trimmed/paired/{sample}/fastq_report.json"
	params:
		min_length=config["min_read_length"],
		quality_threshold=config["quality_threshold"]
	threads: 4
	conda:
		"envs/fastp.yaml"
	shell:
		"""
		mkdir -p trimmed/paired

		fastp \
			--in1 {input.fq1} \
			--in2 {input.fq2} \
			--out1 {output.trimmed_r1} \
			--out2 {output.trimmed_r2} \
			--thread {threads} \
			--length_required {params.min_length} \
			--qualified_quality_phred {params.quality_threshold} \
			--html {output.report_html} \
			--json {output.report_json} \
			--verbose
		"""

rule trim_reads_se:
	input:
		single="samples/single/sra/{sample}/{sample}.fastq"
	output:
		trimmed_single="trimmed/single/{sample}/{sample}_trimmed.fastq",
		report_html="trimmed/single/{sample}/fastq_report.html",
		report_json="trimmed/single/{sample}/fastq_report.json"
	params:
		min_length=config["min_read_length"],
		quality_threshold=config["quality_threshold"]
	threads: 4
	conda:
		"envs/fastp.yaml"
	shell:
		"""
		mkdir -p trimmed/single

		fastp \
			--in1 {input.single} \
			--out1 {output.trimmed_single} \
			--thread {threads} \
			--length_required {params.min_length} \
			--qualified_quality_phred {params.quality_threshold} \
			--html {output.report_html} \
			--json {output.report_json} \
			--verbose
		"""

rule map_to_reference_pe:
	input:
		r1="trimmed/paired/{sample}/{sample}_1_trimmed.fastq",
		r2="trimmed/paired/{sample}/{sample}_2_trimmed.fastq"
	output:
		bam="mapped/{sample}/{sample}_pe_mapped.bam",
		bai="mapped/{sample}/{sample}_pe_mapped.bam.bai"
	threads: 4
	conda:
		"envs/bwa_samtools.yaml"
	shell:
		"""
		mkdir -p mapped/{wildcards.sample}

		bwa-mem2 index {REF}
		bwa-mem2 mem -t {threads} {REF} {input.r1} {input.r2} | \
		samtools view -b -F 4 | \
		samtools sort -@ {threads} -o {output.bam} -
		samtools index {output.bam}
		"""

rule map_to_reference_se:
	input:
		single="trimmed/single/{sample}/{sample}_trimmed.fastq"
	output:
		bam="mapped/{sample}/{sample}_se_mapped.bam",
		bai="mapped/{sample}/{sample}_se_mapped.bam.bai"
	threads:4
	conda:
		"envs/bwa_samtools.yaml"
	shell:
		"""
		mkdir -p mapped/{wildcards.sample}

		bwa-mem2 index {REF}
		bwa-mem2 mem -t {threads} {REF} {input.single} | \
		samtools view -b -F 4 | \
		samtools sort -@ {threads} -o {output.bam} -
		samtools index {output.bam}
		"""

def get_bam(wc):
    if wc.sample in paired_samples:
        return f"mapped/{wc.sample}/{wc.sample}_pe_mapped.bam"
    elif wc.sample in single_samples:
        return f"mapped/{wc.sample}/{wc.sample}_se_mapped.bam"
    else:
        raise ValueError(f"Sample {wc.sample} not found in paired or single lists")

rule vcf_calling:
	input:
		bam=get_bam
	output:
		vcf="vcf/{sample}/{sample}.vcf.gz",
		depth="vcf/{sample}/{sample}.depth.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"""
		mkdir -p vcf/{wildcards.sample}

		bcftools mpileup -f {REF} {input.bam} | \
		bcftools call -mv -Oz -o {output.vcf}

		bcftools index {output.vcf}
		samtools depth -a -H {input.bam} > {output.depth}
		"""

rule analyze_vcfs:
	input:
		vcf="vcf/{sample}/{sample}.vcf.gz",
		depth="vcf/{sample}/{sample}.depth.tsv",
		ref = REF
	output:
		report="assembled_seq/{sample}_report.tsv",
		fasta="assembled_seq/{sample}.fasta",
		entropy_report="assembled_seq/{sample}_entropy_report.tsv"
	params:
		min_depth=config["min_depth"],
		min_alt_freq=config["min_alt_freq"],
		max_freq_for_N = config["max_freq_for_N"]
	conda:
		"envs/vcf_analysis.yaml"
	script:
		"scripts/Combined.py"
