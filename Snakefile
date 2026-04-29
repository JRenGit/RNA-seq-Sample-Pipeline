import pandas as pd
import os
import sys

configfile: "config.yaml"

# Validate sample sheet
REQUIRED_COLS = ["sample", "fq1", "fq2", "condition"]

try:
  samples_df = pd.read_csv(config["samples"], sep="\t")
except FileNotFoundError:
  sys.exit(f"Error: sample sheet '{config['samples']}' not found.")
except pd.errors.EmptyDataError:
  sys.exit(f"Error: sample sheet '{config['samples']}'is empty.")
except pd.errors.ParserError as e:
  sys.exit(f"Error: sample sheet '{config['samples']}' has malformed rows. "
           f"Check for inconsistent column counts.\n{e}")
except PermissionError:
  sys.exit(f"Error: permission denied reading '{config['samples']}'.")
except UnicodeDecodeError:
  sys.exit(f"Error: sample sheet '{config['samples']}' contains non-UTF-8 characters.")

missing_cols = [c for c in REQUIRED_COLS if c not in samples_df.columns]
if missing_cols:
  sys.exit(f"Error: sample sheet missing columns: {missing_cols}")

if samples_df["sample"].duplicated().any():
  dupes = samples_df["sample"][samples_df["sample"].duplicated()].tolist()
  sys.exit(f"Error: duplicate samples names in sample sheet: {dupes}")

# Validate FASTQ paths
for _, row in samples_df.iterrows():
  for fq in ["fq1", "fq2"]:
    if not os.path.exists(row[fq]):
      sys.exit(f"Error: FASTQ not found for sample '{row['sample']}': {row[fq]}")

# Validate reference files
for ref_key in ["fasta", "gtf"]:
    ref_path = config["genome"][ref_key]
    if not os.path.exists(ref_path):
      sys.exit(f"Error: reference file not found: {ref_path}")

SAMPLES = samples_df["sample"].tolist()

# Add helper to look up FASTQ paths
def get_fq(wildcards, read):
  return samples_df.loc[samples_df["sample"] == wildcards.sample, read].iloc[0]

# Target rule
rule all:
  input:
    expand("{results}/fastqc/{sample}_{read}_fastqc.html",
      results=config["results"], samples=SAMPLES, read=["R1", "R2"],

# Rule: FastQC with per-sample read quality assessment
# check quality scores before computationally expensive STAR alignment
rule fastqc:
  input:
    fq=lambda wc: get_fq(wc, wc.read.lower())
  output:
    html="{results}/fastqc/{sample}_{read}_fastqc.html",
    zip="{results}/fastqc/{sample}_{read}_fastqc.zip"
  params:
    outdir="{results}/fastqc"
  log:
    "{logs}/fastqc/{sample}_{read}.log"
  conda:
    "envs/tools.yaml"
  threads: 2
  shell:
    """
    fastqc {input.fq} \
      --outdir {params.dir} \
      --threads {threads} \
      > {log} 2>&1

    if [ ! -s {output.html} ] || [ ! -s {output.zip} ]; then
      echo "Error: FastQC did not produce expected outputs for {wildcards.sample} {wildcards.read}" >> {log}
      exit 1
    fi
    """
# Rule: STAR Alignment
rule star_align:
  input:
    fq1=lambda wc: get_fq(wc, "fq1"),
    fq2=lambda wc: get_fq(wc, "fq2"),
    index=config["genome"]["star_index"],
    qc_r1="{results}/fastqc/{sample}_R1_fastqc.html",
    qc_r2="{results}/fastqc/{sample}_R2_fastqc.html"
  output:
    bam="{results}/star/{sample}/{sample}.Aligned.sortedbyCoord.out.bam",
    log_final="{results}/star/{sample}/{sample}.Log.final.out",
    sj="{results}/star/{sample}/{sample}.SJ.out.tab"
  params:
    out_prefix="{results}/star/{sample}/{sample}.",
    overhang=config["star"]["overhang"]
  log:
    "{logs}/star/{sample}.log"
  conda:
    "envs/tools.yaml"
  threads: config["star"]["threads"]
  shell:
    """
    STAR \
      --runMode alignReads \
      --genomeDir {input.index} \
      --readFilesIn {input.fq1} {input.fq2} \
      --readFilesCommand zcat \
      --outSAMtype BAM SortedByCoordinate \
      --twopassMode Basic \
      --sjdbOverhang {params.overhang} \
      --outFileNamePrefix {params.out_prefix} \
      --runThreadN {threads} \
      > {log} 2>&1

    if [ ! -s {output.bam} ]; then
      echo "Error: STAR produced empty BAM for {wildcards.sample}" >> {log}
      exit 1
    fi

    samtools index {output.bam} >> {log} 2>&1
    """
# Rule: Feature Counts with final output as a count matrix for DESeq2 use
rule featurecounts:
  input:
    bam=expand("{results}/star/{sample}.Aligned.sortedByCoord.out.bam",
               results=config["results"], sample=SAMPLES),
    gtf=config["genome"]["gtf"]
  output:
    counts="{results}/featurecounts/counts.txt",
    summary="{results}/featurecounts/counts.txt.summary"
  params:
    strand=config["featurecounts"]["strand"]
  logs:
    "{logs}/featurecounts/featurecounts.log"
  conda:
    "evns/tools.yaml"
  threads: config["featurecounts"]["threads"]
  shell:
    """
    featureCounts \
      -a {input.gtf} \
      -o {output.counts} \
      -p --countReadPairs
      -s {params.strand} \
      -T {threads} \
      -t exon \
      -g gene_id \
      {input.bams} \
      > {log} 2>&1
    if [ ! -s {output.counts} ]; then
      echo "Error: featureCounts produced empty count matrix" >> {log}
      exit 1
    fi
    """"
