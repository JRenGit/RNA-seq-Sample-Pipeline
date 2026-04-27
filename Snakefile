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
