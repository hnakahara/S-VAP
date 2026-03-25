#!/usr/bin/env python3

##############################################################################################
# example command
# python CreateVcf.py --input /path/to/F1_variant.csv --outdir /path/to/output --prefix test
##############################################################################################

import subprocess
from subprocess import PIPE
import os
import pandas as pd
from io import StringIO
import re
import argparse

# ===============================
# paramerter
# ===============================
parser = argparse.ArgumentParser(description="Create VCF from HGVS using transvar")
parser.add_argument("--input", required=True, help="Input CSV file")
parser.add_argument("--outdir", required=True, help="Output directory")
parser.add_argument("--prefix", required=True, help="Output file prefix")

args = parser.parse_args()

input_csv = args.input
outdir = args.outdir
prefix = args.prefix

# create output directory
os.makedirs(outdir, exist_ok=True)

# output file path
transvar_csv_path = os.path.join(outdir, f"{prefix}_transvar.csv")
vcf_path = os.path.join(outdir, f"{prefix}.vcf")
failed_path = os.path.join(outdir, f"{prefix}_failed.csv")

# ===============================
# environment
# ===============================
# os.environ['TRANSVAR_CFG'] = '/mnt/work/ngs/transvar/hg19/transvar_nas.cfg'

# ===============================
# read input file
# ===============================
variants = pd.read_csv(input_csv)

# ===============================
# empty list for result
# ===============================
all_results = []
vcf_rows = []
failed_rows = []


# ===============================
# run
# ===============================
for _, row in variants.iterrows():
    gene = row["gene"]
    cdna = row["cdna"]
    transcript_id = row["transcript"]

    print(f"analysis for {gene}:c.{cdna}")

    cmd = [
        "transvar", "canno",
        "-i", f"{gene}:c.{cdna}",
        "--refseq",
        "--gseq"
    ]

    # subprocess
    try:
        res = subprocess.run(
            cmd,
            stdout=PIPE,
            stderr=PIPE,
            text=True,
            shell=False
        )
    except FileNotFoundError as e:
        failed_rows.append({
            "gene": gene,
            "cdna": cdna,
            "transcript": transcript_id,
            "reason": "transvar_not_found",
            "error": str(e)
        })
        continue

    # ===============================
    # fail
    # ===============================
    if res.stdout.strip() == "":
        failed_rows.append({
            "gene": gene,
            "cdna": cdna,
            "transcript": transcript_id,
            "reason": "no_transvar_output",
            "stderr": res.stderr.strip()
        })
        continue

    # ===============================
    # read as dataframe
    # ===============================
    try:
        df = pd.read_csv(StringIO(res.stdout), sep="\t")
    except Exception as e:
        failed_rows.append({
            "gene": gene,
            "cdna": cdna,
            "transcript": transcript_id,
            "reason": "parse_error",
            "error": str(e)
        })
        continue

    # ===============================
    # transcript filter
    # ===============================
    mane = df[df["transcript"].str.contains(transcript_id, na=False)]
    if mane.empty:
        failed_rows.append({
            "gene": gene,
            "cdna": cdna,
            "transcript": transcript_id,
            "reason": "transcript_not_found"
        })
        continue

    mane = mane.iloc[0]

    # ===============================
    # parse genomic coordinate
    # ===============================
    try:
        coord = mane["coordinates(gDNA/cDNA/protein)"]
        chr_, g_coord, c_coord, p_coord = re.split(":|/", coord)
    except Exception as e:
        failed_rows.append({
            "gene": gene,
            "cdna": cdna,
            "transcript": transcript_id,
            "reason": "coordinate_parse_error",
            "error": str(e)
        })
        continue

    # ===============================
    # normal processing
    # ===============================
    result = {
        "input": mane["input"],
        "gene": mane["gene"],
        "transcript": mane["transcript"],
        "strand": mane["strand"],
        "chrom": chr_,
        "gDNA": g_coord,
        "cDNA": c_coord,
        "protein": p_coord,
        "region": mane["region"],
        "info": mane["info"],
        "vcf_chr": chr_,
        "vcf_pos": mane["POS"],
        "vcf_ref": mane["REF"],
        "vcf_alt": mane["ALT"],
    }

    all_results.append(result)

    vcf_rows.append(
        (result["vcf_chr"], result["vcf_pos"], result["input"],
         result["vcf_ref"], result["vcf_alt"], ".", "PASS", "NAME=" + result["input"])
    )

# ===============================
# output
# ===============================
pd.DataFrame(all_results).to_csv(transvar_csv_path, index=False)

vcf_df = pd.DataFrame(
    vcf_rows,
    columns=["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
)

# ===============================
# ★ sort
# ===============================

# remove chr （chr1 → 1）
vcf_df["CHROM_clean"] = vcf_df["#CHROM"].astype(str).str.replace("chr", "", regex=False)

# order
chrom_order = {str(i): i for i in range(1, 23)}
chrom_order.update({"X": 23, "Y": 24})

# digitization
vcf_df["CHROM_order"] = vcf_df["CHROM_clean"].map(chrom_order)

# converting POS data into numerical values
vcf_df["POS"] = pd.to_numeric(vcf_df["POS"], errors="coerce")

# sort
vcf_df = vcf_df.sort_values(["CHROM_order", "POS"])

# remove unnecessary
vcf_df = vcf_df.drop(columns=["CHROM_clean", "CHROM_order"])

# ===============================
# output
# ===============================
vcf_df.to_csv(vcf_path, sep="\t", index=False)

pd.DataFrame(failed_rows).to_csv(failed_path, index=False)

print("=== DONE ===")
print(f"VCF: {vcf_path}")
print(f"TRANSVAR: {transvar_csv_path}")
print(f"FAILED: {failed_path}")
