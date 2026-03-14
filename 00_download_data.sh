#!/usr/bin/env bash
# =============================================================================
# Day 01 – Data Download Script
# Dataset: Human PBMC RNA-seq (GSE96960 / SRP100973)
#          Illumina paired-end 150 bp | ~50M reads/sample | 6 samples
#          3 Healthy Controls (HC) + 3 Septic Shock (SS) — good for Day 7 DE too
# =============================================================================

set -euo pipefail

# ── SRA accessions ────────────────────────────────────────────────────────────
# Healthy controls
HC=(SRR5223500 SRR5223501 SRR5223502)
# Septic shock patients
SS=(SRR5223503 SRR5223504 SRR5223505)

ALL_SRRS=("${HC[@]}" "${SS[@]}")

# ── Output directories ─────────────────────────────────────────────────────────
RAW_DIR="data/raw_fastq"
mkdir -p "$RAW_DIR"

echo "================================================================="
echo " Downloading 6 PBMC RNA-seq samples from NCBI SRA (GSE96960)"
echo " Tool required: sra-tools >= 3.0   (conda install -c bioconda sra-tools)"
echo "================================================================="

for SRR in "${ALL_SRRS[@]}"; do
    echo ""
    echo "▶ Downloading $SRR ..."
    fasterq-dump "$SRR" \
        --outdir "$RAW_DIR" \
        --split-files \
        --threads 8 \
        --progress
    # Compress to save disk space
    echo "  Compressing ${SRR} FASTQ files..."
    gzip -f "${RAW_DIR}/${SRR}_1.fastq"
    gzip -f "${RAW_DIR}/${SRR}_2.fastq"
    echo "  ✓ ${SRR} done"
done

echo ""
echo "================================================================="
echo " Download complete. Files in: $RAW_DIR"
ls -lh "$RAW_DIR"
echo "================================================================="

# =============================================================================
# ALTERNATIVE: Download only 500k reads per sample (fast test run ~5 min)
# Uncomment the block below and comment the block above if you want a quick test
# =============================================================================
# for SRR in "${ALL_SRRS[@]}"; do
#     echo "▶ Downloading 500k reads from $SRR ..."
#     fasterq-dump "$SRR" \
#         --outdir "$RAW_DIR" \
#         --split-files \
#         --threads 8 \
#         --maxSpotId 500000   # ← limits to first 500k spots (~1M reads PE)
#     gzip -f "${RAW_DIR}/${SRR}_1.fastq"
#     gzip -f "${RAW_DIR}/${SRR}_2.fastq"
#     echo "  ✓ $SRR subset done"
# done
