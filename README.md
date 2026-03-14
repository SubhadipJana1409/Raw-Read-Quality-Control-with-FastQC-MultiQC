# Day 01 — Raw Read Quality Control with FastQC & MultiQC

**#30DaysOfBioinformatics | SubhadipJana1409**

![FastQC](https://img.shields.io/badge/FastQC-v0.12.1-blue)
![MultiQC](https://img.shields.io/badge/MultiQC-v1.33-blue)
![fastp](https://img.shields.io/badge/fastp-v0.23.4-blue)
![Status](https://img.shields.io/badge/Status-Complete-brightgreen)

> *Establishing NGS baseline quality metrics before any downstream analysis*

---

## Research Context

Quality control is the most critical — and most overlooked — step in any NGS pipeline.
**FastQC** (Andrews, 2010) is the gold standard for per-base quality scores, GC content,
adapter contamination, and duplication rate inspection. **MultiQC** (Ewels et al., 2016)
aggregates reports across all samples into a single interactive dashboard, surfacing
batch-level anomalies invisible in per-sample review.

This day establishes the **preprocessed FASTQ files** that feed directly into:
- **Day 02** — Alignment benchmarking (BWA-MEM2 / STAR / Minimap2)
- **Day 06** — RNA-seq quantification (STAR + Salmon)
- **Day 07** — Differential expression (DESeq2 / edgeR / limma-voom)

---

## Dataset

| Field | Value |
|---|---|
| GEO Accession | [GSE96960](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96960) |
| SRA Project | SRP100973 |
| Organism | *Homo sapiens* |
| Cell Type | PBMC (Peripheral Blood Mononuclear Cells) |
| Library | Paired-end, Illumina HiSeq 2500, 75 bp |
| Condition | Healthy Control (HC) vs Septic Shock (SS) |
| Demo Sample | SRR5223500 — Healthy Control 1 |
| Reference | Liao et al. (2017) *JCI Insight* 2(2) |

---

## Tools & Versions

| Tool | Version | Purpose |
|---|---|---|
| FastQC | 0.12.1 | Per-sample QC report |
| MultiQC | 1.33 | Aggregate QC dashboard |
| fastp | 0.23.4 | Adapter trimming + quality filtering |
| Python | 3.11 | QC parsing, figures, summary tables |

---

## Project Structure

```
Day_01_Raw_Read_QC/
├── data/
│   ├── raw_fastq/                     # Input: SRR5223500_1/2.fastq.gz
│   └── trimmed_fastq/                 # Output: *trimmed.fastq.gz
├── qc/
│   ├── raw/                           # FastQC reports (raw)
│   ├── trimmed/                       # FastQC reports (trimmed)
│   ├── fastp_reports/                 # fastp JSON + HTML
│   ├── multiqc_raw/                   # Pre-trim MultiQC HTML
│   └── multiqc_trimmed/               # Post-trim MultiQC HTML
├── results/
│   ├── qc_summary_table.csv           # Full per-sample QC metrics
│   └── fastqc_module_flags.csv        # PASS/WARN/FAIL per module
├── figures/
│   ├── fig1_fastqc_heatmap.png        # Module status heatmap
│   ├── fig2_before_after_trimming.png # Q30/GC/reads comparison
│   ├── fig3_read_fate.png             # Filtering breakdown
│   └── fig4_library_quality.png      # Duplication + insert size
├── scripts/
│   └── generate_analysis.py          # All figures + tables
├── Day_01_Raw_Read_QC.ipynb           # Full analysis notebook
├── 00_download_data.sh                # SRA download script
├── environment.yml                    # Conda environment
└── README.md
```

---

## Setup

```bash
conda env create -f environment.yml
conda activate day01_qc
```

---

## Quick Start

```bash
# 1. Download data (500k reads subset)
bash 00_download_data.sh

# 2. Run full notebook
jupyter lab Day_01_Raw_Read_QC.ipynb

# 3. Or generate figures/tables only
python scripts/generate_analysis.py
```

---

## Results

### QC Summary — SRR5223500 (Healthy Control 1)

| Metric | Before Trimming | After Trimming | Change |
|---|---|---|---|
| Total Reads | 1,000,000 | 819,902 | −18.0% |
| Q20 Rate | 93.62% | 99.82% | **+6.2%** |
| Q30 Rate | 90.19% | 98.42% | **+8.2%** |
| GC Content | 47.75% | 47.31% | −0.44% |
| Duplication Rate | 1.55% | — | ✅ Excellent |
| Insert Size Peak | 75 bp | — | Expected |

### Read Filtering Breakdown

| Category | Count | % |
|---|---|---|
| ✅ Passed filter | 819,902 | 82.0% |
| ⚠️ Too short (<50bp) | 180,096 | 18.0% |
| ❌ Low quality | 2 | <0.01% |
| Adapter trimmed | 7,916 | 0.8% |

### FastQC Module Flags

| Module | Raw R1 | Raw R2 | Trimmed R1 | Trimmed R2 |
|---|---|---|---|---|
| Basic Statistics | ✅ | ✅ | ✅ | ✅ |
| Per base sequence quality | ✅ | ✅ | ✅ | ✅ |
| Per sequence quality scores | ✅ | ✅ | ✅ | ✅ |
| **Per base sequence content** | ❌ FAIL | ⚠️ WARN | ❌ FAIL | ⚠️ WARN |
| Per sequence GC content | ✅ | ✅ | ✅ | ✅ |
| Adapter Content | ✅ | ✅ | ✅ | ✅ |
| Sequence Duplication Levels | ✅ | ✅ | ✅ | ✅ |
| Sequence Length Distribution | ✅ | ✅ | ⚠️ WARN | ⚠️ WARN |

---

## Key Findings & Biological Interpretation

1. **High quality library** — Q30 rate of 90.19% raw → 98.42% post-trim, well above the ENCODE RNA-seq minimum of 80%.

2. **Per base sequence content FAIL is expected** — The first ~10 bases show nucleotide imbalance, a known artifact of random hexamer priming in RNA-seq library prep. This is not a data quality issue.

3. **Very low duplication rate (1.55%)** — Indicates minimal PCR over-amplification, which preserves accurate expression quantification for downstream DE analysis.

4. **18% reads removed as too short** — Sliding-window trimming cut reads to <50bp. Discarding these is appropriate — short reads map ambiguously to multi-gene loci in RNA-seq.

5. **Minimal adapter contamination (<1%)** — Only 7,916 reads required adapter trimming, indicating high quality library preparation.

6. **GC content 47.75%** — Within ±3% of expected 50% for human PBMC RNA-seq; no GC bias.

> **Overall: SRR5223500 is a high-quality paired-end RNA-seq library, ready for alignment in Day 02.**

---

## Trimming Parameter Justification

| Parameter | Value | Rationale |
|---|---|---|
| `--cut_right` | sliding window 3'→5' | Illumina quality degrades toward 3' end |
| `--cut_window_size 4` | 4 bp | Recommended default; balances sensitivity vs specificity |
| `--cut_mean_quality 20` | Q20 | 99% base call accuracy |
| `--length_required 50` | 50 bp | Sub-50bp reads map unreliably in RNA-seq |
| `--detect_adapter_for_pe` | auto | Overlap-based detection; no manual adapter sequences |
| `--correction` | enabled | Reduces overlap-region base errors |

**Why fastp over Trimmomatic?** fastp is 3–8× faster (SIMD-accelerated), produces native JSON for MultiQC, and auto-detects adapters without manual input (Chen et al., 2018).

---

## Key Papers

1. **Andrews S (2010)** FastQC: A Quality Control Tool for High Throughput Sequence Data. Babraham Bioinformatics.
2. **Ewels P et al. (2016)** MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics* 32(19):3047–3048. doi:10.1093/bioinformatics/btw354
3. **Chen S et al. (2018)** fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics* 34(17):i884–i890. doi:10.1093/bioinformatics/bty560
4. **Liao X et al. (2017)** Concordance of gene expression in human blood and PBMCs. *JCI Insight* 2(2). doi:10.1172/jci.insight.89728
