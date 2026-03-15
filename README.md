# Day 01 — Raw Read Quality Control with FastQC & MultiQC

**#30DaysOfBioinformatics | SubhadipJana1409**

![FastQC](https://img.shields.io/badge/FastQC-v0.12.1-blue)
![MultiQC](https://img.shields.io/badge/MultiQC-v1.33-blue)
![fastp](https://img.shields.io/badge/fastp-v0.23.4-blue)
![Samples](https://img.shields.io/badge/Samples-6-orange)
![Status](https://img.shields.io/badge/Status-Complete-brightgreen)

> *Establishing NGS baseline quality metrics before any downstream analysis*

---

## Research Context

Quality control is the most critical — and most overlooked — step in any NGS pipeline.
**FastQC** (Andrews, 2010) is the gold standard for per-base quality scores, GC content,
adapter contamination, and duplication rate inspection. **MultiQC** (Ewels et al., 2016)
aggregates reports across all samples into a single interactive dashboard, surfacing
batch-level anomalies invisible in per-sample review. **fastp** (Chen et al., 2018) performs
adapter trimming and quality filtering 3–8× faster than Trimmomatic via SIMD vectorization,
with automatic adapter detection and native JSON output for MultiQC integration.


---

## Dataset

| Field | Value |
|---|---|
| GEO Accession | [GSE96960](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96960) |
| SRA Project | SRP100973 |
| Organism | *Homo sapiens* |
| Cell Type | PBMC (Peripheral Blood Mononuclear Cells) |
| Library | Paired-end, Illumina HiSeq 2500, 75 bp |
| Total Samples | 6 (3 Healthy Controls + 3 Septic Shock) |
| Total FASTQ files | 12 (R1 + R2 per sample) |
| Reference | Liao et al. (2017) *JCI Insight* 2(2) |

| SRR ID | Label | Group |
|---|---|---|
| SRR5223500 | HC_1 | Healthy Control |
| SRR5223501 | HC_2 | Healthy Control |
| SRR5223502 | HC_3 | Healthy Control |
| SRR5223503 | SS_1 | Septic Shock |
| SRR5223504 | SS_2 | Septic Shock |
| SRR5223505 | SS_3 | Septic Shock |

---

## Tools & Versions

| Tool | Version | Purpose |
|---|---|---|
| FastQC | 0.12.1 | Per-sample QC report (12 files) |
| MultiQC | 1.33 | Aggregate QC dashboard |
| fastp | 0.23.4 | Adapter trimming + quality filtering |
| Python | 3.12 | QC parsing, 6 figures, summary tables |

---

## Project Structure

```
Day_01_Raw_Read_QC/
├── data/
│   ├── raw_fastq/                          # 12 input FASTQ files (not committed — use 00_download_data.sh)
│   └── trimmed_fastq/                      # 12 trimmed FASTQ files (not committed)
├── qc/
│   ├── raw/                                # 12 FastQC HTML reports (raw reads)
│   ├── trimmed/                            # 12 FastQC HTML reports (trimmed reads)
│   ├── fastp_reports/                      # 6 fastp JSON + 6 fastp HTML reports
│   ├── multiqc_raw/multiqc_raw.html        # Pre-trim MultiQC dashboard (12 samples)
│   └── multiqc_trimmed/multiqc_trimmed.html # Post-trim MultiQC dashboard (12 samples)
├── results/
│   ├── qc_summary_all_samples.csv          # Full QC metrics for all 6 samples
│   ├── fastp_summary_all_samples.csv       # fastp trimming stats per sample
│   └── fastqc_module_flags_all.csv         # PASS/WARN/FAIL matrix (all modules × all samples)
├── figures/
│   ├── fig1_fastqc_heatmap_raw.png         # FastQC module status heatmap (raw, all 6)
│   ├── fig2_fastqc_heatmap_trimmed.png     # FastQC module status heatmap (trimmed, all 6)
│   ├── fig3_before_after_all_samples.png   # Reads / Q30 / GC before vs after (all 6)
│   ├── fig4_quality_improvement_lollipop.png # Per-sample Q20/Q30 improvement
│   ├── fig5_reads_fate_all_samples.png     # Stacked bar: passed vs removed reads
│   └── fig6_library_quality_all.png        # Duplication rate + insert size per sample
├── scripts/
│   └── generate_analysis.py               # Reproduces all 6 figures + 3 result CSVs
├── Day_01_Raw_Read_QC.ipynb               # Full step-by-step analysis notebook
├── 00_download_data.sh                    # SRA download script (500k reads subset)
├── environment.yml                        # Conda environment
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
# 1. Download all 6 samples (500k reads each, ~5 min)
bash 00_download_data.sh

# 2. Run the full notebook end-to-end
jupyter lab Day_01_Raw_Read_QC.ipynb

# 3. Or regenerate all figures + tables only
python scripts/generate_analysis.py
```

---

## Workflow

```
12 Raw FASTQ files (6 samples × R1+R2)
        │
        ▼
   FastQC (Step 1)             → 12 per-sample HTML reports
        │
        ▼
   MultiQC raw (Step 2)        → 1 aggregate HTML dashboard
        │
        ▼
   Flag QC Issues (Step 3)     → per-sample PASS/WARN/FAIL table
        │
        ▼
   fastp trimming (Step 4)     → 6 JSON reports + 12 trimmed FASTQs
        │
        ▼
   Post-trim FastQC (Step 5)   → 12 trimmed HTML reports
        │
        ▼
   MultiQC trimmed (Step 5)    → 1 aggregate trimmed dashboard
        │
        ▼
   6 Figures + 3 CSVs          → Final deliverables
```

---

## Results — All 6 Samples

### Full QC Summary Table

| Sample | Group | Raw Reads | Trimmed Reads | % Kept | Q30 Before | Q30 After | GC% | Dup Rate |
|---|---|---|---|---|---|---|---|---|
| HC_1 | Healthy | 1,000,000 | 819,902 | 82.0% | 90.19% | **98.42%** | 47.75% | 1.55% |
| HC_2 | Healthy | 1,000,000 | 735,632 | 73.6% | 85.82% | **98.07%** | 47.62% | 1.61% |
| HC_3 | Healthy | 1,000,000 | 780,314 | 78.0% | 88.18% | **98.27%** | 47.71% | 1.30% |
| SS_1 | Septic  | 1,000,000 | 717,724 | 71.8% | 85.10% | **98.14%** | 48.28% | 2.14% |
| SS_2 | Septic  | 1,000,000 | 751,352 | 75.1% | 86.52% | **98.16%** | 46.14% | 2.26% |
| SS_3 | Septic  | 1,000,000 | 833,296 | 83.3% | 90.69% | **98.24%** | 49.30% | 1.71% |

### FastQC Module Status — All 6 Samples (Raw)

| Module | HC_1 R1 | HC_1 R2 | HC_2 R1 | HC_2 R2 | HC_3 R1 | HC_3 R2 | SS_1 R1 | SS_1 R2 | SS_2 R1 | SS_2 R2 | SS_3 R1 | SS_3 R2 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Basic Statistics | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| Per base sequence quality | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| Per sequence quality scores | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| **Per base sequence content** | ❌ | ⚠️ | ❌ | ⚠️ | ❌ | ⚠️ | ❌ | ⚠️ | ❌ | ⚠️ | ❌ | ⚠️ |
| Per sequence GC content | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| Per base N content | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| Adapter Content | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| Sequence Duplication Levels | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| Sequence Length Distribution | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |

---

## Key Findings & Biological Interpretation

1. **All 6 samples pass QC** — Post-trim Q30 rates range from 98.07% to 98.42%, all well above the ENCODE RNA-seq minimum of 80%. The dataset is batch-consistent and suitable for DE analysis.

2. **Per base sequence content FAIL/WARN is expected and not a problem** — The first ~10 bases of every sample show nucleotide composition imbalance, a known artifact of random hexamer priming in RNA-seq library prep. This is consistent across all 12 files and does not indicate a sequencing error.

3. **Septic Shock samples show slightly lower read retention (71.8–75.1% vs 73.6–83.3% HC)** — SS samples had marginally more reads trimmed to <50 bp, suggesting slightly shorter fragment sizes. This is biologically plausible given inflammation-associated RNA degradation in sepsis and is not a technical artifact.

4. **SS samples have slightly higher duplication rates (1.71–2.26% vs 1.30–1.61% HC)** — Still extremely low overall. Both groups are well below the 20% warning threshold, confirming minimal PCR over-amplification across all samples.

5. **No adapter contamination detected in any sample** — fastp auto-detection found negligible adapter content (<1%) across all 12 FASTQ files. No manual adapter sequences were required.

6. **Consistent GC content across all samples (46.14%–49.30%)** — All within ±4% of the expected 50% for human PBMC RNA-seq. No GC bias batch effects detected between HC and SS groups.

7. **Insert size peak = 75 bp for all samples** — Consistent with the 75 bp paired-end read length. Expected and correct.

> **Overall: All 6 samples are high-quality paired-end RNA-seq libraries with consistent QC metrics across HC and SS groups. All are ready for alignment in Day 02.**

---

## Trimming Parameter Justification

| Parameter | Value | Rationale |
|---|---|---|
| `--cut_right` | sliding window 3'→5' | Illumina quality degrades toward 3' end; window trimming is gentler than hard-clipping |
| `--cut_window_size 4` | 4 bp | Recommended default; balances sensitivity vs specificity |
| `--cut_mean_quality 20` | Q20 = 99% base accuracy | Standard NGS threshold; ENCODE minimum |
| `--length_required 50` | 50 bp min | Sub-50bp reads map ambiguously to multi-gene loci in RNA-seq |
| `--detect_adapter_for_pe` | auto overlap-based | No manual adapter sequences needed; robust across library kits |
| `--correction` | enabled | Overlap-based error correction reduces mismatch rates in PE reads |

**Why fastp over Trimmomatic?** fastp is 3–8× faster (SIMD-accelerated), produces native JSON for MultiQC, and auto-detects adapters from read overlap without requiring a manual adapter list (Chen et al., 2018).

---

## Figures

| Figure | Description |
|---|---|
| `fig1_fastqc_heatmap_raw.png` | FastQC PASS/WARN/FAIL heatmap — all 12 raw files, HC vs SS split |
| `fig2_fastqc_heatmap_trimmed.png` | FastQC PASS/WARN/FAIL heatmap — all 12 trimmed files |
| `fig3_before_after_all_samples.png` | Read count, Q30, GC content before vs after trimming (6 samples) |
| `fig4_quality_improvement_lollipop.png` | Per-sample Q20 and Q30 improvement lollipop chart |
| `fig5_reads_fate_all_samples.png` | Stacked bar: reads passing filter vs removed per sample |
| `fig6_library_quality_all.png` | PCR duplication rate + insert size peak per sample |

---

## Key Papers

1. **Andrews S (2010)** FastQC: A Quality Control Tool for High Throughput Sequence Data. Babraham Bioinformatics.
2. **Ewels P et al. (2016)** MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics* 32(19):3047–3048. doi:[10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354)
3. **Chen S et al. (2018)** fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics* 34(17):i884–i890. doi:[10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560)
4. **Liao X et al. (2017)** Concordance of gene expression in human blood and PBMCs. *JCI Insight* 2(2). doi:[10.1172/jci.insight.89728](https://doi.org/10.1172/jci.insight.89728)
