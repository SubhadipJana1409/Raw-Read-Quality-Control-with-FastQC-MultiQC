import json
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

BASE    = Path("/home/claude/Day_01_Final")
FIG_DIR = BASE / "figures"
RES_DIR = BASE / "results"

SAMPLES = {
    "SRR5223500": {"label": "HC_1", "group": "HC"},
    "SRR5223501": {"label": "HC_2", "group": "HC"},
    "SRR5223502": {"label": "HC_3", "group": "HC"},
    "SRR5223503": {"label": "SS_1", "group": "SS"},
    "SRR5223504": {"label": "SS_2", "group": "SS"},
    "SRR5223505": {"label": "SS_3", "group": "SS"},
}
HC_C = "#3498db"; SS_C = "#e74c3c"
PASS_C = "#2ecc71"; WARN_C = "#f39c12"; FAIL_C = "#e74c3c"
STATUS_N = {"PASS": 0, "WARN": 1, "FAIL": 2}

# ── Parse fastp JSONs ────────────────────────────────────────────────────────
records = []
for srr, meta in SAMPLES.items():
    fp = json.load(open(BASE / f"qc/fastp_reports/{srr}_fastp.json"))
    records.append({
        "srr": srr, "label": meta["label"], "group": meta["group"],
        "reads_before_M": round(fp["summary"]["before_filtering"]["total_reads"]/1e6, 3),
        "reads_after_M":  round(fp["summary"]["after_filtering"]["total_reads"]/1e6, 3),
        "pct_kept":  round(fp["summary"]["after_filtering"]["total_reads"] /
                           fp["summary"]["before_filtering"]["total_reads"]*100, 1),
        "q20_before": round(fp["summary"]["before_filtering"]["q20_rate"]*100, 2),
        "q20_after":  round(fp["summary"]["after_filtering"]["q20_rate"]*100, 2),
        "q30_before": round(fp["summary"]["before_filtering"]["q30_rate"]*100, 2),
        "q30_after":  round(fp["summary"]["after_filtering"]["q30_rate"]*100, 2),
        "gc_before":  round(fp["summary"]["before_filtering"]["gc_content"]*100, 2),
        "gc_after":   round(fp["summary"]["after_filtering"]["gc_content"]*100, 2),
        "dup_rate":   round(fp["duplication"]["rate"]*100, 2),
        "too_short":  fp["filtering_result"]["too_short_reads"],
        "adapter_trimmed": fp["adapter_cutting"]["adapter_trimmed_reads"],
        "insert_peak": fp["insert_size"]["peak"],
    })
df = pd.DataFrame(records)
df.to_csv(RES_DIR / "fastp_summary_all_samples.csv", index=False)
print("✓ fastp summary loaded for", len(df), "samples")

# ── Parse FastQC summaries ───────────────────────────────────────────────────
def load_summaries(qc_dir, stage):
    rows = []
    for f in sorted(Path(qc_dir).rglob("summary.txt")):
        dname = f.parent.name.replace("_fastqc","").replace(".trimmed","")
        parts = dname.rsplit("_", 1)
        srr = parts[0]; read = f"R{parts[1]}" if len(parts)==2 else "R1"
        meta = SAMPLES.get(srr, {})
        for line in open(f):
            status, module, _ = line.strip().split("\t")
            rows.append({"srr":srr, "label":meta.get("label",srr),
                         "group":meta.get("group","?"),
                         "read":read, "module":module, "status":status, "stage":stage})
    return pd.DataFrame(rows)

df_raw  = load_summaries(BASE/"qc/raw",     "Raw")
df_trim = load_summaries(BASE/"qc/trimmed", "Trimmed")
# add numeric status column
df_raw["status_n"]  = df_raw["status"].map(STATUS_N)
df_trim["status_n"] = df_trim["status"].map(STATUS_N)
print("✓ FastQC summaries loaded:", len(df_raw), "raw rows,", len(df_trim), "trim rows")

# ─────────────────────────────────────────────────────────────────────────────
# FIG 1 — FastQC Heatmap Raw
# ─────────────────────────────────────────────────────────────────────────────
cmap = matplotlib.colors.ListedColormap([PASS_C, WARN_C, FAIL_C])

def fastqc_heatmap(df_stage, title, save_path):
    ds = df_stage.copy()
    ds["col"] = ds["label"] + "\n" + ds["read"]
    pivot_n = ds.pivot_table(index="module", columns="col", values="status_n", aggfunc="first")
    pivot_s = ds.pivot_table(index="module", columns="col", values="status",   aggfunc="first")
    col_order = sorted(pivot_n.columns, key=lambda x: (0 if "HC" in x else 1, x))
    pivot_n = pivot_n[col_order]; pivot_s = pivot_s[col_order]

    fig, ax = plt.subplots(figsize=(max(13, len(col_order)*0.75), 5))
    ax.imshow(pivot_n.values.astype(float), cmap=cmap, vmin=0, vmax=2, aspect="auto")
    ax.set_xticks(range(len(col_order))); ax.set_xticklabels(col_order, fontsize=7.5)
    ax.set_yticks(range(len(pivot_n.index))); ax.set_yticklabels(pivot_n.index, fontsize=8)

    for i in range(len(pivot_n.index)):
        for j in range(len(col_order)):
            val = str(pivot_s.iloc[i,j]) if not pd.isna(pivot_s.iloc[i,j]) else "?"
            ax.text(j, i, val[0], ha="center", va="center", fontsize=7,
                    color="white" if val=="FAIL" else "black", fontweight="bold")

    hc_n = len([c for c in col_order if "HC" in c])
    ax.axvline(hc_n-0.5, color="black", linewidth=2)
    ax.text((hc_n-1)/2, -0.9, "Healthy Controls (HC)", ha="center",
            fontsize=9, fontweight="bold", color=HC_C)
    ax.text(hc_n+(len(col_order)-hc_n-1)/2, -0.9, "Septic Shock (SS)", ha="center",
            fontsize=9, fontweight="bold", color=SS_C)

    patches = [mpatches.Patch(facecolor=PASS_C, label="PASS"),
               mpatches.Patch(facecolor=WARN_C, label="WARN"),
               mpatches.Patch(facecolor=FAIL_C, label="FAIL")]
    ax.legend(handles=patches, loc="upper right", bbox_to_anchor=(1.15,1), fontsize=9)
    ax.set_title(title, fontsize=11, fontweight="bold", pad=22)
    plt.tight_layout()
    fig.savefig(save_path, dpi=150, bbox_inches="tight"); plt.close()
    print(f"✓ {save_path.name}")

fastqc_heatmap(df_raw,  "FastQC Module Status — Raw Reads | GSE96960 PBMC (n=6)",
               FIG_DIR/"fig1_fastqc_heatmap_raw.png")
fastqc_heatmap(df_trim, "FastQC Module Status — Trimmed Reads | GSE96960 PBMC (n=6)",
               FIG_DIR/"fig2_fastqc_heatmap_trimmed.png")

# ─────────────────────────────────────────────────────────────────────────────
# FIG 3 — Before vs After bars (Q30, GC, Reads)
# ─────────────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
x = np.arange(len(df)); w = 0.35
colors = [HC_C if g=="HC" else SS_C for g in df["group"]]

for ax, (b,a,ylabel,title,thresh) in zip(axes, [
    ("reads_before_M","reads_after_M","Total Reads (M)","A. Read Count",None),
    ("q30_before","q30_after","Q30 Rate (%)","B. Q30 Rate",80),
    ("gc_before","gc_after","GC Content (%)","C. GC Content",50),
]):
    ax.bar(x-w/2, df[b], w, color=[c+"88" for c in colors], edgecolor="black", linewidth=0.7, label="Before")
    ax.bar(x+w/2, df[a], w, color=colors,                   edgecolor="black", linewidth=0.7, label="After")
    if thresh:
        ax.axhline(thresh, color="red", linestyle="--", linewidth=1, label=f"{thresh} threshold")
    ax.set_xticks(x); ax.set_xticklabels(df["label"], fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9); ax.set_title(title, fontweight="bold", fontsize=10)
    ax.axvline(2.5, color="grey", linestyle=":", linewidth=1)
    ymax = ax.get_ylim()[1]
    ax.text(1, ymax*0.96, "HC", ha="center", color=HC_C, fontweight="bold", fontsize=9)
    ax.text(4, ymax*0.96, "SS", ha="center", color=SS_C, fontweight="bold", fontsize=9)
    ax.legend(fontsize=8)

fig.suptitle("Before vs After fastp Trimming — All 6 Samples (GSE96960 PBMC)",
             fontsize=12, fontweight="bold", y=1.02)
plt.tight_layout()
fig.savefig(FIG_DIR/"fig3_before_after_all_samples.png", dpi=150, bbox_inches="tight")
plt.close(); print("✓ fig3_before_after_all_samples.png")

# ─────────────────────────────────────────────────────────────────────────────
# FIG 4 — Lollipop Q30/Q20 improvement
# ─────────────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
for ax, (bc, ac, title) in zip(axes, [
    ("q30_before","q30_after","Q30 Rate: Before → After"),
    ("q20_before","q20_after","Q20 Rate: Before → After"),
]):
    for i, row in df.iterrows():
        c = HC_C if row["group"]=="HC" else SS_C
        ax.plot([row[bc], row[ac]], [i,i], color=c, linewidth=2.5, zorder=1)
        ax.scatter(row[bc], i, color=c, s=90, marker="o", zorder=2, alpha=0.55)
        ax.scatter(row[ac], i, color=c, s=130, marker="D", zorder=3)
        ax.text(row[ac]+0.25, i, f"{row[ac]:.1f}%", va="center", fontsize=9)

    ax.set_yticks(range(len(df))); ax.set_yticklabels(df["label"], fontsize=10)
    ax.set_xlabel("Quality Rate (%)", fontsize=9); ax.set_title(title, fontweight="bold", fontsize=10)
    ax.axvline(80 if "30" in bc else 90, color="red", linestyle="--", linewidth=1, alpha=0.5)
    ax.grid(axis="x", alpha=0.3)
    ax.legend(handles=[
        mpatches.Patch(facecolor=HC_C, label="HC"),
        mpatches.Patch(facecolor=SS_C, label="SS"),
        plt.Line2D([0],[0],marker="o",color="grey",label="Before",markersize=8,linestyle="None"),
        plt.Line2D([0],[0],marker="D",color="grey",label="After", markersize=8,linestyle="None"),
    ], fontsize=8, loc="lower right")

fig.suptitle("Per-Sample Quality Score Improvement After fastp",
             fontsize=12, fontweight="bold", y=1.02)
plt.tight_layout()
fig.savefig(FIG_DIR/"fig4_quality_improvement_lollipop.png", dpi=150, bbox_inches="tight")
plt.close(); print("✓ fig4_quality_improvement_lollipop.png")

# ─────────────────────────────────────────────────────────────────────────────
# FIG 5 — Stacked reads fate
# ─────────────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(11, 5))
passed  = df["reads_after_M"].values
removed = df["reads_before_M"].values - df["reads_after_M"].values
ax.bar(df["label"], passed,  color=colors, edgecolor="black", linewidth=0.7, label="Passed Filter")
ax.bar(df["label"], removed, bottom=passed,
       color=["#aed6f1" if g=="HC" else "#f1948a" for g in df["group"]],
       edgecolor="black", linewidth=0.7, alpha=0.85, label="Removed (too short)")
for i, row in df.iterrows():
    ax.text(i, row["reads_before_M"]+0.01, f"{row['pct_kept']}%",
            ha="center", fontsize=9, fontweight="bold")
ax.axvline(2.5, color="grey", linestyle=":", linewidth=1.5)
ymax = ax.get_ylim()[1]
ax.text(1, ymax*0.95, "Healthy Controls", ha="center", color=HC_C, fontweight="bold")
ax.text(4, ymax*0.95, "Septic Shock",     ha="center", color=SS_C, fontweight="bold")
ax.set_ylabel("Reads (millions)", fontsize=10)
ax.set_title("Read Fate After fastp Filtering — All 6 Samples\n(% above bar = reads passing filter)",
             fontsize=11, fontweight="bold")
ax.legend(fontsize=9)
plt.tight_layout()
fig.savefig(FIG_DIR/"fig5_reads_fate_all_samples.png", dpi=150, bbox_inches="tight")
plt.close(); print("✓ fig5_reads_fate_all_samples.png")

# ─────────────────────────────────────────────────────────────────────────────
# FIG 6 — Duplication & Insert Size
# ─────────────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

ax = axes[0]
bars = ax.bar(df["label"], df["dup_rate"], color=colors, edgecolor="black", linewidth=0.7)
ax.bar_label(bars, fmt="%.2f%%", padding=3, fontsize=9, fontweight="bold")
ax.axhline(20, color="orange", linestyle="--", linewidth=1, label="20% warning")
ax.axhline(40, color="red",    linestyle="--", linewidth=1, label="40% fail")
ax.axvline(2.5, color="grey", linestyle=":", linewidth=1)
ax.set_ylabel("Duplication Rate (%)", fontsize=10)
ax.set_title("PCR Duplication Rate per Sample", fontweight="bold", fontsize=10)
ax.legend(fontsize=8)
ymax = ax.get_ylim()[1]
ax.text(1, ymax*0.95, "HC", ha="center", color=HC_C, fontweight="bold")
ax.text(4, ymax*0.95, "SS", ha="center", color=SS_C, fontweight="bold")

ax = axes[1]
bars2 = ax.bar(df["label"], df["insert_peak"], color=colors, edgecolor="black", linewidth=0.7)
ax.bar_label(bars2, fmt="%d bp", padding=3, fontsize=9, fontweight="bold")
ax.axvline(2.5, color="grey", linestyle=":", linewidth=1)
ax.set_ylabel("Insert Size Peak (bp)", fontsize=10)
ax.set_title("Insert Size Peak per Sample", fontweight="bold", fontsize=10)
ymax2 = ax.get_ylim()[1]
ax.text(1, ymax2*0.95, "HC", ha="center", color=HC_C, fontweight="bold")
ax.text(4, ymax2*0.95, "SS", ha="center", color=SS_C, fontweight="bold")

fig.suptitle("Library Quality Metrics — All 6 Samples", fontsize=12, fontweight="bold", y=1.02)
plt.tight_layout()
fig.savefig(FIG_DIR/"fig6_library_quality_all.png", dpi=150, bbox_inches="tight")
plt.close(); print("✓ fig6_library_quality_all.png")

# ─────────────────────────────────────────────────────────────────────────────
# Save final summary table
# ─────────────────────────────────────────────────────────────────────────────
df.to_csv(RES_DIR/"qc_summary_all_samples.csv", index=False)
print("\n" + "="*65)
print("FINAL QC SUMMARY — ALL 6 SAMPLES (GSE96960)")
print("="*65)
print(df[["label","group","reads_before_M","reads_after_M",
          "pct_kept","q30_before","q30_after","dup_rate"]].to_string(index=False))
print("\n✅ All 6 figures + results complete!")
