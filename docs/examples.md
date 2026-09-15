# Examples

This gallery walks through the runnable example scripts in [`examples/`](https://github.com/SalvaaGomez/Pyntegrate/tree/main/examples). Each script reproduces one figure (or figure panel) from a FACT (SSRP1/SUPT16H) knockdown ChIP-seq/ATAC-seq/RNA-seq analysis, wrapping one Pyntegrate pipeline end-to-end — from raw signal/expression files to a saved figure. Every script follows the same shape: a **Configuration** section at the top, a docstring explaining the workflow, one or more pipeline calls, and an `OUTPUT_PATH` (saved to `outputs/figs/` at 300 dpi).

## Figure 1 — Meta-gene Profile and Cross-modality Correlations

!!! example "Pipelines: `GeneScaleRegionsProfilePipeline`, `ChIPCategoryCorrelationPipeline`, `RNAKnockdownCorrelationPipeline`, `ATACKnockdownCorrelationPipeline`, `PeakHeatmapPipeline`"
    Five independent panels: a deepTools-style scaled-gene-body meta-profile (1A), two hexbin correlations between the SSRP1 and SUPT16H knockdowns for ChIP, RNA and ATAC log2 fold-changes (1B–1D), and a ChIP-seq IP vs. INPUT profile on a separate dataset (1E).

=== "1A — Meta-gene profile"
    Scales every annotated gene body to a common length with fixed-length upstream/downstream flanks and averages SSRP1 vs. SUPT16H ChIP-seq signal across all genes.

    ![Meta-gene average-signal profile (SSRP1 vs SUPT16H)](images_examples/fig1a-genescale-profile.png)

    **Script:** [Fig 1A.py](https://github.com/SalvaaGomez/Pyntegrate/tree/main/examples/Fig%201A.py)

=== "1B — ChIP correlation"
    Hexbin density scatter correlating SSRP1 vs. SUPT16H ChIP log2 fold-change at TSS windows classified as FACT-high in at least one component.

    ![Hexbin density scatter: SSRP1 vs SUPT16H ChIP log2FC](images_examples/fig1b-chip-correlation.png)

    **Script:** [Fig 1B.py](https://github.com/SalvaaGomez/Pyntegrate/tree/main/examples/Fig%201B.py)

=== "1C — RNA knockdown correlation"
    Hexbin density scatter correlating RNA-seq log2 fold-change between the SSRP1 and SUPT16H knockdowns, for genes significantly changed in at least one condition.

    ![Hexbin density scatter: SSRP1 KD vs SUPT16H KD RNA log2FC](images_examples/fig1c-rna-knockdown-correlation.png)

    **Script:** [Fig 1C.py](https://github.com/SalvaaGomez/Pyntegrate/tree/main/examples/Fig%201C.py)

=== "1D — ATAC knockdown correlation"
    Hexbin density scatter correlating ATAC-seq log2 fold-change between the SUPT16H and SSRP1 knockdowns, for peaks significantly differential in at least one condition.

    ![Hexbin density scatter: SUPT16H KD vs SSRP1 KD ATAC log2FC](images_examples/fig1d-atac-knockdown-correlation.png)

    **Script:** [Fig 1D.py](https://github.com/SalvaaGomez/Pyntegrate/tree/main/examples/Fig%201D.py)

=== "1E — ChIP IP vs INPUT profile"
    ChIP-seq IP − INPUT average profile around TSS on an independent mouse (mm10) dataset.

    ![ChIP-seq IP vs INPUT average profile](images_examples/fig1e-peak-profile.png)

    **Script:** [Fig 1E.py](https://github.com/SalvaaGomez/Pyntegrate/tree/main/examples/Fig%201E.py)

---

## Figure 2 — ChIP-seq Heatmaps Split by RNA-seq Status

!!! example "Pipeline: `IntegratedHeatmapPipeline`"
    Loads SSRP1 (2A) / SUPT16H (2B) ChIP-seq signal, synchronises it with RNA-seq results, and produces an integrated heatmap + profile with genes overlaid by RNA up/down/unchanged status. No category filtering — all genes are included.

=== "2A — SSRP1"
    ![Integrated ChIP-seq (SSRP1) + RNA heatmap](images_examples/fig2a-integrated-heatmap-rna-ssrp1.png){: style="max-height:600px" }

=== "2B — SUPT16H"
    ![Integrated ChIP-seq (SUPT16H) + RNA heatmap](images_examples/fig2b-integrated-heatmap-rna-supt16h.png){: style="max-height:600px" }

**Script:** [Fig 2.py](https://github.com/SalvaaGomez/Pyntegrate/tree/main/examples/Fig%202.py)

---

## Figure 3 — CSV-only Decile Analysis (ChIP vs. RNA)

!!! example "Pipeline: `DecileScatterPipeline`"
    Works entirely from a pre-merged CSV/TSV. Sorts genes by SSRP1 ChIP log2 fold-change into deciles and plots the corresponding RNA-seq log2 fold-change (SSRP1 in 3A, SUPT16H in 3B) as a dual-axis boxplot + stripplot, with no category filtering.

=== "3A — SSRP1 RNA"
    ![Decile boxplot: ChIP deciles vs SSRP1 RNA log2FC](images_examples/fig3a-decile-rna-ssrp1.png)

=== "3B — SUPT16H RNA"
    ![Decile boxplot: ChIP deciles vs SUPT16H RNA log2FC](images_examples/fig3b-decile-rna-supt16h.png)

**Script:** [Fig 3.py](https://github.com/SalvaaGomez/Pyntegrate/tree/main/examples/Fig%203.py)

---

## Figure 4 — ChIP-seq Signal vs. ATAC-seq Status and Overlay

!!! example "Pipelines: `IntegratedHeatmapPipeline`, `OverlayProfilePipeline`"
    4A/4B integrate SSRP1 and SUPT16H ChIP-seq signal with ATAC-seq up/down/unchanged status (profile overlay only, heatmaps hidden). 4C overlays SSRP1 ChIP-seq against its own ATAC-seq track on twin y-axes, with the CONTROL ATAC-seq track shown for reference.

=== "4A — SSRP1 vs ATAC status"
    ![SSRP1 ChIP-seq signal split by ATAC-seq status](images_examples/fig4a-integrated-heatmap-atac-ssrp1.png)

=== "4B — SUPT16H vs ATAC status"
    ![SUPT16H ChIP-seq signal split by ATAC-seq status](images_examples/fig4b-integrated-heatmap-atac-supt16h.png)

=== "4C — SSRP1 ChIP/ATAC overlay"
    ![SSRP1 ChIP-seq vs ATAC-seq dual-axis overlay profile](images_examples/fig4c-overlay-profile-ssrp1.png)

=== "4D — SUPT16H ChIP/ATAC overlay"
    ![SUPT16H ChIP-seq vs ATAC-seq dual-axis overlay profile](images_examples/fig4d-overlay-profile-supt16h.png)

**Script:** [Fig 4AB.py](https://github.com/SalvaaGomez/Pyntegrate/tree/main/examples/Fig%204AB.py) (4A/4B) · [Fig 4CD.py](https://github.com/SalvaaGomez/Pyntegrate/tree/main/examples/Fig%204CD.py) (4C/4D)

---

## Figure 5 — ChIP-seq Heatmaps vs. RNA-seq Status, FACT-filtered

!!! example "Pipeline: `IntegratedHeatmapPipeline`"
    Same as Figure 2, but restricted to genes classified as FACT `high` or `low` for the corresponding component (SSRP1 FACT category for 5A, SUPT16H FACT category for 5B).

=== "5A — SSRP1"
    ![Integrated ChIP-seq (SSRP1) + RNA heatmap, FACT-filtered](images_examples/fig5a-integrated-heatmap-rna-ssrp1-fact.png){: style="max-height:600px" }

=== "5B — SUPT16H"
    ![Integrated ChIP-seq (SUPT16H) + RNA heatmap, FACT-filtered](images_examples/fig5b-integrated-heatmap-rna-supt16h-fact.png){: style="max-height:600px" }

**Script:** [Fig 5.py](https://github.com/SalvaaGomez/Pyntegrate/tree/main/examples/Fig%205.py)

---

## Figure 6 — RNA Violin Panels

!!! example "Pipeline: `RNAViolinPipeline`"
    Draws four horizontal violin panels of RNA-seq log2 fold-change (SSRP1/SUPT16H × FACT high/low) with significant genes (padj < 0.1) overlaid as points, all in one figure with no subplot divisions.

    ![RNA log2FC violin panels with significance overlay](images_examples/fig6-rna-violin.png)

**Script:** [Fig 6.py](https://github.com/SalvaaGomez/Pyntegrate/tree/main/examples/Fig%206.py)

---

## Figure 7 — CSV-only Decile Analysis (ChIP vs. RNA), FACT-filtered

!!! example "Pipeline: `DecileScatterPipeline`"
    Same as Figure 3, split into FACT `high`/`low` panels for both SSRP1 (7A/7B) and SUPT16H (7C/7D).

=== "7A — SSRP1 FACT high"
    ![Decile boxplot: SSRP1 FACT high](images_examples/fig7a-decile-atac-ssrp1-high.png)

=== "7B — SSRP1 FACT low"
    ![Decile boxplot: SSRP1 FACT low](images_examples/fig7b-decile-atac-ssrp1-low.png)

=== "7C — SUPT16H FACT high"
    ![Decile boxplot: SUPT16H FACT high](images_examples/fig7c-decile-atac-supt16h-high.png)

=== "7D — SUPT16H FACT low"
    ![Decile boxplot: SUPT16H FACT low](images_examples/fig7d-decile-atac-supt16h-low.png)

**Script:** [Fig 7.py](https://github.com/SalvaaGomez/Pyntegrate/tree/main/examples/Fig%207.py)

---

## Figure 8 — ChIP-seq Signal vs. ATAC-seq Status, FACT-filtered

!!! example "Pipeline: `IntegratedHeatmapPipeline`"
    Same as Figure 4A/4B, but restricted to genes classified as FACT `high` or `low` for the corresponding component.

=== "8A — SSRP1"
    ![SSRP1 ChIP-seq signal split by ATAC-seq status, FACT-filtered](images_examples/fig8a-integrated-heatmap-atac-ssrp1-fact.png)

=== "8B — SUPT16H"
    ![SUPT16H ChIP-seq signal split by ATAC-seq status, FACT-filtered](images_examples/fig8b-integrated-heatmap-atac-supt16h-fact.png)

**Script:** [Fig 8.py](https://github.com/SalvaaGomez/Pyntegrate/tree/main/examples/Fig%208.py)

---

## Figure 9 — ATAC-seq Category Heatmaps

!!! example "Pipeline: `CategoryHeatmapPipeline`"
    Loads the CONTROL, SSRP1 and SUPT16H ATAC-seq tracks and splits each into FACT `high`/`low` heatmaps (mean-sorted panel shown).

=== "9A — CONTROL"
    ![ATAC-seq category heatmap, CONTROL, sorted by mean](images_examples/fig9a-category-heatmap-control-atac.png){: style="max-height:600px" }

=== "9B — SSRP1"
    ![ATAC-seq category heatmap, SSRP1, sorted by mean](images_examples/fig9b-category-heatmap-ssrp1-atac.png){: style="max-height:600px" }

=== "9C — SUPT16H"
    ![ATAC-seq category heatmap, SUPT16H, sorted by mean](images_examples/fig9c-category-heatmap-supt16h-atac.png){: style="max-height:600px" }

**Script:** [Fig 9.py](https://github.com/SalvaaGomez/Pyntegrate/tree/main/examples/Fig%209.py)

---

## Figure 10 — Volcano Plots

!!! example "Pipeline: `VolcanoPipeline`"
    Volcano plots of ChIP-seq differential binding for SSRP1 (10A) and SUPT16H (10B), with a fixed gene highlight list (`SERPINE1`, `F3`, `CYR61`, `TERT`, `RPSAP52`, `NEAT1`) and no threshold-based auto-labeling.

=== "10A — SSRP1 ChIP volcano"
    ![Volcano plot: SSRP1 ChIP differential binding](images_examples/fig10a-volcano-ssrp1-chip.png)

=== "10B — SUPT16H ChIP volcano"
    ![Volcano plot: SUPT16H ChIP differential binding](images_examples/fig10b-volcano-supt16h-chip.png)

**Script:** [Fig 10.py](https://github.com/SalvaaGomez/Pyntegrate/tree/main/examples/Fig%2010.py)

---

## Figures 11–12 — Circos Plots Comparing Multi-signal Changes (-log10 padj)

!!! example "Pipeline: `CircPlotPipeline`"
    Generates a circos plot per condition (SSRP1 → Figure 11, SUPT16H → Figure 12) comparing ChIP-seq, ATAC-seq, and RNA-seq `-log10(padj)` (via `value_transform`/`value_label`) for a configurable gene list (defaults to `SERPINE1`, `F3`, `CYR61`, `TERT`, `RPSAP52`, `NEAT1`), with a threshold of `0.1` suited to the transformed scale. Three concentric rings represent the three signal types, with connecting ribbons colored by magnitude. Gene data can optionally be exported to TSV alongside the figure.

=== "11 — SSRP1"
    ![Circos plot for SSRP1 condition, -log10(padj)](images_examples/fig11-circplot-ssrp1-padj.png)

=== "12 — SUPT16H"
    ![Circos plot for SUPT16H condition, -log10(padj)](images_examples/fig12-circplot-supt16h-padj.png)

**Script:** [Fig 11-12.py](https://github.com/SalvaaGomez/Pyntegrate/tree/main/examples/Fig%2011-12.py)

---

## Figures 13–14 — Circos Plots Comparing Multi-signal Changes (log2FC)

!!! example "Pipeline: `CircPlotPipeline`"
    Same gene list and layout as Figures 11–12, but plots log2 fold-change instead of `-log10(padj)` for each signal. SSRP1 → Figure 13, SUPT16H → Figure 14.

=== "13 — SSRP1"
    ![Circos plot for SSRP1 condition](images_examples/fig13-circplot-ssrp1.png)

=== "14 — SUPT16H"
    ![Circos plot for SUPT16H condition](images_examples/fig14-circplot-supt16h.png)

**Script:** [Fig 13-14.py](https://github.com/SalvaaGomez/Pyntegrate/tree/main/examples/Fig%2013-14.py)

---

## Supplementary Figure 1 — Decile Analysis (ChIP vs. ATAC) and Quadrant Profiles

!!! example "Pipelines: `DecileScatterPipeline`, `GeneGroupProfilePipeline`"
    S1A/S1B repeat the Figure 3 decile analysis using ATAC-seq log2 fold-change instead of RNA-seq (SSRP1 / SUPT16H, unfiltered). S1C/S1D classify genes into four ATAC×RNA log2FC quadrants (threshold ±1.0) and plot the average ChIP-seq TSS profile per quadrant.

=== "S1A — SSRP1 ATAC deciles"
    ![Decile boxplot: ChIP deciles vs SSRP1 ATAC log2FC](images_examples/figs1a-decile-atac-ssrp1.png)

=== "S1B — SUPT16H ATAC deciles"
    ![Decile boxplot: ChIP deciles vs SUPT16H ATAC log2FC](images_examples/figs1b-decile-atac-supt16h.png)

=== "S1C — SSRP1 quadrant profile"
    ![ChIP-seq profile by ATAC/RNA quadrant, SSRP1](images_examples/figs1c-quadrant-profile-ssrp1.png)

=== "S1D — SUPT16H quadrant profile"
    ![ChIP-seq profile by ATAC/RNA quadrant, SUPT16H](images_examples/figs1d-quadrant-profile-supt16h.png)

**Script:** [Fig S1AB.py](https://github.com/SalvaaGomez/Pyntegrate/tree/main/examples/Fig%20S1AB.py) (S1A/S1B) · [Fig S1CD.py](https://github.com/SalvaaGomez/Pyntegrate/tree/main/examples/Fig%20S1CD.py) (S1C/S1D)

---

## Supplementary Figure 2 — ChIP/ATAC Overlay Profiles, FACT-filtered

!!! example "Pipeline: `OverlayProfilePipeline`"
    Dual-axis ChIP-seq vs. ATAC-seq overlay profiles (with the CONTROL ATAC-seq track shown for reference), split by FACT `high`/`low` for both SSRP1 (S2A/S2B) and SUPT16H (S2C/S2D).

=== "S2A — SSRP1 FACT high"
    ![SSRP1 ChIP/ATAC overlay, FACT high](images_examples/figs2a-overlay-atac-chip-ssrp1-high.png)

=== "S2B — SSRP1 FACT low"
    ![SSRP1 ChIP/ATAC overlay, FACT low](images_examples/figs2b-overlay-atac-chip-ssrp1-low.png)

=== "S2C — SUPT16H FACT high"
    ![SUPT16H ChIP/ATAC overlay, FACT high](images_examples/figs2c-overlay-atac-chip-supt16h-high.png)

=== "S2D — SUPT16H FACT low"
    ![SUPT16H ChIP/ATAC overlay, FACT low](images_examples/figs2d-overlay-atac-chip-supt16h-low.png)

**Script:** [Fig S2.py](https://github.com/SalvaaGomez/Pyntegrate/tree/main/examples/Fig%20S2.py)

---

## Extra — Multi-signal Comparison (3-way Scatter)

!!! example "Pipeline: `MultiSignalComparisonPipeline`"
    Compares log2 fold-change values across ChIP-seq, ATAC-seq and RNA-seq from a single expression-results CSV and produces a 1×3 panel of 2D scatter plots (RNA vs. ChIP, RNA vs. ATAC, ChIP vs. ATAC). Kept in [`examples/extra/`](https://github.com/SalvaaGomez/Pyntegrate/tree/main/examples/extra) as a standalone reference, saved by default; set `OUTPUT_PATH` to `None` to skip saving.

    ![Multi-signal scatter comparison: RNA vs ChIP vs ATAC log2FC](images_examples/fig-extra-multisignal-comparison.png)

**Script:** [extra/multisignal_comparison.py](https://github.com/SalvaaGomez/Pyntegrate/tree/main/examples/extra/multisignal_comparison.py)

---

## How to Run the Examples

### Prerequisites

1. Install Pyntegrate — see [Installation](installation.md)
2. Have the data files available:
      - A gffutils FeatureDB (`.db` file)
      - Signal files (`.bw`, `.bam`, `.bed`, ...)
      - CSV/TSV tables with gene annotations and expression results

### Basic Usage

Script filenames contain spaces, so quote them:

```bash
cd examples/

# Run Figure 1E (ChIP-seq IP vs INPUT profile)
python "Fig 1E.py"

# Run Figure 10 (volcano plots)
python "Fig 10.py"

# Run the CSV-only decile analysis (Figure 3)
python "Fig 3.py"
```

### Customizing Paths and Parameters

Edit the **Configuration** section at the top of each script to match your data locations and analysis parameters. For example, in `Fig 1E.py`:

```python
DB_PATH = "/path/to/your/annotations.db"
IP_BW = "/path/to/your/ip_signal.bw"
INPUT_BW = "/path/to/your/input_signal.bw"
GENOME = "mm10"  # or "hg19"
```

## Quick Reference Table

| Figure | Script | Pipeline | Output |
|---|---|---|---|
| 1A | `Fig 1A.py` | `GeneScaleRegionsProfilePipeline` | Scaled-gene-body meta-profile |
| 1B | `Fig 1B.py` | `ChIPCategoryCorrelationPipeline` | ChIP log2FC hexbin correlation |
| 1C | `Fig 1C.py` | `RNAKnockdownCorrelationPipeline` | RNA log2FC hexbin correlation |
| 1D | `Fig 1D.py` | `ATACKnockdownCorrelationPipeline` | ATAC log2FC hexbin correlation |
| 1E | `Fig 1E.py` | `PeakHeatmapPipeline` | IP − INPUT profile |
| 2A/2B | `Fig 2.py` | `IntegratedHeatmapPipeline` | Heatmap + profile vs. RNA status |
| 3A/3B | `Fig 3.py` | `DecileScatterPipeline` | Decile boxplot vs. RNA log2FC |
| 4A/4B | `Fig 4AB.py` | `IntegratedHeatmapPipeline` | Profile overlay vs. ATAC status |
| 4C/4D | `Fig 4CD.py` | `OverlayProfilePipeline` | ChIP/ATAC dual-axis overlay |
| 5A/5B | `Fig 5.py` | `IntegratedHeatmapPipeline` | Heatmap + profile vs. RNA status, FACT-filtered |
| 6 | `Fig 6.py` | `RNAViolinPipeline` | Violin panels |
| 7A–7D | `Fig 7.py` | `DecileScatterPipeline` | Decile boxplot vs. RNA log2FC, FACT-filtered |
| 8A/8B | `Fig 8.py` | `IntegratedHeatmapPipeline` | Profile overlay vs. ATAC status, FACT-filtered |
| 9A–9C | `Fig 9.py` | `CategoryHeatmapPipeline` | ATAC-seq category heatmaps |
| 10A/10B | `Fig 10.py` | `VolcanoPipeline` | Volcano plot |
| 11/12 | `Fig 11-12.py` | `CircPlotPipeline` | Circos plots (-log10 padj) |
| 13/14 | `Fig 13-14.py` | `CircPlotPipeline` | Circos plots (log2FC) |
| S1A/S1B | `Fig S1AB.py` | `DecileScatterPipeline` | Decile boxplot vs. ATAC log2FC |
| S1C/S1D | `Fig S1CD.py` | `GeneGroupProfilePipeline` | ChIP profile by ATAC/RNA quadrant |
| S2A–S2D | `Fig S2.py` | `OverlayProfilePipeline` | ChIP/ATAC overlay, FACT-filtered |
| Extra | `extra/multisignal_comparison.py` | `MultiSignalComparisonPipeline` | Multi-signal scatter |

## Troubleshooting

!!! failure "Common issues"
    1. Verify that all data file paths are correct and the files exist
    2. Check that the genome version matches your data
    3. Ensure dependencies (`bedtools`, HOMER if needed) are installed and on your `PATH`
    4. Check that your CSV/TSV files have the columns referenced by the pipeline
    5. For more detail, inspect the full error traceback

For additional help, visit [GitHub Issues](https://github.com/SalvaaGomez/Pyntegrate/issues).
