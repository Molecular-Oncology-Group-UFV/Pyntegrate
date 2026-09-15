# Pyntegrate Examples

This directory contains example scripts reproducing a FACT (SSRP1/SUPT16H) knockdown ChIP-seq/ATAC-seq/RNA-seq analysis, one script per figure (or figure panel), each demonstrating one of the main analysis pipelines available in Pyntegrate.

Most examples follow a consistent pattern:
1. **Configuration section** at the top with file paths and parameters
2. **Docstring** explaining the workflow
3. **One or more pipeline calls** that handle all the complexity
4. **`OUTPUT_PATH`** (or `OUTPUT_PATH_*` for multi-figure examples) — every example saves its figure(s) to `outputs/figs/` at 300 dpi via `fig.savefig(OUTPUT_PATH, dpi=300)` before displaying them. The `outputs/figs/` directory is created automatically if it doesn't exist yet. Set `OUTPUT_PATH` to `None` to skip saving.

Script filenames contain spaces (e.g. `Fig 1A.py`) — quote them when running from a shell.

## Examples

### Fig 1A.py — Figure 1A
**GeneScaleRegionsProfilePipeline** - Meta-gene profile over scaled gene bodies

Scales every annotated gene body to a common length (with fixed-length upstream/downstream flanks), averages SSRP1 vs. SUPT16H ChIP-seq signal across all genes, and plots one profile line per signal file — analogous to deepTools `plotProfile --scaleRegions`.

**Output figures:**
- Meta-gene average-signal profile, one line per input signal file

### Fig 1B.py — Figure 1B
**ChIPCategoryCorrelationPipeline** - ChIP log2FC correlation filtered by FACT category

Loads a merged TSV with ChIP log2FC and FACT category columns, retains rows classified as "high" for at least one of SSRP1/SUPT16H, and produces a hexbin density scatter plot correlating the two ChIP log2FC columns, annotated with the Spearman correlation and sample size.

**Output figures:**
- Hexbin density scatter plot (SSRP1 vs SUPT16H ChIP log2FC), annotated with R and n

### Fig 1C.py — Figure 1C
**RNAKnockdownCorrelationPipeline** - RNA-seq knockdown log2FC correlation for significant genes

Loads a merged TSV with RNA log2FC and adjusted p-value columns for the SSRP1 and SUPT16H knockdowns, retains genes significantly changed in at least one condition, and produces a hexbin density scatter plot correlating the two RNA log2FC columns.

**Output figures:**
- Hexbin density scatter plot (SSRP1 KD vs SUPT16H KD RNA log2FC), annotated with R and n

### Fig 1D.py — Figure 1D
**ATACKnockdownCorrelationPipeline** - ATAC-seq knockdown log2FC correlation for significant peaks

Loads a merged TSV with ATAC log2FC and adjusted p-value columns for the SUPT16H and SSRP1 knockdowns, retains peaks significantly changed in at least one condition, and produces a hexbin density scatter plot correlating the two ATAC log2FC columns.

**Output figures:**
- Hexbin density scatter plot (SUPT16H KD vs SSRP1 KD ATAC log2FC), annotated with R and n

### Fig 1E.py — Figure 1E
**PeakHeatmapPipeline** - ChIP-seq IP vs INPUT profile

Loads IP and INPUT signal files (BigWig) for an independent mouse (mm10) dataset, generates signal arrays around transcription start sites (TSS), computes the IP − INPUT subtracted signal, and produces a profile plot (heatmap panels are also available from the pipeline but not saved by this script).

**Output figures:**
- Average profile (IP vs INPUT)

### Fig 2.py — Figure 2A/2B
**IntegratedHeatmapPipeline** - Signal heatmap integrated with RNA-seq status

Loads SSRP1 (2A) / SUPT16H (2B) ChIP-seq signal, synchronises it with RNA-seq results, and produces an integrated heatmap + profile with genes overlaid by RNA up/down/unchanged status. No category filtering.

**Output figures:**
- Signal heatmap with RNA expression categories overlaid (one per condition)

### Fig 3.py — Figure 3A/3B
**DecileScatterPipeline** - CSV-only decile analysis (ChIP vs RNA)

Works entirely from a pre-merged CSV/TSV. Sorts genes by SSRP1 ChIP log2 fold-change into deciles and plots the corresponding RNA-seq log2 fold-change (SSRP1 in 3A, SUPT16H in 3B) as a dual-axis boxplot + stripplot. No category filtering.

**Output figures:**
- Dual-axis figure: boxplot (ChIP deciles) + stripplot (RNA log2FC), one per condition

### Fig 4AB.py — Figure 4A/4B
**IntegratedHeatmapPipeline** - Signal profile integrated with ATAC-seq status

Same as Figure 2, but against ATAC-seq up/down/unchanged status instead of RNA-seq, with heatmap panels hidden (profile overlay only).

**Output figures:**
- Signal profile with ATAC-seq expression categories overlaid (one per condition)

### Fig 4CD.py — Figure 4C/4D
**OverlayProfilePipeline** - Dual-axis ChIP/ATAC profile overlay

Loads two signal files independently (a condition's ChIP-seq track and its matching ATAC-seq track, with the CONTROL ATAC-seq track shown for reference), finds genes present in both datasets, and produces a dual-y-axis profile figure. 4C covers SSRP1, 4D covers SUPT16H.

**Output figures:**
- Dual-axis profile plot (ChIP vs ATAC), one per condition

### Fig 5.py — Figure 5A/5B
**IntegratedHeatmapPipeline** - Signal heatmap integrated with RNA-seq status, FACT-filtered

Same as Figure 2, but restricted to genes classified as FACT `high` or `low` for the corresponding component (SSRP1 FACT category for 5A, SUPT16H FACT category for 5B).

**Output figures:**
- Signal heatmap with RNA expression categories overlaid, FACT-filtered (one per condition)

### Fig 6.py — Figure 6
**RNAViolinPipeline** - RNA expression distribution with significance overlay

Shows the full distribution of RNA log2 fold-change values as four violin panels (SSRP1/SUPT16H × FACT high/low, via `VIOLIN_SPECS`), all in one horizontal figure with no subplot divisions, with significant genes (padj < 0.1) highlighted as points.

**Output figures:**
- Violin panels with significant-gene overlay, all in one horizontal figure

### Fig 7.py — Figure 7A–7D
**DecileScatterPipeline** - CSV-only decile analysis (ChIP vs RNA), FACT-filtered

Same as Figure 3, split into FACT `high`/`low` panels for both SSRP1 (7A/7B) and SUPT16H (7C/7D).

**Output figures:**
- Dual-axis figure: boxplot (ChIP deciles) + stripplot (RNA log2FC), one per FACT category/condition combination

### Fig 8.py — Figure 8A/8B
**IntegratedHeatmapPipeline** - Signal profile integrated with ATAC-seq status, FACT-filtered

Same as Figure 4A/4B, but restricted to genes classified as FACT `high` or `low` for the corresponding component.

**Output figures:**
- Signal profile with ATAC-seq expression categories overlaid, FACT-filtered (one per condition)

### Fig 9.py — Figure 9A–9C
**CategoryHeatmapPipeline** - ATAC-seq heatmaps split by FACT category

Loads the CONTROL (9A), SSRP1 (9B) and SUPT16H (9C) ATAC-seq tracks, maps genes to FACT category labels (high/low) from a CSV table, and produces heatmaps split by category (only the mean-sorted panel is saved).

**Output figures:**
- Category heatmap sorted by mean, one per track

### Fig 10.py — Figure 10A/10B
**VolcanoPipeline** - Volcano plot from ChIP-seq differential binding results

Loads ChIP-seq differential binding results from a CSV file and produces a volcano plot per condition (SSRP1 for 10A, SUPT16H for 10B) with points coloured by regulation status, intensity-based alpha transparency, and a fixed gene highlight list.

**Output figures:**
- Volcano plot with labelled significant genes, one per condition

### Fig 11-12.py — Figures 11/12
**CircPlotPipeline** - Circos plot comparing multi-signal changes (-log10 adjusted p-value)

Generates circos plots comparing ChIP-seq, ATAC-seq, and RNA-seq `-log10(padj)` values (via `value_transform`/`value_label`) for a configurable gene list (defaults to `SERPINE1`, `F3`, `CYR61`, `TERT`, `RPSAP52`, `NEAT1`) across two conditions (SSRP1 → Figure 11, SUPT16H → Figure 12), with a threshold of `0.1` suited to the transformed scale.

The circular visualization shows three concentric rings for each signal type, with connections highlighting genes that match user-specified combinations of high/low signal changes.

**Output figures:**
- Circos plot for SSRP1 condition
- Circos plot for SUPT16H condition
- (Optional) TSV exports of the underlying gene data

### Fig 13-14.py — Figures 13/14
**CircPlotPipeline** - Circos plot comparing multi-signal changes (log2 fold-change)

Same as Figures 11/12, but plots log2 fold-change instead of `-log10(padj)`. SSRP1 → Figure 13, SUPT16H → Figure 14.

**Output figures:**
- Circos plot for SSRP1 condition
- Circos plot for SUPT16H condition
- (Optional) TSV exports of the underlying gene data

### Fig S1AB.py — Supplementary Figure S1A/S1B
**DecileScatterPipeline** - CSV-only decile analysis (ChIP vs ATAC)

Repeats the Figure 3 decile analysis using ATAC-seq log2 fold-change instead of RNA-seq (SSRP1 for S1A, SUPT16H for S1B), unfiltered by category.

**Output figures:**
- Dual-axis figure: boxplot (ChIP deciles) + stripplot (ATAC log2FC), one per condition

### Fig S1CD.py — Supplementary Figure S1C/S1D
**GeneGroupProfilePipeline** - ChIP-seq profile split by ATAC/RNA log2FC quadrants

Classifies genes into four quadrants based on ATAC-seq and RNA-seq log2 fold-change thresholds (±1.0) read from `final_table_filtered.tsv` (ATAC high/RNA high, ATAC high/RNA low, ATAC low/RNA low, ATAC low/RNA high), then plots the average ChIP-seq signal around the TSS for each quadrant as a separate line. S1C covers SSRP1, S1D covers SUPT16H.

**Output figures:**
- Single-axis ChIP-seq average-signal profile, one line per quadrant (one figure per condition)

### Fig S2.py — Supplementary Figure S2A–S2D
**OverlayProfilePipeline** - Dual-axis ChIP/ATAC profile overlay, FACT-filtered

Same as Figure 4C/4D, split into FACT `high`/`low` panels for both SSRP1 (S2A/S2B) and SUPT16H (S2C/S2D).

**Output figures:**
- Dual-axis profile plot (ChIP vs ATAC), one per FACT category/condition combination

### extra/multisignal_comparison.py — Extra
**MultiSignalComparisonPipeline** - Multi-signal 3-way comparison

Compares log2 fold-change values across multiple modalities (e.g., ChIP-seq, ATAC-seq, RNA-seq) from a single expression-results CSV and produces 2D scatter plots showing the relationships between signal and expression changes. Kept as a standalone reference outside the main figure set (`OUTPUT_PATH` is set, so the figure is saved by default; set it to `None` to skip saving).

**Output figures:**
- 3-way scatter plot panel (RNA vs ChIP, RNA vs ATAC, ChIP vs ATAC)

## Running the Examples

### Prerequisites
1. Install Pyntegrate: `pip install -e ../` (from this directory)
2. Have the data files available (paths configured in each example)
3. Depending on the example, you'll need some combination of:
   - A gffutils FeatureDB (`.db` file)
   - Signal files (`.bw`, `.bam`, `.bed`, ...)
   - CSV/TSV tables with gene annotations and expression results

### Basic Usage
```bash
# Run Figure 1E (ChIP-seq IP vs INPUT profile)
python "Fig 1E.py"

# Run Figure 10 (volcano plots)
python "Fig 10.py"

# Run the CSV-only decile analysis (Figure 3)
python "Fig 3.py"
```

### Customising Paths and Parameters
Edit the **Configuration** section at the top of each script to match your data locations and analysis parameters.

For example, in `Fig 1E.py`:
```python
DB_PATH = "/path/to/your/annotations.db"
IP_BW = "/path/to/your/ip_signal.bw"
INPUT_BW = "/path/to/your/input_signal.bw"
GENOME = "mm10"  # or "hg19"
```

## Quick Reference

| Figure | Script | Pipeline | Input | Output |
|---|---|---|---|---|
| 1A | `Fig 1A.py` | GeneScaleRegionsProfilePipeline | Gene BED + signal files | Meta-gene profile |
| 1B | `Fig 1B.py` | ChIPCategoryCorrelationPipeline | Merged TSV | Hexbin scatter |
| 1C | `Fig 1C.py` | RNAKnockdownCorrelationPipeline | Merged TSV | Hexbin scatter |
| 1D | `Fig 1D.py` | ATACKnockdownCorrelationPipeline | Merged TSV | Hexbin scatter |
| 1E | `Fig 1E.py` | PeakHeatmapPipeline | IP + INPUT signal | Profile |
| 2A/2B | `Fig 2.py` | IntegratedHeatmapPipeline | Signal + RNA CSV | Heatmap + profile |
| 3A/3B | `Fig 3.py` | DecileScatterPipeline | CSV only | Decile boxplot + stripplot |
| 4A/4B | `Fig 4AB.py` | IntegratedHeatmapPipeline | Signal + ATAC CSV | Profile overlay |
| 4C/4D | `Fig 4CD.py` | OverlayProfilePipeline | 2× signal files | Overlay profile |
| 5A/5B | `Fig 5.py` | IntegratedHeatmapPipeline | Signal + RNA CSV, FACT-filtered | Heatmap + profile |
| 6 | `Fig 6.py` | RNAViolinPipeline | RNA CSV only | Violin panels |
| 7A–7D | `Fig 7.py` | DecileScatterPipeline | CSV, FACT-filtered | Decile boxplot + stripplot |
| 8A/8B | `Fig 8.py` | IntegratedHeatmapPipeline | Signal + ATAC CSV, FACT-filtered | Profile overlay |
| 9A–9C | `Fig 9.py` | CategoryHeatmapPipeline | Signal + category CSV | Category heatmaps |
| 10A/10B | `Fig 10.py` | VolcanoPipeline | ChIP CSV | Volcano plot |
| 11/12 | `Fig 11-12.py` | CircPlotPipeline | CSV with 3-signal comparison | Circos plots (-log10 padj) |
| 13/14 | `Fig 13-14.py` | CircPlotPipeline | CSV with 3-signal comparison | Circos plots (log2FC) |
| S1A/S1B | `Fig S1AB.py` | DecileScatterPipeline | CSV only | Decile boxplot + stripplot |
| S1C/S1D | `Fig S1CD.py` | GeneGroupProfilePipeline | Signal + ATAC/RNA log2FC CSV | Quadrant profile |
| S2A–S2D | `Fig S2.py` | OverlayProfilePipeline | 2× signal files, FACT-filtered | Overlay profile |
| Extra | `extra/multisignal_comparison.py` | MultiSignalComparisonPipeline | CSV with 3 LFC columns | 3-way scatter plot |

## For More Information
- See the docstrings in each script for details on the workflow
- Check `Pyntegrate/pipelines.py` for the pipeline class implementations
- Read the Pyntegrate documentation for full API reference
