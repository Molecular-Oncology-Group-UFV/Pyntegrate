"""
Fig 7: CSV-only decile analysis (unified example), FACT-filtered.

This example demonstrates the DecileScatterPipeline which works entirely
from CSV data without extracting signal arrays from BigWig files. It builds
four panels, one per subunit x FACT-class combination:
- 7A/7B: SSRP1 occupancy deciles, FACT-high / FACT-low
- 7C/7D: SUPT16H occupancy deciles, FACT-high / FACT-low

Each panel bins promoters into deciles of that panel's own subunit occupancy
(`CHIP_COL_AB` / `CHIP_COL_CD`) and overlays the matching subunit's RNA
log2 fold-change on the right axis.

Optionally filters data by category (e.g., FACT category).

Workflow:
1. Load CSV with ChIP and RNA columns
2. Optionally filter by category
3. Sort by selected column (log2FC or baseMean) descending
4. Create decile groups
5. Generate dual-axis boxplot figures

"""

import os

import matplotlib.pyplot as plt
from Pyntegrate.pipelines import DecileScatterPipeline

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")

# ============================================================================
# Configuration (adjust these paths to your data)
# ============================================================================
CSV_PATH = os.path.join(DATA_DIR, "final_table_filtered.tsv")

CHIP_YLABEL = "ChIP (log2 fold change)"
N_DECILES = 10

# The occupancy column must follow the subunit of each panel pair. A single
# global CHIP_COL was previously passed to all four pipelines, so panels C and
# D binned and plotted SSRP1 occupancy under SUPT16H titles while their RNA
# overlay correctly switched to SUPT16H (fixed 2026-08-20).
CHIP_COL_AB = "ssrp1 ChIP log2 fold change"
CHIP_COL_CD = "suptt16h ChIP log2 fold change"   # header spelling as in the input table

##############################################
################# FIGURE 7AB #################
##############################################
RNA_COL_AB = "ssrp1 RNA log2 fold change"
RNA_YLABEL_AB = "RNA (log2 fold change)"
CATEGORY_COL_AB = "SSRP1 FACT category"

### 7A ###
OUTPUT_PATH_A = "examples/outputs/figs/fig7A.svg" # Output file (set to None to skip saving)
TITLE_A = "SSRP1 FACT High"
CATEGORIES_A = "high"  # Filter to only show "high" category for SSRP1

### 7B ###
OUTPUT_PATH_B = "examples/outputs/figs/fig7B.svg" # Output file (set to None to skip saving)
TITLE_B = "SSRP1 FACT Low"
CATEGORIES_B = "low"  # Filter to only show "low" category for SSRP1

##############################################
################# FIGURE 7CD #################
##############################################

RNA_COL_CD = "supt16h RNA log2 fold change"
RNA_YLABEL_CD = "RNA (log2 fold change)"
CATEGORY_COL_CD = "SUPT16H FACT category"

### 7C ###
OUTPUT_PATH_C = "examples/outputs/figs/fig7C.svg" # Output file (set to None to skip saving)
TITLE_C = "SUPT16H FACT High"
CATEGORIES_C = "high"  # Filter to only show "high" category for SUPT16H

### 7D ###
OUTPUT_PATH_D = "examples/outputs/figs/fig7D.svg" # Output file (set to None to skip saving)
TITLE_D = "SUPT16H FACT Low"
CATEGORIES_D = "low"  # Filter to only show "low" category for SUPT16H

# ============================================================================
# Run Pipeline
# ============================================================================
if __name__ == "__main__":
    print(f"Starting CSV-based decile analysis...")
    
    pipeA = DecileScatterPipeline(
        csv_path=CSV_PATH,
        chip_col=CHIP_COL_AB,
        rna_col=RNA_COL_AB,
        n_deciles=N_DECILES,
        chip_ylabel=CHIP_YLABEL,
        rna_ylabel=RNA_YLABEL_AB,
        label_fontsize=25,
        title= TITLE_A,
        category_col=CATEGORY_COL_AB,
        categories=CATEGORIES_A,
    )

    fig = pipeA.run()

    if OUTPUT_PATH_A:
        os.makedirs(os.path.dirname(OUTPUT_PATH_A), exist_ok=True)
        fig.savefig(OUTPUT_PATH_A, dpi=300)
        print(f"  - Figure A saved to: {OUTPUT_PATH_A}")

    pipeB = DecileScatterPipeline(
        csv_path=CSV_PATH,
        chip_col=CHIP_COL_AB,
        rna_col=RNA_COL_AB,
        n_deciles=N_DECILES,
        chip_ylabel=CHIP_YLABEL,
        rna_ylabel=RNA_YLABEL_AB,
        label_fontsize=25,
        title= TITLE_B,
        category_col=CATEGORY_COL_AB,
        categories=CATEGORIES_B,
    )

    fig = pipeB.run()

    if OUTPUT_PATH_B:
        os.makedirs(os.path.dirname(OUTPUT_PATH_B), exist_ok=True)
        fig.savefig(OUTPUT_PATH_B, dpi=300)
        print(f"  - Figure B saved to: {OUTPUT_PATH_B}")

    pipeC = DecileScatterPipeline(
        csv_path=CSV_PATH,
        chip_col=CHIP_COL_CD,
        rna_col=RNA_COL_CD,
        n_deciles=N_DECILES,
        chip_ylabel=CHIP_YLABEL,
        rna_ylabel=RNA_YLABEL_CD,
        label_fontsize=25,
        title= TITLE_C,
        category_col=CATEGORY_COL_CD,
        categories=CATEGORIES_C,
    )
    
    fig = pipeC.run()
    
    if OUTPUT_PATH_C:
        os.makedirs(os.path.dirname(OUTPUT_PATH_C), exist_ok=True)
        fig.savefig(OUTPUT_PATH_C, dpi=300)
        print(f"  - Figure C saved to: {OUTPUT_PATH_C}")

    pipeD = DecileScatterPipeline(
        csv_path=CSV_PATH,
        chip_col=CHIP_COL_CD,
        rna_col=RNA_COL_CD,
        n_deciles=N_DECILES,
        chip_ylabel=CHIP_YLABEL,
        rna_ylabel=RNA_YLABEL_CD,
        label_fontsize=25,
        title= TITLE_D,
        category_col=CATEGORY_COL_CD,
        categories=CATEGORIES_D,
    )

    fig = pipeD.run()
    
    if OUTPUT_PATH_D:
        os.makedirs(os.path.dirname(OUTPUT_PATH_D), exist_ok=True)
        fig.savefig(OUTPUT_PATH_D, dpi=300)
        print(f"  - Figure D saved to: {OUTPUT_PATH_D}")

    plt.show()

