from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent

INPUT_FILE = BASE_DIR / "data" / "final_table.tsv"
OUTPUT_FILE = BASE_DIR / "data" / "final_table_filtered_test.tsv"

REQUIRED_COLUMNS = [
    "EnsemblGeneID",
    "ssrp1 ChIP log2 fold change",
    "suptt16h ChIP log2 fold change",
    "ssrp1 RNA log2 fold change",
    "supt16h RNA log2 fold change",
    "ssrp1 ATAC log2 Fold Change",
    "supt16h ATAC log2 Fold Change",
    "Mean_TPM_Control", "Mean_TPM_SSRP1", "Mean_TPM_SUPT16H"]

# Values considered "empty". Add "NA", "NaN", "-", etc. if your table
# encodes missing values that way and you also want them removed.
NULL_VALUES = {""}

# If the file raises a UnicodeDecodeError, switch to "latin-1".
ENCODING = "utf-8"

# ---------------------------------------------------------------------------


def main() -> None:
    with INPUT_FILE.open("r", encoding=ENCODING, newline="") as fin, \
         OUTPUT_FILE.open("w", encoding=ENCODING, newline="") as fout:

        header = fin.readline()
        if not header:
            raise ValueError("The input file is empty.")

        names = header.rstrip("\r\n").split("\t")
        missing = [c for c in REQUIRED_COLUMNS if c not in names]
        if missing:
            raise ValueError(f"Columns not found in header: {missing}")
        indices = [names.index(c) for c in REQUIRED_COLUMNS]
        gene_idx = names.index("EnsemblGeneID")

        fout.write(header)  # the header is copied untouched

        seen = set()
        read_count = kept_count = duplicate_count = 0
        for line in fin:
            if not line.strip():           # ignore stray blank lines
                continue
            read_count += 1

            fields = line.rstrip("\r\n").split("\t")
            is_empty = any(
                i >= len(fields) or fields[i].strip() in NULL_VALUES
                for i in indices
            )
            if is_empty:
                continue

            gene = fields[gene_idx].strip() if gene_idx < len(fields) else ""
            if gene in seen:                # a row with this ID was already kept
                duplicate_count += 1
                continue
            seen.add(gene)

            fout.write(line)                # write the original line as-is
            kept_count += 1

    print(f"Data rows read:    {read_count}")
    print(f"Rows kept:         {kept_count}")
    print(f"Rows removed (empty values):      {read_count - kept_count - duplicate_count}")
    print(f"Rows removed (duplicate EnsemblGeneID): {duplicate_count}")
    print(f"File saved to:     {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
