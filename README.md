# duck-aligner

Lazy FASTQ sequence alignment and targeted variant detection powered by DuckDB and SIMD-accelerated alignment libraries.

## Overview

Traditional genomic alignment pipelines (e.g., BWA-MEM / Bowtie2 followed by SAMtools) require aligning entire sequencing datasets before querying specific loci or variants. **duck-aligner** enables targeted, on-demand sequence alignment directly via SQL queries in DuckDB:

- **On-the-fly Ingestion:** Stream and parse raw FASTQ files (`.fastq` and `.fastq.gz`) directly into relational tables via `read_fastq(filepath)`.
- **Targeted Alignment UDFs:** Fast SIMD-accelerated Smith-Waterman (`local`), Needleman-Wunsch (`global`), and semi-global (`semi_global`) alignment kernels via `align_sequence` and `align_score`.
- **Push-Down Filtering & Early Termination:** Exploit SQL `WHERE` clauses and `LIMIT` early exit to stop reading files as soon as target criteria or match counts are satisfied.

---

## SQL Schema & Functions

### 1. `read_fastq(path VARCHAR)`
Table macro that reads raw `.fastq` or compressed `.fastq.gz` files and projects the following schema:
- `id VARCHAR`: Read header identifier line (e.g. `@READ_1`).
- `sequence VARCHAR`: Raw nucleotide sequence (A, C, G, T, N).
- `quality VARCHAR`: Phred quality scores.
- `length INTEGER`: Read sequence length.

### 2. `align_score(read_seq, ref_seq, mode := 'local', match_score := 2, mismatch_penalty := -3, gap_open := -5, gap_extend := -2) -> DOUBLE`
Calculates a normalized similarity score (0.0 to 1.0) between a read sequence and a target reference sequence.

### 3. `align_sequence(read_seq, ref_seq, mode := 'local', match_score := 2, mismatch_penalty := -3, gap_open := -5, gap_extend := -2) -> STRUCT`
Performs alignment and returns a structured record:
- `score INTEGER`: Raw alignment score.
- `ref_start INTEGER`: 0-based start coordinate on the reference.
- `ref_end INTEGER`: 0-based end coordinate on the reference.
- `read_start INTEGER`: 0-based start coordinate on the read.
- `read_end INTEGER`: 0-based end coordinate on the read.
- `cigar VARCHAR`: Standard SAM-format CIGAR string (e.g., `50M`, `10M1I39M`).
- `identity_pct DOUBLE`: Percentage identity over the aligned region.

Supported modes: `'local'`, `'global'`, and `'semi_global'`.

---

## Python API Usage

```python
import duck_aligner

# Connect to DuckDB with duck_aligner functions automatically registered
con = duck_aligner.connect()

# Rapid Variant Verification
query = """
SELECT
    count(*) AS matching_reads,
    avg((align_result).identity_pct) AS avg_identity
FROM (
    SELECT
        align_sequence(sequence, 'ACTGCAACTAGCTAGCTAGCTAGTCGATCGATCGATCGATCGATCGATCG', 'local') AS align_result
    FROM read_fastq('sample.fastq.gz')
)
WHERE (align_result).score > 80;
"""
print(con.execute(query).fetchall())

# Early-Termination Screening with LIMIT
query_lazy = """
SELECT
    id AS sequence_id,
    (align_res).score AS score,
    (align_res).cigar AS cigar
FROM (
    SELECT
        id,
        align_sequence(sequence, 'TGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG', 'semi_global') AS align_res
    FROM read_fastq('sample.fastq')
)
WHERE (align_res).identity_pct >= 0.95
LIMIT 10;
"""
print(con.execute(query_lazy).fetchall())
```

---

## CLI Usage

```bash
# Scan a FASTQ file for a target sequence with a match threshold
duck-aligner scan sample.fastq.gz ACTGCAACTAGCTAGCTAGCTAGTCGATCGATCGATCGATCGATCGATCG --min-identity 0.90 --limit 10

# Execute custom SQL queries
duck-aligner query "SELECT count(*) FROM read_fastq('sample.fastq') WHERE align_score(sequence, 'ACTGCA') > 0.8"
```

---

## Running Tests

```bash
poetry run pytest -v
```
