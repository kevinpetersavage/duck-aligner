"""
FASTQ reader utilities and DuckDB integration for lazy sequence ingestion.
"""

import gzip
import os
from typing import Generator, Tuple


def is_gzip_file(filepath: str) -> bool:
    """Check if file is gzip/bgzip compressed by examining magic bytes."""
    if not os.path.exists(filepath):
        return False
    try:
        with open(filepath, "rb") as f:
            magic = f.read(2)
            return magic == b"\x1f\x8b"
    except Exception:
        return False


def parse_fastq_records(filepath: str) -> Generator[Tuple[str, str, str, int], None, None]:
    """
    Parse a FASTQ file (plain or gzip-compressed) yielding:
    (id, sequence, quality, length)
    """
    open_fn = gzip.open if (filepath.endswith(".gz") or is_gzip_file(filepath)) else open
    mode = "rt" if (filepath.endswith(".gz") or is_gzip_file(filepath)) else "r"

    with open_fn(filepath, mode, encoding="utf-8", errors="replace") as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline()
            plus = f.readline()
            qual = f.readline()

            if not seq or not plus or not qual:
                break

            h_str = header.strip()
            s_str = seq.strip()
            q_str = qual.strip()

            yield (h_str, s_str, q_str, len(s_str))


# SQL Macro to register read_fastq directly into DuckDB
# Note: DuckDB's read_csv natively decompresses .gz files, reads line by line,
# and by using aggregation over line groups (rn // 4), avoids materializing
# intermediate tables before query optimization.
READ_FASTQ_MACRO_SQL = """
CREATE OR REPLACE MACRO read_fastq(path) AS TABLE
    WITH raw_lines AS (
        SELECT line, row_number() OVER () - 1 AS rn
        FROM read_csv(path, header=false, columns={'line': 'VARCHAR'}, quote='', escape='', delim=NULL, auto_detect=false)
    )
    SELECT
        MAX(CASE WHEN rn % 4 = 0 THEN line END) AS id,
        MAX(CASE WHEN rn % 4 = 1 THEN line END) AS sequence,
        MAX(CASE WHEN rn % 4 = 3 THEN line END) AS quality,
        CAST(length(MAX(CASE WHEN rn % 4 = 1 THEN line END)) AS INTEGER) AS length
    FROM raw_lines
    GROUP BY rn // 4
    ORDER BY rn // 4;
"""


def register_fastq_reader(con) -> None:
    """Register the read_fastq table macro into a DuckDB connection."""
    con.execute(READ_FASTQ_MACRO_SQL)
