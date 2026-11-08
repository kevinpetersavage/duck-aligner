"""
DuckDB UDF registration for Duck Aligner.
"""

import duckdb
from duck_aligner.aligner import align, compute_align_score
from duck_aligner.fastq import register_fastq_reader


# Define DuckDB STRUCT return type matching specification:
# STRUCT(score INTEGER, ref_start INTEGER, ref_end INTEGER, read_start INTEGER, read_end INTEGER, cigar VARCHAR, identity_pct DOUBLE)
ALIGN_STRUCT_TYPE = duckdb.struct_type({
    "score": "INTEGER",
    "ref_start": "INTEGER",
    "ref_end": "INTEGER",
    "read_start": "INTEGER",
    "read_end": "INTEGER",
    "cigar": "VARCHAR",
    "identity_pct": "DOUBLE",
})


def _align_sequence_py(
    read_seq,
    ref_seq,
    mode="local",
    match_score=2,
    mismatch_penalty=-3,
    gap_open=-5,
    gap_extend=-2,
):
    return align(
        read_seq,
        ref_seq,
        mode=mode,
        match_score=match_score,
        mismatch_penalty=mismatch_penalty,
        gap_open=gap_open,
        gap_extend=gap_extend,
    )


def _align_score_py(
    read_seq,
    ref_seq,
    mode="local",
    match_score=2,
    mismatch_penalty=-3,
    gap_open=-5,
    gap_extend=-2,
):
    return compute_align_score(
        read_seq,
        ref_seq,
        mode=mode,
        match_score=match_score,
        mismatch_penalty=mismatch_penalty,
        gap_open=gap_open,
        gap_extend=gap_extend,
    )


def register_alignment_functions(con) -> None:
    """
    Register align_sequence and align_score UDFs and Macros with DuckDB.
    Supports default argument values through SQL macros wrapping the native Python UDFs.
    """
    # 1. Register underlying python functions with full signatures
    con.create_function(
        "_align_sequence_internal",
        _align_sequence_py,
        parameters=[
            "VARCHAR",
            "VARCHAR",
            "VARCHAR",
            "INTEGER",
            "INTEGER",
            "INTEGER",
            "INTEGER",
        ],
        return_type=ALIGN_STRUCT_TYPE,
    )

    con.create_function(
        "_align_score_internal",
        _align_score_py,
        parameters=[
            "VARCHAR",
            "VARCHAR",
            "VARCHAR",
            "INTEGER",
            "INTEGER",
            "INTEGER",
            "INTEGER",
        ],
        return_type="DOUBLE",
    )

    # 2. Register user-facing macros with default parameters
    con.execute("""
    CREATE OR REPLACE MACRO align_sequence(
        read_seq,
        ref_seq,
        mode := 'local',
        match_score := 2,
        mismatch_penalty := -3,
        gap_open := -5,
        gap_extend := -2
    ) AS
        _align_sequence_internal(read_seq, ref_seq, mode, match_score, mismatch_penalty, gap_open, gap_extend);
    """)

    con.execute("""
    CREATE OR REPLACE MACRO align_score(
        read_seq,
        ref_seq,
        mode := 'local',
        match_score := 2,
        mismatch_penalty := -3,
        gap_open := -5,
        gap_extend := -2
    ) AS
        _align_score_internal(read_seq, ref_seq, mode, match_score, mismatch_penalty, gap_open, gap_extend);
    """)


def register(con=None):
    """
    Register all Duck Aligner functions (read_fastq, align_score, align_sequence)
    into the provided DuckDB connection (or default connection if None).
    """
    if con is None:
        con = duckdb.default_connection()
    register_fastq_reader(con)
    register_alignment_functions(con)
    return con


def connect(database=":memory:", **kwargs):
    """
    Create a new DuckDB connection with Duck Aligner functions registered.
    """
    con = duckdb.connect(database=database, **kwargs)
    register(con)
    return con
