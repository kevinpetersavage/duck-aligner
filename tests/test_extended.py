import pytest
import duck_aligner
from duck_aligner.aligner import align, compute_align_score


def test_align_various_modes_and_penalties():
    read = "ACGTACGT"
    ref = "ACGTTCGT"

    # Global with mismatch
    res_nw = align(read, ref, mode="global", match_score=2, mismatch_penalty=-3)
    assert res_nw["score"] == 7 * 2 - 3  # 11
    assert res_nw["identity_pct"] == 7 / 8
    assert res_nw["read_start"] == 0
    assert res_nw["read_end"] == 7
    assert res_nw["ref_start"] == 0
    assert res_nw["ref_end"] == 7

    # Semi-global with flanking extensions
    flanked_read = "NNN" + read + "NNN"
    res_sg = align(flanked_read, ref, mode="semi_global", match_score=2, mismatch_penalty=-3)
    assert res_sg["score"] > 0
    assert "M" in res_sg["cigar"]

    # Local alignment of embedded target
    long_read = "TTTTTT" + read + "GGGGGG"
    res_sw = align(long_read, read, mode="local")
    assert res_sw["score"] == 16
    assert res_sw["identity_pct"] == 1.0
    assert res_sw["read_start"] == 6
    assert res_sw["read_end"] == 13
    assert res_sw["ref_start"] == 0
    assert res_sw["ref_end"] == 7


def test_align_indels():
    # Insertion in read relative to ref
    read = "ACGTTTACGT"
    ref = "ACGTACGT"
    res = align(read, ref, mode="global", match_score=2, mismatch_penalty=-3, gap_open=-5, gap_extend=-2)
    assert "I" in res["cigar"]
    assert res["score"] > 0

    # Deletion in read relative to ref
    read_del = "ACGTACGT"
    ref_del = "ACGTTTACGT"
    res_del = align(read_del, ref_del, mode="global", match_score=2, mismatch_penalty=-3, gap_open=-5, gap_extend=-2)
    assert "D" in res_del["cigar"]
    assert res_del["score"] > 0


def test_duckdb_default_connection_registration():
    con = duck_aligner.register()
    res = con.execute("SELECT align_score('ACGT', 'ACGT')").fetchone()[0]
    assert res == 1.0


def test_duckdb_early_termination_limit(tmp_path):
    fastq_file = tmp_path / "large.fastq"
    with open(fastq_file, "w") as f:
        for i in range(100):
            seq = "ACGTACGTACGT" if i < 10 else "TGCATGCATGCA"
            f.write(f"@READ_{i}\n{seq}\n+\n{'I'*len(seq)}\n")

    con = duck_aligner.connect()
    # Query with LIMIT 5 should terminate early
    query = f"""
    SELECT id, (res).score, (res).cigar
    FROM (
        SELECT id, align_sequence(sequence, 'ACGTACGTACGT', 'local') AS res
        FROM read_fastq('{fastq_file}')
    )
    WHERE (res).score = 24
    LIMIT 5;
    """
    rows = con.execute(query).fetchall()
    assert len(rows) == 5
    for r in rows:
        assert r[1] == 24
        assert r[2] == "12M"
