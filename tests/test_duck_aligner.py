import gzip
import pytest
import duck_aligner
from duck_aligner.aligner import align, compute_align_score
from duck_aligner.fastq import parse_fastq_records


@pytest.fixture
def sample_fastq_plain(tmp_path):
    fastq_file = tmp_path / "sample.fastq"
    content = (
        "@READ_1\n"
        "ACTGCAACTAGCTAGCTAGCTAGTCGATCGATCGATCGATCGATCGATCG\n"
        "+\n"
        "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
        "@READ_2\n"
        "TGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG\n"
        "+\n"
        "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
        "@READ_3\n"
        "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
        "+\n"
        "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
    )
    fastq_file.write_text(content)
    return str(fastq_file)


@pytest.fixture
def sample_fastq_gz(tmp_path):
    fastq_file = tmp_path / "sample.fastq.gz"
    content = (
        b"@READ_1\n"
        b"ACTGCAACTAGCTAGCTAGCTAGTCGATCGATCGATCGATCGATCGATCG\n"
        b"+\n"
        b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
        b"@READ_2\n"
        b"TGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG\n"
        b"+\n"
        b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
        b"@READ_3\n"
        b"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
        b"+\n"
        b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
    )
    with gzip.open(fastq_file, "wb") as f:
        f.write(content)
    return str(fastq_file)


def test_align_basic():
    res = align("ACTG", "ACTG", mode="global")
    assert res["score"] == 8
    assert res["cigar"] == "4M"
    assert res["identity_pct"] == 1.0
    assert res["ref_start"] == 0
    assert res["ref_end"] == 3
    assert res["read_start"] == 0
    assert res["read_end"] == 3


def test_align_modes():
    # Local alignment
    res_local = align("AAAAACTGCCCC", "ACTG", mode="local")
    assert res_local["score"] == 8
    assert res_local["identity_pct"] == 1.0

    # Semi-global alignment
    res_sg = align("ACTGCA", "ACTG", mode="semi_global")
    assert res_sg["score"] > 0

    # Global alignment with mismatch
    res_mismatch = align("ACTG", "ACCG", mode="global")
    assert res_mismatch["identity_pct"] == 0.75


def test_align_empty():
    res = align("", "ACTG")
    assert res["score"] == 0
    assert res["identity_pct"] == 0.0

    res2 = align("ACTG", None)
    assert res2["score"] == 0
    assert res2["identity_pct"] == 0.0


def test_compute_align_score():
    score = compute_align_score("ACTG", "ACTG")
    assert score == 1.0

    score_mismatch = compute_align_score("ACCGT", "ACTGT", mode="global")
    assert 0.0 < score_mismatch < 1.0

    # Non-matching sequence should have low or 0 score
    score_low = compute_align_score("GGGGGGGGGG", "ACTGCAACTA")
    assert score_low < 0.2


def test_fastq_parser(sample_fastq_plain, sample_fastq_gz):
    records = list(parse_fastq_records(sample_fastq_plain))
    assert len(records) == 3
    assert records[0][0] == "@READ_1"
    assert records[0][3] == len(records[0][1])

    records_gz = list(parse_fastq_records(sample_fastq_gz))
    assert len(records_gz) == 3
    assert records_gz[0][0] == "@READ_1"


def test_duckdb_read_fastq(sample_fastq_plain, sample_fastq_gz):
    con = duck_aligner.connect()

    rows = con.execute(f"SELECT id, sequence, quality, length FROM read_fastq('{sample_fastq_plain}')").fetchall()
    assert len(rows) == 3
    assert rows[0][0] == "@READ_1"
    assert rows[0][3] == 50

    rows_gz = con.execute(f"SELECT id, sequence, quality, length FROM read_fastq('{sample_fastq_gz}')").fetchall()
    assert len(rows_gz) == 3
    assert rows_gz[0][0] == "@READ_1"


def test_duckdb_align_score(sample_fastq_plain):
    con = duck_aligner.connect()
    query = f"""
    SELECT count(*)
    FROM read_fastq('{sample_fastq_plain}')
    WHERE align_score(sequence, 'ACTGCAACTAGCTAGCTAGCTAGTCGATCGATCGATCGATCGATCGATCG') > 0.90;
    """
    res = con.execute(query).fetchone()[0]
    assert res == 1


def test_duckdb_example_a_from_notebook(sample_fastq_plain):
    con = duck_aligner.connect()
    query = f"""
    SELECT
        count(*) AS matching_reads,
        avg((align_result).identity_pct) AS avg_identity
    FROM (
        SELECT
            align_sequence(sequence, 'ACTGCAACTAGCTAGCTAGCTAGTCGATCGATCGATCGATCGATCGATCG', 'local') AS align_result
        FROM read_fastq('{sample_fastq_plain}')
    )
    WHERE (align_result).score > 80;
    """
    res = con.execute(query).fetchall()
    assert len(res) == 1
    count_matches, avg_id = res[0]
    assert count_matches == 1
    assert avg_id == 1.0


def test_duckdb_example_b_from_notebook(sample_fastq_plain):
    con = duck_aligner.connect()
    query = f"""
    SELECT
        sequence_id,
        (align_res).score AS score,
        (align_res).cigar AS cigar
    FROM (
        SELECT
            id AS sequence_id,
            align_sequence(sequence, 'TGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG', 'semi_global') AS align_res
        FROM read_fastq('{sample_fastq_plain}')
    )
    WHERE (align_res).identity_pct >= 0.95
    LIMIT 1;
    """
    res = con.execute(query).fetchall()
    assert len(res) == 1
    seq_id, score, cigar = res[0]
    assert seq_id == "@READ_2"
    assert score > 0
    assert "M" in cigar
