#!/usr/bin/env python3
"""Run duck-aligner against reads from a public 1000 Genomes CRAM region."""

import argparse
import gzip
import re
import shutil
import tempfile
from pathlib import Path
from urllib.request import urlopen

import duck_aligner


CRAM_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    "1000_genomes_project/data/CHS/HG00419/high_cov_alignment/"
    "HG00419.alt_bwamem_GRCh38DH.20150917.CHS.high_coverage.cram"
)
DEFAULT_REGION = "chr20:10000000-10002000"


def parse_region(value):
    try:
        contig, coordinates = value.rsplit(":", 1)
        start, end = (int(part.replace(",", "")) for part in coordinates.split("-", 1))
    except (ValueError, TypeError):
        raise argparse.ArgumentTypeError("region must be CONTIG:START-END (1-based, inclusive)")

    if not re.fullmatch(r"[A-Za-z0-9_.-]+", contig) or start < 1 or end < start:
        raise argparse.ArgumentTypeError("region must be CONTIG:START-END (1-based, inclusive)")

    return contig, start - 1, end


def positive_int(value):
    try:
        result = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("value must be a positive integer")
    if result < 1:
        raise argparse.ArgumentTypeError("value must be a positive integer")
    return result


def download_reference(contig, directory, pysam):
    url = f"https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/{contig}.fa.gz"
    compressed_path = directory / f"{contig}.fa.gz"
    reference_path = directory / f"{contig}.fa"

    with urlopen(url, timeout=60) as response, compressed_path.open("wb") as output:
        shutil.copyfileobj(response, output)
    with gzip.open(compressed_path, "rb") as source, reference_path.open("wb") as output:
        shutil.copyfileobj(source, output)

    pysam.faidx(str(reference_path))
    return reference_path


def write_fastq(
    cram_path,
    index_path,
    reference_path,
    contig,
    start,
    end,
    output_path,
    max_reads,
    pysam,
):
    written = 0
    with pysam.AlignmentFile(
        cram_path,
        "rc",
        index_filename=index_path,
        reference_filename=str(reference_path),
        require_index=True,
    ) as cram, output_path.open("w", encoding="utf-8") as fastq:
        if cram.get_tid(contig) < 0:
            raise ValueError(f"CRAM does not contain reference sequence {contig!r}")

        for read in cram.fetch(contig, start, end):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            sequence = read.query_sequence
            if not sequence:
                continue
            quality = read.qual or ("I" * len(sequence))
            fastq.write(f"@{read.query_name}\n{sequence}\n+\n{quality}\n")
            written += 1
            if written >= max_reads:
                break

    return written


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cram", default=CRAM_URL, help="HTTP(S) URL or local path to an indexed CRAM")
    parser.add_argument(
        "--crai",
        help="HTTP(S) URL or local path to its CRAI index (defaults to CRAM path plus .crai)",
    )
    parser.add_argument(
        "--region",
        type=parse_region,
        default=parse_region(DEFAULT_REGION),
        metavar="CONTIG:START-END",
        help=f"1-based inclusive region (default: {DEFAULT_REGION})",
    )
    parser.add_argument("--max-reads", type=positive_int, default=100)
    parser.add_argument("--min-identity", type=float, default=0.90)
    parser.add_argument("--limit", type=positive_int, default=20)
    parser.add_argument("--reference", type=Path, help="Indexed GRCh38 FASTA; downloads the region's chromosome if omitted")
    args = parser.parse_args()

    if not 0 <= args.min_identity <= 1:
        parser.error("--min-identity must be between 0 and 1")

    try:
        import pysam
    except ImportError:
        parser.error("CRAM support requires pysam; install it with `poetry install --extras cram`")

    contig, start, end = args.region
    index_path = args.crai or f"{args.cram}.crai"

    with tempfile.TemporaryDirectory(prefix="duck-aligner-1000g-") as temp_dir:
        temp_path = Path(temp_dir)
        reference_path = args.reference or download_reference(contig, temp_path, pysam)
        if not Path(f"{reference_path}.fai").exists():
            pysam.faidx(str(reference_path))

        with pysam.FastaFile(str(reference_path)) as reference:
            target_sequence = reference.fetch(contig, start, end).upper()
        if not target_sequence:
            parser.error(f"reference {reference_path} has no sequence for {contig}:{start + 1}-{end}")

        fastq_path = temp_path / "reads.fastq"
        read_count = write_fastq(
            args.cram,
            index_path,
            reference_path,
            contig,
            start,
            end,
            fastq_path,
            args.max_reads,
            pysam,
        )
        if read_count == 0:
            parser.error(f"no reads found in {contig}:{start + 1}-{end}")

        query = """
            SELECT id, (result).score, (result).identity_pct, (result).cigar
            FROM (
                SELECT id, align_sequence(sequence, ?, 'local') AS result
                FROM read_fastq(?)
            )
            WHERE (result).identity_pct >= ?
            ORDER BY (result).score DESC
            LIMIT ?
        """
        rows = duck_aligner.connect().execute(
            query,
            [target_sequence, str(fastq_path), args.min_identity, args.limit],
        ).fetchall()

    sample = "HG00419 (1000 Genomes)" if args.cram == CRAM_URL else args.cram.rsplit("/", 1)[-1]
    print(f"Sample: {sample}")
    print(f"Region: {contig}:{start + 1}-{end}")
    print(f"Reads aligned: {read_count}; matches >= {args.min_identity:.0%}: {len(rows)}")
    print("read_id\tscore\tidentity\tcigar")
    for read_id, score, identity, cigar in rows:
        print(f"{read_id}\t{score}\t{identity:.1%}\t{cigar}")


if __name__ == "__main__":
    main()
