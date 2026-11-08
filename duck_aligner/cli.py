"""
Command line interface for Duck Aligner.
"""

import argparse
import sys
import duck_aligner


def main():
    parser = argparse.ArgumentParser(
        prog="duck-aligner",
        description="Lazy FASTQ sequence alignment and variant detection via DuckDB SQL",
    )
    subparsers = parser.add_subparsers(dest="command", help="Available subcommands")

    # 'query' subcommand
    query_parser = subparsers.add_parser("query", help="Execute an arbitrary SQL query over FASTQ files")
    query_parser.add_argument("sql", help="SQL query to execute (use read_fastq, align_score, align_sequence)")

    # 'scan' subcommand
    scan_parser = subparsers.add_parser("scan", help="Scan a FASTQ file for a target sequence")
    scan_parser.add_argument("fastq", help="Path to FASTQ file (.fastq or .fastq.gz)")
    scan_parser.add_argument("target", help="Target reference sequence")
    scan_parser.add_argument("--mode", default="local", choices=["local", "global", "semi_global"], help="Alignment mode")
    scan_parser.add_argument("--min-identity", type=float, default=0.90, help="Minimum identity percentage (0.0 to 1.0)")
    scan_parser.add_argument("--min-score", type=int, default=0, help="Minimum raw score")
    scan_parser.add_argument("--limit", type=int, default=10, help="Maximum matching reads to return")

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    con = duck_aligner.connect()

    if args.command == "query":
        res = con.execute(args.sql)
        df = res.fetchdf()
        print(df.to_string(index=False))

    elif args.command == "scan":
        sql = f"""
        SELECT
            id,
            (res).score AS score,
            (res).cigar AS cigar,
            (res).identity_pct AS identity_pct,
            (res).read_start AS read_start,
            (res).read_end AS read_end,
            (res).ref_start AS ref_start,
            (res).ref_end AS ref_end
        FROM (
            SELECT
                id,
                align_sequence(sequence, '{args.target}', '{args.mode}') AS res
            FROM read_fastq('{args.fastq}')
        )
        WHERE (res).identity_pct >= {args.min_identity}
          AND (res).score >= {args.min_score}
        LIMIT {args.limit};
        """
        res = con.execute(sql)
        df = res.fetchdf()
        print(df.to_string(index=False))


if __name__ == "__main__":
    main()
