"""
Core alignment logic using SIMD-accelerated libraries (Parasail and Edlib).
"""

from typing import Dict, Any, Optional
import edlib
import parasail

# Cache for custom substitution matrices (match, mismatch)
_MATRIX_CACHE: Dict[tuple, Any] = {}


def get_matrix(match_score: int, mismatch_penalty: int):
    """Retrieve or create a parasail substitution matrix for DNA sequences (ACGTN)."""
    key = (match_score, mismatch_penalty)
    if key not in _MATRIX_CACHE:
        _MATRIX_CACHE[key] = parasail.matrix_create("ACGTNacgtn", match_score, mismatch_penalty)
    return _MATRIX_CACHE[key]


def _format_cigar(ops) -> str:
    """Format cigar operations into a standard SAM CIGAR string."""
    if not ops:
        return ""
    combined = []
    for op, length in ops:
        # Map '=' (match) and 'X' (mismatch) to 'M' for standard SAM format
        sam_op = "M" if op in ("=", "X") else op
        if combined and combined[-1][0] == sam_op:
            combined[-1] = (sam_op, combined[-1][1] + length)
        else:
            combined.append((sam_op, length))
    return "".join(f"{length}{op}" for op, length in combined)


def align(
    read_seq: Optional[str],
    ref_seq: Optional[str],
    mode: str = "local",
    match_score: int = 2,
    mismatch_penalty: int = -3,
    gap_open: int = -5,
    gap_extend: int = -2,
) -> Dict[str, Any]:
    """
    Perform sequence alignment between read_seq and ref_seq.

    Returns dict matching DuckDB STRUCT:
        score: int
        ref_start: int
        ref_end: int
        read_start: int
        read_end: int
        cigar: str
        identity_pct: float
    """
    if not read_seq or not ref_seq:
        return {
            "score": 0,
            "ref_start": 0,
            "ref_end": 0,
            "read_start": 0,
            "read_end": 0,
            "cigar": "",
            "identity_pct": 0.0,
        }

    read_seq = str(read_seq)
    ref_seq = str(ref_seq)
    mode_lower = (mode or "local").lower().strip()

    open_cost = abs(gap_open)
    extend_cost = abs(gap_extend)
    matrix = get_matrix(match_score, mismatch_penalty)

    # In alignment: query is read_seq, target/reference is ref_seq
    if mode_lower in ("local", "sw", "smith_waterman"):
        align_fn = parasail.sw_trace_striped_32
    elif mode_lower in ("global", "nw", "needleman_wunsch"):
        align_fn = parasail.nw_trace_striped_32
    elif mode_lower in ("semi_global", "semiglobal", "sg"):
        align_fn = parasail.sg_trace_striped_32
    else:
        align_fn = parasail.sw_trace_striped_32

    res = align_fn(read_seq, ref_seq, open_cost, extend_cost, matrix)

    score = int(res.score) if res.score > -1000000000 else 0
    cigar_obj = res.cigar

    if cigar_obj and cigar_obj.len > 0 and score > 0:
        raw_ops = [
            (cigar_obj.decode_op(cigar_obj.seq[i]).decode("ascii"), cigar_obj.decode_len(cigar_obj.seq[i]))
            for i in range(cigar_obj.len)
        ]

        read_start = int(cigar_obj.beg_query)
        read_end = int(res.end_query)
        ref_start = int(cigar_obj.beg_ref)
        ref_end = int(res.end_ref)

        ops = list(raw_ops)
        if mode_lower in ("local", "sw", "smith_waterman"):
            while ops and ops[0][0] in ("I", "D"):
                op, l = ops.pop(0)
                if op == "I":
                    read_start += l
                elif op == "D":
                    ref_start += l
            while ops and ops[-1][0] in ("I", "D"):
                ops.pop()

        cigar_str = _format_cigar(ops)

        # Calculate percentage identity over the aligned region
        matches = sum(l for op, l in ops if op == "=")
        aligned_len = sum(l for op, l in ops)
        identity_pct = float(matches / aligned_len) if aligned_len > 0 else 0.0
    else:
        read_start = 0
        read_end = 0
        ref_start = 0
        ref_end = 0
        cigar_str = ""
        identity_pct = 0.0

    return {
        "score": score,
        "ref_start": ref_start,
        "ref_end": ref_end,
        "read_start": read_start,
        "read_end": read_end,
        "cigar": cigar_str,
        "identity_pct": identity_pct,
    }


def compute_align_score(
    read_seq: Optional[str],
    ref_seq: Optional[str],
    mode: str = "local",
    match_score: int = 2,
    mismatch_penalty: int = -3,
    gap_open: int = -5,
    gap_extend: int = -2,
) -> float:
    """
    Calculate a normalized similarity score (0.0 to 1.0) between a read sequence
    and a reference target sequence.

    Normalized relative to the maximum theoretical score achievable on the target
    reference sequence (len(ref_seq) * match_score).
    """
    if not read_seq or not ref_seq:
        return 0.0

    ref_seq = str(ref_seq)
    max_score = len(ref_seq) * max(1, match_score)
    res = align(
        read_seq,
        ref_seq,
        mode=mode,
        match_score=match_score,
        mismatch_penalty=mismatch_penalty,
        gap_open=gap_open,
        gap_extend=gap_extend,
    )
    if res["score"] <= 0:
        return 0.0
    return min(1.0, max(0.0, float(res["score"]) / float(max_score)))
