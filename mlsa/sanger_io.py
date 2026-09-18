"""AB1 trace parsing, quality trimming, and locus classification for the
lab's routine Sanger sequencing filenames (hsp65 + 16S NTM workup).
"""
from __future__ import annotations

import re
from pathlib import Path

from Bio import SeqIO

TNR_RE = re.compile(r"(\d{10})")


def is_hsp65(filename: str) -> bool:
    """Ported from the immensekansasii Rust classifier's `is_hsp65`
    (`!is_fasta() && (contains("hsp65") || contains("65kda"))`). Callers only
    ever pass already-.ab1-filtered filenames (see
    scripts/link_sanger_kansasii.sh), so the `!is_fasta()` guard is always
    true here and is omitted."""
    lower = filename.lower()
    return "hsp65" in lower or "65kda" in lower


def is_16s(filename: str) -> bool:
    """Ported from the Rust classifier's `is_16s`
    (`contains("mbak14") || contains("mbak-14")`), restricted to .ab1
    filenames as for is_hsp65 above."""
    lower = filename.lower()
    return "mbak14" in lower or "mbak-14" in lower


def classify_locus(filename: str) -> str | None:
    """Return "hsp65", "16S", or None for a Sanger filename. hsp65 is checked
    first since a small number of filenames coincidentally contain both
    tokens (e.g. batch/run identifiers); hsp65 is the more specific match."""
    if is_hsp65(filename):
        return "hsp65"
    if is_16s(filename):
        return "16S"
    return None


def extract_tnr(filename: str, known_tnrs: set[str]) -> str | None:
    """Return the TNR embedded in filename, validated against the known TNR
    set from screening_map.csv (guards against any other coincidental
    10-digit run, e.g. part of a timestamp)."""
    for m in TNR_RE.finditer(filename):
        if m.group(1) in known_tnrs:
            return m.group(1)
    return None


def mott_trim(qualities: list[int], error_threshold: float = 0.05) -> tuple[int, int]:
    """Modified-Mott quality trimming: find the contiguous window maximizing
    sum(error_threshold - 10**(-q/10)) over the read (Kadane's algorithm, the
    standard formulation used by e.g. Staden/BioPerl trimming). Returns
    (start, end) half-open indices; (0, 0) if no usable window exists."""
    scores = [error_threshold - 10 ** (-q / 10) for q in qualities]
    best_sum = 0.0
    cur_sum = 0.0
    best_start = best_end = cur_start = 0
    for i, val in enumerate(scores):
        if cur_sum <= 0:
            cur_start = i
            cur_sum = val
        else:
            cur_sum += val
        if cur_sum > best_sum:
            best_sum = cur_sum
            best_start = cur_start
            best_end = i + 1
    return best_start, best_end


def load_trimmed_ab1(path: Path, min_length: int = 100) -> tuple[str, float] | None:
    """Parse an .ab1 trace and quality-trim it. Returns (trimmed_sequence,
    mean_post_trim_phred_quality), or None if the trimmed region is shorter
    than min_length (unusable read), the trace has no quality track, or the
    file can't be parsed at all (some of these traces are 20 years old and a
    handful are truncated/corrupted; skip rather than abort the whole run)."""
    try:
        record = SeqIO.read(path, "abi")
    except Exception:
        return None
    qualities = record.letter_annotations.get("phred_quality")
    if not qualities:
        return None
    start, end = mott_trim(qualities)
    if end - start < min_length:
        return None
    seq = str(record.seq[start:end])
    mean_q = sum(qualities[start:end]) / (end - start)
    return seq, mean_q
