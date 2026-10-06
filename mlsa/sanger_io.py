"""AB1 trace parsing, quality trimming, and locus classification for the
lab's routine Sanger sequencing filenames (hsp65 + 16S NTM workup).
"""
from __future__ import annotations

import re
from collections.abc import Container
from pathlib import Path
from typing import Any, cast

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

TNR_RE = re.compile(r"(\d{10})")


# hsp65 reads: the gene name in its many historical spellings ("hsp65",
# "hsp_11", "hsp 65", "65 hsp", typos "hps65"/"hso65", "65kDa") or only the
# TB-11/TB-12(w) primer name ("Tb-11", "TB12w", a standalone "12w"). Broader
# than the immensekansasii Rust classifier's `is_hsp65`, which only knows
# "hsp65"/"65kda" and misses these reads.
_HSP65_RE = re.compile(
    r"hsp|hps65|hso65|65\s*kda|tb[-_ ]?1[12]w?|(?<![a-z0-9])12w(?![a-z0-9])",
    re.IGNORECASE,
)
# 16S reads: the lab's 16S primers MBAK-14 (as in the Rust classifier's
# `is_16s`), Mbak259r and Mbak264r. The latter two also appear as just
# "bak259r", "259r", "r259", "264r" or "264" between separators.
_16S_RE = re.compile(
    r"mbak[-_ ]?14(?!\d)|m?bak[-_ ]?(?:259|264)(?!\d)"
    r"|(?<![a-z0-9])r?(?:259|264)r?(?![a-z0-9])",
    re.IGNORECASE,
)


def is_hsp65(filename: str) -> bool:
    """Callers only pass .ab1 filenames (see scripts/link_sanger_kansasii.sh)."""
    return bool(_HSP65_RE.search(filename))


def is_16s(filename: str) -> bool:
    return bool(_16S_RE.search(filename))


def classify_locus(filename: str) -> str | None:
    """Return "hsp65", "16S", or None for a Sanger filename. hsp65 is checked
    first since a small number of filenames coincidentally contain both
    tokens (e.g. batch/run identifiers); hsp65 is the more specific match."""
    if is_hsp65(filename):
        return "hsp65"
    if is_16s(filename):
        return "16S"
    return None


def extract_tnr(filename: str, known_tnrs: Container[str]) -> str | None:
    """Return the TNR embedded in filename, validated against the known TNRs
    from screening_map_link.csv (guards against any other coincidental
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


def mixed_peak_fraction(record: SeqRecord, start: int, end: int, min_ratio: float = 0.25) -> float | None:
    """Fraction of base calls in [start, end) whose second-highest trace channel
    reaches min_ratio of the highest one at the base-call position (PLOC2). A
    second template in the PCR (mixed culture, cross-contamination) shows as
    such secondary peaks, which the instrument writes as Y/K/M/R/S/W or N. None
    if the trace lacks the channels or call positions."""
    raw = cast("dict[str, Any]", record.annotations.get("abif_raw", {}))
    order: str | bytes | None = raw.get("FWO_1")
    locs = raw.get("PLOC2")
    if order is None or locs is None:
        return None
    if isinstance(order, bytes):
        order = order.decode()
    try:
        channels = {base: raw[f"DATA{i}"] for base, i in zip(order, (9, 10, 11, 12))}
    except KeyError:
        return None
    n_mixed = 0
    for i in range(start, end):
        pos = min(locs[i], len(channels["A"]) - 1)
        top, second = sorted((channels[b][pos] for b in "ACGT"), reverse=True)[:2]
        n_mixed += top > 0 and second / top >= min_ratio
    return n_mixed / (end - start)


def load_trimmed_ab1_mixed(path: Path, min_length: int = 100) -> tuple[str, float, float | None] | None:
    """Like load_trimmed_ab1 but returns (trimmed_sequence, mean_quality,
    mixed_fraction); mixed_fraction is mixed_peak_fraction over the trimmed
    region (None when the trace has no peak data). None under the same
    conditions as load_trimmed_ab1."""
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
    return seq, mean_q, mixed_peak_fraction(record, start, end)


def load_trimmed_ab1(path: Path, min_length: int = 100) -> tuple[str, float] | None:
    """Parse an .ab1 trace and quality-trim it. Returns (trimmed_sequence,
    mean_post_trim_phred_quality), or None if the trimmed region is shorter
    than min_length (unusable read), the trace has no quality track, or the
    file can't be parsed at all (some of these traces are 20 years old and a
    handful are truncated/corrupted; skip rather than abort the whole run)."""
    loaded = load_trimmed_ab1_mixed(path, min_length)
    return None if loaded is None else loaded[:2]
