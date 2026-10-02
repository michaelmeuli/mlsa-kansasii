#!/usr/bin/env python
"""Main program 3: species identification of Sanger reads by plain alignment to
one reference genome per species (the GTDB representatives of the 7 species in
data/gtdb_genomes/Mycobacteriaceae/kansasii_complex_gtdb_representatives).

No locus extraction and no MLSA: every read is aligned to every genome, and the
reference with the fewest differences over the read is the closest species. See
mlsa/refalign.py for the method and README.md for when this is (not) enough.

Usage:
  python main3_reference_alignment/run_reference_alignment.py READ [READ ...]
  READ = .ab1 trace (quality-trimmed here), .fasta/.fa read (used as is), or a directory of those.

Writes to OUTPUT (default output/mlsa/main3_reference_alignment/):
  reference_alignment.tsv   one row per read
  <read>_ref_alignment.pdf  difference matrix + alignment of the read and all 7 references
About 3 s per read on one core, so fine on the login node for a handful of reads; use a
job for hundreds.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

import pandas as pd  # noqa: E402
from Bio import SeqIO  # noqa: E402

from mlsa import DATA_ROOT, GTDB_REPRESENTATIVES, SPECIES  # noqa: E402
from mlsa.refalign import identify_read, load_references, write_pdf  # noqa: E402
from mlsa.sanger_io import classify_locus, load_trimmed_ab1  # noqa: E402

OUTPUT = DATA_ROOT.parent / "output" / "mlsa" / "main3_reference_alignment"


def iter_reads(paths: list[Path]):
    """Yield (name, locus, sequence); .ab1 reads are Mott-trimmed, reads < 100 bp are skipped."""
    files = []
    for p in paths:
        files += sorted(f for f in p.rglob("*") if f.suffix in {".ab1", ".fasta", ".fa"}) if p.is_dir() else [p]
    for f in files:
        locus = classify_locus(f.name)
        if f.suffix == ".ab1":
            r = load_trimmed_ab1(f)
            if r is None:
                print(f"skipped (unparsable, no qualities or < 100 bp after trimming): {f.name}")
                continue
            yield f.stem, locus, r[0].upper()
        else:
            for rec in SeqIO.parse(f, "fasta"):
                yield rec.id, locus, str(rec.seq).upper()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reads", nargs="+", type=Path)
    ap.add_argument("--ref-dir", type=Path, default=GTDB_REPRESENTATIVES)
    ap.add_argument("--out", type=Path, default=OUTPUT)
    ap.add_argument("--min-identity", type=float, default=0.99,
                    help="closest reference must be at least this identical to call (default 0.99)")
    ap.add_argument("--min-margin", type=int, default=2,
                    help="closest reference must have at least this many fewer differences than the next (default 2)")
    ap.add_argument("--no-pdf", action="store_true")
    args = ap.parse_args()

    refs = load_references(args.ref_dir, SPECIES)
    print(f"{len(refs)} references: " + ", ".join(f"{r.species} ({r.accession})" for r in refs))
    args.out.mkdir(parents=True, exist_ok=True)
    rows = []
    for name, locus, seq in iter_reads(args.reads):
        res = identify_read(refs, name, locus, seq, args.min_identity, args.min_margin)
        row = {"read": name, "locus": locus, "read_bp": len(seq), "status": res.status,
               "closest_species": res.closest_species}
        if res.hits:
            b, s = res.best, res.runner_up
            row.update(closest_accession=res.hits[b].ref.accession, closest_identity=round(100 * res.identity(b), 2),
                       closest_diffs=res.diffs[b], compared_bp=res.compared[b])
            if s is not None:
                row.update(second_species=res.hits[s].ref.species, second_identity=round(100 * res.identity(s), 2),
                           second_diffs=res.diffs[s], margin=res.diffs[s] - res.diffs[b])
            for i, h in enumerate(res.hits):
                row[f"identity_{h.ref.species}"] = round(100 * res.identity(i), 2)
            if not args.no_pdf:
                row["pdf"] = write_pdf(res, args.out / f"{name.replace(' ', '_')}_ref_alignment.pdf").name
        print(f"{name}: {res.status}, closest {res.closest_species}"
              + (f" {row['closest_identity']}% ({row['closest_diffs']} diffs)" if res.hits else ""), flush=True)
        rows.append(row)
    if rows:
        pd.DataFrame(rows).to_csv(args.out / "reference_alignment.tsv", sep="\t", index=False)
        print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
