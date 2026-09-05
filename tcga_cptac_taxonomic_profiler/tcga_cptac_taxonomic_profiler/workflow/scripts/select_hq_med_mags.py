#!/usr/bin/env python3
"""Select medium- and high-quality Binette bins and stage them as ALL_MAGS.

This is the hand-off between binning and annotation. Nothing previously
connected them: get_hq_bins populated ALL_MAGS from the older per-sample
CheckM2/VAMB path, so with Binette as the consensus binner the annotation
workflow had no input at all.

Tiers follow the usual completeness/contamination convention, which is also what
MIMAG uses for its quantitative criteria:

    high    completeness >  90  and contamination <  5
    medium  completeness >= 50  and contamination < 10

Only the quantitative half of MIMAG is applied here - the rRNA and tRNA criteria
need bakta, which runs downstream, and mimag.py then assembles the full call.
Naming a bin "high quality" on completeness alone would overstate it.

Standard library only: no conda environment, so it cannot break on a re-solve
and needs no network on a compute node.

Usage:
    select_hq_med_mags.py --report final_bins_quality_reports.tsv \
        --bin-dir final_bins --out-dir ALL_MAGS \
        --min-completeness 50 --max-contamination 10
"""

import argparse
import os
import shutil
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", required=True)
    ap.add_argument("--bin-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--min-completeness", type=float, default=50.0)
    ap.add_argument("--max-contamination", type=float, default=10.0)
    a = ap.parse_args()

    with open(a.report) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for col in ("name", "completeness", "contamination"):
            if col not in header:
                sys.exit("ERROR: %s lacks a %r column; header was %s"
                         % (a.report, col, header))
        i_n, i_comp, i_cont = (header.index(c) for c in
                               ("name", "completeness", "contamination"))
        rows = [l.rstrip("\n").split("\t") for l in fh if l.strip()]

    if not rows:
        sys.exit("ERROR: %s has no data rows" % a.report)

    # Binette writes final_bins/bin_<name>.fa; fall back to any matching stem.
    def find_fasta(name):
        for cand in ("bin_%s.fa" % name, "%s.fa" % name, "%s.fna" % name,
                     "bin_%s.fna" % name):
            p = os.path.join(a.bin_dir, cand)
            if os.path.isfile(p):
                return p
        return None

    os.makedirs(a.out_dir, exist_ok=True)
    kept = high = medium = 0
    missing = []
    for r in rows:
        try:
            comp, cont = float(r[i_comp]), float(r[i_cont])
        except (ValueError, IndexError):
            continue
        if not (comp >= a.min_completeness and cont < a.max_contamination):
            continue
        name = r[i_n]
        src = find_fasta(name)
        if src is None:
            missing.append(name)
            continue
        # ALL_MAGS is globbed for *.fna by the annotation workflow.
        shutil.copyfile(src, os.path.join(a.out_dir, "%s.fna" % name))
        kept += 1
        if comp > 90 and cont < 5:
            high += 1
        else:
            medium += 1

    if missing:
        sys.exit("ERROR: %d selected bins have no FASTA in %s (e.g. %s)"
                 % (len(missing), a.bin_dir, ", ".join(missing[:3])))
    if kept == 0:
        sys.exit("ERROR: no bin met completeness >= %g and contamination < %g; "
                 "refusing to write an empty ALL_MAGS"
                 % (a.min_completeness, a.max_contamination))

    sys.stderr.write(
        "staged %d MAGs to %s (%d high: >90/<5, %d medium: >=%g/<%g) from %d bins\n"
        % (kept, a.out_dir, high, medium, a.min_completeness,
           a.max_contamination, len(rows)))


if __name__ == "__main__":
    main()
