#!/usr/bin/env python3
"""Aggregate per-sample Bracken reports into one taxon x sample count matrix.

Replaces the kraken-biom rule, which produced a valid but completely empty table
({"shape": [0, 0], "data": [], "rows": [], "columns": []}) and let the pipeline
record the stage as complete. Two faults, either of which was fatal:

  shell: kraken-biom {input[0]} -o {output[0]}

  1. {input[0]} is the FIRST file of the expanded list, not the list. The rule
     built the space-joined file list in params and then never referenced it, so
     kraken-biom saw exactly one sample.
  2. kraken-biom parses Kraken2 *report* format. Bracken's tabular output is a
     different shape, so it matched nothing and emitted an empty table rather
     than failing.

Written against the standard library only, deliberately: it needs no conda
environment, so it cannot be broken by a re-solve and needs no network on a
compute node.

Input  (Bracken tabular, tab-separated):
    name  taxonomy_id  taxonomy_lvl  kraken_assigned_reads  added_reads
    new_est_reads  fraction_total_reads

Output (tab-separated, taxa as rows, samples as columns):
    taxonomy_id  taxon  <sample_1>  <sample_2>  ...

Counts are new_est_reads, i.e. Bracken's re-estimated abundance.

Usage:
    bracken_to_matrix.py --level species --out matrix.tsv FILE [FILE ...]
"""

import argparse
import os
import sys


def sample_name(path, level):
    """Strip the directory and the .kraken_bracken_<level>.txt suffix."""
    base = os.path.basename(path)
    suffix = ".kraken_bracken_%s.txt" % level
    return base[: -len(suffix)] if base.endswith(suffix) else base


def read_one(path):
    """Return {taxon: (taxonomy_id, reads)} for one Bracken report."""
    counts = {}
    with open(path) as fh:
        header = fh.readline()
        if not header:
            raise ValueError("empty file: %s" % path)
        cols = header.rstrip("\n").split("\t")
        try:
            i_name = cols.index("name")
            i_taxid = cols.index("taxonomy_id")
            i_reads = cols.index("new_est_reads")
        except ValueError:
            raise ValueError(
                "%s does not look like Bracken output; header was: %s"
                % (path, cols)
            )
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(i_name, i_taxid, i_reads):
                raise ValueError("short line in %s: %r" % (path, line[:120]))
            counts[f[i_name]] = (f[i_taxid], int(float(f[i_reads])))
    return counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", required=True, choices=("species", "genus"))
    ap.add_argument("--out", required=True)
    ap.add_argument("files", nargs="+")
    a = ap.parse_args()

    if not a.files:
        sys.exit("ERROR: no input files given")

    per_sample, taxid_of, order = {}, {}, []
    for path in a.files:
        s = sample_name(path, a.level)
        if s in per_sample:
            sys.exit("ERROR: duplicate sample name %r (from %s)" % (s, path))
        c = read_one(path)
        per_sample[s] = c
        order.append(s)
        for taxon, (taxid, _) in c.items():
            taxid_of.setdefault(taxon, taxid)

    taxa = sorted(taxid_of)

    # Refuse to write an empty or truncated matrix - the exact failure mode of
    # the rule this replaces.
    if not taxa:
        sys.exit("ERROR: no taxa parsed from %d files; refusing to write an "
                 "empty matrix" % len(a.files))
    if len(order) != len(a.files):
        sys.exit("ERROR: %d files in, %d samples out" % (len(a.files), len(order)))

    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    total = 0
    nonzero = 0
    with open(a.out, "w") as out:
        out.write("taxonomy_id\ttaxon\t" + "\t".join(order) + "\n")
        for taxon in taxa:
            row = []
            for s in order:
                n = per_sample[s].get(taxon, (None, 0))[1]
                total += n
                if n:
                    nonzero += 1
                row.append(str(n))
            out.write("%s\t%s\t%s\n" % (taxid_of[taxon], taxon, "\t".join(row)))

    cells = len(taxa) * len(order)
    sys.stderr.write(
        "%s matrix: %d taxa x %d samples, %d reads, %.1f%% zeros -> %s\n"
        % (a.level, len(taxa), len(order), total,
           100.0 * (cells - nonzero) / cells, a.out)
    )


if __name__ == "__main__":
    main()
