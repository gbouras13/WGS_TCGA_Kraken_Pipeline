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

TAXONOMY. kraken-biom attached a full lineage (Rank1..Rank7), and the
downstream decontamination depends on it - it prunes eukaryotes with
Rank1 == "k__Eukaryota", which matters because Homo sapiens is the second
largest taxon in this cohort. Bracken's tabular output carries no lineage, so
--tax-out reconstructs one from the Kraken2-format reports (.rep) that sit
alongside: rank code in column 4, taxid in column 5, and two spaces of
indentation per level in column 6.

    taxonomy_id  taxon  domain  phylum  class  order  family  genus  species

Usage:
    bracken_to_matrix.py --level species --out matrix.tsv FILE [FILE ...]
    bracken_to_matrix.py --level species --out matrix.tsv \
        --tax-out taxonomy.tsv --report-dir RESULTS/KRAKEN FILE [FILE ...]
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


# Kraken rank codes we keep, in order. Sub-ranks (D1, G1, ...) are skipped:
# kraken-biom's Rank1..Rank7 held the primary ranks only.
PRIMARY_RANKS = [("D", "domain"), ("P", "phylum"), ("C", "class"),
                 ("O", "order"), ("F", "family"), ("G", "genus"),
                 ("S", "species")]
RANK_INDEX = {code: i for i, (code, _) in enumerate(PRIMARY_RANKS)}


def parse_report(path, lineage_of):
    """Accumulate taxid -> lineage from one Kraken2-format report.

    Columns: pct, clade_reads, taxon_reads, rank_code, taxid, indented name.
    Depth comes from the leading spaces of the name, two per level.
    """
    stack = []  # (depth, rank_index, name)
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            rank, taxid, raw = f[3], f[4], f[5]
            depth = (len(raw) - len(raw.lstrip(" "))) // 2
            name = raw.strip()
            while stack and stack[-1][0] >= depth:
                stack.pop()
            ri = RANK_INDEX.get(rank)
            stack.append((depth, ri, name))
            if ri is None:
                continue  # sub-rank or root: contributes no lineage column
            lin = [""] * len(PRIMARY_RANKS)
            for _, sri, sname in stack:
                if sri is not None:
                    lin[sri] = sname
            prev = lineage_of.get(taxid)
            # keep the most complete lineage seen for this taxid
            if prev is None or sum(bool(x) for x in lin) > sum(bool(x) for x in prev[1]):
                lineage_of[taxid] = (name, lin)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", required=True, choices=("species", "genus"))
    ap.add_argument("--out", required=True)
    ap.add_argument("--tax-out", help="write a taxonomy table here")
    ap.add_argument("--report-dir",
                    help="directory of Kraken2-format .rep files (needed for --tax-out)")
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

    if a.tax_out:
        if not a.report_dir:
            sys.exit("ERROR: --tax-out requires --report-dir")
        suffix = ".kraken_bracken_%s.rep" % (
            "genuses" if a.level == "genus" else a.level)
        reports = sorted(
            os.path.join(a.report_dir, f)
            for f in os.listdir(a.report_dir) if f.endswith(suffix))
        if not reports:
            sys.exit("ERROR: no *%s under %s" % (suffix, a.report_dir))
        lineage_of = {}
        for r in reports:
            parse_report(r, lineage_of)
        by_name = {n: lin for n, lin in lineage_of.values()}
        missing = [t for t in taxa if t not in by_name]
        if missing:
            sys.exit("ERROR: no lineage for %d of %d taxa (e.g. %s); "
                     "reports and counts disagree"
                     % (len(missing), len(taxa), ", ".join(missing[:3])))
        with open(a.tax_out, "w") as out:
            out.write("taxonomy_id\ttaxon\t"
                      + "\t".join(n for _, n in PRIMARY_RANKS) + "\n")
            for t in taxa:
                out.write("%s\t%s\t%s\n"
                          % (taxid_of[t], t, "\t".join(by_name[t])))
        sys.stderr.write("%s taxonomy: %d taxa from %d reports -> %s\n"
                         % (a.level, len(taxa), len(reports), a.tax_out))

    cells = len(taxa) * len(order)
    sys.stderr.write(
        "%s matrix: %d taxa x %d samples, %d reads, %.1f%% zeros -> %s\n"
        % (a.level, len(taxa), len(order), total,
           100.0 * (cells - nonzero) / cells, a.out)
    )


if __name__ == "__main__":
    main()
