#!/usr/bin/env bash
# Count read pairs in each unmapped FASTQ, for the pre/post host-depletion
# comparison and for correcting the bacterial-load denominator.
#
# Usage: scripts/count_input_reads.sh <FASTQ_DIR> > input_read_counts.tsv
#
# Output: sample <TAB> pairs
# Note the R1 record count IS the pair count for properly paired files.

set -euo pipefail
dir="${1:?usage: $0 <FASTQ_DIR>}"

printf "sample\tpairs\n"
for f in "$dir"/*_R1.fastq.gz; do
    [ -e "$f" ] || continue
    s=$(basename "$f" _R1.fastq.gz)
    n=$(( $(zcat -f "$f" | wc -l) / 4 ))
    printf "%s\t%s\n" "$s" "$n"
done
