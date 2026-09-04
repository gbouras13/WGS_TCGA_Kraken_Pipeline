"""Aggregate the per-sample Bracken reports into taxon x sample count matrices.

REPLACES kraken-biom, which produced a valid but completely EMPTY table:

    {"generated_by": "kraken-biom v1.2.0", "shape": [0, 0],
     "data": [], "rows": [], "columns": []}

342 bytes, and the pipeline recorded the stage as complete. The old rule was:

    params: ' '.join(expand(... genus ...)), ' '.join(expand(... species ...))
    shell:  kraken-biom {input[0]} -o {output[0]} --fmt json
            kraken-biom {input[1]} -o {output[1]} --fmt json

Two independent faults:

  1. `{input[0]}` is the FIRST FILE of the expanded input list, not the list.
     The rule built the correct space-joined file lists in params and then never
     referenced them, so kraken-biom was handed exactly one sample.
  2. kraken-biom parses Kraken2 *report* format. Bracken's tabular output has a
     different shape, so it matched nothing and wrote an empty table instead of
     failing.

The consequence was not confined to this rule: the downstream decontamination
(W1a in the analysis repo) reads bracken_species.biom, so an empty file silently
blocked the whole contaminant-filtering step.

The replacement writes plain TSV via workflow/scripts/bracken_to_matrix.py,
which uses the python standard library only - no conda environment, so it cannot
be broken by a re-solve and needs no network on a compute node. It refuses to
write an empty or truncated matrix and exits non-zero on any parse failure.

Verified on the 336-sample v2 run: 1,346 species x 336 and 453 genera x 336,
516,388,104 reads, matching an independent R implementation cell for cell.
"""


rule bracken_matrix:
    input:
        genus = expand(os.path.join(BRACKEN, "{sample}.kraken_bracken_genus.txt"), sample=SAMPLES),
        species = expand(os.path.join(BRACKEN, "{sample}.kraken_bracken_species.txt"), sample=SAMPLES)
    output:
        genus = os.path.join(BIOM, "bracken_genus_counts.tsv"),
        species = os.path.join(BIOM, "bracken_species_counts.tsv")
    params:
        script = os.path.join(workflow.basedir, 'scripts', 'bracken_to_matrix.py')
    log:
        os.path.join(LOGS, "biom", "bracken_matrix.log")
    benchmark:
        os.path.join(BENCHMARKS, "biom", "bracken_matrix.log")
    resources:
        mem_mb = config.resources.sml.mem,
        time = config.resources.med.time
    threads:
        1
    shell:
        """
        python3 {params.script} --level genus   --out {output.genus}   {input.genus}   >> {log} 2>&1
        python3 {params.script} --level species --out {output.species} {input.species} >> {log} 2>&1

        # Belt and braces: the rule this replaces wrote a file that existed and
        # was empty, so existence is not the check. Require a header plus at
        # least one data row in each.
        for f in {output.genus} {output.species}; do
            n=$(wc -l < "$f")
            echo "$f: $n lines" >> {log}
            if [ "$n" -lt 2 ]; then
                echo "ERROR: $f has no data rows" >> {log}
                exit 1
            fi
        done
        """


rule aggr_biom:
    """Aggregate the Bracken count matrices."""
    input:
        os.path.join(BIOM, "bracken_genus_counts.tsv"),
        os.path.join(BIOM, "bracken_species_counts.tsv")
    output:
        os.path.join(FLAGS, "aggr_biom.flag")
    resources:
        mem_mb = config.resources.sml.mem,
        time = config.resources.sml.time
    threads:
        1
    shell:
        """
        touch {output[0]}
        """
