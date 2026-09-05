rule run_bakta:
    """
    generate mash db
    """
    input:
        mag =  os.path.join(ALL_MAGS, '{mag}.fna')
    output:
        out_tsv = os.path.join(BAKTA, '{mag}', '{mag}.tsv') ,
        gbk = os.path.join(BAKTA, '{mag}', '{mag}.gbff'),
        # baktfold consumes bakta's JSON. bakta writes it by default; declaring
        # it here is what makes it a tracked output rather than an incidental
        # file, so a partial bakta run cannot silently satisfy baktfold.
        json = os.path.join(BAKTA, '{mag}', '{mag}.json')
    threads:
        config.resources.med.cpu
    resources:
        mem_mb = config.resources.med.mem,
        time = config.resources.med.time
    conda:
        os.path.join("..", "envs", "bakta.yaml")
    benchmark:
        os.path.join(BENCHMARKS, 'bakta', "{mag}.txt")
    log:
        os.path.join(LOGS, 'bakta', "{mag}.log")
    params:
        outdir = os.path.join(BAKTA, '{mag}'),
        db=config.databases.bakta,
    shell:
        """
        # Redirect to {log}. Without this the rule declared a log it never
        # wrote to, and bakta's own errors landed only in the scheduler's
        # stderr - which is how a pyhmmer incompatibility looked like silent
        # failure across 127 of 128 MAGs.
        bakta --db {params.db} --output {params.outdir} -f -t {threads} {input.mag} >> {log} 2>&1

        # Verify rather than trust the exit code.
        if [ ! -s {output.json} ]; then
            echo "ERROR: bakta produced no JSON for {wildcards.mag}" >> {log}
            exit 1
        fi
        """

rule run_baktfold:
    """Structural annotation of the proteins bakta could only call hypothetical.

    baktfold converts amino acid sequences to 3Di tokens with ProstT5 and
    searches SwissProt, the AlphaFold database, PDB and CATH with foldseek, so it
    assigns function where sequence similarity alone found none. It runs AFTER
    bakta and takes bakta's JSON.

    Applied to the same HQ_MED_MAGS set as the rest of the annotation: there is
    no value in structurally annotating a bin that is too incomplete or too
    contaminated to interpret.
    """
    input:
        json = os.path.join(BAKTA, '{mag}', '{mag}.json')
    output:
        inference = os.path.join(BAKTFOLD, '{mag}', '{mag}.inference.tsv')
    threads:
        config.resources.med.cpu
    resources:
        mem_mb = config.resources.med.mem,
        time = config.resources.med.time,
        # Sends this rule to a GPU partition; every other rule takes the profile
        # default, which is a no-op flag. Verified 2026-09-05 on p2-gpu-29: the
        # environment reports cuda available True and completes a GPU matmul.
        slurm_extra = config.resources.gpu.slurm_extra
    conda:
        os.path.join("..", "envs", "baktfold.yaml")
    benchmark:
        os.path.join(BENCHMARKS, 'baktfold', "{mag}.txt")
    log:
        os.path.join(LOGS, 'baktfold', "{mag}.log")
    params:
        outdir = os.path.join(BAKTFOLD, '{mag}'),
        db = config.databases.baktfold
    shell:
        """
        baktfold run -i {input.json} -o {params.outdir} -f -t {threads} -d {params.db} \
            --foldseek-gpu >> {log} 2>&1

        # Verify rather than trust the exit code - four tools in this workflow
        # return 0 on failure (see RERUN.md).
        if [ ! -s {output.inference} ]; then
            echo "ERROR: baktfold produced no inference table for {wildcards.mag}" >> {log}
            exit 1
        fi
        """


rule run_bakta_mimag:
    """
    generate counts of tRNA, 5S, 16S and 23S genes for each isolate
    """
    input:
        tsvs = expand(os.path.join(BAKTA, '{mag}', '{mag}.tsv') , mag=HQ_MED_MAGS),
    output:
        out_tsv = os.path.join(BAKTA, 'mimag.tsv')
    threads:
        config.resources.sml.cpu
    resources:
        mem_mb = config.resources.sml.mem,
        time = config.resources.sml.time
    conda:
        os.path.join("..", "envs", "biopython.yaml")
    benchmark:
        os.path.join(BENCHMARKS, 'bakta_mimag', "bakta_mimag.txt")
    log:
        os.path.join(LOGS, 'bakta_mimag', "bakta_mimag.log")
    params:
        bakta_dir=BAKTA
    script:
        '../scripts/mimag.py'

rule aggr_bakta:
    """
    generate mash db
    """
    input:
        mimag = os.path.join(BAKTA, 'mimag.tsv'),
        baktfold = expand(os.path.join(BAKTFOLD, '{mag}', '{mag}.inference.tsv'), mag=HQ_MED_MAGS)
    output:
        outtouch = os.path.join(FLAGS, 'bakta.flag'),
    threads:
        config.resources.sml.cpu
    resources:
        mem_mb = config.resources.sml.mem,
        time = config.resources.sml.time
    shell:
        """
        touch {output.outtouch}

        """
