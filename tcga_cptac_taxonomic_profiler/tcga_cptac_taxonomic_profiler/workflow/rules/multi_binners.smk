# Additional binners + consensus refinement.
#
# The pipeline previously ran VAMB alone. The v1 MAGs in the manuscript were
# produced by Aviary using CONCOCT (42 bins), MaxBin2 (23) and MetaBAT (16) -
# none of them VAMB - so VAMB-only was a regression on what actually generated
# the published genomes.
#
# All binners here reuse the multi-sample catalogue and per-sample sorted BAMs
# built for VAMB. Multi-sample is the best-performing binning mode in recent
# benchmarks (Nat Commun 2025; Brief Bioinform 2025).
#
# Three algorithm families, so the consensus diversifies the failure mode rather
# than just the tool count:
#   composition + coverage models : MetaBAT2, CONCOCT
#   deep learning embeddings      : VAMB, SemiBin2, COMEBin
#   consensus / refinement        : Binette (scores with CheckM2)

rule metabat2_depth:
    """Contig depth from the sorted BAMs, for MetaBAT2 and CONCOCT."""
    input:
        bams = expand(os.path.join(VAMB_BAMS, '{sample}_sorted.bam'), sample=SAMPLES)
    output:
        depth = os.path.join(METABAT2_RESULTS, 'depth.txt')
    conda:
        os.path.join('..', 'envs', 'metabat2.yaml')
    log:
        os.path.join(LOGS, 'binners', 'metabat2_depth.log')
    resources:
        mem_mb = config.resources.big.mem,
        time = config.resources.med.time,
        tmpdir = config.tmpdir
    threads:
        config.resources.med.cpu
    shell:
        """
        mkdir -p $(dirname {output.depth})
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bams} 2> {log}
        """


rule run_metabat2:
    """MetaBAT2: tetranucleotide frequency + coverage correlation."""
    input:
        catalogue = os.path.join(VAMB_CATALOGUE, 'catalogue.fna.gz'),
        depth = os.path.join(METABAT2_RESULTS, 'depth.txt')
    output:
        flag = os.path.join(FLAGS, 'metabat2.flag')
    params:
        outprefix = os.path.join(METABAT2_RESULTS, 'bins', 'bin'),
        min_contig = config.binning.min_contig_length
    conda:
        os.path.join('..', 'envs', 'metabat2.yaml')
    log:
        os.path.join(LOGS, 'binners', 'metabat2.log')
    benchmark:
        os.path.join(BENCHMARKS, 'binners', 'metabat2.txt')
    resources:
        mem_mb = config.resources.big.mem,
        time = config.resources.big.time,
        tmpdir = config.tmpdir
    threads:
        config.resources.big.cpu
    shell:
        """
        mkdir -p $(dirname {params.outprefix})
        metabat2 --inFile {input.catalogue} --abdFile {input.depth} \
            --outFile {params.outprefix} --minContig {params.min_contig} \
            --numThreads {threads} > {log} 2>&1
        touch {output.flag}
        """


rule run_concoct:
    """CONCOCT: Gaussian mixture model on composition + coverage.

    Produced 42 of the 81 bins in the v1 manuscript - the single largest
    contributor - so it is included for continuity as well as diversity.
    """
    input:
        catalogue = os.path.join(VAMB_CATALOGUE, 'catalogue.fna.gz'),
        # Reuse the jgi depth table rather than concoct_coverage_table.py, which
        # is broken on modern python (see the shell block).
        depth = os.path.join(METABAT2_RESULTS, 'depth.txt')
    output:
        flag = os.path.join(FLAGS, 'concoct.flag')
    params:
        outdir = CONCOCT_RESULTS,
        min_contig = config.binning.min_contig_length
    conda:
        os.path.join('..', 'envs', 'concoct.yaml')
    log:
        os.path.join(LOGS, 'binners', 'concoct.log')
    benchmark:
        os.path.join(BENCHMARKS, 'binners', 'concoct.txt')
    resources:
        mem_mb = config.resources.big.mem,
        time = config.resources.big.time,
        tmpdir = config.tmpdir
    threads:
        config.resources.big.cpu
    shell:
        """
        mkdir -p {params.outdir}/bins
        cd {params.outdir}
        zcat {input.catalogue} > catalogue.fna

        # concoct_coverage_table.py is broken under modern python:
        #   TypeError: write() argument must be str, not bytes
        # (it writes subprocess bytes straight to sys.stderr). CONCOCT 1.1.0 is
        # unmaintained, so rather than patch a shipped script the coverage table
        # is derived from the jgi depth file already built for MetaBAT2.
        #
        # depth.txt columns: 1=contigName 2=contigLen 3=totalAvgDepth then
        # per-BAM (depth, var) pairs. CONCOCT wants contigName + the depth
        # columns, i.e. 1,4,6,8,...
        #
        # Built with cut rather than awk on purpose: snakemake passes shell
        # blocks through python .format(), so a backslash-t or backslash-n in an
        # awk program is interpreted THERE and awk receives a literal tab or
        # newline - which silently breaks the program mid-string. cut defaults to
        # tab-delimited, so no escape sequences are needed anywhere.
        NCOL=$(head -1 {input.depth} | awk '{{print NF}}')
        COLS=1
        for i in $(seq 4 2 "$NCOL"); do COLS="$COLS,$i"; done
        cut -f"$COLS" {input.depth} > coverage_table.tsv
        echo "coverage table: $(wc -l < coverage_table.tsv) rows from $NCOL depth columns; cols=$COLS" >> {log}

        # Run on whole contigs. The cut-up step is optional and would require
        # regenerating coverage against the chunk names.
        concoct --composition_file catalogue.fna --coverage_file coverage_table.tsv \
            -l {params.min_contig} \
            -b concoct_out --threads {threads} >> {log} 2>&1
        merge_cutup_clustering.py concoct_out_clustering_gt{params.min_contig}.csv > clustering_merged.csv 2>> {log} || \
          cp concoct_out_clustering_gt{params.min_contig}.csv clustering_merged.csv
        extract_fasta_bins.py catalogue.fna clustering_merged.csv --output_path bins >> {log} 2>&1
        rm -f catalogue.fna
        touch {output.flag}
        """


rule run_semibin2:
    """SemiBin2 in co-assembly mode. Among the best performers in recent benchmarks.

    Uses single_easy_bin, NOT multi_easy_bin. multi_easy_bin splits the
    concatenated catalogue back into per-sample sub-assemblies on the separator,
    then trains per sample -- and in this cohort some samples retain too few
    contigs after the 1500 bp filter, which crashes training with

        ValueError: operands could not be broadcast together with shapes (0,30)

    That will only get worse on the full 336, where many blood normals carry
    almost no microbial sequence. single_easy_bin treats the catalogue as one
    co-assembly with N BAMs, which is exactly how VAMB, MetaBAT2 and CONCOCT
    treat it -- so all four now bin the same object as well as the same contigs.
    """
    input:
        catalogue = os.path.join(VAMB_CATALOGUE, 'catalogue.fna.gz'),
        bams = expand(os.path.join(VAMB_BAMS, '{sample}_sorted.bam'), sample=SAMPLES)
    output:
        flag = os.path.join(FLAGS, 'semibin2.flag')
    params:
        outdir = SEMIBIN2_RESULTS,
        separator = config.binning.separator,
        min_contig = config.binning.min_contig_length
    conda:
        os.path.join('..', 'envs', 'semibin2.yaml')
    log:
        os.path.join(LOGS, 'binners', 'semibin2.log')
    benchmark:
        os.path.join(BENCHMARKS, 'binners', 'semibin2.txt')
    resources:
        mem_mb = config.resources.big.mem,
        time = config.resources.big.time,
        tmpdir = config.tmpdir
    threads:
        config.resources.big.cpu
    shell:
        """
        SemiBin2 single_easy_bin \
            --input-fasta {input.catalogue} \
            --input-bam {input.bams} \
            --output {params.outdir} \
            --min-len {params.min_contig} \
            --threads {threads} > {log} 2>&1

        # NORMALISE: multi_easy_bin writes per-sample bins to
        # OUTDIR/samples/SAMPLE/output_bins/, not a flat bins/ directory.
        # Binette takes explicit --bin_dirs, so every binner is normalised to
        # OUTDIR/bins to remove any dependence on a tool's internal layout.
        mkdir -p {params.outdir}/bins
        find {params.outdir} -path "{params.outdir}/bins" -prune -o \
             -name "*.fa" -print -o -name "*.fna" -print -o -name "*.fasta" -print 2>/dev/null \
          | while read f; do ln -sf "$f" "{params.outdir}/bins/$(echo "$f" | md5sum | cut -c1-8)_$(basename "$f")"; done
        # Verify before flagging. run_comebin.sh prints "Data augmentation
        # exited with status 1" and then EXITS 0, so bash strict mode does not
        # catch it: the rule ran to completion, normalised nothing, touched the
        # flag, and snakemake recorded COMPLETED with zero bins produced. Never
        # trust this wrapper's exit code - check for output instead.
        n=$(find {params.outdir}/bins -maxdepth 1 -name "*.fa" -print -o -maxdepth 1 -name "*.fna" -print -o -maxdepth 1 -name "*.fasta" -print 2>/dev/null | wc -l)
        echo "normalised bins: $n" >> {log}
        if [ "$n" -eq 0 ]; then
            echo "ERROR: COMEBin produced no bins; refusing to write comebin.flag" >> {log}
            exit 1
        fi
        touch {output.flag}
        """


rule run_comebin:
    """COMEBin: contrastive multiple-view learning.

    Best-performing deep learning binner in recent benchmarks but also the
    slowest, which is acceptable here since throughput is not the constraint.

    NOTE: run_comebin.sh defaults to -d cuda. Compute nodes are CPU, so the
    device is set explicitly from config; leaving it default would fail.
    """
    input:
        catalogue = os.path.join(VAMB_CATALOGUE, 'catalogue.fna.gz'),
        bams = expand(os.path.join(VAMB_BAMS, '{sample}_sorted.bam'), sample=SAMPLES)
    output:
        flag = os.path.join(FLAGS, 'comebin.flag')
    params:
        outdir = COMEBIN_RESULTS,
        bamdir = VAMB_BAMS,
        device = config.binning.comebin_device,
        min_contig = config.binning.min_contig_length
    conda:
        os.path.join('..', 'envs', 'comebin.yaml')
    log:
        os.path.join(LOGS, 'binners', 'comebin.log')
    benchmark:
        os.path.join(BENCHMARKS, 'binners', 'comebin.txt')
    resources:
        mem_mb = config.resources.big.mem,
        time = config.resources.big.time,
        tmpdir = config.tmpdir
    threads:
        config.resources.big.cpu
    shell:
        """
        mkdir -p {params.outdir}
        # COMEBin has NO minimum-contig-length option -- its -l is the loss
        # function temperature, not a length. So the catalogue is filtered here
        # instead, ensuring COMEBin bins the same contig set as the other four.
        # Length-filter the catalogue so COMEBin bins the same contig set as the
        # others. Uses ORS (default newline) rather than a backslash-n literal:
        # snakemake formats shell blocks through python, so backslash escapes in
        # an awk program are interpreted there and awk receives a real newline,
        # which breaks the program mid-string. Same trap that broke CONCOCT.
        zcat {input.catalogue} \
          | awk -v m={params.min_contig} '
                /^>/ {{ if (h != "" && length(s) >= m) print h ORS s; h=$0; s=""; next }}
                     {{ s = s $0 }}
                END  {{ if (h != "" && length(s) >= m) print h ORS s }}' \
          > {params.outdir}/catalogue.fna
        echo "comebin input contigs >= {params.min_contig}: $(grep -c '^>' {params.outdir}/catalogue.fna)" >> {log}

        # COMEBin's -p takes a DIRECTORY and reads every .bam in it. VAMB_BAMS
        # holds both {{sample}}.bam and {{sample}}_sorted.bam - 60 files for 30
        # samples - so pointing at it directly makes COMEBin see each sample
        # twice and abort with
        #     ValueError: BAM coverage contains repeated groups for retained
        #     contig S3:NODE_1564_...
        # Stage exactly the declared sorted inputs into a private directory so
        # the rule cannot pick up intermediates that happen to sit alongside.
        rm -rf {params.outdir}/bams
        mkdir -p {params.outdir}/bams
        for b in {input.bams}; do ln -sf "$b" "{params.outdir}/bams/$(basename "$b")"; done
        echo "comebin bam inputs: $(ls {params.outdir}/bams | wc -l)" >> {log}

        # Append, do not truncate: > would wipe the two lines just written.
        run_comebin.sh -a {params.outdir}/catalogue.fna \
            -o {params.outdir} \
            -p {params.outdir}/bams \
            -t {threads} \
            -d {params.device} >> {log} 2>&1
        rm -f {params.outdir}/catalogue.fna

        # NORMALISE to OUTDIR/bins (see the SemiBin2 rule). COMEBin writes to
        # comebin_res/comebin_res_bins/, but that layout is version-dependent.
        mkdir -p {params.outdir}/bins
        find {params.outdir} -path "{params.outdir}/bins" -prune -o \
             -name "*.fa" -print -o -name "*.fna" -print -o -name "*.fasta" -print 2>/dev/null \
          | while read f; do ln -sf "$f" "{params.outdir}/bins/$(echo "$f" | md5sum | cut -c1-8)_$(basename "$f")"; done
        echo "normalised bins: $(ls {params.outdir}/bins 2>/dev/null | wc -l)" >> {log}
        touch {output.flag}
        """


rule run_binette:
    """Consensus refinement across all five bin sets, scored with CheckM2.

    Binette selects the best-scoring bin per genome across the input bin
    directories, so a diverse candidate set matters more than any single
    binner's performance.
    """
    input:
        vamb = os.path.join(FLAGS, 'vamb.flag'),
        metabat2 = os.path.join(FLAGS, 'metabat2.flag'),
        concoct = os.path.join(FLAGS, 'concoct.flag'),
        semibin2 = os.path.join(FLAGS, 'semibin2.flag'),
        comebin = os.path.join(FLAGS, 'comebin.flag'),
        catalogue = os.path.join(VAMB_CATALOGUE, 'catalogue.fna.gz')
    output:
        flag = os.path.join(FLAGS, 'binette.flag')
    params:
        outdir = BINETTE_RESULTS,
        checkm2_db = config.databases.checkm2,
        # Every binner normalises its output to {RESULTS}/bins, so these paths
        # do not depend on any tool's internal layout (SemiBin2 in particular
        # writes samples/{sample}/output_bins/, and COMEBin's nesting is
        # version-dependent).
        bin_dirs = " ".join([
            os.path.join(VAMB_RESULTS, 'bins_flat'),
            os.path.join(METABAT2_RESULTS, 'bins'),
            os.path.join(CONCOCT_RESULTS, 'bins'),
            os.path.join(SEMIBIN2_RESULTS, 'bins'),
            os.path.join(COMEBIN_RESULTS, 'bins'),
        ]),
        contamination_weight = config.binning.binette_contamination_weight,
        # Must go through params: a shell block can only interpolate input,
        # output, params, wildcards, threads, resources, config and log. A bare
        # global such as VAMB_RESULTS raises NameError at runtime.
        vamb_results = VAMB_RESULTS
    conda:
        os.path.join('..', 'envs', 'binette.yaml')
    log:
        os.path.join(LOGS, 'binners', 'binette.log')
    benchmark:
        os.path.join(BENCHMARKS, 'binners', 'binette.txt')
    resources:
        mem_mb = config.resources.big.mem,
        time = config.resources.big.time,
        tmpdir = config.tmpdir
    threads:
        config.resources.big.cpu
    shell:
        """
        mkdir -p {params.outdir}

        # VAMB is the one binner not normalised in its own rule (it predates
        # these), and it is also the only one whose native layout is NOT flat:
        # it writes bins/<sample>/*.fna, because the samples_with_bins
        # checkpoint globs bins/*/ for per-sample directories.
        #
        # The previous version of this block pruned bins/ and searched the rest
        # of VAMB_RESULTS - i.e. it excluded the only directory containing any
        # fasta, and silently linked nothing. Binette then received a directory
        # holding 21 sub-directories and no genomes.
        #
        # Flatten into a SEPARATE bins_flat/ so the checkpoint's bins/*/ layout
        # is left intact. The path-hash prefix matches what every other binner
        # here does; VAMB's names happened to be unique across all 30 samples in
        # the subset test (145 files, 145 distinct basenames), but flattening a
        # per-sample tree on that assumption is not worth the risk.
        mkdir -p {params.vamb_results}/bins_flat
        find {params.vamb_results}/bins -mindepth 2 -name "*.fna" -print -o -name "*.fa" -print 2>/dev/null \
          | while read f; do ln -sf "$f" "{params.vamb_results}/bins_flat/$(echo "$f" | md5sum | cut -c1-8)_$(basename "$f")"; done

        # Fail loudly on an empty bin set. Binette accepts an empty --bin_dirs
        # entry without complaint, which would silently reduce the consensus
        # from five binners to fewer and be very hard to notice afterwards.
        # Count FASTA FILES, not directory entries. `ls | wc -l` counted
        # VAMB's 21 per-sample sub-directories as 21 bins and let an
        # effectively empty bin set through the check it exists to prevent.
        for d in {params.bin_dirs}; do
            n=$(find "$d" -maxdepth 1 -name "*.fa" -print -o -maxdepth 1 -name "*.fna" -print -o -maxdepth 1 -name "*.fasta" -print 2>/dev/null | wc -l)
            echo "bin_dir $d -> $n bins" >> {log}
            if [ "$n" -eq 0 ]; then
                echo "ERROR: $d contains no bins; refusing to run a degraded consensus" >> {log}
                exit 1
            fi
        done

        zcat {input.catalogue} > {params.outdir}/catalogue.fna
        binette --bin_dirs {params.bin_dirs} \
            --contigs {params.outdir}/catalogue.fna \
            --outdir {params.outdir} \
            --checkm2_db {params.checkm2_db} \
            --contamination_weight {params.contamination_weight} \
            --threads {threads} >> {log} 2>&1
        rm -f {params.outdir}/catalogue.fna
        touch {output.flag}
        """
