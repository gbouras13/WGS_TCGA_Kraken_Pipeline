# Full re-run: TCGA/CPTAC HNSC taxonomic profiling (v2)

Operational record for the complete re-profiling of the HNSC cohort with
in-workflow Deacon host depletion. Everything downstream of host depletion is
invalidated by this change and is regenerated from scratch.

Branch: `deacon-t2t`. Pipeline package: `tcga_cptac_taxonomic_profiler` v0.1.0.

## What changed from the original run

| | Original | v2 |
|---|---|---|
| Host depletion | Performed outside the workflow; parameters not recorded | `host-removal` stage, Deacon vs `panhuman-1`, recorded in config |
| Host reference | `human-t2t-hla-phix174.fa` (T2T-CHM13 + HLA + phiX174) | `panhuman-1` = HPRC year-1 assemblies + CHM13v2.0 + GRCh38.p14, bacterial/viral removed |
| Residual human QC | Not measured | Per-sample `--summary` JSON from every sample |
| Singletons | An `_RS` file existed from the manual step | Deacon keeps/drops pairs together; empty `_RS` emitted, spades `-s` now conditional |
| Classifiers | Kraken2 only (plus SingleM for fraction) | Kraken2 **and** Metabuli, run on identical input and compared |

SingleM is deliberately unaffected: `singlem.smk` reads the **raw** unaligned
FASTQs, not the host-depleted ones, because its microbial-fraction calculation
needs the whole metagenome. This means SingleM remains an entirely independent
measurement of the same quantity — worth stating explicitly in the Methods,
because it makes SingleM a genuine orthogonal control rather than a
re-derivation of the Kraken2 result.

## Starting point: FASTQs, not BAMs

The source BAMs are no longer held. What remains is the output of the original
`extract` stage: the hg19-unmapped reads as `{sample}_R{1,2}.fastq.gz`
(`samtools view -f 12 -F 256`, i.e. both mates unmapped against
HG19_Broad_Variant), plus the per-sample read counts in `Data/READCOUNT/`.

This is sufficient. **Stage 1 (`extract`) and the `read_counts` rule are not
re-run** — their outputs are already in hand and are deterministic. The re-run
begins at stage 2.

Verified: all 308 metadata rows match their saved `_readcount_all.txt` value
exactly, so the bacterial-load denominator is fully preserved. There are 336
readcount files against 308 de-duplicated metadata rows, which is the source of
the "336 WGS samples" figure in the 2025 Methods draft.

Two consequences to state in the Methods as limitations, both conservative:

- Reads that mapped to hg19 are gone. hg19 contains known bacterial sequence in
  unplaced scaffolds, so a small number of genuinely microbial reads were
  discarded at source. This loses signal; it cannot create false positives.
- `-f 12` requires *both* mates unmapped, so a microbial read whose mate
  spuriously mapped to hg19 was also discarded. Again conservative.

## Before starting

1. **The denominator is confirmed wrong by a factor of exactly 2.**
   `kraken.smk` runs `kraken2 --paired`, classifying each *pair* as one unit,
   while `read_counts.smk` counted *records* with `samtools view -c`.

   Verified on `TCGA-IQ-7632-01A-11D-2317` (51 bp):

   | | |
   |---|---|
   | `samtools view -c` (the denominator used) | 406,415,890 records |
   | Ge et al. `tot_rds` for the same aliquot | 203,207,945 fragments |
   | ratio | **2.0000** |

   So every reported bacterial-load percentage is exactly half its true value.
   (Per-sample figures are held in the private analysis repository.)

   Hazard ratios are unaffected — a constant factor on a logged covariate shifts
   the intercept only. Affected: the Figure 1 y-axis, every absolute percentage
   in the text, the SingleM concordance figure, and any cross-study comparison.

   Fix in the rebuild by using pair counts on both sides. `scripts/count_input_reads.sh`
   produces per-sample pair counts from the FASTQs.

   (Across 112 aliquots matched to Ge et al. the ratio ranged 1.84-2.32 rather
   than sitting at exactly 2. That spread is most likely run selection: the
   Methods note that where multiple sequencing runs existed, the BAM with more
   bases was used, so the matched aliquot is not always the same run Ge
   profiled. The exact 2.0000 on a directly comparable sample is the signal.)

2. **Decide on the Kraken2 database.** Currently `k2_pluspf_20230605` (June
   2023), whose human reference is GRCh38 — so human reads get a second chance
   to be misclassified even after clean depletion. A build including CHM13v2.0
   closes that gap. Since a full re-run is happening anyway, updating is the
   right moment, but note it changes taxonomy assignments and the decontamination
   thresholds must be re-derived, not carried over.
3. **Metabuli database: GTDB R226, fetched directly from S3.**

   `metabuli databases GTDB` **cannot be used**. metabuli 1.2.0 ships stale
   URLs, the bucket has been reorganised (the old names now sit under
   `metabuli/archive/`), and every advertised URL 404s. Worse, the command
   **exits 0 on failure**, leaving a corrupt tarball — so it must never be
   wrapped in a Snakemake rule as-is.

   Current databases, listed from the bucket rather than the index page:

   | Database | Size |
   |---|---|
   | `metabuli/gtdb226.tar.gz` | 298 GiB |
   | `metabuli/gtdb232.tar.gz` | 634 GB |
   | `metabuli/refseq_standard.tar.gz` | 80 GB |

   `scripts/`-style helper: `metagenomics_dbs/metabuli/fetch_gtdb226.sh` (login
   node only, single-threaded, niced, resumable with `wget -c`). Peak RAM during
   download is ~13 MB.

   GTDB was chosen over RefSeq. The comparison of interest is bacterial *load*
   — "is this read bacterial" — which is insensitive to taxonomy version, so
   the earlier argument for matching Kraken2's NCBI taxonomy carries little
   weight. GTDB instead aligns Metabuli with SingleM and GTDB-Tk.

4. **GTDB releases: R226 for Metabuli and GTDB-Tk, R232 for SingleM.**

   Full alignment on R226 turned out to be impossible. SingleM 0.20.3 reads the
   R226 metapackage but is broken: it treats any bytes on DIAMOND's stderr as a
   fatal error, and DIAMOND >180 days old prints a version nag, so every sample
   failed while DIAMOND was in fact exiting 0. The fix is in singlem 0.21.4 —
   which hard-rejects the R226 metapackage:

       Exception: The metapackage defined by the SINGLEM_METAPACKAGE_PATH
       environment variable is either malformed or does not match the version
       encoded in the version of SingleM being used.

   The changelog wording ("default metapackage updated to GTDB R232") reads like
   a default; it is a requirement. There is no SingleM version that both reads
   R226 and handles DIAMOND correctly.

   | Tool | Was | Now | Reference data |
   |---|---|---|---|
   | Metabuli | — | `gtdb226` | `GTDB226/gtdb+human+virus`, 378 GB |
   | SingleM | 0.18.3, r220 | **`>=0.21.4`** | **`S6.5.0.GTDB_r232.metapackage_20260319.smpkg.zb`**, 2.0 GB (Zenodo 20150069) |
   | GTDB-Tk | 2.3, release214 | `=2.6.1` | `GTDB/release226`, 147 GB |

   **Why the mismatch is acceptable.** SingleM contributes a microbial
   *fraction* — a per-sample scalar — as an independent check on bacterial load.
   It is not used for taxon-by-taxon comparison against Kraken2 or Metabuli, so
   the GTDB release it was built from does not affect the quantity being
   compared. Metabuli and GTDB-Tk, which *do* produce taxonomy that is compared
   and quoted, remain aligned on R226.

   State this explicitly in the Methods rather than leaving it to be noticed.

   The rejected alternative was moving everything to R232: it would restore
   alignment but costs a 634 GB Metabuli re-download, re-running Metabuli, a
   ~150 GB GTDB-Tk re-download and a bump to GTDB-Tk 2.7+, to regain a property
   that mainly matters for MAG taxonomy.

   **MAG taxonomy must be regenerated, not carried over.** The existing calls
   are GTDB r214 assignments, and GTDB renames and reclassifies between
   releases. Any MAG taxonomy quoted in the manuscript has to come from the
   R226 re-run.

   Both bumped environments need rebuilding on the login node before submitting
   anything that uses them:

   ```bash
   tcga_cptac_taxonomic_profiler singlem  --input <dir> --output <dir> --conda-create-envs-only
   tcga_cptac_taxonomic_profiler binning  --input <dir> --output <dir> --conda-create-envs-only
   ```

5. **Install the Deacon index** (~3.3 GB), alongside the existing host reference:

```bash
tcga_cptac_taxonomic_profiler install-host --database /hpcfs/users/a1667917/metagenomics_dbs/host_genome_db
```

## One sample needs --phred-offset 33 (assembled manually)

`TCGA-CN-4738-01A-02D-1509_120420_SN1120_0135_BD0T5FACXX_s_7_rg.sorted` was the
only one of 336 that would not assemble. SPAdes exits after 0.05 s, before
reading any sequence:

    ERROR General (main.cpp:90) Failed to determine offset!
                                Specify it manually and restart
    SPAdes finished with the following error code: -1

**Cause: the sample has no usable quality scores.** Every quality character is
`B`, with no variation across 500,000 reads sampled. SPAdes infers the PHRED
offset from the range of observed quality characters, and a single character
carries no information to separate Phred+33 from Phred+64.

For comparison, working samples show a normal spread:

| Sample | Quality characters (500k reads) |
|---|---|
| TCGA-CN-4738 (failed) | `B` only |
| G23062.TCGA_BA_4077_01B | `"#%&'()*+,-./:;<=>?@$0123456789ABCDEFG` |

The working range starts at `"` (ASCII 34), which would be negative under
Phred+64, so the cohort is unambiguously **Phred+33**. Uniform `B` (ASCII 66)
is therefore Q33 — a placeholder written at some earlier stage, not measured
quality. Note that under Phred+64 the same character would be Q2, the Illumina
"read segment quality control indicator", so the encoding assumption is not
cosmetic; it is what makes these reads usable rather than instrument-flagged.

The reads themselves are fine (3,959,689 pairs) — only the quality track is
uninformative.

**Handled manually rather than in the rule.** `run_phred_fix.sh` in the analysis
directory reruns just this sample with `--phred-offset 33` and otherwise
identical parameters. Adding the flag to `sample_assembly.smk` would change the
rule for all 336 and risk re-running a completed stage, for the sake of one
sample.

Worth stating in the Methods: one sample carried placeholder quality scores and
was assembled with the offset specified explicitly. Any future quality-aware
step (trimming, quality filtering) will behave oddly on this sample for the same
reason.

## Gotcha: the CLI caches config in the OUTPUT directory

`tcga_cptac_taxonomic_profiler <stage> --help` documents
`--configfile [default: (outputDir)/config.yaml]`. On the first run of a stage
the merged config is written into the **output** directory, and every later run
against that same output directory reuses it — the repo's `config/config.yaml`
is ignored.

This silently defeats config edits. After switching SingleM to the R232
metapackage, the repo config was correct, the dry-run still emitted the R226
path, and every job died with:

    ZenodoBackpackMalformedException: Version in CONTENTS.json: 5.4.0
    does not match version provided: 6.5.0

The stale `v2_singlem/config.yaml` still held the R226 path. Deleting the output
directory is what makes a config change take effect:

```bash
rm -rf <outputDir>            # or just <outputDir>/config.yaml
```

Then confirm the change actually reached the rule before submitting:

```bash
tcga_cptac_taxonomic_profiler singlem --input <dir> --output <dir> --dry-run \
  | grep -oE "SINGLEM_METAPACKAGE_PATH=[^ ]*" | sort -u
```

The same applies to any database or parameter change for the kraken, metabuli
and binning stages. Two related traps already cost a submission each:

- **Pull on the HPC before submitting.** A config change committed and pushed
  locally does nothing until the cluster checkout is updated.
- **Conda env hashes drift.** Environments present in `workflow/conda` from the
  original run no longer matched the hashes snakemake computes, so kraken tried
  to build them on a compute node with no network. Always run
  `--conda-create-envs-only` on the login node first, even for stages that have
  run before.

## Phoenix: compute nodes have no internet

Compute nodes cannot reach the network, so anything that downloads must happen
on the **login node** first. Two consequences:

**1. Pre-build every conda environment on the login node.** With `--use-conda`,
snakemake materialises environments at job runtime — which on a compute node
will fail with no network. Build them ahead of time:

```bash
tcga_cptac_taxonomic_profiler host-removal --input <FASTQ_DIR> --output v2.out/hostrm --conda-create-envs-only
```

This populates `workflow/conda/` (already ~42 GB from the original run) with the
new `deacon` environment. Only then submit the real job. The same applies to any
stage whose environment has not been built before.

**2. The index must be fetched on the login node.** The `get_deacon_index` rule
added to `run_install_host.smk` calls `deacon index fetch`, which needs network.
Run `install-host` on the login node, or download the index directly:

```bash
wget -c -O /hpcfs/users/a1667917/metagenomics_dbs/host_genome_db/panhuman-1.k31w15.idx \
  "https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/deacon/3/panhuman-1.k31w15.idx"
```

Expected size: 3,279,263,603 bytes. Verify before use — a truncated index will
silently under-deplete rather than error.

## Cluster configuration has rotted since the original run

Two things in the older `run_*.sh` scripts no longer work:

| Reference | Status | Use instead |
|---|---|---|
| `#SBATCH -p skylake` | Partition no longer exists | `icelake` (93 nodes), or `sacgf` |
| `/hpcfs/users/a1667917/snakemake_slurm_profile/bact_assembly` | Directory does not exist | `<repo>/tcga_cptac_taxonomic_profiler/profile/slurm` (shipped in-repo) |

`run_host_removal.sh` in the analysis directory is written against the working
values. The in-repo profile sets no partition, so snakemake-submitted jobs take
the cluster default; only the outer submit script needs `-p icelake`.

## Three classifiers, three different failure modes

Metabuli is added alongside Kraken2 rather than replacing it. The point is not
to pick a winner but to report where two independent methods agree, which is
what the 2026 *Cancer Cell* consensus asks for and what the retracted
literature never did.

| Tool | Basis | Fails when |
|---|---|---|
| Kraken2 | exact DNA k-mer matching against a reference | reads are short or the true organism is divergent from anything in the database; spurious exact matches from unremoved host |
| Metabuli | 'metamers' spanning amino acid **and** DNA space — AA conservation for sensitivity, DNA substitutions to separate close relatives | shares database bias with Kraken2, but not the exact-DNA-match failure mode |
| SingleM | universal single-copy marker genes, run on the **raw** reads | low-biomass samples where no marker read is captured; cannot see non-marker taxa |

This matters most for the 51 bp samples — 60% of the primary tumours — where
exact DNA k-mer matching is weakest and where residual host is highest. That is
precisely the regime that produced the false positives behind the retractions.

**Metabuli also supersedes the MMseqs2 easy-taxonomy stage.** That step was the
existing amino-acid-based classifier; Metabuli does the same job purpose-built
and far faster, and unlike MMseqs2 it emits a Kraken2-format report, so it feeds
the existing `kraken-biom` path unchanged. Keep the MMseqs2 stage only if a
reviewer asks for it.

**Metabuli gives per-read assignments, so it takes the full decontamination.**
This is the key practical difference from SingleM. `{sample}_classifications.tsv`
is one line per read with eight columns — `is_classified`, read name, `taxID`,
`query_length`, `score`, `e_value`, `rank`, and `taxID:match_count` pairs — a
superset of Kraken2's per-read output. The Methods currently note that SingleM
could not take the rigorous decontamination *because* it emits no per-read rank;
that objection does not apply to Metabuli.

So Metabuli goes through the identical cascade — prevalence prefilter, then
ALDEx2 centre filter — and can be compared to Kraken2 like-for-like at every
stage: raw counts, post-prefilter, post-decontamination, and final bacterial
load. SingleM remains the orthogonal check on microbial *fraction*; Metabuli is
the like-for-like check on *composition*. The report is Kraken2-format, so
`kraken-biom` consumes it unchanged, and `metabuli extract` replaces
`Kraken_Tools/extract_kraken_reads.py` for pulling reads by taxon.

The per-read `score` and `e_value` also allow confidence thresholding via
`--min-score` (currently 0), a finer-grained equivalent of Kraken2's
`--confidence 0.15`. Worth running as a sensitivity analysis rather than fixing
one value.

**What to report.** Per-sample bacterial load from both classifiers, their
concordance (as done for SingleM), and — most usefully — the taxa called by one
and not the other. A genus that only Kraken2 sees, in only the 51 bp samples,
is a contamination signature. A genus both call, in both read lengths, and that
also assembles into a MAG, is real.

## Host depletion: read-length dependence (336/336 complete)

726 GB in, 261 GB out.

**Depletion is strongly read-length dependent, and not in the direction
predicted** — longer reads are depleted substantially harder. Per-sample and
per-stratum figures are held in the private analysis repository; what belongs
here is the mechanism, because it is a property of the tool and will recur on
any mixed-read-length cohort.

The Ge et al. data predicted the opposite — that 51 bp samples carried *more*
recoverable human. Two lines of evidence say this is a property of Deacon
rather than of the samples:

1. **Both read-length groups start with comparable unmapped fractions**, and
   the 101 bp group — which began with slightly *more* unmapped reads — lost
   far more of them. So the difference is introduced by depletion, not inherited
   from the input.
2. **Minimizer count scales with read length.** At k=31, w=15 a 51 bp read
   carries ~1.4 minimizers and a 101 bp read ~4.7 — a 3.4x difference in
   opportunities to match the index. With `abs_threshold: 1`, a longer read is
   far more likely to hit the threshold regardless of its true origin.

**Consequence: the 51 bp samples are systematically under-depleted**, and they
are 60% of the primary tumours. Residual human is therefore likely *higher* in
exactly the samples where Kraken2's exact k-mer matching is weakest.

**The decisive test is already available.** Both the Kraken2 `k2_pluspf`
database and the Metabuli `gtdb+human+virus` database contain human references,
so each will report directly how many post-depletion reads are still human.
Tabulate that by read length before interpreting any load figure — it converts
this from inference to measurement. Until then, treat read length as a live
confounder of everything downstream, and re-check that the load effect remains
independent of it — as it was on the pre-Deacon data. Those models and their
estimates live in the private analysis repository.

The `abs_threshold: 2` sensitivity run is less informative than planned: it
would deplete less everywhere without addressing the read-length differential.
A read-length-aware threshold, or explicit adjustment for read length in the
analysis, is the better response.

## Planned but NOT implemented: Phables (phage resolution)

Recorded as a plan only — deliberately not wired in.

**What it would add.** Phables (v1.5.0, bioconda) resolves phage genomes from
assembly *graphs* rather than contigs, recovering complete phage genomes that
contig-based tools miss. The pipeline already has geNomad (MGE) and PhiSpy
(prophage), so this extends rather than duplicates: geNomad classifies contigs,
PhiSpy finds prophage in bacterial genomes, Phables reconstructs free phage
genomes from graph bubbles.

**Prerequisite that must be fixed first.** `sample_assembly.smk` declares only
`contigs.fasta` as an output. metaSPAdes *writes*
`assembly_graph_with_scaffolds.gfa` into the assembly directory, but because it
is undeclared, snakemake does not track it and it cannot be consumed by a
downstream rule. Phables requires that GFA. So the assembly rule needs the graph
added as a declared output **before** any Phables work — and that change alone
is worth making, since the graph is currently produced and discarded.

**Scope judgement.** Phage detection is outside the current scope. It opens a
new results section and a new set of reviewer questions on top of the two
existing blockers (decontamination, MAGs). Recommend running it only
after the MAG chain completes, then deciding whether it belongs in this paper or
a follow-up.

**Sketch of the implementation, when wanted:**

1. Add `gfa = os.path.join(SAMPLE_ASSEMBLIES, '{sample}', 'assembly_graph_with_scaffolds.gfa')`
   to `individual_sample_assembly` outputs.
2. `workflow/envs/phables.yaml` — `phables >=1.5.0`.
3. `workflow/rules/phables.smk` — per-sample `phables run` against the GFA,
   with its own database (Phables ships a `phables install` step, which is a
   download and therefore login-node only).
4. `workflow/run_phables.smk` entry point + `phables` CLI subcommand,
   mirroring the sylph/metabuli stages.
5. Pre-build the env on the login node before submitting.

## Traps that cost a submission on the full run

Both of these passed every test on the 30-sample subset and failed on the full
336-sample run. Neither is specific to this project.

**Pre-build environments against the target set you will actually run.**
`--conda-create-envs-only` only creates environments for rules in the DAG it is
given. The subset runs had targeted individual flags and so never reached the
CheckM2 or GTDB-Tk rules; their environments were never created. The full run's
driver then reached CheckM2 on a compute node and died:

    CondaHTTPError: HTTP 000 CONNECTION FAILED for repo.anaconda.com
    Creating conda environment .../envs/checkm2.yaml

Compute nodes have no internet. Worse, SLURM recorded the driver as
`COMPLETED 0:0` — the scheduler's state said nothing was wrong, and only the
job's stderr showed the failure.

**Declare a rule's real product as an output, not just its flag.** A rule that
ends in `touch flag` and declares only that flag hides its actual outputs from
the DAG. When a downstream rule takes one of those files as input, snakemake
sees a file no rule produces and the entire DAG build fails:

    MissingInputException in rule stage_hq_med_mags
        affected files: .../final_bins_quality_reports.tsv

This is invisible while a previous run has left the file on disk, which is why
it survived every subset test. It appears only against a clean output directory
— that is, on the real run.

**Conda resolves GPU packages against the BUILD host.** An environment needing
CUDA pytorch, solved on a login node whose driver is older than the compute
nodes', silently gets a CPU build: no error, no warning, just a job that is
inexplicably slow later. Environments cannot be built on compute nodes (no
internet), so set `CONDA_OVERRIDE_CUDA` to the version the job's hardware
supports, and verify afterwards rather than assuming:

    <env>/bin/python -c "import torch; print(torch.version.cuda, torch.cuda.is_available())"

Related: on this cluster `-p a100cpu` is the CPU-only side of the A100 nodes
(its QOS sets `gres/gpu=0`), so a GPU request there is rejected with
`QOSMaxGRESPerNode`. GPU work belongs on `a100`, `newa100` or `singlegpu`.

## Stage sequence

Each stage takes the previous stage's output directory as `--input`. Add
`--profile tcga_cptac_taxonomic_profiler/profile/slurm` for cluster execution,
and `--dry-run` first on every stage.

Stage 1 (`extract`) is skipped: the BAMs are gone and their output is already
held. `<FASTQ_DIR>` below is the directory of existing
`{sample}_R{1,2}.fastq.gz` files.

```bash
# 2. Host depletion (NEW - this is the stage that was previously manual)
tcga_cptac_taxonomic_profiler host-removal --input <FASTQ_DIR> --output v2.out/hostrm
```

```bash
# 3. Kraken2 + Bracken + biom
tcga_cptac_taxonomic_profiler kraken --input v2.out/hostrm/RESULTS/HOST_RM_FASTQ --output v2.out/kraken
```

```bash
# 3b. Metabuli - same input as kraken, algorithmically independent
tcga_cptac_taxonomic_profiler metabuli --input v2.out/hostrm/RESULTS/HOST_RM_FASTQ --output v2.out/metabuli
```

```bash
# 4. SingleM - note: reads the RAW fastqs, not the depleted ones
tcga_cptac_taxonomic_profiler singlem --input <FASTQ_DIR> --output v2.out/singlem
```

```bash
# 5. Per-sample assembly
tcga_cptac_taxonomic_profiler assembly --input v2.out/hostrm/RESULTS/HOST_RM_FASTQ --output v2.out/assembly
```

```bash
# 6. Binning, CheckM2, GTDB-Tk
tcga_cptac_taxonomic_profiler binning --input v2.out/hostrm/RESULTS/HOST_RM_FASTQ --output v2.out/binning
```

```bash
# 7. MAG annotation (Bakta, PhiSpy)
tcga_cptac_taxonomic_profiler annotate --input v2.out/binning --output v2.out/annotate
```

## Sensitivity analysis to run alongside

The depletion threshold is a judgement call, so run it both ways on the full
cohort and report the difference rather than defending one setting:

```bash
tcga_cptac_taxonomic_profiler host-removal --input <FASTQ_DIR> --output v2.out/hostrm_a2 --config deacon.abs_threshold=2
```

`abs_threshold: 1` (the config default here) is more aggressive than Deacon's
own default of 2. For tumour WGS a retained human read costs more than a
discarded microbial one, but the difference between the two settings is a
number the Methods should contain.

## Stratify everything by read length

The cohort is a roughly 60/40 mix of 51 bp and 101 bp reads (94 vs 63 primary
tumours), and read length strongly determines how much residual human sequence
survives host depletion. In Ge et al.'s HNSC data, moving from GRCh38 to
GRCh38+CHM13 reclaims a median **21.7%** of unmapped reads as human in 51 bp
samples but only **0.066%** in 101 bp samples — a near-perfect split
(221/228 high-gain samples are 51 bp; all 107 low-gain samples are 101 bp).

Short reads carry far less information, so they fail to place uniquely against
a reference far more often. Since the input here is hg19-unmapped — a worse
reference still — the 51 bp majority is where Deacon will do most of its work.

**Quantified on a real sample.** For `TCGA-IQ-7632-01A-11D-2317` (51 bp), the
hg19-unmapped FASTQ holds 7,221,628 pairs. Ge et al. report 4,794,777 reads
unmapped against GRCh38 for the same aliquot, and 3,028,240 against
GRCh38+CHM13. So of the reads that actually entered Kraken2:

| Reference | Reads surviving | % of the input used |
|---|---|---|
| hg19 (what was used) | 7,221,628 | 100% |
| GRCh38 | 4,794,777 | 66.4% |
| GRCh38 + CHM13 | 3,028,240 | 41.9% |

**58% of that sample's classifier input is human sequence a modern reference
would have removed.** This is a 51 bp sample, so it is the bad end of the range
rather than typical — but 51 bp samples are 60% of the primary tumours. Run
`scripts/count_input_reads.sh` across the cohort and tabulate this properly
before and after depletion; it is both a QC figure and a Methods table.

Four checks were run on the pre-Deacon data and all were reassuring; all are
worth repeating after the re-run. They are, in order: whether raw microbial
fraction differs by read length; whether any difference survives
decontamination; whether read length predicts survival on its own; and whether
the load effect holds once adjusted for read length. The estimates are held in
the private analysis repository.

Deacon buys cleaner input, recorded parameters and a per-sample QC metric rather
than a rescue — but the stratified comparison belongs in the paper either way.

## To capture for the Methods

Collect these as the run proceeds — they are the reporting items a
post-retraction reviewer will look for, and they are far easier to record now
than to reconstruct later:

- Residual human fraction per sample, from `HOST_RM_FASTQ/{sample}.deacon.json`,
  reported separately for 51 bp and 101 bp samples
- Reads surviving depletion per sample, at both thresholds
- Kraken2 database name and build date
- Deacon version and index version (`panhuman-1`, April 2025)
- SingleM metapackage version (`S4.3.0.GTDB_r220.metapackage_20240523`)
- GTDB release used for MAG taxonomy (currently r214)
- Per-sample MAG counts with CheckM2 completeness/contamination

## After the pipeline completes

The R analysis is rebuilt on the new outputs in the separate paper repository.
Do not carry forward any `.rds` intermediate from the previous run: the
decontamination thresholds, the taxa retained, and every bacterial load value
change with the new depletion. The five defects recorded in the
`archive/pre-revival-2026-08` tag of `TCGA_Post_Analysis` are fixed during that
rebuild, not ported.
