/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Check already if long reads are provided
include { hasExtension } from '../modules/local/functions'
def hybrid = false
if(hasExtension(params.input, "csv")){
    Channel
        .from(file(params.input))
        .splitCsv(header: true)
        .map { row ->
                if (row.size() == 5) {
                    if (row.long_reads) hybrid = true
                } else {
                    log.error "Input samplesheet contains row with ${row.size()} column(s). Expects 5."
                    System.exit(1)
                }
            }
}

// Validate input parameters
WorkflowMag.initialise(params, log, hybrid)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.phix_reference, params.host_fasta, params.centrifuge_db, params.kraken2_db, params.cat_db, params.gtdb, params.lambda_reference, params.busco_reference ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()
def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS                               } from '../modules/local/get_software_versions'       addParams( options: [publish_files : ['csv':'']]          )
include { BOWTIE2_REMOVAL_BUILD as BOWTIE2_HOST_REMOVAL_BUILD } from '../modules/local/bowtie2_removal_build'
include { BOWTIE2_REMOVAL_ALIGN as BOWTIE2_HOST_REMOVAL_ALIGN } from '../modules/local/bowtie2_removal_align'       addParams( options: modules['bowtie2_host_removal_align'] )
include { BOWTIE2_REMOVAL_BUILD as BOWTIE2_PHIX_REMOVAL_BUILD } from '../modules/local/bowtie2_removal_build'
include { BOWTIE2_REMOVAL_ALIGN as BOWTIE2_PHIX_REMOVAL_ALIGN } from '../modules/local/bowtie2_removal_align'       addParams( options: modules['bowtie2_phix_removal_align'] )
include { PORECHOP                                            } from '../modules/local/porechop'
include { NANOLYSE                                            } from '../modules/local/nanolyse'                    addParams( options: modules['nanolyse']                   )
include { FILTLONG                                            } from '../modules/local/filtlong'
include { NANOPLOT as NANOPLOT_RAW                            } from '../modules/local/nanoplot'                    addParams( options: modules['nanoplot_raw']               )
include { NANOPLOT as NANOPLOT_FILTERED                       } from '../modules/local/nanoplot'                    addParams( options: modules['nanoplot_filtered']          )
include { CENTRIFUGE_DB_PREPARATION                           } from '../modules/local/centrifuge_db_preparation'
include { CENTRIFUGE                                          } from '../modules/local/centrifuge'                  addParams( options: modules['centrifuge']                 )
include { KRAKEN2_DB_PREPARATION                              } from '../modules/local/kraken2_db_preparation'
include { KRAKEN2                                             } from '../modules/local/kraken2'                     addParams( options: modules['kraken2']                    )
include { KRONA_DB                                            } from '../modules/local/krona_db'
include { KRONA                                               } from '../modules/local/krona'                       addParams( options: modules['krona']                      )
include { POOL_SINGLE_READS                                   } from '../modules/local/pool_single_reads'
include { POOL_PAIRED_READS                                   } from '../modules/local/pool_paired_reads'
include { POOL_SINGLE_READS as POOL_LONG_READS                } from '../modules/local/pool_single_reads'
include { MEGAHIT                                             } from '../modules/local/megahit'                     addParams( options: modules['megahit']                    )
include { SPADES                                              } from '../modules/local/spades'                      addParams( options: modules['spades']                     )
include { SPADESHYBRID                                        } from '../modules/local/spadeshybrid'                addParams( options: modules['spadeshybrid']               )
include { QUAST                                               } from '../modules/local/quast'                       addParams( options: modules['quast']                      )
include { QUAST_BINS                                          } from '../modules/local/quast_bins'                  addParams( options: modules['quast_bins']                 )
include { QUAST_BINS_SUMMARY                                  } from '../modules/local/quast_bins_summary'          addParams( options: modules['quast_bins_summary']         )
include { CAT_DB                                              } from '../modules/local/cat_db'
include { CAT_DB_GENERATE                                     } from '../modules/local/cat_db_generate'             addParams( options: modules['cat_db_generate']            )
include { CAT                                                 } from '../modules/local/cat'                         addParams( options: modules['cat']                        )
include { BIN_SUMMARY                                         } from '../modules/local/bin_summary'                 addParams( options: modules['bin_summary']                )
include { MULTIQC                                             } from '../modules/local/multiqc'                     addParams( options: multiqc_options                       )

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK         } from '../subworkflows/local/input_check'
include { METABAT2_BINNING    } from '../subworkflows/local/metabat2_binning'      addParams( bowtie2_align_options: modules['bowtie2_assembly_align'], metabat2_options: modules['metabat2'], mag_depths_options: modules['mag_depths'], mag_depths_plot_options: modules['mag_depths_plot'], mag_depths_summary_options: modules['mag_depths_summary'])
include { BUSCO_QC            } from '../subworkflows/local/busco_qc'              addParams( busco_db_options: modules['busco_db_preparation'], busco_options: modules['busco'], busco_save_download_options: modules['busco_save_download'], busco_plot_options: modules['busco_plot'], busco_summary_options: modules['busco_summary'])
include { GTDBTK              } from '../subworkflows/local/gtdbtk'                addParams( gtdbtk_classify_options: modules['gtdbtk_classify'], gtdbtk_summary_options: modules['gtdbtk_summary'])

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC as FASTQC_RAW     } from '../modules/nf-core/modules/fastqc/main'              addParams( options: modules['fastqc_raw']            )
include { FASTQC as FASTQC_TRIMMED } from '../modules/nf-core/modules/fastqc/main'              addParams( options: modules['fastqc_trimmed']        )
include { FASTP                    } from '../modules/nf-core/modules/fastp/main'               addParams( options: modules['fastp']                 )
include { PRODIGAL                 } from '../modules/nf-core/modules/prodigal/main'            addParams( options: modules['prodigal']              )
include { PROKKA                   } from '../modules/nf-core/modules/prokka/main'              addParams( options: modules['prokka']                )

////////////////////////////////////////////////////
/* --  Create channel for reference databases  -- */
////////////////////////////////////////////////////

if ( params.host_genome ) {
    host_fasta = params.genomes[params.host_genome].fasta ?: false
    Channel
        .value(file( "${host_fasta}" ))
        .set { ch_host_fasta }

    host_bowtie2index = params.genomes[params.host_genome].bowtie2 ?: false
    Channel
        .value(file( "${host_bowtie2index}/*" ))
        .set { ch_host_bowtie2index }
} else if ( params.host_fasta ) {
    Channel
        .value(file( "${params.host_fasta}" ))
        .set { ch_host_fasta }
} else {
    ch_host_fasta = Channel.empty()
}

if(params.busco_reference){
    Channel
        .value(file( "${params.busco_reference}" ))
        .set { ch_busco_db_file }
} else {
    ch_busco_db_file = Channel.empty()
}
if (params.busco_download_path) {
    Channel
        .value(file( "${params.busco_download_path}" ))
        .set { ch_busco_download_folder }
} else {
    ch_busco_download_folder = Channel.empty()
}

if(params.centrifuge_db){
    Channel
        .value(file( "${params.centrifuge_db}" ))
        .set { ch_centrifuge_db_file }
} else {
    ch_centrifuge_db_file = Channel.empty()
}

if(params.kraken2_db){
    Channel
        .value(file( "${params.kraken2_db}" ))
        .set { ch_kraken2_db_file }
} else {
    ch_kraken2_db_file = Channel.empty()
}

if(params.cat_db){
    Channel
        .value(file( "${params.cat_db}" ))
        .set { ch_cat_db_file }
} else {
    ch_cat_db_file = Channel.empty()
}

if(!params.keep_phix) {
    Channel
        .value(file( "${params.phix_reference}" ))
        .set { ch_phix_db_file }
}

if (!params.keep_lambda) {
    Channel
        .value(file( "${params.lambda_reference}" ))
        .set { ch_nanolyse_db }
}

gtdb = params.skip_busco ? false : params.gtdb
if (gtdb) {
    Channel
        .value(file( "${gtdb}" ))
        .set { ch_gtdb }
} else {
    ch_gtdb = Channel.empty()
}

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report    = []
def busco_failed_bins = [:]

workflow MAG {

    ch_software_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK ()
    ch_raw_short_reads = INPUT_CHECK.out.raw_short_reads
    ch_raw_long_reads  = INPUT_CHECK.out.raw_long_reads

    /*
    ================================================================================
                                    Preprocessing and QC for short reads
    ================================================================================
    */

    FASTQC_RAW (
        ch_raw_short_reads
    )
    ch_software_versions = ch_software_versions.mix(FASTQC_RAW.out.version.first().ifEmpty(null))

    FASTP (
        ch_raw_short_reads
    )
    ch_short_reads = FASTP.out.reads
    ch_software_versions = ch_software_versions.mix(FASTP.out.version.first().ifEmpty(null))

    if (params.host_fasta){
        BOWTIE2_HOST_REMOVAL_BUILD (
            ch_host_fasta
        )
        ch_host_bowtie2index = BOWTIE2_HOST_REMOVAL_BUILD.out.index
    }
    ch_bowtie2_removal_host_multiqc = Channel.empty()
    if (params.host_fasta || params.host_genome){
        BOWTIE2_HOST_REMOVAL_ALIGN (
            ch_short_reads,
            ch_host_bowtie2index
        )
        ch_short_reads = BOWTIE2_HOST_REMOVAL_ALIGN.out.reads
        ch_bowtie2_removal_host_multiqc = BOWTIE2_HOST_REMOVAL_ALIGN.out.log
        ch_software_versions = ch_software_versions.mix(BOWTIE2_HOST_REMOVAL_ALIGN.out.version.first().ifEmpty(null))
    }

    if(!params.keep_phix) {
        BOWTIE2_PHIX_REMOVAL_BUILD (
            ch_phix_db_file
        )
        BOWTIE2_PHIX_REMOVAL_ALIGN (
            ch_short_reads,
            BOWTIE2_PHIX_REMOVAL_BUILD.out.index
        )
        ch_short_reads = BOWTIE2_PHIX_REMOVAL_ALIGN.out.reads
    }

    FASTQC_TRIMMED (
        ch_short_reads
    )

    /*
    ================================================================================
                                    Preprocessing and QC for long reads
    ================================================================================
    */
    NANOPLOT_RAW (
        ch_raw_long_reads
    )
    ch_software_versions = ch_software_versions.mix(NANOPLOT_RAW.out.version.first().ifEmpty(null))

    ch_long_reads = ch_raw_long_reads
    if (!params.skip_adapter_trimming) {
        PORECHOP (
            ch_raw_long_reads
        )
        ch_long_reads = PORECHOP.out.reads
        ch_software_versions = ch_software_versions.mix(PORECHOP.out.version.first().ifEmpty(null))
    }

    if (!params.keep_lambda) {
        NANOLYSE (
            ch_long_reads,
            ch_nanolyse_db
        )
        ch_long_reads = NANOLYSE.out.reads
        ch_software_versions = ch_software_versions.mix(NANOLYSE.out.version.first().ifEmpty(null))
    }

    // join long and short reads by sample name
    ch_short_reads
        .map { meta, sr -> [ meta.id, meta, sr ] }
        .set {ch_short_reads_tmp}

    ch_long_reads
        .map { meta, lr -> [ meta.id, meta, lr ] }
        .join(ch_short_reads_tmp, by: 0)
        .map { id, meta_lr, lr, meta_sr, sr -> [ meta_lr, lr, sr[0], sr[1] ] }  // should not occur for single-end, since SPAdes (hybrid) does not support single-end
        .set{ ch_short_and_long_reads }

    FILTLONG (
        ch_short_and_long_reads
    )
    ch_long_reads = FILTLONG.out.reads
    ch_software_versions = ch_software_versions.mix(FILTLONG.out.version.first().ifEmpty(null))

    NANOPLOT_FILTERED (
        ch_long_reads
    )

    /*
    ================================================================================
                                    Taxonomic information
    ================================================================================
    */
    CENTRIFUGE_DB_PREPARATION ( ch_centrifuge_db_file )
    CENTRIFUGE (
        ch_short_reads,
        CENTRIFUGE_DB_PREPARATION.out.db
    )
    ch_software_versions = ch_software_versions.mix(CENTRIFUGE.out.version.first().ifEmpty(null))

    KRAKEN2_DB_PREPARATION (
        ch_kraken2_db_file
    )
    KRAKEN2 (
        ch_short_reads,
        KRAKEN2_DB_PREPARATION.out.db
    )
    ch_software_versions = ch_software_versions.mix(KRAKEN2.out.version.first().ifEmpty(null))

    if (( params.centrifuge_db || params.kraken2_db ) && !params.skip_krona){
        KRONA_DB ()
        CENTRIFUGE.out.results_for_krona.mix(KRAKEN2.out.results_for_krona)
            . map { classifier, meta, report ->
                def meta_new = meta.clone()
                meta_new.classifier  = classifier
                [ meta_new, report ]
            }
            .set { ch_tax_classifications }
        KRONA (
            ch_tax_classifications,
            KRONA_DB.out.db.collect()
        )
        ch_software_versions = ch_software_versions.mix(KRONA.out.version.first().ifEmpty(null))
    }

    /*
    ================================================================================
                                    Assembly
    ================================================================================
    */

    // Co-assembly: prepare grouping for MEGAHIT and for pooling for SPAdes
    if (params.coassemble_group) {
        // short reads
        // group and set group as new id
        ch_short_reads
            .map { meta, reads -> [ meta.group, meta, reads ] }
            .groupTuple(by: 0)
            .map { group, metas, reads ->
                    def meta = [:]
                    meta.id          = "group-$group"
                    meta.group       = group
                    meta.single_end  = params.single_end
                    if (!params.single_end) [ meta, reads.collect { it[0] }, reads.collect { it[1] } ]
                    else [ meta, reads.collect { it }, [] ]
            }
            .set { ch_short_reads_grouped }
        // long reads
        // group and set group as new id
        ch_long_reads
            .map { meta, reads -> [ meta.group, meta, reads ] }
            .groupTuple(by: 0)
            .map { group, metas, reads ->
                def meta = [:]
                meta.id          = "group-$group"
                meta.group       = group
                [ meta, reads.collect { it } ]
            }
            .set { ch_long_reads_grouped }
    } else {
        ch_short_reads
            .map { meta, reads ->
                    if (!params.single_end){ [ meta, [reads[0]], [reads[1]] ] }
                    else [ meta, [reads], [] ] }
            .set { ch_short_reads_grouped }
    }

    ch_assemblies = Channel.empty()
    if (!params.skip_megahit){
        MEGAHIT ( ch_short_reads_grouped )
        MEGAHIT.out.assembly
            .map { meta, assembly ->
                def meta_new = meta.clone()
                meta_new.assembler  = "MEGAHIT"
                [ meta_new, assembly ]
            }
            .set { ch_megahit_assemblies }
        ch_assemblies = ch_assemblies.mix(ch_megahit_assemblies)
        ch_software_versions = ch_software_versions.mix(MEGAHIT.out.version.first().ifEmpty(null))
    }

    // Co-assembly: pool reads for SPAdes
    if (params.coassemble_group) {
        // short reads
        if (!params.single_end && (!params.skip_spades || !params.skip_spadeshybrid)){
            if (params.single_end){
                POOL_SINGLE_READS ( ch_short_reads_grouped )
                ch_short_reads_spades = POOL_SINGLE_READS.out.reads
            } else {
                POOL_PAIRED_READS ( ch_short_reads_grouped )
                ch_short_reads_spades = POOL_PAIRED_READS.out.reads
            }
        }
        // long reads
        if (!params.single_end && !params.skip_spadeshybrid){
            POOL_LONG_READS ( ch_long_reads_grouped )
            ch_long_reads_spades = POOL_LONG_READS.out.reads
        }
    } else {
        ch_short_reads_spades = ch_short_reads
        ch_long_reads
            .map { meta, reads -> [ meta, [reads] ] }
            .set { ch_long_reads_spades }
    }

    if (!params.single_end && !params.skip_spades){
        SPADES ( ch_short_reads_spades )
        SPADES.out.assembly
            .map { meta, assembly ->
                def meta_new = meta.clone()
                meta_new.assembler  = "SPAdes"
                [ meta_new, assembly ]
            }
            .set { ch_spades_assemblies }
        ch_assemblies = ch_assemblies.mix(ch_spades_assemblies)
        ch_software_versions = ch_software_versions.mix(SPADES.out.version.first().ifEmpty(null))
    }

    if (!params.single_end && !params.skip_spadeshybrid){
        ch_short_reads_spades
            .map { meta, reads -> [ meta.id, meta, reads ] }
            .set {ch_short_reads_spades_tmp}
        ch_long_reads_spades
            .map { meta, reads -> [ meta.id, meta, reads ] }
            .combine(ch_short_reads_spades_tmp, by: 0)
            .map { id, meta_long, long_reads, meta_short, short_reads -> [ meta_short, long_reads, short_reads ] }
            .set { ch_reads_spadeshybrid }
        SPADESHYBRID ( ch_reads_spadeshybrid )
        SPADESHYBRID.out.assembly
            .map { meta, assembly ->
                def meta_new = meta.clone()
                meta_new.assembler  = "SPAdesHybrid"
                [ meta_new, assembly ]
            }
            .set { ch_spadeshybrid_assemblies }
        ch_assemblies = ch_assemblies.mix(ch_spadeshybrid_assemblies)
        ch_software_versions = ch_software_versions.mix(SPADESHYBRID.out.version.first().ifEmpty(null))
    }

    ch_quast_multiqc = Channel.empty()
    if (!params.skip_quast){
        QUAST ( ch_assemblies )
        ch_quast_multiqc = QUAST.out.qc
        ch_software_versions = ch_software_versions.mix(QUAST.out.version.first().ifEmpty(null))
    }

    /*
    ================================================================================
                                    Predict proteins
    ================================================================================
    */

    if (!params.skip_prodigal){
        PRODIGAL (
            ch_assemblies,
            modules['prodigal']['output_format']
        )
        ch_software_versions = ch_software_versions.mix(PRODIGAL.out.versions.first().ifEmpty(null))
    }

    /*
    ================================================================================
                                    Binning
    ================================================================================
    */

    ch_bowtie2_assembly_multiqc = Channel.empty()
    ch_busco_summary            = Channel.empty()
    ch_busco_multiqc            = Channel.empty()
    if (!params.skip_binning){

        METABAT2_BINNING (
            ch_assemblies,
            ch_short_reads
        )
        ch_bowtie2_assembly_multiqc = METABAT2_BINNING.out.bowtie2_assembly_multiqc
        ch_software_versions = ch_software_versions.mix(METABAT2_BINNING.out.bowtie2_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(METABAT2_BINNING.out.metabat2_version.first().ifEmpty(null))

        if (!params.skip_busco){
            /*
            * BUSCO subworkflow: Quantitative measures for the assessment of genome assembly
            */
            BUSCO_QC (
                ch_busco_db_file,
                ch_busco_download_folder,
                METABAT2_BINNING.out.bins.transpose()
            )
            ch_busco_summary = BUSCO_QC.out.summary
            ch_busco_multiqc = BUSCO_QC.out.multiqc
            ch_software_versions = ch_software_versions.mix(BUSCO_QC.out.version.first().ifEmpty(null))
            // process information if BUSCO analysis failed for individual bins due to no matching genes
            BUSCO_QC.out
                .failed_bin
                .splitCsv(sep: '\t')
                .map { bin, error -> if (!bin.contains(".unbinned.")) busco_failed_bins[bin] = error }
        }

        ch_quast_bins_summary = Channel.empty()
        if (!params.skip_quast){
            QUAST_BINS ( METABAT2_BINNING.out.bins )
            ch_software_versions = ch_software_versions.mix(QUAST_BINS.out.version.first().ifEmpty(null))
            QUAST_BINS_SUMMARY ( QUAST_BINS.out.quast_bin_summaries.collect() )
            ch_quast_bins_summary = QUAST_BINS_SUMMARY.out.summary
        }

        /*
         * CAT: Bin Annotation Tool (BAT) are pipelines for the taxonomic classification of long DNA sequences and metagenome assembled genomes (MAGs/bins)
         */
        ch_cat_db = Channel.empty()
        if (params.cat_db){
            CAT_DB ( ch_cat_db_file )
            ch_cat_db = CAT_DB.out.db
        } else if (params.cat_db_generate){
            CAT_DB_GENERATE ()
            ch_cat_db = CAT_DB_GENERATE.out.db
        }
        CAT (
            METABAT2_BINNING.out.bins,
            ch_cat_db
        )
        ch_software_versions = ch_software_versions.mix(CAT.out.version.first().ifEmpty(null))

        /*
         * GTDB-tk: taxonomic classifications using GTDB reference
         */
        ch_gtdbtk_summary = Channel.empty()
        if ( gtdb ){
            GTDBTK (
                METABAT2_BINNING.out.bins,
                ch_busco_summary,
                ch_gtdb
            )
            ch_software_versions = ch_software_versions.mix(GTDBTK.out.version.first().ifEmpty(null))
            ch_gtdbtk_summary = GTDBTK.out.summary
        }

        if (!params.skip_busco || !params.skip_quast || gtdb){
            BIN_SUMMARY (
                METABAT2_BINNING.out.depths_summary,
                ch_busco_summary.ifEmpty([]),
                ch_quast_bins_summary.ifEmpty([]),
                ch_gtdbtk_summary.ifEmpty([])
            )
        }

        /*
         * Prokka: Genome annotation
         */
        METABAT2_BINNING.out.bins.transpose()
            .map { meta, bin ->
                def meta_new = meta.clone()
                meta_new.id  = bin.getBaseName()
                [ meta_new, bin ]
            }
            .set { ch_bins_for_prokka }

        if (!params.skip_prokka){
            PROKKA (
                ch_bins_for_prokka,
                [],
                []
            )
            ch_software_versions = ch_software_versions.mix(PROKKA.out.versions.first().ifEmpty(null))
        }
    }

    //
    // MODULE: Pipeline reporting
    //
    ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }

    GET_SOFTWARE_VERSIONS (
        ch_software_versions.map { it }.collect()
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowMag.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_custom_config.collect().ifEmpty([]),
        FASTQC_RAW.out.zip.collect{it[1]}.ifEmpty([]),
        FASTP.out.json.collect{it[1]}.ifEmpty([]),
        FASTQC_TRIMMED.out.zip.collect{it[1]}.ifEmpty([]),
        ch_bowtie2_removal_host_multiqc.collect{it[1]}.ifEmpty([]),
        ch_quast_multiqc.collect().ifEmpty([]),
        ch_bowtie2_assembly_multiqc.collect().ifEmpty([]),
        ch_busco_multiqc.collect().ifEmpty([])
    )
    multiqc_report       = MULTIQC.out.report.toList()
    ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report, busco_failed_bins)
    }
    NfcoreTemplate.summary(workflow, params, log, busco_failed_bins)
}

/*
========================================================================================
    THE END
========================================================================================
*/
