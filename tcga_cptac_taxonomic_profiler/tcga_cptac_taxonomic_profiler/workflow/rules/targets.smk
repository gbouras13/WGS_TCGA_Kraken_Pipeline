"""
All target output files are declared here
"""


# Preprocessing files
KrakenTargets = [
    os.path.join(FLAGS, "aggr_kraken.flag"),
    os.path.join(FLAGS, "aggr_biom.flag"),
    os.path.join(FLAGS, "aggr_bracken.flag"),
]

# Targets for extract
ExtractFiles = [
    os.path.join(FLAGS, "aggr_fastq.flag"),
    os.path.join(FLAGS, "aggr_read_count.flag")
]

MMseqsTargets =[

    os.path.join(FLAGS, "aggr_mmseqs2.flag")    

]


CoAssemblyTargets = [
    os.path.join(FLAGS, "aggr_coassembly.flag")
]

SampleAssemblyTargets = [
    os.path.join(FLAGS, "aggr_sample_assembly.flag")
]


SampleAssemblyBins = [

    os.path.join(FLAGS, 'vamb.flag'),
    os.path.join(FLAGS, 'checkm2.flag'),
    os.path.join(CHECKM2_RESULTS, "combined_check2_quality_report_hq.tsv"),
    os.path.join(FLAGS, 'gtdb.flag')
]

AnnotationTargets = [os.path.join(FLAGS, 'bakta.flag')]

MGETargets = [os.path.join(FLAGS, 'genomad.flag')]