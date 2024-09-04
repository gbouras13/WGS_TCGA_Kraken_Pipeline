#!/usr/bin/env python

"""
taken and adapted from
https://github.com/ac1513/MAGqual/blob/main/workflow/scripts/python/qual_parse.py
"""

import pandas as pd
import os

def parse_mags(samples, bakta_dir, outfile):
    # opens the catalogue

    fiveS_dict = {}
    trna_dict = {}
    sixteenS_dict = {}
    twentythreeS_dict = {}

    print(samples)
    for mag in samples:
        os.path.join(bakta_dir, mag, f'{mag}.tsv') 

        bakta_file = os.path.join(bakta_dir, mag, f'{mag}.tsv')  

        with open(bakta_file, 'r') as bakta_in:
            trna_set = set()
            fiveS = False
            sixteenS = False
            twentythreeS = False
            for line in bakta_in:
                if "tRNA-Ala" in line:
                    trna_set.add("tRNA-Ala")
                if "tRNA-Arg" in line:
                    trna_set.add("tRNA-Arg")
                if "tRNA-Asn" in line:
                    trna_set.add("tRNA-Asn")
                if "tRNA-Asp" in line:
                    trna_set.add("tRNA-Asp")
                if "tRNA-Cys" in line:
                    trna_set.add("tRNA-Cys")
                if "tRNA-Gln" in line:
                    trna_set.add("tRNA-Gln")
                if "tRNA-Glu" in line:
                    trna_set.add("tRNA-Glu")
                if "tRNA-Gly" in line:
                    trna_set.add("tRNA-Gly")
                if "tRNA-His" in line:
                    trna_set.add("tRNA-His")
                if "tRNA-Ile" in line:
                    trna_set.add("tRNA-Ile")
                if "tRNA-Leu" in line:
                    trna_set.add("tRNA-Leu")
                if "tRNA-Lys" in line:
                    trna_set.add("tRNA-Lys")
                if "tRNA-Met" in line:
                    trna_set.add("tRNA-Met")
                if "tRNA-Phe" in line:
                    trna_set.add("tRNA-Phe")
                if "tRNA-Pro" in line:
                    trna_set.add("tRNA-Pro")
                if "tRNA-Ser" in line:
                    trna_set.add("tRNA-Ser")
                if "tRNA-Thr" in line:
                    trna_set.add("tRNA-Thr")
                if "tRNA-Trp" in line:
                    trna_set.add("tRNA-Trp")
                if "tRNA-Tyr" in line:
                    trna_set.add("tRNA-Tyr")
                if "tRNA-Val" in line:
                    trna_set.add("tRNA-Val")
                if ("5S ribosomal RNA" in line) and ("partial" not in line):
                    fiveS = True
                if ("16S ribosomal RNA" in line) and ("partial" not in line):
                    sixteenS = True
                if ("23S ribosomal RNA" in line) and ("partial" not in line):
                    twentythreeS = True
            #trna count
            trna_count = len(trna_set)

            fiveS_dict[mag] = fiveS
            sixteenS_dict[mag] = sixteenS
            twentythreeS_dict[mag] = twentythreeS
            trna_dict[mag] = trna_count

    print(fiveS_dict)
    combined_df = pd.DataFrame({
    'sample': list(fiveS_dict.keys()),
    'fiveS': list(fiveS_dict.values()),
    'sixteenS': list(sixteenS_dict.values()),
    'twentythreeS': list(twentythreeS_dict.values()),
    'tRNA_count': list(trna_dict.values())
})
    combined_df.to_csv(
            outfile,
            sep="\t",
            index=False,
        )




# to actually run the script
parse_mags(snakemake.params.samples, snakemake.params.bakta_dir, snakemake.output.out_tsv)
