#!/usr/bin/env python

import pandas as pd
import os
import shutil

def parse_mags(samples, checkm2_directory, vamb_bin_dir, combined_mag_directory):
    # opens the catalogue

    high_qual_dfs =[]
    
    for sample in samples:
        # Read the TSV file
        file_path = os.path.join(checkm2_directory, sample, "quality_report.tsv")
        data = pd.read_csv(file_path, sep='\t')

        # Filter rows based on criteria
        high_qual_df = data[(data['Completeness'] > 90) & (data['Contamination'] < 5)]
        data['Sample'] = sample
        high_qual_dfs.append(high_qual_df)

        # copy bin to combined directory

        for name in high_qual_df['Name']:

            input_fasta = os.path.join(vamb_bin_dir, sample, f"{name}.fna")
            out_fasta = os.path.join(combined_mag_directory, f"{sample}_{name}.fna")

            shutil.copy(input_fasta, out_fasta)

    combined_df = pd.concat(high_qual_dfs, ignore_index=True)

    combined_df.to_csv(
            os.path.join(checkm2_directory, "combined_check2_quality_report.tsv"),
            sep="\t",
            index=False,
        )


# to actually run the script
parse_mags(snakemake.params.samples, snakemake.params.checkm2_directory, snakemake.params.vamb_bin_dir, snakemake.params.combined_mag_directory)
