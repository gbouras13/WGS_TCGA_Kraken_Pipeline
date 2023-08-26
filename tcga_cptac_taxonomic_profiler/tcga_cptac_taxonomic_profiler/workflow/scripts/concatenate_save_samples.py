#!/usr/bin/env python


from Bio import SeqIO
import gzip as gz
import pandas as pd

def concatenate_contigs(catalogue, fastas, samples, csv_path, min_length, separator ):
    # opens the catalogue
    headers = []
    
    with gz.open(catalogue, "wt") as fout:

        i = 1
        
        for sample, input_fasta in zip(samples, fastas):

            for record in SeqIO.parse(input_fasta, 'fasta'):
                if len(record.seq) > min_length:
                    header = f"S{i}{separator}{record.id}"
                    record.id = header
                    record.description = ""
                    SeqIO.write(record, fout, 'fasta')

            headers.append(f"S{i}")
            i += 1

    
    sample_vamb_df = pd.DataFrame(
        {
            "sample" : samples,
            "Binner_ID" : headers
    }
    )

    sample_vamb_df.to_csv(csv_path, sep = "\t", index = False, header = True)




# to actually run the script
concatenate_contigs(snakemake.output.catalogue, snakemake.input.fastas, snakemake.params.samples, snakemake.output.csv_path, snakemake.params.min_contig_length, snakemake.params.separator)


