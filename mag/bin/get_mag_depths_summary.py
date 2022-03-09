#!/usr/bin/env python

import sys
import argparse
import pandas as pd


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--depths'       , required=True, nargs="+", metavar='FILE'                             , help="TSV file for each assembly containing bin depths for samples: bin, sample1, ....")
    parser.add_argument('-o', "--out"          , required=True           , metavar='FILE', type=argparse.FileType('w'), help="Output file containing depths for all assemblies and all samples.")
    return parser.parse_args(args)

def main(args=None):
    args = parse_args(args)

    results = pd.DataFrame()
    for assembly_depths_file in args.depths:
        assembly_results = pd.read_csv(assembly_depths_file, index_col="bin", sep="\t")
        results          = results.append(assembly_results, sort=True, verify_integrity=True)

    results.to_csv(args.out, sep='\t')

if __name__ == "__main__":
    sys.exit(main())
