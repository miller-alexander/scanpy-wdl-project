#!/usr/bin/python

import argparse
import pandas as pd
import numpy as np
import scanpy as sc

sc.settings.verbosity = 3

parser = argparse.ArgumentParser(description="Perform Leiden clustering and embed via UMAP")

parser.add_argument("-i", "--infile", help="Input file path")
parser.add_argument("-o", "--outfile", help="Output file name")

args = parser.parse_args()

results_file = str(args.outfile)

adata = sc.read_h5ad(str(args.infile))

#statistically rank genes via Mann-Whitney U-test
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save=True)

adata.write(results_file)