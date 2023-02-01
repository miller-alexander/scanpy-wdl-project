#!/usr/bin/python

import argparse
import pandas as pd
import numpy as np
import scanpy as sc

sc.settings.verbosity = 3

parser = argparse.ArgumentParser(description="Select the top variable genes")

parser.add_argument("-i", "--infile", help="Input file path")
parser.add_argument("-o", "--outfile", help="Output file name")
parser.add_argument("-n", "--ntopgenes", help="Number of highly variable genes")

args = parser.parse_args()

results_file = str(args.outfile)

adata = sc.read_h5ad(str(args.infile))

#Calculate the highly variable genes
sc.pp.highly_variable_genes(adata, n_top_genes=2000)

#Snapshot the log-normalized expression in the .raw attribute and filter
adata.raw = adata
adata = adata[:, adata.var.highly_variable]

#Regression and scaling
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)

adata.write(results_file)