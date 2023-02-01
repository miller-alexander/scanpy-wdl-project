#!/usr/bin/python

import argparse
import pandas as pd
import numpy as np
import scanpy as sc

sc.settings.verbosity = 3

parser = argparse.ArgumentParser(description="Compute PCA decomposition")

parser.add_argument("-i", "--infile", help="Input file path")
parser.add_argument("-o", "--outfile", help="Output file name")
parser.add_argument("-c", "--components", help="Number of components used for PCA")

args = parser.parse_args()

results_file = str(args.outfile)

adata = sc.read_h5ad(str(args.infile))

#calculate PCA
sc.tl.pca(adata, n_comps=args.components, svd_solver='arpack')

adata.write(results_file)