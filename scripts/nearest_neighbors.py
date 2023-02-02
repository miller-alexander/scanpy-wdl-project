#!/usr/bin/python

import argparse
import pandas as pd
import numpy as np
import scanpy as sc

sc.settings.verbosity = 3

parser = argparse.ArgumentParser(description="Compute nearest neighbors graph")

parser.add_argument("-i", "--infile", help="Input file path")
parser.add_argument("-o", "--outfile", help="Output file name")
parser.add_argument("-c", "--components", help="Number of components used for PCA")
parser.add_argument("-n", "--neighbors", help="Number of neighbors used to construct graph")

args = parser.parse_args()

results_file = str(args.outfile)

adata = sc.read_h5ad(str(args.infile))
adata.uns['log1p']['base'] = None

#compute nearest neighbors graph
sc.pp.neighbors(adata, n_neighbors=int(args.neighbors), n_pcs=int(args.components))

adata.write(results_file)