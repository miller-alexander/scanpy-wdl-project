#!/usr/bin/python

import argparse
import pandas as pd
import numpy as np
import scanpy as sc

sc.settings.verbosity = 3

parser = argparse.ArgumentParser(description="Perform Leiden clustering and embed via UMAP")

parser.add_argument("-i", "--infile", help="Input file path")
parser.add_argument("-o", "--outfile", help="Output file name")
parser.add_argument("-r", "--resolution", help="Resolution for the Leiden clustering")

args = parser.parse_args()

results_file = str(args.outfile)

adata = sc.read_h5ad(str(args.infile))

#embed the graph via UMAP
sc.tl.paga(adata)
sc.pl.paga(adata, plot=False)
sc.tl.umap(adata, init_pos='paga')

#perform Leiden clustering
sc.tl.leiden(adata, resolution=args.resolution)

#plot UMAP with leiden clusters labeled
sc.pl.umap(adata, color='leiden', save=True)

adata.write(results_file)
