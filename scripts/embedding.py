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
adata.uns['log1p']['base'] = None

#embed the graph via UMAP
sc.tl.umap(adata)

#perform Leiden clustering
sc.tl.leiden(adata, resolution=float(args.resolution))

#Strip .h5ad off outfile to retrieve name of channel
channel_name = "_" + results_file.split(".")[0].split("/")[-1] + ".png"

#set plotting parameters and plot UMAP with leiden clusters labeled
sc.settings.figdir = "."
sc.settings.set_figure_params(format="png")
sc.pl.umap(adata, color='leiden', save=channel_name)

adata.write(results_file)
