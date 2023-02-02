#!/usr/bin/python

import argparse
import pandas as pd
import numpy as np
import scanpy as sc

sc.settings.verbosity = 3

parser = argparse.ArgumentParser(description="Apply quality control filters to single cell RNA-sequencing data")
parser.add_argument("-i", "--infile", help="Input file path")
parser.add_argument("-o", "--outfile", help="Output file name")
parser.add_argument("-n", "--mingenes", help="Minimum number of genes to pass cells")
parser.add_argument("-x", "--maxgenes", help="Maximum number of genes to pass cells")
parser.add_argument("-m", "--mtpct", help="Maximum percentage of mitochondrial genes to pass cells")
parser.add_argument("-c", "--cellpct", help="Minimum percentage of cells to pass genes")
parser.add_argument("-r", "--reads", help="Number of reads per cell to normalize the data")

args = parser.parse_args()

results_file = str(args.outfile)

adata = sc.read_10x_h5(str(args.infile))
adata.var_names_make_unique()

#Filter out the cells with below args.mingenes number of genes and above args.maxgenes number of genes
sc.pp.filter_cells(adata, min_genes=int(args.mingenes))
sc.pp.filter_cells(adata, max_genes=int(args.maxgenes))

#Filter out the genes that are not expressed in a sufficient number of cells
num_cells = np.shape(adata)[0]
min_cells = int(args.cellpct)/100 * num_cells
sc.pp.filter_genes(adata, min_cells=min_cells)

#Filter out the cells which contain over args.mtpct percent of genes as mitochondrial
adata.var['mt'] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs.pct_counts_mt < int(args.mtpct), :]

#Log-normalize the data
sc.pp.normalize_total(adata, target_sum=int(args.reads))
sc.pp.log1p(adata)

adata.write(results_file)