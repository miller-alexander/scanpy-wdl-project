FROM  python:3.10
MAINTAINER Alex Miller <alexander.miller326@gmail.com>

LABEL Image for simple scanpy pre-analysis of single cell RNA-sequencing data

ADD qc_filtering.py highly_variable.py pca.py nearest_neighbors.py embedding.py rank_genes.py .

RUN pip install --no-cache-dir numpy pandas scanpy leidenalg
