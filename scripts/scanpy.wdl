version 1.0

workflow scanpy{
	input{
		Array[File] channels
		Array[String] names
		Int mingenes
		Int maxgenes
		Int mtpct
		Int cellpct
		Int reads
		Int ntopgenes
		Int components
		Int neighbors
		Float resolution
	}

	scatter(pair in zip(channels, names)){
		call qc_filtering{
			input: h5 = pair.left,
			name = pair.right,
			mingenes = mingenes,
			maxgenes = maxgenes,
			mtpct = mtpct,
			cellpct = cellpct,
			reads = reads
		}

		call highly_variable{
			input: h5ad = qc_filtering.post_filtering,
			name = pair.right,
			ntopgenes = ntopgenes
		}

		call pca{
			input: h5ad = highly_variable.highly_variable,
			name = pair.right,
			components = components
		}

		call nearest_neighbors{
			input: h5ad = pca.pca,
			name = pair.right,
			components = components,
			neighbors = neighbors
		}

		call embedding{
			input: h5ad = nearest_neighbors.nearest_neighbors,
			name = pair.right,
			resolution = resolution
		}

		call rank_genes{
			input: h5ad = embedding.embedding,
			name = pair.right
		}
	}

	output{
		Array[File] count_matrix_h5ad = rank_genes.ranked_genes
		Array[File] umap_png = embedding.umap_png
		Array[File] gene_rank_png = rank_genes.ranked_genes_png
	}
}

task qc_filtering{
	input{
		String h5
		String name
		Int mingenes
		Int maxgenes
		Int mtpct
		Int cellpct
		Int reads
	}

	command<<<
		python3 qc_filtering.py -i ~{h5} -o "~{name}.h5ad" -n ~{mingenes} -x ~{maxgenes} -m ~{mtpct} -c ~{cellpct} -r ~{reads}
	>>>

	output{
		File post_filtering = "~{name}.h5ad"
	}

	runtime{
		cpu: "1"
		memory: "4 G"
		docker: "atex91/scanpy-project"
	}
}

task highly_variable{
	input{
		String h5ad
		String name
		Int ntopgenes
	}

	command<<<
		python3 highly_variable.py -i ~{h5ad} -o "~{name}.h5ad" -n ~{ntopgenes}
	>>>

	output{
		File highly_variable = "~{name}.h5ad"
	}

	runtime{
		cpu: "1"
		memory: "8 G"
		docker: "atex91/scanpy-project"
	}
}

task pca{
	input{
		String h5ad
		String name
		Int components
	}

	command<<<
		python3 pca.py -i ~{h5ad} -o "~{name}.h5ad" -c ~{components}
	>>>

	output{
		File pca = "~{name}.h5ad"
	}

	runtime{
		cpu: "1"
		memory: "8 G"
		docker: "atex91/scanpy-project"
	}
}

task nearest_neighbors{
	input{
		String h5ad
		String name
		Int components
		Int neighbors
	}

	command<<<
		python3 nearest_neighbors.py -i ~{h5ad} -o "~{name}.h5ad" -c ~{components} -n ~{neighbors}
	>>>

	output{
		File nearest_neighbors = "~{name}.h5ad"
	}

	runtime{
		cpu: "1"
		memory: "8"
		docker: "atex91/scanpy-project"
	}
}

task embedding{
	input{
		String h5ad
		String name
		Float resolution
	}

	command<<<
		python3 embedding.py -i ~{h5ad} -o "~{name}.h5ad" -r ~{resolution}
	>>>

	output{
		File embedding = "~{name}.h5ad"
		File umap_png = "umap_~{name}.png"
	}

	runtime{
		cpu: "1"
		memory: "8"
		docker: "atex91/scanpy-project"
	}
}

task rank_genes{
	input{
		String h5ad
		String name
	}

	command<<<
		python3 rank_genes -i ~{h5ad} -o "~{name}.h5ad"
	>>>

	output{
		File ranked_genes = "~{name}.h5ad"
		File ranked_genes_png = "rank_genes_groups_~{name}.png"
	}

	runtime{
		cpu: "1"
		memory: "8"
		docker: "atex91/scanpy-project"
	}
}