#!/usr/bin/env python3
import anndata
import pandas as pd
from argparse import ArgumentParser
from pathlib import Path
from cross_dataset_common import get_pval_dfs, make_quant_df, create_minimal_dataset, precompute_dataset_percentages, make_minimal_adata
from concurrent.futures import ThreadPoolExecutor

def annotate_genes(adata):
    df = adata.to_df()
    df['dataset'] = adata.obs['dataset']
    for gene in df.columns:
        if gene != 'dataset':
            gene_df = df[df[gene] > 0]
            adata.var.at[gene, 'num_cells'] = len(gene_df.index)
            adata.var.at[gene, 'num_datasets'] = len(gene_df['dataset'].unique())

    return adata.var.copy()

def main(concatenated_annotated_file: Path, old_cluster_file:Path):
    adata = anndata.read_h5ad(concatenated_annotated_file)
    adata.var_names_make_unique()

    cell_df = adata.obs.copy()
    cell_df["clusters"] = cell_df['leiden']
    cell_df = cell_df[['cell_id', 'barcode', 'dataset', 'organ', 'modality', 'clusters']]

    gene_df = annotate_genes(adata)

    organ_df, cluster_df = get_pval_dfs(adata)

    dataset_adatas = [adata[adata.obs['dataset'] == dataset] for dataset in adata.obs['dataset'].unique()]
    with ThreadPoolExecutor(max_workers=len(dataset_adatas)) as e:
         percentage_dfs = e.map(precompute_dataset_percentages, dataset_adatas)

    percentage_df = pd.concat(percentage_dfs)

    with pd.HDFStore(old_cluster_file) as store:
        cluster_df = store.get('cluster')

    minimal_adata = make_minimal_adata(adata, "atac")
    minimal_adata.write('atac.h5ad')

    with pd.HDFStore('atac.hdf5') as store:
        store.put('cell', cell_df, format='t')
        store.put('organ', organ_df)
        store.put('cluster', cluster_df)
        store.put('gene', gene_df)

    with pd.HDFStore('atac_precompute.hdf5') as store:
        store.put('percentages', percentage_df)

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('concatenated_annotated_file', type=Path)
    p.add_argument('old_cluster_file', type=Path)
    args = p.parse_args()

    main(args.concatenated_annotated_file, args.old_cluster_file)
