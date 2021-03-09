#!/usr/bin/env python3
import anndata
import pandas as pd
from argparse import ArgumentParser
from pathlib import Path
from cross_dataset_common import get_pval_dfs, make_quant_df, create_minimal_dataset, tar_zip_scp


def main(concatenated_annotated_file: Path, old_cluster_file:Path, known_hosts_file: Path):
    adata = anndata.read_h5ad(concatenated_annotated_file)

    cell_df = adata.obs.copy()
#    cell_df["clusters"] = pd.Series([[cluster] for cluster in cell_df["cluster"]])
    cell_df = cell_df[['cell_id', 'barcode', 'dataset', 'organ', 'modality', 'clusters']]

    quant_df = make_quant_df(adata)
    quant_df.to_csv('atac.csv')

    organ_df, cluster_df = get_pval_dfs(adata, "atac")

    with pd.HDFStore(old_cluster_file) as store:
        cluster_df = store.get('cluster')

    with pd.HDFStore('atac.hdf5') as store:
        store.put('cell', cell_df, format='t')
        store.put('organ', organ_df)
        store.put('cluster', cluster_df)

    create_minimal_dataset(cell_df, quant_df, organ_df, cluster_df, 'atac')

    tar_zip_scp("atac", known_hosts_file)


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('concatenated_annotated_file', type=Path)
    p.add_argument('old_cluster_file', type=Path)
    p.add_argument('known_hosts_file', type=Path)
    args = p.parse_args()

    main(args.concatenated_annotated_file, args.old_cluster_file)
