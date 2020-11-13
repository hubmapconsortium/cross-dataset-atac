#!/usr/bin/env python3
import anndata
import pandas as pd
from argparse import ArgumentParser
from pathlib import Path
from cross_dataset_common import get_pval_dfs, make_quant_csv


def main(concatenated_annotated_file: Path):
    adata = anndata.read_h5ad(concatenated_annotated_file)

    adata.obs['cell_id'] = adata.obs.index

    cell_df = adata.obs.copy()

    make_quant_csv(adata, 'atac')

    organ_df, cluster_df = get_pval_dfs(adata)

    with pd.HDFStore('atac.hdf5') as store:
        store.put('cell', cell_df, format='t')
        store.put('organ', organ_df)
        store.put('cluster', cluster_df)


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('concatenated_annotated_file', type=Path)
    args = p.parse_args()

    main(args.concatenated_annotated_file)
