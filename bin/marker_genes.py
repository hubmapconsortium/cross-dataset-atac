#!/usr/bin/env python3
import anndata
import pandas as pd
from argparse import ArgumentParser
from pathlib import Path
from cross_dataset_common import get_pval_dfs, flatten_quant_df

def get_quant_df(adata: anndata.AnnData)->pd.DataFrame:
    return pd.DataFrame(adata.X, columns=adata.var.index, index=adata.obs.index)


def main(concatenated_annotated_file: Path):
    adata = anndata.read_h5ad(concatenated_annotated_file)

    adata.obs['cell_id'] = adata.obs.index

    cell_df = adata.obs.copy()

    quant_df = get_quant_df(adata)

    long_df = flatten_quant_df(quant_df)
    long_df.to_csv('atac.csv')

    pval_df = get_pval_dfs(adata)

    with pd.HDFStore('atac.hdf5') as store:
        store.put('cell', cell_df, format='t')
        store.put('p_values', pval_df)


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('concatenated_annotated_file', type=Path)
    args = p.parse_args()

    main(args.concatenated_annotated_file)
