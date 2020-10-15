#!/usr/bin/env python3
import anndata
import pandas as pd
from argparse import ArgumentParser
from pathlib import Path
from cross_dataset_common import get_rows, add_quant_columns

def get_quant_df(adata: anndata.AnnData)->pd.DataFrame:
    return pd.DataFrame(adata.X, columns=adata.var.index, index=adata.obs.index)

def make_long_df(quant_df: pd.DataFrame):
    dict_list = [{'cell_id': i, 'gene_id': column, 'value':quant_df.at[i, column], 'modality':'atac'} for i in quant_df.index for column in quant_df.columns]
    return pd.DataFrame(dict_list)

def main(concatenated_annotated_file: Path):
    adata = anndata.read_h5ad(concatenated_annotated_file)

    groupings = ['cluster', 'dataset', 'tissue_type']

    adata.obs['cell_id'] = adata.obs.index

    group_rows = get_rows(adata, groupings)

    cell_df = adata.obs.copy()

    quant_df = get_quant_df(adata)

    long_df = make_long_df(quant_df)
    long_df.to_csv('long_atac_quant.csv')

    group_df = pd.DataFrame(group_rows, dtype=object)

    with pd.HDFStore('atac.hdf5') as store:
        store.put('cell', cell_df, format='t')
        store.put('group', group_df)
        store.put('quant', quant_df)


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('concatenated_annotated_file', type=Path)
    args = p.parse_args()

    main(args.concatenated_annotated_file)
