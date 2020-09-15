#!/usr/bin/env python3
import anndata
import pandas as pd
from argparse import ArgumentParser
from pathlib import Path
from cross_dataset_common import get_rows


def main(concatenated_annotated_file: Path):
    adata = anndata.read_h5ad(concatenated_annotated_file)

    groupings = ['cluster', 'dataset', 'tissue_type']

    group_rows = get_rows(adata, groupings)

    cell_df = adata.obs.copy()

    group_df = pd.DataFrame(group_rows)

    cell_df.to_csv('atac.csv')
    group_df.to_csv('atac_group.csv')


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('concatenated_annotated_file', type=Path)
    args = p.parse_args()

    main(args.concatenated_annotated_file)
