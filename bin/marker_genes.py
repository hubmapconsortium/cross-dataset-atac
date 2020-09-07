#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from typing import List, Dict

import anndata
import numpy as np
import pandas as pd
import scanpy as sc


def get_rows(adata: anndata.AnnData, groupings: List[str]) -> List[Dict]:
    group_rows = []

    cutoff = 0.9
    marker_cutoff = .001

    num_genes = len(adata.var_names)

    cell_df = adata.obs.copy()

    for group_by in groupings:
        # for each thing we want to group by

        sc.tl.rank_genes_groups(adata, group_by, method='t-test', rankby_abs=True, n_geness=num_genes)

        # get the group_ids and then the gene_names and scores for each
        for group_id in cell_df[group_by].unique():

            if type(group_id) == float and np.isnan(group_id):
                continue

            gene_names = adata.uns['rank_genes_groups']['names'][group_id]
            pvals = adata.uns['rank_genes_groups']['pvals'][group_id]
            names_and_pvals = zip(gene_names, pvals)

            genes = [n_p[0] for n_p in names_and_pvals if n_p[1] < cutoff]
            marker_genes = [n_p[0] for n_p in names_and_pvals if n_p[1] < marker_cutoff]

            group_rows.append(
                {'group_type': group_by, 'group_id': group_id, 'genes': genes, 'marker_genes': marker_genes})

        return group_rows


def main(concatenated_annotated_file: Path):
    adata = anndata.read_h5ad(concatenated_annotated_file)

    groupings = ['cluster', 'dataset']

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
