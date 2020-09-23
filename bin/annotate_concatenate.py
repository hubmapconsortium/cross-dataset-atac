#!/usr/bin/env python3
import json
from argparse import ArgumentParser
from pathlib import Path
from typing import Dict, List, Tuple

import anndata
import h5py
import numpy as np
import pandas as pd
from cross_dataset_common import get_tissue_type

CELL_GY_GENE_FILENAME = 'cell_by_gene.hdf5'
CELL_CLUSTER_FILENAME = 'umap_coords_clusters.csv'

GENE_MAPPING_DIRECTORIES = [
    Path(__file__).parent.parent / 'data',
    Path('/opt/data'),
]

def read_gene_mapping() -> Dict[str, str]:
    """
    Try to find the Ensembl to HUGO symbol mapping, with paths suitable
    for running this script inside and outside a Docker container.
    :return:
    """
    for directory in GENE_MAPPING_DIRECTORIES:
        mapping_file = directory / 'ensembl_to_symbol.json'
        if mapping_file.is_file():
            with open(mapping_file) as f:
                return json.load(f)
    message_pieces = ["Couldn't find Ensembl → HUGO mapping file. Tried:"]
    message_pieces.extend(f'\t{path}' for path in GENE_MAPPING_DIRECTORIES)
    raise ValueError('\n'.join(message_pieces))

def get_cell_by_gene_data(cell_by_gene_file: Path) -> Tuple[np.ndarray, List[str], List[str]]:
    with h5py.File(cell_by_gene_file, 'r') as f:
        cell_by_gene = np.array(f['cell_by_gene'])
        cells = [row.decode('utf-8') for row in np.array(f['row_names'])]
        genes = [col.decode('utf-8') for col in np.array(f['col_names'])]

    return cell_by_gene, cells, genes

def read_cell_by_gene(directory: Path, nexus_token: str) -> anndata.AnnData:
    dataset = directory.stem
    tissue_type = get_tissue_type(dataset, nexus_token)

    cell_by_gene, cells, genes = get_cell_by_gene_data(directory / CELL_GY_GENE_FILENAME)

    cluster_df = pd.read_csv(directory / CELL_CLUSTER_FILENAME, index_col=0)
    cluster_list = [f'{dataset}-{cluster}' for cluster in cluster_df['cluster']]
    cluster_series = pd.Series(cluster_list, index=cluster_df.index)

    data_for_obs_df = {
        'cluster': cluster_series.loc[cells],
        'dataset': dataset,
        'tissue_type': tissue_type,
        'modality': 'atac',
    }
    obs_df = pd.DataFrame(data_for_obs_df, index=cells)

    adata = anndata.AnnData(
        X=cell_by_gene,
        obs=obs_df,
        var=pd.DataFrame(index=genes),
    )
    return adata

def main(nexus_token: str, output_directories: List[Path]):
    adatas = [read_cell_by_gene(directory, nexus_token) for directory in output_directories]
    first, *rest = adatas
    concatenated = first.concatenate(rest)
    concatenated.uns['omic'] = 'ATAC'

    gene_mapping = read_gene_mapping()
    keep_vars = [gene in gene_mapping for gene in concatenated.var.index]
    concatenated = concatenated[:, keep_vars]
    concatenated.var.index = [gene_mapping[var] for var in concatenated.var.index]
    # This introduces duplicate gene names, use Pandas for aggregation
    # since anndata doesn't have that functionality
    temp_df = pd.DataFrame(concatenated.X, index=concatenated.obs.index, columns=concatenated.var.index)
    aggregated = temp_df.groupby(level=0, axis=1).sum()

    adata = anndata.AnnData(aggregated, obs=concatenated.obs)
    adata.uns['omic'] = 'ATAC'

    adata.write('concatenated_annotated.h5ad')

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('nexus_token')
    p.add_argument('data_directory', type=Path, nargs='+')
    args = p.parse_args()

    main(args.nexus_token, args.data_directory)
