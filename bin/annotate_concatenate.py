#!/usr/bin/env python3
import json
from argparse import ArgumentParser
from pathlib import Path
from typing import Dict, List, Tuple

import anndata
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
from cross_dataset_common import get_tissue_type, hash_cell_id, get_cluster_df, precompute_dataset_percentages, precomptue_percentages

from hubmap_cell_id_gen_py import get_sequencing_cell_id
from concurrent.futures import ThreadPoolExecutor

CELL_GY_GENE_FILENAME = 'cell_by_gene.hdf5'
CELL_BY_GENE_H5AD_FILENAME = 'cell_by_gene.h5ad'
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

def read_cell_by_gene_hdf5(directory: Path, nexus_token: str) -> anndata.AnnData:
    dataset = directory.stem
    tissue_type = get_tissue_type(dataset, nexus_token)

    cell_by_gene, cells, genes = get_cell_by_gene_data(directory / CELL_GY_GENE_FILENAME)

    cluster_df = pd.read_csv(directory / CELL_CLUSTER_FILENAME, index_col=0)
    cluster_list = [f'leiden-UMAP-{dataset}-{cluster}' for cluster in cluster_df['cluster']]
    cluster_series = pd.Series(cluster_list, index=cluster_df.index)

    barcodes = [cell_id for cell_id in cells]
    semantic_cell_ids = [get_sequencing_cell_id(dataset, barcode) for barcode in barcodes]

    data_for_obs_df = {
        'cell_id': semantic_cell_ids,
        'barcode': barcodes,
        'leiden': cluster_series.loc[cells],
        'dataset': dataset,
        'organ': tissue_type,
        'modality':'atac',
    }
    obs_df = pd.DataFrame(data_for_obs_df, index=cells)

    adata = anndata.AnnData(
        X=cell_by_gene,
        obs=obs_df,
        var=pd.DataFrame(index=genes),
    )
    return adata

def read_cell_by_gene_h5ad(directory, nexus_token):
    dataset = directory.stem
    tissue_type = get_tissue_type(dataset, nexus_token)

    adata = anndata.read(directory / CELL_BY_GENE_H5AD_FILENAME)

    cluster_list = [f'leiden-UMAP-{dataset}-{cluster}' for cluster in adata.obs['cluster']]
    cluster_series = pd.Series(cluster_list, index=adata.obs.index)
    adata.obs['leiden'] = cluster_series

    barcodes = [cell_id for cell_id in adata.obs.index]
    semantic_cell_ids = [get_sequencing_cell_id(dataset, barcode) for barcode in barcodes]

    adata.obs['cell_id'] = pd.Series(semantic_cell_ids, index=adata.obs.index)
    adata.obs['barcode'] = pd.Series(barcodes, index=adata.obs.index)
    adata.obs['dataset'] = dataset
    adata.obs['organ'] = tissue_type
    adata.obs['modality'] = 'atac'

    return adata

def read_cell_by_gene(directory, nexus_token):
    try:
        return read_cell_by_gene_h5ad(directory, nexus_token)
    except:
        print(directory)
        return read_cell_by_gene_hdf5(directory, nexus_token)

def map_gene_ids(adata):
    gene_mapping = read_gene_mapping()
    keep_vars = [gene in gene_mapping for gene in adata.var.index]
    adata = adata[:, keep_vars]
    if isinstance(adata.X, np.ndarray):
        temp_df = pd.DataFrame(adata.X, index=adata.obs.index, columns=adata.var.index)
    else:
        temp_df = pd.DataFrame(adata.X.todense(), index=adata.obs.index, columns=adata.var.index)
    aggregated = temp_df.groupby(level=0, axis=1).sum()
    adata = anndata.AnnData(aggregated, obs=adata.obs)
    adata.var.index = [gene_mapping[var] for var in adata.var.index]
    # This introduces duplicate gene names, use Pandas for aggregation
    # since anndata doesn't have that functionality
    adata.var_names_make_unique()
    return adata

def main(nexus_token:str, output_directories: List[Path]=[]):

    if nexus_token == "None":
        nexus_token = None

    adatas = [read_cell_by_gene(directory, nexus_token) for directory in output_directories]
    for adata in adatas:
        print(adata.obs["dataset"].unique())
        print(type(adata.X))
        if isinstance(adata.X, np.ndarray):
            print(adata.X.shape)
        else:
            print(adata.X.todense().shape)
    gene_mapped_adatas = [map_gene_ids(adata) for adata in adatas]
    for adata in gene_mapped_adatas:
        sc.tl.rank_genes_groups(adata, 'leiden', method='t-test', rankby_abs=True, n_genes=len(adata.var.index))

    with ThreadPoolExecutor(max_workers=len(output_directories)) as e:
        percentage_dfs = e.map(precompute_dataset_percentages, gene_mapped_adatas)

    percentage_df = pd.concat(percentage_dfs)

    cluster_dfs = [get_cluster_df(adata) for adata in gene_mapped_adatas]
    cluster_df = pd.concat(cluster_dfs)
    with pd.HDFStore('cluster.hdf5') as store:
        store.put('cluster', cluster_df)

    first, *rest = adatas
    concatenated = first.concatenate(rest)
    concatenated = map_gene_ids(concatenated)

    concatenated.uns['percentages'] = percentage_df

    concatenated.write('concatenated_annotated.h5ad')

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('nexus_token', type=str)
    p.add_argument('data_directories', type=Path, nargs='+')
    p.add_argument("--enable-manhole", action="store_true")
    args = p.parse_args()

    if args.enable_manhole:
        import manhole

        manhole.install(activate_on="USR1")

    main(args.nexus_token, args.data_directories)
