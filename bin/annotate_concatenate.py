#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from typing import Iterable, Tuple
import h5py
import pandas as pd
import numpy as np
import anndata
import scanpy as sc
import json
from typing import List

from cross_dataset_common import get_tissue_type, get_gene_dicts

def get_dataset(cell_by_gene_file: Path) -> str:
    return cell_by_gene_file.parent.stem

def get_cell_by_gene_df(cell_by_gene_file: Path) -> pd.DataFrame:
    f = h5py.File(cell_by_gene_file, 'r')
    cell_by_gene = np.array(f['cell_by_gene'])
    genes = np.array(f['col_names'])
    cells = np.array(f['row_names'])
    cell_by_gene_df = pd.DataFrame(cell_by_gene, columns=genes, index=cells)
    return cell_by_gene_df


def get_output_files(directories: List[Path]) -> Iterable[Tuple[Path, Path, Path]]:
    relative_paths = (Path('cell_by_gene.hdf5'), Path('cellMotif.csv'), Path('cellClusterAssignment.csv'))
    output_files = [(directory / relative_paths[0], directory / relative_paths[1], directory / relative_paths[2]) for
                    directory in directories]
    return output_files


def merge_dfs(cell_by_gene_file: Path, cell_motif_file: Path, cell_cluster_file: Path,
              nexus_token: str) -> pd.DataFrame:
    dataset = get_dataset(cell_by_gene_file)
    tissue_type = get_tissue_type(dataset, nexus_token)

    cell_by_gene_df = get_cell_by_gene_df(cell_by_gene_file)
    cell_by_gene_df['cell_id'] = cell_by_gene_df.index

    cluster_df = pd.read_csv(cell_cluster_file)
    cluster_list = [dataset + '-' + str(cluster) for cluster in cluster_df['Cluster']]
    cluster_df['cluster'] = pd.Series(cluster_list, dtype=str)
    cluster_df['cluster'] = cluster_df['cluster'].astype(str)
    print(cluster_df['cluster'].dtype)
    cluster_df['cell_id'] = cluster_df['BarcodeID']
    cluster_df = cluster_df[['cell_id', 'cluster']].copy()
    print(cluster_df['cluster'].dtype)

    merge_df = cell_by_gene_df.merge(cluster_df, how='outer')

    merge_df['dataset'] = dataset
    merge_df['tissue_type'] = tissue_type
    merge_df['modality'] = 'atac'

    print(merge_df['cluster'].dtype)

    return merge_df


def make_adata(modality_df: pd.DataFrame, ensembl_to_symbol_path: Path, symbol_to_ensembl_path: Path):
    obs_columns = ['cell_id', 'cluster', 'dataset', 'tissue_type', 'modality']
    obs = modality_df[obs_columns].copy()
    var_columns = modality_df.drop(obs_columns, axis=1).columns
    var_columns = [column.decode('UTF-8') for column in var_columns]

    symbol_to_ensembl_dict = {}
    ensembl_to_symbol_dict = {}

    if ensembl_to_symbol_path.exists():
        with open(ensembl_to_symbol_path, 'r') as json_file:
            ensembl_to_symbol_dict = json.load(json_file)
        with open(symbol_to_ensembl_path, 'r') as json_file:
            symbol_to_ensembl_dict = json.load(json_file)

    else:
        ensembl_to_symbol_dict, symbol_to_ensembl_dict = get_gene_dicts(var_columns)
        with open('ensembl_to_symbol.json', 'w') as json_file:
            json.dump(ensembl_to_symbol_dict, json_file)
        with open('symbol_to_ensembl.json', 'w') as json_file:
            json.dump(symbol_to_ensembl_dict, json_file)

    keep_vars = [key for key in ensembl_to_symbol_dict.keys()]

    var = pd.DataFrame(index=var_columns)
    x = modality_df.drop(obs_columns, axis=1).to_numpy()

    uns = {'omic': 'ATAC'}

    adata = anndata.AnnData(X=x, var=var, obs=obs, uns=uns)

    adata = adata[:,keep_vars]
    adata.var.index = [ensembl_to_symbol_dict[ensembl_id] for ensembl_id in adata.var.index]

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    return adata


def main(nexus_token: str, output_directories: List[Path],
         ensembl_to_symbol_path: Path = Path('/opt/ensembl_to_symbol.json'),
         symbol_to_ensembl_path: Path = Path('/opt/symbol_to_ensembl.json')):
    dataset_dfs = [merge_dfs(cell_by_gene_file, cell_motif_file, cell_cluster_file, nexus_token) for
                   cell_by_gene_file, cell_motif_file, cell_cluster_file in get_output_files(output_directories)]

    modality_df = pd.concat(dataset_dfs)

    adata = make_adata(modality_df, ensembl_to_symbol_path, symbol_to_ensembl_path)

    adata.write('concatenated_annotated.h5ad')


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('nexus_token', type=str)
    p.add_argument('data_directories', type=Path, nargs='+')
    args = p.parse_args()

    main(args.nexus_token, args.data_directories)
