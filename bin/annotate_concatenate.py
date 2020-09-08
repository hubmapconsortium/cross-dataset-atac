#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from typing import Iterable, Tuple, List
import h5py
import pandas as pd
import numpy as np
import anndata
import scanpy as sc
import yaml
import requests
import json


def get_dataset(cell_by_gene_file: Path) -> str:
    return cell_by_gene_file.parent.stem


def ensembl_to_symbol(ensembl_id: str) -> str:
    ensembl_id = ensembl_id.split('.')[0]
    request_url = 'https://mygene.info/v3/gene/' + ensembl_id.split('.')[0] + '?fields=symbol&dotfield=True'
    r = requests.get(request_url)
    return r.json()['symbol']


def get_tissue_type(dataset: str, token: str) -> str:
    organ_dict = yaml.load(open('/opt/organ_types.yaml'), Loader=yaml.BaseLoader)

    dataset_query_dict = {
        "query": {
            "bool": {
                "must": [],
                "filter": [
                    {
                        "match_all": {}
                    },
                    {
                        "exists": {
                            "field": "files.rel_path"
                        }
                    },
                    {
                        "match_phrase": {
                            "uuid": {
                                "query": dataset
                            },
                        }

                    }
                ],
                "should": [],
                "must_not": [
                    {
                        "match_phrase": {
                            "status": {
                                "query": "Error"
                            }
                        }
                    }
                ]
            }
        }
    }

    dataset_response = requests.post(
        'https://search-api.hubmapconsortium.org/search',
        json=dataset_query_dict,
        headers={'Authorization': 'Bearer ' + token})
    hits = dataset_response.json()['hits']['hits']

    for hit in hits:
        for ancestor in hit['_source']['ancestors']:
            if 'organ' in ancestor.keys():
                return organ_dict[ancestor['organ']]['description']


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


def make_adata(modality_df: pd.DataFrame):
    obs_columns = ['cell_id', 'cluster', 'dataset', 'tissue_type', 'modality']
    obs = modality_df[obs_columns].copy()
    var_columns = modality_df.drop(obs_columns, axis=1).columns
    symbol_ensembl_dict = {ensembl_to_symbol(ensembl_id): ensembl_id for ensembl_id in var_columns}
    with open("symbol_to_ensemble.json", "w") as dictionary_file:
        json.dump(symbol_ensembl_dict, dictionary_file)
    var_columns = [ensembl_to_symbol(ensembl_id) for ensembl_id in var_columns]
    var = pd.DataFrame(index=var_columns)
    x = modality_df.drop(obs_columns, axis=1).to_numpy()

    uns = {'omic': 'ATAC'}

    adata = anndata.AnnData(X=x, var=var, obs=obs, uns=uns)

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    return adata


def main(nexus_token: str, output_directories: List[Path]):
    dataset_dfs = [merge_dfs(cell_by_gene_file, cell_motif_file, cell_cluster_file, nexus_token) for
                   cell_by_gene_file, cell_motif_file, cell_cluster_file in get_output_files(output_directories)]

    modality_df = pd.concat(dataset_dfs)

    adata = make_adata(modality_df)

    adata.write('concatenated_annotated.h5ad')


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('nexus_token', type=str)
    p.add_argument('data_directories', type=Path, nargs='+')
    args = p.parse_args()

    main(args.nexus_token, args.data_directories)
