#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from typing import Iterable, Tuple, List, Dict
import h5py
import pandas as pd
import numpy as np
import anndata
from os import fspath
import scanpy as sc

def get_dataset(cell_by_gene_file: Path)->str:
    return cell_by_gene_file.parent.stem

def ensemble_to_symbol(ensemble_id:str)->str:
    request_url = 'https://mygene.info/v3/gene/' + ensemble_id + '?fields=symbol&dotfield=True'
    r = requests.get(request_url)
    return r.json()['symbol']

def get_tissue_type(dataset:str, token:str)->str:

    organ_dict = yaml.load(open('organ_types.yaml'), Loader=yaml.BaseLoader)

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
    'https://search-api.dev.hubmapconsortium.org/search',
    json = dataset_query_dict,
    headers = {'Authorization': 'Bearer ' + token})
    hits = dataset_response.json()['hits']['hits']

    for hit in hits:
        for ancestor in hit['_source']['ancestors']:
            if 'organ' in ancestor.keys():
                return organ_dict[ancestor['organ']]['description']

def get_cell_by_gene_df(cell_by_gene_file: Path)->pd.DataFrame:
  f = h5py.File(cell_by_gene_file, 'r')
  cell_by_gene = np.array(f['cell_by_gene'])
  genes = np.array(f['col_names'])
  genes = [ensemble_to_symbol(gene) for gene in genes]
  cells = np.array(f['row_names'])
  cell_by_gene_df = pd.DataFrame(cell_by_gene, columns=genes, index=cells)
  return cell_by_gene_df

def get_cell_by_gene_files(directory: Path) -> Iterable[Path]:
    pattern = '**/cell_by_gene.hdf5'
    yield from directory.glob(pattern)

def get_output_files(directory:Path) -> Iterable[Tuple[Path, Path, Path]]:
    cell_by_gene_files = get_cell_by_gene_files(directory)
    output_files = [(cell_by_gene_file, cell_by_gene_file.parent / Path('cellMotif.csv'), cell_by_gene_file.parent / Path('cellClusterAssignment.csv')) for cell_by_gene_file in cell_by_gene_files]
    return output_files

def merge_dfs(cell_by_gene_file:Path, cell_motif_file:Path, cell_cluster_file:Path, nexus_token:str)->pd.DataFrame:

    dataset = get_dataset(cell_by_gene_file)
    tissue_type = get_tissue_type(cell_by_gene_file.parent, nexus_token)

    cell_by_gene_df = get_cell_by_gene_df(cell_by_gene_file)
    cell_by_gene_df['cell_id'] = cell_by_gene_df.index

    cluster_df = pd.read_csv(cell_cluster_file)
    cluster_df['cluster'] = dataset + "-" + cluster_df['Cluster'].astype(str)
    cluster_df['cell_id'] = cluster_df['BarcodeID']
    cluster_df = cluster_df[['cell_id', 'cluster']].copy()

    cell_by_gene_df = cell_by_gene_df.head(1000).copy()
    cluster_df = cluster_df.head(1000).copy()

    merge_df = cell_by_gene_df.merge(cluster_df, how='outer')

    merge_df['dataset'] = dataset
    merge_df['tissue_type'] = tissue_type
    merge_df['modality'] = 'atac'

    return merge_df

def make_adata(modality_df: pd.DataFrame):

    obs_columns = ['cell_id', 'cluster', 'dataset', 'tissue_type', 'modality']
    obs = modality_df[obs_columns].copy()
    var_columns = modality_df.drop(obs_columns, axis=1).columns
    var = pd.DataFrame(index=var_columns)
    x = modality_df.drop(obs_columns, axis=1).to_numpy()

    uns = {'omic':'ATAC'}

    adata = anndata.AnnData(X=x, var=var, obs=obs, uns=uns)

    return adata


def main(output_directory: Path):

    dataset_dfs = [merge_dfs(cell_by_gene_file, cell_motif_file, cell_cluster_file) for cell_by_gene_file, cell_motif_file, cell_cluster_file in get_output_files(output_directory)]

    modality_df = pd.concat(dataset_dfs)

    adata = make_adata(modality_df)

    adata.write('concatenated_annotated.h5ad')


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('output_directory', type=Path)
    args = p.parse_args()

    main(args.output_directory)
