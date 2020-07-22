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

def get_tissue_type(dataset_directory: Path) -> str:

    data_set_spreadsheet = "/home/sean/Documents/code/cross-dataset-diffexpr/bin/spreadsheet.csv"

    data_set_dir = fspath(dataset_directory)

    spreadsheet_df = pd.read_csv(data_set_spreadsheet)

    for i in range(len(spreadsheet_df.index)):
        if data_set_dir in str(spreadsheet_df['localPath'][i]):
            return spreadsheet_df['Organ/Tissue'][i]

def get_cell_by_gene_df(cell_by_gene_file: Path)->pd.DataFrame:
  f = h5py.File(cell_by_gene_file, 'r')
  cell_by_gene = np.array(f['cell_by_gene'])
  genes = np.array(f['col_names'])
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

def merge_dfs(cell_by_gene_file:Path, cell_motif_file:Path, cell_cluster_file:Path)->pd.DataFrame:

    dataset = get_dataset(cell_by_gene_file)
    tissue_type = get_tissue_type(cell_by_gene_file.parent)

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

def get_gene_rows(gene_groupings:Dict[str, List[str]], marker_gene_groupings:Dict[str, List[str]]):

    gene_rows = []

    for gene in gene_groupings.keys():
        gene_id = gene
        gene_list = gene_groupings[gene]
        if gene in marker_gene_groupings.keys():
            marker_gene_list = marker_gene_groupings[gene]
        else:
            marker_gene_list = []
        gene_rows.append({'gene_id':gene_id, 'groups':gene_list, 'marker_groups':marker_gene_list})

    return gene_rows

def get_rows(adata:anndata.AnnData, groupings:List[str])->List[Dict]:

    group_rows = []

    gene_groupings = {}
    marker_gene_groupings = {}

    cutoff = 0.9
    marker_cutoff = .001

    num_genes = len(adata.var_names)

    cell_df = adata.obs.copy()

    for group_by in groupings:
    #for each thing we want to group by

        sc.tl.rank_genes_groups(adata, group_by, method='t-test', rankby_abs=True, n_geness=num_genes)

        #get the group_ids and then the gene_names and scores for each
        for group_id in cell_df[group_by].unique():

            if type(group_id) == float and np.isnan(group_id):
                continue

            condition = group_by + "==" + str(group_id)

            gene_names = adata.uns['rank_genes_groups']['names'][group_id]
            pvals = adata.uns['rank_genes_groups']['pvals'][group_id]
            names_and_pvals = zip(gene_names, pvals)

            for n_p in names_and_pvals:

                if n_p[1] < cutoff:
                    if n_p[0] not in gene_groupings.keys():
                        gene_groupings[n_p[0]] = []
                    gene_groupings[n_p[0]].append(condition)

                if n_p[1] < marker_cutoff:
                    if n_p[0] not in marker_gene_groupings.keys():
                        marker_gene_groupings[n_p[0]] = []
                    marker_gene_groupings[n_p[0]].append(condition)

            genes = [n_p[0] for n_p in names_and_pvals if np[1] < cutoff]
            marker_genes = [n_p[0] for n_p in names_and_pvals if np[1] < marker_cutoff]

            group_rows.append({'condition':condition, 'genes':genes, 'marker_genes':marker_genes})

        gene_rows = get_gene_rows(gene_groupings, marker_gene_groupings)

        return group_rows, gene_rows

def main(output_directory: Path):

    dataset_dfs = [merge_dfs(cell_by_gene_file, cell_motif_file, cell_cluster_file) for cell_by_gene_file, cell_motif_file, cell_cluster_file in get_output_files(output_directory)]

    modality_df = pd.concat(dataset_dfs)

    adata = make_adata(modality_df)

    groupings = ['cluster', 'dataset']

    group_rows, gene_rows = get_rows(adata, groupings)

    cell_df = adata.obs.copy()

    gene_df = pd.DataFrame(gene_rows)
    group_df = pd.DataFrame(group_rows)

    cell_df.to_csv('atac.csv')
    group_df.to_csv('atac_group.csv')
    gene_df.to_csv('atac_gene.csv')

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('output_directory', type=Path)
    args = p.parse_args()

    main(args.output_directory)
