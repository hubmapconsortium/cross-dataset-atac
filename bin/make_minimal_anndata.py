#!/usr/bin/env python3
import anndata
import pandas as pd
from argparse import ArgumentParser
from pathlib import Path
from cross_dataset_common import get_pval_dfs, make_quant_df, create_minimal_dataset, precompute_dataset_percentages
from hubmap_cell_id_gen_py import get_sequencing_cell_id
from concurrent.futures import ThreadPoolExecutor
import json
from typing import Dict
import numpy as np
import scipy
from scipy.sparse import csr_matrix

def annotate_single_gene(param_tuple):
    df = param_tuple[0]
    var_id = param_tuple[1]
    return_dict = {'gene_symbol':var_id}
    gene_df = df[df[var_id] > 0]
    return_dict['num_cells'] = len(gene_df.index)
    return_dict['num_datasets'] = len(gene_df['dataset'].unique())
    return return_dict

def annotate_genes(adata):
    df = adata.to_df()
    df['dataset'] = adata.obs['dataset']
    params_tuples = [(df[['dataset', var_id]], var_id) for var_id in df.columns if var_id != 'dataset']

    with ThreadPoolExecutor(max_workers=20) as e:
        dict_list = e.map(annotate_single_gene, params_tuples)

    return pd.DataFrame(dict_list)

def make_minimal_adata(csv_file: Path):
    quant_df = pd.read_csv(csv_file)
    var_ids = list(quant_df['q_var_id'].unique())
    var = pd.DataFrame(index=var_ids)
    obs_ids = list(quant_df['q_cell_id'].unique())
    obs = pd.DataFrame(index=obs_ids)
    X = np.zeros((len(obs_ids, len(var_ids))))
    adata = anndata.AnnData(X=X, var=var, obs=obs)
    for i in quant_df.index:
        obs_id = quant_df["q_cell_id"][i]
        var_id = quant_df["q_gene_id"][i]
        value = quant_df["value"][i]
        adata[obs_id, var_id] = value

    return adata

def make_pval_adata(hdf_file: Path):
    organ_df = pd.read_hdf(hdf_file, 'organ')
    cluster_df = pd.read_hdf(hdf_file, 'cluster')
    var_ids = list(organ_df['gene_id'].unique())
    var = pd.DataFrame(index=var_ids)
    grouping_ids = list(cluster_df['grouping_name'].unique()) + list(organ_df['grouping_name'].unique())
    obs = pd.DataFrame(index=grouping_ids)
    X = np.zeros((len(grouping_ids, len(var_ids))))
    adata = anndata.AnnData(X=X, var=var, obs=obs)
    for i in organ_df.index:
        grouping = organ_df["grouping_name"][i]
        var_id = organ_df["gene_id"][i]
        value = organ_df["value"][i]
        adata[grouping, var_id] = value
    for i in cluster_df.index:
        grouping = cluster_df["grouping_name"][i]
        var_id = cluster_df["gene_id"][i]
        value = cluster_df["value"][i]
        adata[grouping, var_id] = value
    return adata

def main(csv_file: Path, hdf_file: Path):
    min_adata = make_minimal_adata(csv_file)
    min_adata.write_h5ad("atac.h5ad")

    pval_adata = make_pval_adata(hdf_file)
    pval_adata.write("atac_pvals.h5ad")

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('h5ad_file', type=Path)
    p.add_argument('hdf_file', type=Path)
    args = p.parse_args()

    main(args.h5ad_file, args.hdf_file)
