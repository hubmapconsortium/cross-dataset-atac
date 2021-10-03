#!/usr/bin/env python3
import anndata
import pandas as pd
from argparse import ArgumentParser
from pathlib import Path
from cross_dataset_common import get_pval_dfs, make_quant_df, create_minimal_dataset, precompute_dataset_percentages
from concurrent.futures import ThreadPoolExecutor

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

def main(concatenated_annotated_file: Path):
    adata = anndata.read_h5ad(concatenated_annotated_file)
    adata.var_names_make_unique()

    gene_df = annotate_genes(adata)

    dataset_adatas = [adata[adata.obs['dataset'] == dataset] for dataset in adata.obs['dataset'].unique()]
    with ThreadPoolExecutor(max_workers=len(dataset_adatas)) as e:
         percentage_dfs = e.map(precompute_dataset_percentages, dataset_adatas)

    percentage_df = pd.concat(percentage_dfs)

    with pd.HDFStore('atac_precompute.hdf5') as store:
        store.put('percentages', percentage_df)
        store.put('gene', gene_df)

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('concatenated_annotated_file', type=Path)
    args = p.parse_args()

    main(args.concatenated_annotated_file)
