#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from typing import Iterable, Tuple
import h5py
import pandas as pd
import numpy as np

def get_datasets():
    return

def get_metadata():
    return

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

def get_output_files(directory:Path) -> Iterable[Tuple[Path, Path]]:
    cell_by_gene_files = get_cell_by_gene_files(output_directory)
    output_files = [(cell_by_gene_file, Path(cell_by_gene_file.name.replace('cell_by_gene.hdf5', 'cellMotif.csv')) for cell_by_gene_file in cell_by_gene_files]
    return output_files

def main(output_directory: Path):

    conn = sqlite3.connect(database_file)
    cur = conn.cursor()

    non_motif_columns = ['cell_id', 'motif']
    non_gene_columns = ['cell_id', 'expressed_genes', 'overexpressed_genes']

    for cell_by_gene_file, cell_motif_file in get_output_files(directory):

        df1 = get_cell_by_gene_df(cell_by_gene_file)

        df1['expressed_genes'] = ''
        df1['overexpressed_genes'] = ''

        for i, row in df1.iterrows():
            gene_columns = [column for column in df1.columns if column not in non_gene_columns]
            genes_list = [column for column in gene_columns if row[column] > 0.0]
            genes_string = ", ".join(genes_list)
            df1.at[i, 'expressed_genes'] = genes_string
            df1.at[i, 'overexpressed_genes'] = genes_string

        df1 = df1[['cell_id', 'expressed_genes', 'overexpressed_genes']].copy()

        df2 = pd.read_csv(cell_motif_file)
        df2['cell_id'] = df2['Unnamed: 0']
        df2 = df2.drop('Unnamed: 0', axis=1)

        df2['motif'] = ''

        for i, row in df2.iterrows():
            motif_columns = [column for column in df2.columns if column not in non_motif_columns]
            motifs_list = [column for column in motif_columns if row[column] > 0.0]
            motifs_string = ", ".join(motifs_list)
            df2.at[i, 'motif'] = motifs_string

        df2 = df2[['cell_id', 'motif']].copy()

        merge_df = df1.merge(df2)

        merge_df.to_sql('atac', conn, if_exists='replace', index=True)

    conn.close()

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('output_directory', type=Path)
    args = p.parse_args()

    main(args.output_directory)
