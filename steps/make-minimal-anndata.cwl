cwlVersion: v1.0
class: CommandLineTool
label: Finds marker genes associated with different groupings in atac-seq data

hints:
  DockerRequirement:
    dockerPull: hubmap/cross-dataset-scanpy:latest
baseCommand: /opt/make_minimal_anndata.py

inputs:
  csv_file:
    type: File
    doc: csv file containing cell by gene ATAC seq data
    inputBinding:
      position: 1

  hdf_file:
    type: File
    doc: hdf file containing p value ATAC seq data
    inputBinding:
      position: 2

outputs:
  h5ad_file:
    type: File
    outputBinding:
      glob: "atac.h5ad"
    doc: h5ad containing minimal anndata representation of numeric data

  pval_h5ad_file:
    type: File
    outputBinding:
      glob: "atac_pvals.h5ad"
    doc: h5ad containing anndata representation of pval data
