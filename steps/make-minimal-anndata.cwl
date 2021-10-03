cwlVersion: v1.0
class: CommandLineTool
label: Finds marker genes associated with different groupings in atac-seq data

hints:
  DockerRequirement:
    dockerPull: hubmap/cross-dataset-scanpy:latest
baseCommand: /opt/make_minimal_anndata.py

inputs:
  concatenated_file:
    type: File
    doc: h5ad file containing batch corrected RNA seq data
    inputBinding:
      position: 1

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
