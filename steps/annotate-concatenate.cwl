cwlVersion: v1.0
class: CommandLineTool
label: Annotates each h5ad file with dataset and tissue type, then concatenates

hints:
  DockerRequirement:
    dockerPull: hubmap/cross-dataset-atac:latest
baseCommand: /opt/annotate_concatenate.py

inputs:

  nexus_token:
    type: string?
    doc: Valid nexus token for search-api
    inputBinding:
      position: 1
    default: "None"

  data_directories:
    type: Directory[]
    doc: List of paths to processed dataset directories
    inputBinding:
      position: 2

outputs:
  concatenated_annotated_file:
    type: File
    doc: h5ad file containing annotated and concatenated atac-seq data
    outputBinding:
      glob: 'concatenated_annotated.h5ad'

  gene_dictionaries:
    type: File[]
    doc: json files mapping from gene symbols to versioned ensembl_ids and back
    outputBinding:
      glob: '*.json'

  old_cluster_file:
    type: File
    outputBinding:
      glob: "cluster.hdf5"
    doc: Hdf file containing old cluster p values