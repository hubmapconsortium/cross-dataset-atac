cwlVersion: v1.0
class: CommandLineTool
label: Finds marker genes associated with different groupings in atac-seq data

hints:
  DockerRequirement:
    dockerPull: hubmap/cross-dataset-atac:latest
baseCommand: /opt/marker_genes.py

inputs:
  concatenated_annotated_file:
    type: File
    doc: h5ad file containing annotated and concatenated atac-seq data
    inputBinding:
      position: 1

outputs:
  csv_files:
    type: File[]
    outputBinding:
      glob: "*.csv"
    doc: csvs containing cell level, group level, and gene level data from across atac-seq datasets
