cwlVersion: v1.0
class: CommandLineTool
label: Annotates each h5ad file with dataset and tissue type, then concatenates

hints:
  DockerRequirement:
    dockerPull: hubmap/cross-dataset-atac:latest
baseCommand: /opt/parse_sprm_output.py

inputs:
  data_dir:
    type: Directory
    doc: Base directory to be recursively searched for h5ad files to be annotated and concatenated
    inputBinding:
      position: 1

outputs:
  csv_files:
    type: File
    outputBinding:
      glob: "*.csv"
    doc: csvs containing cell level, group level, and gene level data from across atac-seq datasets
