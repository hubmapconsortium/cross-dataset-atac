cwlVersion: v1.0
class: CommandLineTool
label: Annotates each h5ad file with dataset and tissue type, then concatenates

hints:
  DockerRequirement:
    dockerPull: hubmap/cross-dataset-atac:latest
baseCommand: /opt/annotate_concatenate.py

inputs:
  data_dir:
    type: Directory
    doc: Base directory to be recursively searched for h5ad files to be annotated and concatenated
    inputBinding:
      position: 1

outputs:
  concatenated_annotated_file:
    type: File
    doc: h5ad file containing annotated and concatenated atac-seq data
    outputBinding:
      glob: 'concatenated_annotated.h5ad'
