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
  hdf5_file:
    type: File
    outputBinding:
      glob: "atac.hdf5"
    doc: hdf5 file with layers containing dataframes for cell, group, and quant data

  csv_file:
    type: File
    outputBinding:
      glob: "long_atac_quant.csv"
    doc: csv file containing long narrow expression data

  minimal_files:
    type: File[]
    outputBinding:
      glob: "mini*"
