#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.0
label: Pipeline for parsing and aggregating atac output across datasets

inputs:
  data_dir:
    label: "Directory containing h5ad data files"
    type: Directory

outputs:
  csv_files:
    outputSource: marker-genes/csv_files
    type: File[]

steps:
  - id: annotate-concatenate
    in:
      - id: data_dir
        source: data_dir

    out:
      - concatenated_annotated_file

    run: steps/annotate-concatenate.cwl
    label: "Annotates and concatenates csv and hdf5 files, writes out to h5ad"

  - id: marker-genes
    in:
      - id: concatenated_annotated_file
        source: annotate-concatenate/concatenated_annotated_file
    out:
      - csv_files

    run: steps/marker-genes.cwl
    label: "Finds marker genes associated with different groupings in atac-seq data"
