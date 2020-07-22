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
    outputSource: join-annotate/csv_files
    type: File[]

steps:
  - id: join-annotate
    in:
      - id: data_dir
        source: data_dir

    out:
      - csv_files

    run: steps/join-annotate.cwl
    label: "Annotates and concatenates csv and hdf5 files, writes out csv_files"
