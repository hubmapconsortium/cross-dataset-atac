#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.0
label: Pipeline for parsing and aggregating atac output across datasets

inputs:

  data_directories:
    label: "List of containing h5ad data files"
    type: Directory[]

  nexus_token:
    label: "Valid nexus token for search-api"
    type: String

outputs:

  csv_files:
    outputSource: marker-genes/csv_files
    type: File[]

  concatenated_file:
    outputSource: annotate-concatenate/concatenated_file
    type: File

steps:

  - id: annotate-concatenate
    in:
      - id: data_directories
        source: data_directories
      - id: nexus_token
        source: nexus_token

    out:
      - concatenated_file

    run: steps/annotate-concatenate.cwl
    label: "Annotates and concatenates h5ad data files in directory"

  - id: marker-genes
    in:
      - id: concatenated_annotated_file
        source: annotate-concatenate/concatenated_annotated_file
    out:
      - csv_files

    run: steps/marker-genes.cwl
    label: "Finds marker genes associated with different groupings in atac-seq data"
