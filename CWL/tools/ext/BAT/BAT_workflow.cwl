#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
inputs:
  reference: File
  query: File
  mate_pair: File
  prefix_db: string
  prefix_location: Directory
  path_outfiles: string

outputs:
  bam:
    type: File
    outputSource: mapping/bam
  vcf_file:
    type: File
    outputSource: calling/vcf

steps:
  mapping:
    run: BAT_mapping.cwl
    in:
      reference: reference
      query: query
      mate_pair: mate_pair
      prefix_db: prefix_db
      prefix_location: prefix_location
      path_outfiles: path_outfiles
    out: [bam]

  calling:
    run: BAT_calling.cwl
    in: 
      reference: reference
      query: mapping/bam
    out: [vcf]