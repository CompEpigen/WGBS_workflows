#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [] #
requirements:
  DockerRequirement:
    dockerPull: gemBS_main 
      # ./docker_ext/gemBS_main
    dockerOutputDirectory: /opt
  InitialWorkDirRequirement:
    listing:
      - $(inputs.fastq)
      - $(inputs.ref)
      - $(inputs.conf)
      - $(inputs.csv)

inputs: 
  conf:
    type: File
    doc: configuration
  csv:
    type: File
    doc: list of inputs datasets
  ref:
    type: Directory  
    doc: directory of reference
  fastq:
    type: Directory  
    doc: directory of fastq
outputs:           
  report:
    type:
      type: array
      items: Directory
    outputBinding:
      glob: "report"
  indexes:
    type:
      type: array
      items: Directory
    outputBinding:
      glob: "indexes"
  mapping:
    type:
      type: array
      items: Directory
    outputBinding:
      glob: "mapping"         
  calls:
    type:
      type: array
      items: Directory
    outputBinding:
      glob: "calls" 
  exract:
    type:
      type: array
      items: Directory
    outputBinding:
      glob: "extract"                       