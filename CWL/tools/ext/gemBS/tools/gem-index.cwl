#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: gem-indexer
#hints:
#  DockerRequirement:
#    dockerPull: yylin/bs_call:1.0

requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    ramMin: 15000
    
## Mandatory INPUTs
inputs: 
  input: 
    doc: Multi-FASTA
    type: File
    inputBinding:
      position: 1
      prefix: -i
  output:
    doc: output file name
    type: string
    inputBinding:
      position: 2
      prefix: -o

## Optional Index args
  bisulfite-index:
    doc: bisulfite-index
    type: boolean?
    inputBinding:
      position: 3
      prefix: -b

## Optional System args
  threads:
    doc: set threads 
    type: int?
    inputBinding:
      position: 4
      prefix: -t

## OUTPUT PART
outputs: 
  gem:
     type: File
     outputBinding:
        glob: "*.gem"