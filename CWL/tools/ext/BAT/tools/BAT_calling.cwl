#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [BAT_calling] 
requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: bat 
  InitialWorkDirRequirement:
    listing:
      - $(inputs.query)

inputs: 
  reference:
    type: File
    doc: path/filename of reference genome fasta
    inputBinding:
      prefix: -d
  query:
    type: File
    doc: path/filename of query sequences
    inputBinding:
      prefix: -q   

## OUTPUT PART      
outputs:           
  log:
    type: File
    outputBinding:
      glob: "*.calling.log"
  vcf:
    type: File
    outputBinding:
      glob: "*.vcf.gz"   
