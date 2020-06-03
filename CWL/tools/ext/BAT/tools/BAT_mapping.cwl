#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [BAT_mapping] 
requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: bat 
    dockerOutputDirectory: /opt
  InitialWorkDirRequirement:
    listing:
      - $(inputs.prefix_location)

inputs: 
  reference:
    type: File
    doc: path/filename of reference genome fasta
    inputBinding:
      prefix: -g
      position: 1
  query:
    type: File
    doc: path/filename of query sequences
    inputBinding:
      prefix: -q
      position: 2
  mate_pair:
    type: File?
    doc: path/filename of mate pair sequences
    inputBinding:
      prefix: -p   
      position: 3  
  prefix_db:
    type: string  
    doc: path/prefix of database indices
    inputBinding:
      prefix: -i  
      position: 4  
  prefix_location:
    type: Directory
  path_outfiles:
    type: string  
    doc: path/prefix of outfiles   
    inputBinding:
      prefix: -o 
      position: 5

## OUTPUT PART      
outputs:           
  bam:
    type: File
    outputBinding:
      glob: $(inputs.path_outfiles + ".bam")
  bam.bai:
    type: File
    outputBinding:
      glob: $(inputs.path_outfiles + ".bam.bai")      
  excluded_bam:
    type: File
    outputBinding:
      glob: $(inputs.path_outfiles + ".excluded.bam") 
  excluded_bam.bai:
    type: File
    outputBinding:
      glob: $(inputs.path_outfiles + ".excluded.bam.bai")
  log:
    type: File
    outputBinding:
      glob: "*.log" 
  unmapped.gz:
    type: File
    outputBinding:
      glob: $(inputs.path_outfiles + ".unmapped.gz")             
           
