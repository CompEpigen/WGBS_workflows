#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: ["gsnap", "-A", "sam", "--gunzip"] 
requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: gsnap
    dockerOutputDirectory: /data

inputs: 
  ref:
    doc: name of the db index. (built in image with hg19 and hg38)
    type: string
    inputBinding:
      prefix: -d
      position: 1
  r1:
    doc: read1.fa.gz
    type: File
    inputBinding:
      position: 2
  r2:
    doc: read2.fa.gz
    type: File
    inputBinding: 
      position: 3  
  output_name:
    doc: prefix of the output name 
    type: string?
    default: "gsnap_aligned"

## OUTPUT PART      
outputs: 
  sam:
    type: stdout
  log:
    type: stderr

stdout: $(inputs.output_name).sam      
stderr: $(inputs.output_name).log

           