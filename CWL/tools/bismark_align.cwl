
cwlVersion: v1.0
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
hints:
  ResourceRequirement:
    coresMin: $(inputs.threads)
    ramMin: ${return(Math.ceil(inputs.threads/5)*14000)}
    tmpdirMin: 30000
  DockerRequirement:
    dockerPull: kerstenbreuer/bismark:0.22.3

baseCommand: bismark
arguments:
  - valueFrom: --bam

inputs:
  # main input
  genome:
    type: Directory
    inputBinding:
      prefix: --genome
      position: 10
  read1:
    type: File
    inputBinding:
      prefix: "-1"
      position: 11
  read2:
    type: File
    inputBinding:
      prefix: "-2"
      position: 12
  pbat:
    type: boolean
    inputBinding:
      prefix: --pbat
  threads:
    type: int
    default: 1
    inputBinding:
      prefix: --multicore
      valueFrom: ${return(Math.ceil(self/5))}
      position: 1
  bismark_local:
    type: boolean
    inputBinding:
      prefix: --local
      position: 1
  non_directional:
    type: boolean
    inputBinding:
      prefix: --non_directional
      position: 1
  dovetail:
    type: boolean
    inputBinding:
      prefix: --dovetail
      position: 1

outputs:
  aligned_reads:
    type: File
    outputBinding:
      glob: "*.bam"
  log:
    type: File
    outputBinding:
      glob: "*.txt"
