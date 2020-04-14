
cwlVersion: v1.0
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
hints:
  ResourceRequirement:
    coresMin: $(inputs.threads)
    ramMin: ${return(Math.ceil(inputs.threads/5)*14000)}
  DockerRequirement:
    dockerPull: kerstenbreuer/bismark:0.22.3

baseCommand: bismark
arguments:
  - valueFrom: --bam
    position: 1
  - valueFrom: --non_directional
    position: 1

inputs:
  # main input
  genome:
    type: Directory
    inputBinding:
      prefix: --genome
      position: 10
  reads:
    type: File
    inputBinding:
      position: 11
  threads:
    type: int
    default: 1
    inputBinding:
      prefix: --multicore
      valueFrom: ${return(Math.ceil(self/5))}
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
