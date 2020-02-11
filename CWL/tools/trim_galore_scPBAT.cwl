
cwlVersion: v1.0
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
hints:
  ResourceRequirement:
    coresMin: $(Math.min(inputs.threads, 8))
    ramMin: 7000
    tmpdirMin: 10000
  DockerRequirement:
    dockerPull: kerstenbreuer/trim_galore:0.6.4_2.6_0.11.8_scPBAT

baseCommand: trim_galore

arguments:
  - valueFrom: --gzip
    position: 1

inputs:
  # main input
  reads:
    type: File
    inputBinding:
      position: 10
  quality:
    type: int
    default: 20
    inputBinding:
      position: 1
      prefix: --quality
  clip_r1:
    type: int?
    default: 6
    inputBinding:
      position: 1
      prefix: --clip_r1
  three_prime_clip_r1:
    type: int?
    default: 6
    inputBinding:
      position: 1
      prefix: --three_prime_clip_r1
  threads:
    type: int?
    inputBinding:
      valueFrom: $(Math.min(self, 8))
      prefix: --cores
      position: 1

outputs:
  reads_trimmed:
    type: File
    outputBinding:
      glob: "*trimmed.fq.gz"
  log:
    type: File[]
    outputBinding:
      glob:  "*trimming_report.txt"
