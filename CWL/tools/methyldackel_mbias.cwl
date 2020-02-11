cwlVersion: v1.0
class: CommandLineTool

requirements:
  ShellCommandRequirement: {}
  InlineJavascriptRequirement: {}
hints:
  ResourceRequirement:
    coresMin: $( inputs.threads )
    ramMin: 28000
    tmpdirMin: 30000
  DockerRequirement:
    dockerPull: nfcore/methylseq

baseCommand: ["MethylDackel", "mbias"]
arguments:
  - valueFrom: --txt
    position: 1
  
stdout: $(inputs.output_basename + ".mbias.txt")

inputs:
  output_basename:
    type: string
    inputBinding:
      valueFrom: $(self + ".mbias_plot")
      position: 13
  reference:
    type: File
    secondaryFiles:
      - .fai
      - .bwameth.c2t
      - .bwameth.c2t.amb
      - .bwameth.c2t.ann
      - .bwameth.c2t.bwt
      - .bwameth.c2t.pac
      - .bwameth.c2t.sa
    inputBinding:
      position: 11
  bam:
    type: File
    secondaryFiles: .bai
    inputBinding:
      position: 12

  nOT:
    type: string?
    inputBinding:
      position: 1
      prefix: --nOT
  nOB:
    type: string?
    inputBinding:
      position: 1
      prefix: --nOB
  nCTOT:
    type: string?
    inputBinding:
      position: 1
      prefix: --nCTOT
  nCTOB:
    type: string?
    inputBinding:
      position: 1
      prefix: --nCTOB
  noCG: 
    type: boolean?
    inputBinding:
      position: 1
      prefix: --noCpG

  threads:
    type: int
    default: 1
    inputBinding:
      position: 1
      prefix: "-@"


outputs:
  mbias_output:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*$(inputs.output_basename)*"