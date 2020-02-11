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

baseCommand: ["MethylDackel", "extract"]
inputs:
  output_basename:
    type: string
    inputBinding:
      valueFrom: $(self + ".mcall")
      position: 10
      prefix: -o
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


  min_mapq:
    type: int
    default: 0
    inputBinding:
      position: 1
      prefix: -q
  min_phred:
    type: int
    default: 1
    inputBinding:
      position: 1
      prefix: -p
  min_depth:
    type: int
    default: 1
    inputBinding:
      position: 1
      prefix: --minDepth
  OT:
    type: string?
    inputBinding:
      position: 1
      prefix: --OT
  OB:
    type: string?
    inputBinding:
      position: 1
      prefix: --OB
  CTOT:
    type: string?
    inputBinding:
      position: 1
      prefix: --CTOT
  CTOB:
    type: string?
    inputBinding:
      position: 1
      prefix: --CTOB
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
    type: boolean
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
  mcall_bedgraph:
    type: File
    outputBinding:
      glob: "*.bedGraph"