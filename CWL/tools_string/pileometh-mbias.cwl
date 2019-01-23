cwlVersion: v1.0
class: CommandLineTool
requirements:
  ShellCommandRequirement: {}
inputs:
  bam_file:
    doc: "the input bam file\n"
    type: string
    inputBinding:
      position: 5
  reference:
    doc: "FASTA file with the reference genome\n"
    type: string
    inputBinding:
      position: 4
  mbiasfile_name:
    doc: "FASTA file \n"
    type: string
    inputBinding:
      position: 1000000
      separate: true
  min_mapq:
    doc: "min_mapq\n"
    type: int
    default: 0
    inputBinding:
      position: 100000
      prefix: -q
      separate: true
  min_phred:
    doc: "min_phred\n"
    type: int
    default: 1
    inputBinding:
      position: 100000
      prefix: -p
      separate: true
  min_depth:
    doc: "min_depth\n"
    type: int
    default: 1
    inputBinding:
      position: 100000
      prefix: -D
      separate: true
baseCommand: /ngs_share/tools/PileOMeth/PileOMeth
arguments:
- position: 2
  valueFrom: mbias
- shellQuote: false
  position: 3
  valueFrom: --txt
stdout: $(inputs.mbiasfile_name).txt
outputs:
  mbias_file:
    type:
      items: File
      type: array
    outputBinding:
      glob: '*$(inputs.mbiasfile_name)*'

