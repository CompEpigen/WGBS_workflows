cwlVersion: v1.0
class: CommandLineTool
inputs:
  bam_file:
    doc: "the input bam file\n"
    type: File
    secondaryFiles:
    - ^.bai
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

