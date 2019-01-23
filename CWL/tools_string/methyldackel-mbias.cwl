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
  OT:
    doc: "OT \n"
    type: string
    inputBinding:
      position: 100000
      prefix: --OT
      separate: true
  OB:
    doc: "OB \n"
    type: string
    inputBinding:
      position: 100000
      prefix: --OB
      separate: true
  CTOT:
    doc: "CTOT \n"
    type: string
    inputBinding:
      position: 100000
      prefix: --CTOT
      separate: true
  CTOB:
    doc: "CTOB \n"
    type: string
    inputBinding:
      position: 100000
      prefix: --CTOB
      separate: true
  nOT:
    doc: "nOT \n"
    type: string
    inputBinding:
      position: 100000
      prefix: --nOT
      separate: true
  nOB:
    doc: "nOB \n"
    type: string
    inputBinding:
      position: 100000
      prefix: --nOB
      separate: true
  nCTOT:
    doc: "nCTOT \n"
    type: string
    inputBinding:
      position: 100000
      prefix: --nCTOT
      separate: true
  nCTOB:
    doc: "nCTOB \n"
    type: string
    inputBinding:
      position: 100000
      prefix: --nCTOB
      separate: true
baseCommand: /ngs_share/tools/MethylDackel_dev/MethylDackel
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

