cwlVersion: v1.0
class: CommandLineTool
requirements:
  ShellCommandRequirement: {}
hints:
  ResourceRequirement:
    coresMin: 1
    diskMin: 15000
    ramMin: 15000
inputs:
  file:
    doc: "the input fastq file\n"
    type: File
    inputBinding:
      position: 2
  size:
    doc: "size of a chunck in lines\n"
    type: int
    inputBinding:
      position: 10004
      prefix: --lines=
      separate: false
  suffix:
    doc: "a suffix for chuncks\n"
    type: string
    inputBinding:
      position: 10003
      prefix: --additional-suffix=
      separate: false
baseCommand: zcat
arguments:
- shellQuote: false
  position: 10000
  valueFrom: '|'
- position: 10001
  valueFrom: split
- shellQuote: false
  position: 10002
  valueFrom: '-'
outputs:
  output:
    type:
      items: File
      type: array
    outputBinding:
      glob: '*$(inputs.suffix)*'

