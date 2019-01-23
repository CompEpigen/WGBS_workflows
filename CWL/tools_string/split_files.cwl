cwlVersion: v1.0
class: CommandLineTool
inputs:
  file:
    doc: "the input fastq file\n"
    type: File
    inputBinding:
      position: 3
  size:
    doc: "size of a chunck in lines\n"
    type: int
    inputBinding:
      position: 1
      prefix: --lines=
      separate: false
  suffix:
    doc: "a suffix for chuncks\n"
    type: string
    inputBinding:
      position: 2
      prefix: --additional-suffix=
      separate: false
baseCommand: split
outputs:
  output:
    type:
      items: File
      type: array
    outputBinding:
      glob: '*$(inputs.suffix)*'

