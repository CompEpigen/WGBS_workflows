cwlVersion: v1.0
class: CommandLineTool
inputs:
  file_dir: string
  file:
    doc: "the input log file\n"
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
baseCommand: cd
arguments:
- position: 1
  valueFrom: $(inputs.file_dir)
- position: 2
  valueFrom: ;
- position: 100
  valueFrom: cat
- position: 1001
  valueFrom: '>'
- position: 1002
  valueFrom: $(inputs.file_dir + "/" + inputs['alignment_filename']).sam
outputs:
  output:
    type:
      items: File
      type: array
    outputBinding:
      glob: '*$(inputs.suffix)*'

