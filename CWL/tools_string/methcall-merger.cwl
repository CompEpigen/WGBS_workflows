cwlVersion: v1.0
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
inputs:
  input_bed_files:
    items: string
    type: array
  output_file_name: string
baseCommand: tail
arguments:
- position: 1
  valueFrom: '+2'
  prefix: -n
- position: 2
  valueFrom: -q
- position: 100000
  valueFrom: $(inputs.input_bed_files)
stdout: $(inputs['output_file_name'])
outputs:
  merged_bed_file:
    type: File
    outputBinding:
      glob: $(inputs['output_file_name'])
    streamable: true

