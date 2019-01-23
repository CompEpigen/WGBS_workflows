cwlVersion: v1.0
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
inputs:
  input_bed_file:
    type: string
    inputBinding:
      position: 100000
  output_file_name: string
baseCommand: awk
arguments:
- shellQuote: false
  position: 3
  valueFrom: ${ return "{SUM1 += $4/100; SUM2 +=1} END {print 1-SUM1/SUM2}" }
stdout: $(inputs['output_file_name'])
outputs:
  bisulfite_conversion_file:
    type: File
    outputBinding:
      glob: $(inputs['output_file_name'])
    streamable: true

