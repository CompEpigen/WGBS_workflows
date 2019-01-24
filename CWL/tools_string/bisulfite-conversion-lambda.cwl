cwlVersion: v1.0
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 4000
  DockerRequirement:
    dockerPull: nfcore/methylseq

baseCommand: awk
arguments:
- shellQuote: false
  position: 1
  valueFrom: ${ return "{SUM1 += $4/100; SUM2 +=1} END {print 1-SUM1/SUM2}" }

stdout: $(inputs.output_file_name)

inputs:
  input_bed_file:
    type: string
    inputBinding:
      position: 10
  output_file_name: string

outputs:
  bisulfite_conversion_file:
    type: File
    outputBinding:
      glob: $(inputs.output_file_name)
    streamable: true

