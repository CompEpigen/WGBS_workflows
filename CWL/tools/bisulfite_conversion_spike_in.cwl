cwlVersion: v1.0
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}

hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 4000
    tmpdirMin: 10000
  DockerRequirement:
    dockerPull: kerstenbreuer/samtools:1.7

baseCommand: awk
arguments:
- shellQuote: false
  position: 1
  valueFrom: ${ return "{SUM1 += $4/100; SUM2 +=1} END {print 1-SUM1/SUM2}" }

stdout: $(inputs.output_basename + ".bs_convers.txt")

inputs:
  spike_in_mcall_bedgraph:
    type: File
    inputBinding:
      position: 10
  output_basename: string

outputs:
  bisulfite_conversion_file:
    type: File
    outputBinding:
      glob: $(inputs.output_basename + ".bs_convers.txt")

