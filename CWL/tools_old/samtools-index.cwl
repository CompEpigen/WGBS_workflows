cwlVersion: v1.0
class: CommandLineTool
requirements:
  InlineJavascriptRequirement:
    expressionLib:
    - var new_ext = function() { var ext=inputs.bai?'.bai':inputs.csi?'.csi':'.bai';
      return inputs.input.path.split('/').slice(-1)[0]+ext; };
inputs:
  input:
    doc: "Input bam file.\n"
    type: File
    inputBinding:
      position: 2
  bai:
    doc: "Generate BAI-format index for BAM files [default]\n"
    type: boolean
    default: false
  csi:
    doc: "Generate CSI-format index for BAM files\n"
    type: boolean
    default: false
  interval:
    doc: "Set minimum interval size for CSI indices to 2^INT [14]\n"
    type: int?
    inputBinding:
      position: 1
      prefix: -m
baseCommand:
- samtools
- index
arguments:
- position: 1
  valueFrom: $(inputs.bai?'-b':inputs.csi?'-c':[])
- position: 3
  valueFrom: $(new_ext())
outputs:
  index:
    doc: The index file
    type: File
    outputBinding:
      glob: $(new_ext())

