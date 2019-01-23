cwlVersion: v1.0
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
inputs:
  inputFileName_mergedSam:
    type:
      items: File
      type: array
    inputBinding:
      position: 100
  outputFileName_mergedSam:
    doc: "SAM or BAM file to write merged result to Required\n"
    type: string
baseCommand: cat
stdout: $(inputs.outputFileName_mergedSam)
outputs:
  mergeSam_output:
    type: File
    outputBinding:
      glob: $(inputs.outputFileName_mergedSam)

