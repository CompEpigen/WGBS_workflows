cwlVersion: v1.0
class: CommandLineTool
doc: Extract mapped reads from BAM file using Samtools flagstat command
requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}
hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 10000
inputs:
  input_bam_file:
    doc: Aligned BAM file to filter
    type: string
    inputBinding:
      position: 1
  output_suffix:
    type: string
    default: .flagStat
baseCommand:
- samtools
- flagstat
stdout: $(inputs.input_bam_file.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/,
  "") + inputs.output_suffix)
outputs:
  output:
    doc: Samtools Flagstat report file
    type: File
    outputBinding:
      glob: $(inputs.input_bam_file.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/,
        "") + inputs.output_suffix)

