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
  file_dir: string
  input_bam_file:
    doc: Aligned BAM file to filter
    type: string
    inputBinding:
      position: 10
  output_suffix:
    type: string
    default: .flagStat
baseCommand: cd
arguments:
- position: 1
  valueFrom: $(inputs.file_dir)
- position: 2
  valueFrom: ;
- position: 3
  valueFrom: samtools
- position: 4
  valueFrom: flagstat
- position: 5
  valueFrom: '>'
- position: 6
  valueFrom: $(inputs.file_dir + "/" + inputs.input_bam_file.replace(/^.*[\\\/]/,
    '').replace(/\.[^/.]+$/, "") + inputs.output_suffix)
outputs:
  output:
    doc: Samtools Flagstat report file
    type: string
    outputBinding:
      outputEval: $(inputs.file_dir + "/" + inputs.input_bam_file.replace(/^.*[\\\/]/,
        '').replace(/\.[^/.]+$/, "") + inputs.output_suffix)

