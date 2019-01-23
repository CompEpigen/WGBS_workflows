inputs:
- inputBinding:
    position: 1
  type: File
  id: '#input_bam_file'
  description: Aligned BAM file to filter
- type: string
  id: '#output_suffix'
requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement
description: Extract mapped reads from BAM file using Samtools flagstat command
stdout: $(inputs.input_bam_file.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/,
  "") + inputs.output_suffix)
outputs:
- outputBinding:
    glob: $(inputs.input_bam_file.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/,
      "") + inputs.output_suffix)
  type: File
  id: '#output_read_count'
  description: Samtools Flagstat report file
baseCommand:
- samtools
- flagstat
arguments:
- shellQuote: false
  position: 10000
  valueFrom: " | head -n1 | cut -f 1 -d ' '"
class: CommandLineTool
hints:
- dockerImageId: dukegcb/samtools
  class: DockerRequirement

