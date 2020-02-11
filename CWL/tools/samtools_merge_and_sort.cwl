cwlVersion: v1.0
class: CommandLineTool
requirements:
  ShellCommandRequirement: {}
hints:
  ResourceRequirement:
    coresMin: $(inputs.threads)
    ramMin: 20000
    tmpdirMin: 30000
  DockerRequirement:
    dockerPull: kerstenbreuer/samtools:1.7

baseCommand: ["samtools", "merge"]
arguments:
  - valueFrom: "-"
    position: 1
  - valueFrom: "|"
    position: 3
    shellQuote: false
  - valueFrom: samtools
    position: 4
  - valueFrom: sort
    position: 5
  - prefix: "-@"
    valueFrom: $(inputs.threads)
    position: 6
  - prefix: "-o"
    valueFrom: $(inputs.output_name)
    position: 7
  - valueFrom: "-"
    position: 8


inputs:
  output_name:
    type: string
    default: "merged_reads.bam"
  bams:
    doc: bam files to be merged
    type:
      type: array
      items: File
    inputBinding:
      position: 2
  name_sort:
    type: boolean
    default: false
    inputBinding:
      prefix: -n
      position: 7
  threads:
    type: int
    default: 1

outputs:
  - id: bam_merged
    type: File
    outputBinding:
      glob: $(inputs.output_name)
    