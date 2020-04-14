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
    valueFrom: $(inputs.bams_doublex1.nameroot)_$(inputs.bams_doublex2.nameroot).bam
    position: 7
  - valueFrom: "-"
    position: 8


inputs:
  bams_doublex1:
    type: File
    inputBinding:
      position: 2
  bams_doublex2:
    type: File
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
      glob: $(inputs.bams_doublex1.nameroot)_$(inputs.bams_doublex2.nameroot).bam
    