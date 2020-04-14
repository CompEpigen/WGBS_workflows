
cwlVersion: v1.0
class: CommandLineTool
hints:
  ResourceRequirement:
    coresMin: 4
    ramMin: 20000
    tmpdirMin: 10000
  DockerRequirement:
    dockerPull: kerstenbreuer/bismark:0.22.3

baseCommand: bismark2report

inputs:
  alignment_report:
    type: File?
    inputBinding:
      prefix: --alignment_report
      position: 1
  dedup_report:
    type: File?
    inputBinding:
      prefix: --dedup_report
      position: 1
  splitting_report:
    type: File?
    inputBinding:
      prefix: --splitting_report
      position: 1
  mbias_report:
    type: File?
    inputBinding:
      prefix: --mbias_report
      position: 1
outputs:
  report:
    type: File
    outputBinding:
      glob: "*.html"
