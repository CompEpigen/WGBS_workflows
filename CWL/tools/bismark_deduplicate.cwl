
cwlVersion: v1.0
class: CommandLineTool
hints:
  ResourceRequirement:
    coresMin: 4
    ramMin: 20000
    tmpdirMin: 30000
  DockerRequirement:
    dockerPull: kerstenbreuer/bismark:0.22.3

baseCommand: deduplicate_bismark
arguments:
  - valueFrom: |
      ${
        if (inputs.paired_end){
          return("-p")
        }
        else {
          return("-s")
        }
      }

inputs:
  aligned_reads:
    type: File
    inputBinding:
      prefix: --bam
      position: 10
  paired_end:
    type: boolean
    default: true
outputs:
  dedup_reads:
    type: File
    outputBinding:
      glob: "*.deduplicated.bam"
  log:
    type: File
    outputBinding:
      glob: "*.deduplication_report.txt"
