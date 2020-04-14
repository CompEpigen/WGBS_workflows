
cwlVersion: v1.0
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
hints:
  ResourceRequirement:
    coresMin: $(inputs.threads)
    ramMin: ${return(inputs.threads*30000)}
    tmpdirMin: 30000
  DockerRequirement:
    dockerPull: kerstenbreuer/bismark:0.22.3

baseCommand: bismark_methylation_extractor
arguments:
  - valueFrom: --bedGraph
    position: 1
  - valueFrom: --gzip
    position: 1
  - valueFrom: --counts
    position: 1
  - valueFrom: |
      ${
        if (inputs.paired_end){
          return("-p")
        }
        else {
          return("-s")
        }
      }
    position: 1
  - valueFrom: --report
    position: 1
  - valueFrom: --cytosine_report
    position: 1

inputs:
  # main input
  aligned_reads:
    type: File
    inputBinding:
      position: 10
  genome:
    type: Directory
    inputBinding:
      prefix: --genome
      position: 10
  no_overlap:
    type: boolean
    default: true
    inputBinding:
      prefix: --no_overlap
      position: 2
  paired_end:
    type: boolean
    default: true
  ignore:
    type: int?
    inputBinding:
      prefix: --ignore
      position: 3
  ignore_r2:
    type: int?
    inputBinding:
      prefix: --ignore_r2
      position: 3
  ignore_3prime:
    type: int?
    inputBinding:
      prefix: --ignore_3prime
      position: 3
  ignore_3prime_r2:
    type: int?
    inputBinding:
      prefix: --ignore_3prime_r2
      position: 3
  threads:
    type: int
    default: 1
    inputBinding:
      prefix: --multicore
      valueFrom: ${return(Math.ceil(self/5))}
      position: 1

outputs:
  methylation_calls_bedgraph:
    type: File
    outputBinding:
      glob: "*bedGraph.gz"
  methylation_calls_bismark:
    type: File
    outputBinding:
      glob: "*bismark.cov.gz"
  mbias_report:
    type: File
    outputBinding:
      glob: "*.M-bias.txt*"
  splitting_report:
    type: File
    outputBinding:
      glob: "*splitting_report.txt*"
  genome_wide_methylation_report:
    type: File
    outputBinding:
      glob: "*CpG_report.txt*"
  context_specific_methylation_reports:
    type: File[]
    outputBinding:
      glob: "C*_O*.txt*"
  
