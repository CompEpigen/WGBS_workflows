cwlVersion: v1.0
class: CommandLineTool

hints:
  ResourceRequirement:
    coresMin: 8
    ramMin: 28000
  DockerRequirement:
    dockerPull: kerstenbreuer/fame:0.2_8cores_125bp_0minpdist_1000maxpdist_22chr
requirements:
  InlineJavascriptRequirement: {}  
  InitialWorkDirRequirement:
    listing: |
      ${
        var lines = []
        var line
        for (var i=0; i < inputs.sc_id.length; i++){
          line = inputs.sc_id[i] + " " + inputs.fastq1[i].path
          if (inputs.fastq2){
            line += " " + inputs.fastq2[i].path
          }
          lines.push(line)
        }
        var file = [{
          class: "File",
          basename: "meta.tsv",
          contents: lines.join('\u000a')
        }]

        return(file)
      } 

baseCommand: ["FAME"]
arguments:
  - valueFrom: |
      ${
        if(inputs.fastq2){
          return("--paired")
        }
        else{
          return(null)
        }
      }
  - prefix: --sc_input
    valueFrom: meta.tsv
  - valueFrom: "--unord_reads"

inputs:
  load_index:
    type: File
    secondaryFiles:
      - _strands
    inputBinding:
      prefix: --load_index 
  sc_id:
    type: string[]
  fastq1:
    type: File[]
  fastq2:
    type: 
      - "null"
      - type: array
        items: File
  gzip_reads:
    type: boolean
    default: true
    inputBinding:
      prefix: --gzip_reads
  sc_output:
    type: string
    default: sc_methylation_calls
    inputBinding:
      prefix: --sc_output
  out_basename:
    type: string
    default: bulk_methylation_calls
    inputBinding:
      prefix: --out_basename
    

outputs:
  sc_methylation_calls:
    type: File
    outputBinding:
      glob: $(inputs.sc_output)*
  bulk_methylation_calls:
    type: File
    outputBinding:
      glob: $(inputs.out_basename)*
  meta_tsv:
    type: File
    outputBinding:
      glob: meta.tsv