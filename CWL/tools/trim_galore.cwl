
cwlVersion: v1.0
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
hints:
  ResourceRequirement:
    coresMin: $(Math.min(inputs.threads, 8))
    ramMin: 7000
    tmpdirMin: 10000
  DockerRequirement:
    dockerPull: kerstenbreuer/trim_galore:0.6.4_2.6_0.11.8

baseCommand: trim_galore

arguments:
  - valueFrom: --gzip
    position: 1
    # gzip output fastq
  - valueFrom: --paired
    position: 1
  
  ## variable arguments:
  - valueFrom: |
      ${
        if ( inputs.adapter1 == "illumina" ){ return "--illumina" }
        else if ( inputs.adapter1 == "nextera" ){ return "--nextera" }
        else if ( inputs.adapter1 == "small_rna" ){ return "--small_rna" }
        else { return null }
      }
    position: 1
  - prefix: --adapter
    valueFrom: |
      ${
        if ( inputs.apdater1 != null && inputs.adapter1 != "illumina" && inputs.adapter1 != "nextera" && inputs.adapter1 != "small_rna" ){
          return inputs.adapter1
        } else {
          return null
        }
      }
    position: 1
  - prefix: --adapter2
    valueFrom: |
      ${
        if (inputs.apdater2 != null && inputs.adapter1 != "illumina" && inputs.adapter1 != "nextera" && inputs.adapter1 != "small_rna" ){
          return inputs.adapter2
        } else {
          return null
        }
      }
    position: 1

inputs:
  # main input
  read1:
    type: File
    inputBinding:
      position: 10
  read2:
    type: File
    inputBinding:
      position: 11
  adapter1:
    doc: |
      Adapter sequence for first reads.
      if not specified, trim_galore will try to autodetect whether ...
      - Illumina universal adapter (AGATCGGAAGAGC)
      - Nextera adapter (CTGTCTCTTATA)
      - Illumina Small RNA 3' Adapter (TGGAATTCTCGG)
      ... was used.
      You can directly choose one of the above configurations
      by setting the string to "illumina", "nextera", or "small_rna".
    type: string?
  adapter2:
    doc: |
      Adapter sequence for second reads.
      if not specified, trim_galore will try to autodetect whether ...
      - Illumina universal adapter (AGATCGGAAGAGC)
      - Nextera adapter (CTGTCTCTTATA)
      - Illumina Small RNA 3' Adapter (TGGAATTCTCGG)
      ... was used.
      You can directly choose one of the above configurations
      by setting the adapter1 string to "illumina", "nextera", or "small_rna".
    type: string?
  quality:
    type: int
    default: 20
    inputBinding:
      position: 1
      prefix: --quality
  rrbs:
    type: boolean
    inputBinding:
      position: 1
      prefix: --rrbs
  clip_r1:
    type: int?
    inputBinding:
      position: 1
      prefix: --clip_r1
  clip_r2:
    type: int?
    inputBinding:
      position: 1
      prefix: --clip_r2
  three_prime_clip_r1:
    type: int?
    inputBinding:
      position: 1
      prefix: --three_prime_clip_r1
  three_prime_clip_r2:
    type: int?
    inputBinding:
      position: 1
      prefix: --three_prime_clip_r2
  threads:
    type: int
    default: 1
    inputBinding:
      valueFrom: $(Math.min(self, 8))
      prefix: --cores
      position: 1

outputs:
  read1_trimmed:
    type: File
    outputBinding:
      glob: "*val_1.fq.gz"
  read2_trimmed:    
    type: File
    outputBinding:
      glob: "*val_2.fq.gz"
  log:
    type: File[]
    outputBinding:
      glob:  "*trimming_report.txt"
