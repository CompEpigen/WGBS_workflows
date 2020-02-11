cwlVersion: v1.0
class: CommandLineTool

requirements:
  ShellCommandRequirement: {}
  InlineJavascriptRequirement: {}
hints:
  ResourceRequirement:
    coresMin: $( inputs.threads )
    ramMin: 28000
    tmpdirMin: 30000
  DockerRequirement:
    dockerPull: kerstenbreuer/bwameth:latest

baseCommand: ["bwameth.py"]

stdout: $(inputs.fastq1.nameroot).sam

inputs:
  reference:
    doc: the reference fasta file location
    type: File
    secondaryFiles:
      - .fai
      - .bwameth.c2t
      - .bwameth.c2t.amb
      - .bwameth.c2t.ann
      - .bwameth.c2t.bwt
      - .bwameth.c2t.pac
      - .bwameth.c2t.sa
    inputBinding:
      position: 10
      prefix: --reference
      separate: true
  fastq1:
    doc: the input fastq file with the first mate
    type: File
    inputBinding:
      position: 11
  fastq2:
    doc: the input fastq file with the second mate
    type: File?
    inputBinding:
      position: 12
  threads:
    type: int?
    default: 1
    inputBinding:
      position: 1
      prefix: --threads
  is_non_directional:
    doc: Is library type non-directional
    type: boolean
    default: false
    inputBinding:
      position: 1
      prefix: --non-directional

outputs:
  sam:
    type: stdout

