cwlVersion: v1.0
class: CommandLineTool
inputs:
  bam_file:
    doc: "the input bam file\n"
    type: File
    inputBinding:
      position: 3
  reference:
    doc: "FASTA file with the reference genome\n"
    type: string
    inputBinding:
      position: 2
  bedfile_name:
    doc: "FASTA file \n"
    type: string
    inputBinding:
      position: 3
      prefix: -o
      separate: true
  min_mapq:
    doc: "min_mapq \n"
    type: int
    default: 0
    inputBinding:
      position: 100000
      prefix: -q
      separate: true
  min_phred:
    doc: "min_phred \n"
    type: int
    default: 1
    inputBinding:
      position: 100000
      prefix: -p
      separate: true
  min_depth:
    doc: "min_depth \n"
    type: int
    default: 1
    inputBinding:
      position: 100000
      prefix: --minDepth
      separate: true
  OT:
    doc: "OT \n"
    type: string
    inputBinding:
      position: 100000
      prefix: --OT
      separate: true
  OB:
    doc: "OB \n"
    type: string
    inputBinding:
      position: 100000
      prefix: --OB
      separate: true
  CTOT:
    doc: "CTOT \n"
    type: string
    inputBinding:
      position: 100000
      prefix: --CTOT
      separate: true
  CTOB:
    doc: "CTOB \n"
    type: string
    inputBinding:
      position: 100000
      prefix: --CTOB
      separate: true
  nOT:
    doc: "nOT \n"
    type: string
    inputBinding:
      position: 100000
      prefix: --nOT
      separate: true
  nOB:
    doc: "nOB \n"
    type: string
    inputBinding:
      position: 100000
      prefix: --nOB
      separate: true
  nCTOT:
    doc: "nCTOT \n"
    type: string
    inputBinding:
      position: 100000
      prefix: --nCTOT
      separate: true
  nCTOB:
    doc: "nCTOB \n"
    type: string
    inputBinding:
      position: 100000
      prefix: --nCTOB
      separate: true
  noCG: boolean
baseCommand: /ngs_share/tools/methyldackel_dev/methyldackel
arguments:
- position: 2
  valueFrom: extract
- position: 10000
  valueFrom: ${ if(inputs.noCG){ return "--noCpG"; }else{ return null; } }
- position: 10001
  valueFrom: ${ if(inputs.noCG){ return "--CHH"; }else{ return null; } }
outputs:
  methcall_bed:
    type: File
    outputBinding:
      glob: '*$(inputs.bedfile_name)*'

