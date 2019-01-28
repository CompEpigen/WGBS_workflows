cwlVersion: v1.0
class: CommandLineTool
requirements:
  ShellCommandRequirement: {}
inputs:
  file_dir: string
  bam_file:
    doc: "the input bam file\n"
    type: string
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
baseCommand: cd
arguments:
- position: -2
  valueFrom: $(inputs.file_dir)
- position: -1
  valueFrom: ;
- position: 1
  valueFrom: /ngs_share/tools/methyldackel_dev/methyldackel
- position: 2
  valueFrom: extract
- position: 10000
  valueFrom: ${ if(inputs.noCG){ return "--noCpG"; }else{ return null; } }
- position: 10001
  valueFrom: ${ if(inputs.noCG){ return "--CHH"; }else{ return null; } }
- position: 1000000
  valueFrom: ;
- position: 1000001
  valueFrom: ls
- shellQuote: false
  position: 1000002
  valueFrom: '*bedGraph'
- position: 1000003
  valueFrom: '|'
- position: 1000004
  valueFrom: grep
- position: 1000005
  valueFrom: $(inputs.bedfile_name)
stdout: meth_calls
outputs:
  methcall_bed:
    type: string
    outputBinding:
      glob: meth_calls
      loadContents: true
      outputEval: $(self[0].contents.split("\n")[0])

