cwlVersion: v1.0
inputs:
  isbam:
    default: false
    doc: "output in BAM format\n"
    inputBinding:
      position: 2
      prefix: -b
    type: boolean
  readswithbits:
    doc: "only include reads with all bits set in INT set in FLAG [0]\n"
    inputBinding:
      position: 1
      prefix: -f
    type: int?
  cigar:
    doc: "only include reads with number of CIGAR operations\nconsuming query sequence\
      \ >= INT [0]\n"
    inputBinding:
      position: 1
      prefix: -m
    type: int?
  readswithoutbits:
    doc: "only include reads with none of the bits set in INT set in FLAG [0]\n"
    inputBinding:
      position: 1
      prefix: -F
    type: int?
  fastcompression:
    default: false
    doc: "use fast BAM compression (implies -b)\n"
    inputBinding:
      position: 1
      prefix: '-1'
    type: boolean
  iscram:
    default: false
    doc: "output in CRAM format\n"
    inputBinding:
      position: 2
      prefix: -C
    type: boolean
  collapsecigar:
    default: false
    doc: "collapse the backward CIGAR operation\n"
    inputBinding:
      position: 1
      prefix: -B
    type: boolean
  readsingroup:
    doc: "only include reads in read group STR [null]\n"
    inputBinding:
      position: 1
      prefix: -r
    type: string?
  randomseed:
    doc: "integer part sets seed of random number generator [0];\nrest sets fraction\
      \ of templates to subsample [no subsampling]\n"
    inputBinding:
      position: 1
      prefix: -s
    type: float?
  samheader:
    default: false
    doc: "include header in SAM output\n"
    inputBinding:
      position: 1
      prefix: -h
    type: boolean
  count:
    default: false
    doc: "print only the count of matching records\n"
    inputBinding:
      position: 1
      prefix: -c
    type: boolean
  threads:
    doc: "number of BAM compression threads [0]\n"
    inputBinding:
      position: 1
      prefix: -@
    type: int?
  referencefasta:
    doc: "reference sequence FASTA FILE [null]\n"
    inputBinding:
      position: 1
      prefix: -T
    type: File?
  region:
    doc: "[region ...]\n"
    inputBinding:
      position: 5
    type: string?
  bedoverlap:
    doc: "only include reads overlapping this BED FILE [null]\n"
    inputBinding:
      position: 1
      prefix: -L
    type: File?
  readsingroupfile:
    doc: "only include reads with read group listed in FILE [null]\n"
    inputBinding:
      position: 1
      prefix: -R
    type: File?
  uncompressed:
    default: false
    doc: "uncompressed BAM output (implies -b)\n"
    inputBinding:
      position: 1
      prefix: -u
    type: boolean
  readtagtostrip:
    doc: "read tag to strip (repeatable) [null]\n"
    inputBinding:
      position: 1
    type: string[]?
  input:
    doc: "Input bam file.\n"
    inputBinding:
      position: 4
    type: string
  output_name:
    inputBinding:
      position: 2
      prefix: -o
    type: string
  readsquality:
    doc: "only include reads with mapping quality >= INT [0]\n"
    inputBinding:
      position: 1
      prefix: -q
    type: int?
  readsinlibrary:
    doc: "only include reads in library STR [null]\n"
    inputBinding:
      position: 1
      prefix: -l
    type: string?
requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement
outputs:
  output:
    outputBinding:
      outputEval: $(inputs.output_name)
    type: string
baseCommand:
- samtools
- view
class: CommandLineTool

