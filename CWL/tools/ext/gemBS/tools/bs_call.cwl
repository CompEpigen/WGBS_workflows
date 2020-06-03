#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: bs_call
#hints:
#  DockerRequirement:
#    dockerPull: yylin/bs_call:1.0

requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    ramMin: 15000

## Mandatory INPUTs
inputs: 
  sample:
    doc: sample name
    type: string
    inputBinding:
      position: 1
      prefix: -n
  reference:
    doc: MultiFASTA/FASTA
    type: File
    inputBinding:
      position: 2
      prefix: -r
  bam:
    doc: input bam file
    type: File
    inputBinding:
      position: 0
    secondaryFiles:
     - .bai  
  output:
    doc: output file name
    type: string
    inputBinding:
      position: 3
      prefix: -o

  #--contig-bed|-C <file> (BED)
  #--contig-sizes|-s <file>
  #--dbsnp|-D <file> (dbSNP processed file)

## Optional I/O args
  output-type: 
    doc: b - compressed BCF, u - uncompressed BCF, z - compressed VCF, v - uncompressed VCF
    type: string?
    inputBinding:
      prefix: -O
  report-file:
    doc: output JSON file name
    type: string?
    inputBinding:
      prefix: --report-file
  all-positions:
    type: boolean?
    inputBinding:
      prefix: -A
  benchmark-mode:
    type: boolean?
    inputBinding:
      prefix: --benchmark-mode
  

## Optional Operation args
  keep-duplicates:
    doc: Don't merge duplicate reads
    type: boolean?
    inputBinding:
      prefix: -d
  ignore-duplicates:
    doc: Ignore duplicate flag from SAM/BAM files
    type: boolean?
    inputBinding:
      prefix: -e
  keep-unmatched:
    doc: Don't discard reads that don't form proper pairs
    type: boolean?
    inputBinding:
      prefix: -k
  right-trim: 
    doc: Bases to trim from right of read pair
    type: boolean?
    inputBinding:
      prefix: -R
  left-trim: 
    doc: Bases to trim from left of read pair
    type: boolean?
    inputBinding:
      prefix: -L
  blank-trim: 
    doc: Don't use trimmed bases for genotype estimation  
    type: boolean?
    inputBinding:
      prefix: -B
  mapq-threshold:
    doc: Set MAPQ threshold for selecting reads (default 20)
    type: int?
    inputBinding:
      prefix: -q
  bq-threshold:
    doc: Set base quality threshold for calling (default 20)
    type: int?
    inputBinding:
      prefix: -Q
  max-template-length:
    doc: Set maximum template length for a pair (default 1000)
    type: int?
    inputBinding:
      prefix: -l

## Optional Model args
  conversion:
    doc: Set under and over conversion rates (default 0.01,0.05)
    type: float[]?
    inputBinding:
      prefix: -c
      itemSeparator: ","
  reference-bias:
    doc: Set bias to reference homozygote (default 2)
    type: float?
    inputBinding:
      prefix: --reference-bias

## Optional Misc args
  threads:
    doc: Set threads
    type: int?
    default: 1
    inputBinding:
      prefix: -t
      valueFrom: ${return(Math.ceil(self/5))}

## OUTPUT PART
outputs: 
  sam_output:
    type: File
    outputBinding:
      glob: $(inputs.output)  
  stderr_log:
    type: stderr

stderr: ${return inputs.output.replace(/\.[^/.]+$/, "") + ".log"}

