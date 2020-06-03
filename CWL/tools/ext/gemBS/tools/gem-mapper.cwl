#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: gem-mapper
arguments: 
  - valueFrom: $(runtime.cores) # set the number of threads
    prefix: "-t"

#hints:
#  DockerRequirement:
#    dockerPull: yylin/bs_call:1.0

requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    ramMin: 15000

## Mandatory INPUTs
## implement paired-end only.
inputs: 
  index:
    doc: GEM index file
    type: File
    inputBinding:
      position: 1
      prefix: -I
  fasta1:
    doc: paired-end, end-1
    type: File
    inputBinding:
      position: 2
      prefix: --i1
  fasta2:
    doc: paired-end, end-2
    type: File
    inputBinding:
      position: 3
      prefix: --i2
  output:
    type: string
    inputBinding:
      position: 4
      prefix: -o

## Optional I/O args
  gzip-input:
    doc: gzip input
    type: boolean?
    inputBinding:
      prefix: -z
  bzip-input:  
    doc: bzip input
    type: boolean?
    inputBinding:
      prefix: -j
  gzip-output:
    doc: gzip output
    type: boolean?
    inputBinding:
      prefix: --gzip-output
  bzip-output:
    doc: bzip output
    type: boolean?
    inputBinding:
      prefix: --bzip-output
  report-file:
    type: boolean?
    inputBinding:
      prefix: --report-file

## Optional Paired-end Alignment args
  paired-end-alignment:
    type: boolean?
    inputBinding:
      prefix: -p
  min-template-length:
    type: int?
    inputBinding:
      prefix: -l
  max-template-length:
    type: int?
    inputBinding:
      prefix: -L
  discordant-pair-search:
    #doc: 'always'|'if-no-concordant'|'never' (default=if-no-concordant)
    type: string?
    inputBinding:
      prefix: --discordant-pair-search

## Optional Bisulfite Alignment args
  bisulfite-read:
    #doc: 'inferred','1','2','interleaved','non-stranded' (default=inferred)
    type: string?
    inputBinding:
      prefix: --bisulfite-read
  underconversion_sequence:
    doc: underconversion_sequence <sequence name> (default=NC_001416.1)
    type: string?
    inputBinding:
      prefix: --underconversion_sequence
  overconversion_sequence:
    doc: overconversion_sequence <sequence name> (default=NC_001604.1)
    type: string?
    inputBinding:
      prefix: --overconversion_sequence
  control_sequence:
    doc: control_sequence <sequence name> (default=NC_001422.1)
    type: string?
    inputBinding:
      prefix: --control_sequence      
  
## Optional Alignment Score args
  gap-affine-penalties:
    #doc: A,B,O,X (default=1,4,6,1)
    type: int[]?
    inputBinding:
      prefix: -gap-affine-penalties
      itemSeparator: ","      

## Optional Reporting args
  # max-reported-matches:

## Optional Output-format args
  output_format:
    #doc: 'MAP'|'SAM' (default=SAM)
    type: string?
    inputBinding:
      prefix: -F
      valueFrom: |
        ${
          if(inputs.output_format == 'MAP'){
            return "MAP";
          }else{
            return "SAM";
          }
        }

  sam_compact: 
    #doc: 'true'|'false' (default=true)
    type: boolean?
    inputBinding:
      prefix: --sam-compact 
      valueFrom: |
        ${
          if(inputs.sam_compact){
            return "true";
          }else{
            return "false";
          }
        }

  sam-read-group-header:
    #doc: i.e. '@RG\tID:xx\tSM:yy' (default=NULL)
    type: string?
    inputBinding:
      prefix: -r 

## Optional System args
  # threads:

outputs: 
  sam_output:
     type: File
     outputBinding:
        glob: "*.sam"  

