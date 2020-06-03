#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: ["java", "-jar", "/usr/local/BisSNP-0.71.jar", "-T", "BisulfiteGenotyper", "-vfn1", "cpg.raw.vcf", "-vfn2", "snp.raw.vcf"] 
requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: bissnp
    # dockerOutputDirectory: /data
  InitialWorkDirRequirement:
    listing: 
      - $(inputs.bam) 
      - $(inputs.bai)     
      - $(inputs.reference) 

inputs: 
  reference:
    type: File
    inputBinding:
      prefix: -R
      position: 1
    secondaryFiles:
      - ^.dict
      - .fai
  dbsnp:
    doc: dbSNP file
    type: File
    inputBinding:
      prefix: -D
      position: 2
  bam:
    type: File
    inputBinding:
      prefix: -I   
      position: 3  
    #secondaryFiles:
    #    - .bai 
  bai:
    doc: substitute "secondaryFiles" for flow control.
    type: File
  call_conf:
    type: int  
    inputBinding:
      prefix: -stand_call_conf 
      position: 4  
  emit_conf:
    type: int  
    inputBinding:
      prefix: -stand_emit_conf 
      position: 5 
  bed:
    doc: the genome regions for calling.
    type: File?  
    inputBinding:
      prefix: -L 
      position: 6

## OUTPUT PART      
outputs: 
  cpg_vcf:
    type: File
    outputBinding:
      glob: cpg.raw.vcf  
  snp_vcf:
    type: File
    outputBinding:
      glob: snp.raw.vcf             
           