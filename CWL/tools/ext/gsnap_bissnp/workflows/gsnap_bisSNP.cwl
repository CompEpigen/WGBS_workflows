#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
inputs:
  #GSNAP
  ref: string
  r1: File
  r2: File

  #BisSNP
  reference: 
    type: File
    secondaryFiles:
      - ^.dict
      - .fai    
  dbsnp: File
  bed: File
  call_conf: int
  emit_conf: int 

  #trim_galore
  adapter1:
    type: string?
  adapter2:
    type: string?
  threads:
    type: int
    default: 16    
  trim_galore_quality:
    type: int
    default: 20
  trim_galore_rrbs:
    type: boolean
    default: false
  trim_galore_clip_r1:
    type: int?
  trim_galore_clip_r2:
    type: int?
  trim_galore_three_prime_clip_r1:
    type: int?
  trim_galore_three_prime_clip_r2:
    type: int?  
  
outputs: 
  sam:
    type: File
    outputSource: mapping/sam
  bam_unsorted:
    type: File
    outputSource: samtools_view_sam2bam/bam_unsorted
  bam_sorted:
    type: File
    outputSource: samtools_sort/bam_sorted
  bam_sorted_indexed:
    type: File
    outputSource: samtools_index/bam_sorted_indexed  
  bam_withRG:
    type: File
    outputSource: addRG/bam_withRG   
  bam_duprem:
    type: File
    outputSource: remove_duplicates/bam_duprem   
  cpg_vcf:
    type: File
    outputSource: calling/cpg_vcf
  trim_log:
    type: File[]
    outputSource: trimming/log 

steps:
  trimming:
    run: trim_galore.cwl
    in:
      read1: r1
      read2: r2
      adapter1: adapter1
      adapter2: adapter2
      quality: trim_galore_quality
      rrbs: trim_galore_rrbs
      clip_r1: trim_galore_clip_r1
      clip_r2: trim_galore_clip_r2
      three_prime_clip_r1: trim_galore_three_prime_clip_r1
      three_prime_clip_r2: trim_galore_three_prime_clip_r2
      threads: threads      
    out: [read1_trimmed, read2_trimmed, log]

  mapping:
    run: gsnap.cwl
    in:
      ref: ref
      r1: trimming/read1_trimmed
      r2: trimming/read2_trimmed
    out: [sam]

  samtools_view_sam2bam:
    run: samtools_view_sam2bam.cwl  
    in:
        sam: mapping/sam
    out: [bam_unsorted]

  addRG:
    run: picard_addRG.cwl
    in:
      bam_withoutRG: samtools_view_sam2bam/bam_unsorted
    out: [bam_withRG]

  samtools_sort:
    run: samtools_sort.cwl  
    in:
        bam_unsorted: addRG/bam_withRG
    out: [bam_sorted]

  remove_duplicates:
    run: picard_markdup.cwl
    in:
        bam_sorted: samtools_sort/bam_sorted
    out: [bam_duprem]

  samtools_index:
    run: samtools_index.cwl  
    in:
        bam_sorted: remove_duplicates/bam_duprem
    out: [bam_sorted_indexed]    

  calling:
    run: bissnp_bisulfite_genotyper.cwl
    in: 
      reference: reference
      dbsnp: dbsnp
      bed: bed
      call_conf: call_conf
      emit_conf: emit_conf   
      bam: remove_duplicates/bam_duprem
      bai: samtools_index/bam_sorted_indexed #workflow control
    out: [cpg_vcf, snp_vcf]