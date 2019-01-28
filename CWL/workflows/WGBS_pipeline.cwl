cwlVersion: v1.0
class: Workflow

requirements:
  ScatterFeatureRequirement: {}
  MultipleInputFeatureRequirement: {}
  SubworkflowFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}

inputs:
  sample_id: string
  fastq1:
    type:
      type: array
      items: File
  fastq2:
    type:
      type: array
      items: File
  reference: 
    type: File
    secondaryFiles:
      - .fai
      - .bwameth.c2t
      - .bwameth.c2t.amb
      - .bwameth.c2t.ann
      - .bwameth.c2t.bwt
      - .bwameth.c2t.pac
      - .bwameth.c2t.sa
  trimmomatic_adapters_file: File
  spikein_chr_name:
    doc: chromosome name in the reference genome which corresponds to the unmethylated spikein DNA
    type: string
    default: "Lamda"

  # # optional parameter that likely need revision:
  is_non_directional:
    type: boolean
    default: false
  spike_in_chr_name:
    doc: chromosome name in the reference genome which corresponds to the unmethylated spikein DNA
    type: string
    default: "Lamda"
  max_threads:
    type: int
    default: 10

  # optional parameters:
  trimmomatic_phred: 
    type: string
    default: "64"
  trimmomatic_leading: 
    type: int
    default: 0
  trimmomatic_trailing: 
    type: int
    default: 0
  trimmomatic_crop: 
    type: int
    default: 1000
  trimmomatic_headcrop: 
    type: int
    default: 10
  trimmomatic_tailcrop: 
    type: int
    default: 10
  trimmomatic_minlen:
    type: int
    default: 0
  trimmomatic_avgqual:
    type: int
    default: 1
  trimmomatic_illuminaclip: 
    type: string
    default: "2:30:10:8:true"

  methyldackel_min_phred: 
    type: int
    default: 0
  methyldackel_min_depth: 
    type: int
    default: 1
  methyldackel_min_mapq: 
    type: int
    default: 0
  methyldackel_ot: 
    type: string
    default: "0,0,0,0"
  methyldackel_ob: 
    type: string
    default: "0,0,0,0"
  methyldackel_ctot: 
    type: string
    default: "0,0,0,0"
  methyldackel_ctob: 
    type: string
    default: "0,0,0,0"
  methyldackel_not: 
    type: string
    default: "10,10,10,10"
  methyldackel_nob: 
    type: string
    default: "10,10,10,10"
  methyldackel_nctot: 
    type: string
    default: "10,10,10,10"
  methyldackel_nctob: 
    type: string
    default: "10,10,10,10"
  methyldackel_noCG:
    type: boolean
    default: false

steps:
  trim_map_duprem:
    scatter: [fastq1, fastq2]
    scatterMethod: 'dotproduct'
    run: "../workflow_modules/trim_map_duprem.cwl"
    in:
      fastq1: fastq1
      fastq2: fastq2
      trimmomatic_adapters_file: trimmomatic_adapters_file
      trimmomatic_phred: trimmomatic_phred
      trimmomatic_illuminaclip: trimmomatic_illuminaclip
      trimmomatic_leading: trimmomatic_leading
      trimmomatic_trailing: trimmomatic_trailing
      trimmomatic_crop: trimmomatic_crop
      # trimmomatic_headcrop: trimmomatic_headcrop
      # trimmomatic_tailcrop: trimmomatic_tailcrop
      trimmomatic_minlen: trimmomatic_minlen
      trimmomatic_avgqual: trimmomatic_avgqual
      max_threads: max_threads
      reference: reference
      is_non_directional: is_non_directional
    out:
      - trimmomatic_log
      - picard_markdup_stdout
      - bam
      - pre_trim_fastqc_zip
      - pre_trim_fastqc_html
      - post_trim_fastqc_zip
      - post_trim_fastqc_html

  lane_replicate_merging:
    doc: samtools merge - merging bam files of lane replicates
    run: "../tools/samtools_merge.cwl"
    in:
      bams:
        source: trim_map_duprem/bam
      output_name:
        source: sample_id
        valueFrom: $(self + ".bam")
    out:
       - bam_merged

  sorting_merged_bam:
    doc: samtools sort - sorting of merged bam
    run: "../tools/samtools_sort.cwl"
    in:
      bam_unsorted:
        source: lane_replicate_merging/bam_merged
    out:
       - bam_sorted

  index_bam:
    doc: |
      samtools index - indexes sorted bam
    run: "../tools/samtools_index_hack.cwl"
    in:
      bam_sorted: sorting_merged_bam/bam_sorted
    out:
       - bam_sorted_indexed
  
  qc_post_mapping:
    doc: fastqc - quality control for reads after mapping and duplication removal
    run: "../tools/fastqc.cwl"
    in:
      bam:
        source: index_bam/bam_sorted_indexed
    out:
      - fastqc_zip
      - fastqc_html
  
  mbias_calculation:
    run: "../tools/methyldackel_mbias.cwl"
    in:
      bam: index_bam/bam_sorted_indexed
      output_basename: sample_id
      reference: reference
      nOT: methyldackel_not
      nOB: methyldackel_nob
      nCTOT: methyldackel_nctot
      nCTOB: methyldackel_nctob
      noCG: methyldackel_noCG
      threads: max_threads
    out:
    - mbias_output
  
  methylation_calling:
    run: "../tools/methyldackel_extract.cwl"
    in:
      bam: index_bam/bam_sorted_indexed
      output_basename: sample_id
      reference: reference
      OT: methyldackel_ot
      OB: methyldackel_ob
      CTOT: methyldackel_ctot
      CTOB: methyldackel_ctob
      nOT: methyldackel_not
      nOB: methyldackel_nob
      nCTOT: methyldackel_nctot
      nCTOB: methyldackel_nctob
      min_phred: methyldackel_min_phred
      min_depth: methyldackel_min_depth
      min_mapq: methyldackel_min_mapq
      noCG: methyldackel_noCG
      threads: max_threads
    out:
    - mcall_bedgraph

  get_spike_in_reads:
    doc: extracts reads mapping to the unmethylated spike in
    run: "../tools/samtools_view_extract_spike_in.cwl"
    in:
      bam: index_bam/bam_sorted_indexed
      region: spike_in_chr_name
    out:
      - bam_spike_in

  index_spike_in_bam:
    doc: |
      samtools index - indexes sorted bam
    run: "../tools/samtools_index_hack.cwl"
    in:
      bam_sorted: get_spike_in_reads/bam_spike_in
    out:
       - bam_sorted_indexed

  methylation_calling_spike_in:
    run: "../tools/methyldackel_extract.cwl"
    in:
      bam: index_spike_in_bam/bam_sorted_indexed
      output_basename: 
        source: sample_id
        valueFrom: $(self + "_spike_in")
      reference: reference
      OT: methyldackel_ot
      OB: methyldackel_ob
      CTOT: methyldackel_ctot
      CTOB: methyldackel_ctob
      nOT: methyldackel_not
      nOB: methyldackel_nob
      nCTOT: methyldackel_nctot
      nCTOB: methyldackel_nctob
      min_phred: methyldackel_min_phred
      min_depth: methyldackel_min_depth
      min_mapq: methyldackel_min_mapq
      noCG: methyldackel_noCG
      threads: max_threads
    out:
    - mcall_bedgraph

  conversion_estimation_spike_in:
    run: "../tools/bisulfite_conversion_spike_in.cwl"
    in:
      spike_in_mcall_bedgraph: methylation_calling_spike_in/mcall_bedgraph
      output_basename: sample_id
    out:
      - bisulfite_conversion_file

  create_summary_qc_report:
    doc: |
      multiqc summarizes the qc results from fastqc 
      and other tools
    run: "../tools/multiqc_hack.cwl"
    in:
      qc_files_array_of_array:
        source:
          - trim_map_duprem/pre_trim_fastqc_zip
          - trim_map_duprem/pre_trim_fastqc_html
          - trim_map_duprem/post_trim_fastqc_html
          - trim_map_duprem/post_trim_fastqc_zip
        linkMerge: merge_flattened
      qc_files_array:
        source:
          - trim_map_duprem/picard_markdup_stdout
          - trim_map_duprem/trimmomatic_log
          - qc_post_mapping/fastqc_zip
          - qc_post_mapping/fastqc_html
        linkMerge: merge_flattened
      report_name:
        source: sample_id
    out:
      - multiqc_zip
      - multiqc_html

outputs:
  trimmomatic_log:
    type:
      type: array
      items: File
    outputSource: trim_map_duprem/trimmomatic_log

  picard_markdup_stdout:
    type:
      type: array
      items: File
    outputSource: trim_map_duprem/picard_markdup_stdout

  bam:
    type: File
    outputSource: index_bam/bam_sorted_indexed
  bam_spike_in:
    type: File
    outputSource: index_spike_in_bam/bam_sorted_indexed

  mbias_output:
    type:
      type: array
      items: File
    outputSource: mbias_calculation/mbias_output

  mcall_bedgraph:
    type: File
    outputSource: methylation_calling/mcall_bedgraph
  mcall_bedgraph_spike_in:
    type: File
    outputSource: methylation_calling_spike_in/mcall_bedgraph

  bisulfite_conversion_file:
    type: File
    outputSource: conversion_estimation_spike_in/bisulfite_conversion_file

  
  pre_trim_fastqc_zip:
    type:
      type: array
      items: 
        type: array
        items: File
    outputSource: trim_map_duprem/pre_trim_fastqc_zip
  pre_trim_fastqc_html:
    type:
      type: array
      items: 
        type: array
        items: File
    outputSource: trim_map_duprem/pre_trim_fastqc_html
  post_trim_fastqc_zip:
    type:
      type: array
      items: 
        type: array
        items: File
    outputSource: trim_map_duprem/post_trim_fastqc_zip
  post_trim_fastqc_html:
    type:
      type: array
      items: 
        type: array
        items: File
    outputSource: trim_map_duprem/post_trim_fastqc_html
  post_mapping_fastqc_zip:
    type:
      type: array
      items: File
    outputSource: qc_post_mapping/fastqc_zip
  post_mapping_fastqc_html:
    type:
      type: array
      items: File
    outputSource: qc_post_mapping/fastqc_html

  multiqc_zip:
    type: File
    outputSource: create_summary_qc_report/multiqc_zip
  multiqc_html:
    type: File
    outputSource: create_summary_qc_report/multiqc_html