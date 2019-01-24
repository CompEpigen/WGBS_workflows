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
    type: File
    # type: array
    # items: File
  fastq2:
    type: File
    # type: array
    # items: File
  reference_fasta: 
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

  pileometh_min_phred: 
    type: int
    default: 0
  pileometh_min_depth: 
    type: int
    default: 1
  pileometh_min_mapq: 
    type: int
    default: 0
  pileometh_ot: 
    type: string
    default: "0,0,0,0"
  pileometh_ob: 
    type: string
    default: "0,0,0,0"
  pileometh_ctot: 
    type: string
    default: "0,0,0,0"
  pileometh_ctob: 
    type: string
    default: "0,0,0,0"
  pileometh_not: 
    type: string
    default: "10,10,10,10"
  pileometh_nob: 
    type: string
    default: "10,10,10,10"
  pileometh_nctot: 
    type: string
    default: "10,10,10,10"
  pileometh_nctob: 
    type: string
    default: "10,10,10,10"
  pilometh_noCG:
    type: boolean
    default: false


# conversion_chr_name: string

steps:
  adaptor_trimming:
    run: "../tools/trimmomatic.cwl"
    # scatterMethod: dotproduct
    # scatter: [fastq1, fastq2] 
    # scatterMethod: 'dotproduct'
    in:
      fastq1: fastq1
      fastq2: fastq2
      adapters_file: trimmomatic_adapters_file
      phred: trimmomatic_phred
      illuminaclip: trimmomatic_illuminaclip
      leading: trimmomatic_leading
      trailing: trimmomatic_trailing
      crop: trimmomatic_crop
      # headcrop: trimmomatic_headcrop
      # tailcrop: trimmomatic_tailcrop
      minlen: trimmomatic_minlen
      avgqual: trimmomatic_avgqual
      threads: max_threads
    out:
    - fastq1_trimmed
    - fastq2_trimmed
    - fastq1_trimmed_unpaired
    - fastq2_trimmed_unpaired
    - trimmomatic_log

  mapping:
    run: "../tools/bwameth.cwl"
    # scatterMethod: dotproduct
    # scatter: [] #!
    in:
      fastq1: adaptor_trimming/fastq1_trimmed
      fastq2: adaptor_trimming/fastq2_trimmed
      reference: reference_fasta
      output_basename: sample_id
      threads: max_threads
      is_non_directional: is_non_directional
    out:
    - sam

  sam2bam:
    doc: samtools view - convert sam to bam
    run: "../tools/samtools_view_sam2bam.cwl"
    in:
      sam: mapping/sam
    out:
      - bam_unsorted

  sort_bam: 
    doc: samtools sort - sorts unsorted bam file by coordinates.
    run: "../tools/samtools_sort.cwl"
    in:
      bam_unsorted: sam2bam/bam_unsorted
    out:
      - bam_sorted

  remove_duplicates:
    doc: picard markdup - emoves duplicates from a single sorted bam file.
    run: "../tools/picard_markdup.cwl"
    in:
      bam_sorted: sort_bam/bam_sorted
    out:
      - bam_duprem
      - picard_markdup_stdout

  index_bam:
    doc: |
      samtools index - indexes sorted bam
    run: "../tools/samtools_index_hack.cwl"
    in:
      bam_sorted: remove_duplicates/bam_duprem
    out:
       - bam_sorted_indexed
  
  methylation_calling:
    # scatter: '#methylation_calling/bam_file'
    run: "../tools/methyldackel-extract.cwl"
    in:
      bam: index_bam/bam_sorted_indexed
      output_basename: sample_id
      reference: reference_fasta
      noCG: pilometh_noCG
      OT: pileometh_ot
      OB: pileometh_ob
      CTOT: pileometh_ctot
      CTOB: pileometh_ctob
      nOT: pileometh_not
      nOB: pileometh_nob
      nCTOT: pileometh_nctot
      nCTOB: pileometh_nctob
      min_phred: pileometh_min_phred
      min_depth: pileometh_min_depth
      min_mapq: pileometh_min_mapq
      noCG: pilometh_noCG
    out:
    - methylcall_bedgraph

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
    # scatter: '#methylation_calling/bam_file'
    run: "../tools/methyldackel-extract.cwl"
    in:
      bam: index_spike_in_bam/bam_sorted_indexed
      output_basename: 
        source: sample_id
        valueFrom: $(self + "_spike_in")
      reference: reference_fasta
      noCG: pilometh_noCG
      OT: pileometh_ot
      OB: pileometh_ob
      CTOT: pileometh_ctot
      CTOB: pileometh_ctob
      nOT: pileometh_not
      nOB: pileometh_nob
      nCTOT: pileometh_nctot
      nCTOB: pileometh_nctob
      min_phred: pileometh_min_phred
      min_depth: pileometh_min_depth
      min_mapq: pileometh_min_mapq
      noCG: pilometh_noCG
    out:
    - methylcall_bedgraph

  conversion_estimation_spike_in:
    run: "../tools/bisulfite_conversion_spike_in.cwl"
    in:
      spike_in_methylcall_bedgraph: methylation_calling_spike_in/methylcall_bedgraph
      output_basename: sample_id
    out:
      - bisulfite_conversion_file

outputs:
  fastq1_trimmed:
    type: File
    outputSource: adaptor_trimming/fastq1_trimmed
  fastq2_trimmed:
    type: File
    outputSource: adaptor_trimming/fastq2_trimmed
  fastq1_trimmed_unpaired:
    type: File
    outputSource: adaptor_trimming/fastq1_trimmed_unpaired
  fastq2_trimmed_unpaired:
    type: File
    outputSource: adaptor_trimming/fastq2_trimmed_unpaired
  trimmomatic_log:
    type: File
    outputSource: adaptor_trimming/trimmomatic_log

  picard_markdup_stdout:
    type: File
    outputSource: remove_duplicates/picard_markdup_stdout

  bam:
    type: File
    outputSource: index_bam/bam_sorted_indexed
  bam_spike_in:
    type: File
    outputSource: index_spike_in_bam/bam_sorted_indexed

  methylcall_bedgraph:
    type: File
    outputSource: methylation_calling/methylcall_bedgraph
  methylcall_bedgraph_spike_in:
    type: File
    outputSource: methylation_calling_spike_in/methylcall_bedgraph

  bisulfite_conversion_file:
    type: File
    outputSource: conversion_estimation_spike_in/bisulfite_conversion_file

