cwlVersion: v1.0
class: Workflow

requirements:
  ScatterFeatureRequirement: {}
  MultipleInputFeatureRequirement: {}
  SubworkflowFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}

inputs:
  read1:
    type:
      type: array
      items: File
  read2:
    type:
      type: array
      items: File
  adapter1:
    type: string?
  adapter2:
    type: string?
  genome: 
    type: Directory
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
  bismark_pbat:
    type: boolean
    default: false
  bismark_ignore:
    type: int
    default: 9
  bismark_ignore_r2:
    type: int
    default: 12
  bismark_ignore_3prime:
    type: int
    default: 9
  bismark_ignore_3prime_r2:
    type: int
    default: 2
  bismark_no_overlap:
    type: boolean
    default: true
  bismark_local:
    type: boolean
  non_directional:
    type: boolean
  dovetail:
    type: boolean


steps:
  qc_pretrim:
    scatter: [read1, read2]
    scatterMethod: 'dotproduct'
    run: "../tools/fastqc.cwl"
    in:
      read1: read1
      read2: read2
    out:
      - fastqc_zip
      - fastqc_html

  trim:
    scatter: [read1, read2]
    scatterMethod: 'dotproduct'
    run: "../tools/trim_galore.cwl"
    in:
      read1: read1
      read2: read2
      adapter1: adapter1
      adapter2: adapter2
      quality: trim_galore_quality
      rrbs: trim_galore_rrbs
      clip_r1: trim_galore_clip_r1
      clip_r2: trim_galore_clip_r2
      three_prime_clip_r1: trim_galore_three_prime_clip_r1
      three_prime_clip_r2: trim_galore_three_prime_clip_r2
      threads: threads
    out:
      - log
      - read1_trimmed
      - read2_trimmed
  
  qc_posttrim:
    scatter: [read1, read2]
    scatterMethod: 'dotproduct'
    run: "../tools/fastqc.cwl"
    in:
      read1: trim/read1_trimmed
      read2: trim/read2_trimmed
    out:
      - fastqc_zip
      - fastqc_html
      
  align:
    scatter: [read1, read2]
    scatterMethod: 'dotproduct'
    run: "../tools/bismark_align.cwl"
    in:
      read1: trim/read1_trimmed
      read2: trim/read2_trimmed
      pbat: bismark_pbat
      bismark_local: bismark_local
      non_directional: non_directional
      dovetail: dovetail
      threads: threads
      genome: genome
    out:
      - aligned_reads
      - log

  remove_duplicates:
    run: "../tools/bismark_deduplicate.cwl"
    scatter: [aligned_reads]
    scatterMethod: 'dotproduct'
    in:
      aligned_reads: align/aligned_reads
    out:
      - dedup_reads
      - log

  qc_post_mapping:
    doc: |
      samtools flagstat
    run: "../tools/samtools_flagstat.cwl"
    scatter: [bam]
    scatterMethod: 'dotproduct'
    in:
      bam: remove_duplicates/dedup_reads
    out:
       - flagstat_output

  merge_and_sort:
    run: "../tools/samtools_merge_and_sort.cwl"
    in:
      bams:
        source: remove_duplicates/dedup_reads
      name_sort:
        valueFrom: $(true)
      threads: threads
    out:
       - bam_merged

  extract_methylation:
    run: "../tools/bismark_methylation_extractor.cwl"
    in:
      aligned_reads: merge_and_sort/bam_merged
      no_overlap: bismark_no_overlap
      ignore: bismark_ignore
      ignore_r2: bismark_ignore_r2
      ignore_3prime: bismark_ignore_3prime
      ignore_3prime_r2: bismark_ignore_3prime_r2
      threads: threads
      genome: genome
    out:
      - methylation_calls_bedgraph
      - methylation_calls_bismark
      - mbias_report
      - splitting_report
      - genome_wide_methylation_report
      - context_specific_methylation_reports
  
  bismark_report:
    scatter: [alignment_report, dedup_report]
    scatterMethod: 'dotproduct'
    run: "../tools/bismark_report.cwl"
    in:
      alignment_report: align/log
      dedup_report: remove_duplicates/log
      splitting_report: extract_methylation/splitting_report
      mbias_report: extract_methylation/mbias_report
    out:
      - report

outputs:
  qc_pretrim_fastqc_zip:
    type:
      type: array
      items:
        type: array
        items: File
    outputSource: qc_pretrim/fastqc_zip
  qc_pretrim_fastqc_html:
    type:
      type: array
      items:
        type: array
        items: File
    outputSource: qc_pretrim/fastqc_html
  trim_log:
    type:
      type: array
      items:
        type: array
        items: File
    outputSource: trim/log
  read1_trimmed:
    type: File[]
    outputSource: trim/read1_trimmed
  read2_trimmed:
    type: File[]
    outputSource: trim/read2_trimmed
  qc_posttrim_fastqc_zip:
    type:
      type: array
      items:
        type: array
        items: File
    outputSource: qc_posttrim/fastqc_zip
  qc_posttrim_fastqc_html:
    type:
      type: array
      items:
        type: array
        items: File
    outputSource: qc_posttrim/fastqc_html
  align_log:
    type: File[]
    outputSource: align/log
  dedup_reads:
    type: File[]
    outputSource: remove_duplicates/dedup_reads
  dedup_log:
    type: File[]
    outputSource: remove_duplicates/log
  methylation_calls_bedgraph:
    type: File
    outputSource: extract_methylation/methylation_calls_bedgraph
  methylation_calls_bismark:
    type: File
    outputSource: extract_methylation/methylation_calls_bismark
  mbias_report:
    type: File
    outputSource: extract_methylation/mbias_report
  splitting_report:
    type: File
    outputSource: extract_methylation/splitting_report
  genome_wide_methylation_report:
    type: File
    outputSource: extract_methylation/genome_wide_methylation_report
  context_specific_methylation_reports:
    type: File[]
    outputSource: extract_methylation/context_specific_methylation_reports
  bismark_report_html:
    type: File[]
    outputSource: bismark_report/report
  flagstats_post_mapping:
    type: File[]
    outputSource: qc_post_mapping/flagstat_output