cwlVersion: v1.0
class: Workflow

requirements:
  ScatterFeatureRequirement: {}
  MultipleInputFeatureRequirement: {}
  SubworkflowFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}

inputs:
  reads:
    doc: |
      Even though seqeuncing has been done in PE mode, trimming and alignment are performed 
      in SE-mode. Please specify a list of all fastq files belonging to one sample 
      (first and second read, multiple lanes).
    type:
      type: array
      items: File
  genome: 
    type: Directory
  threads:
    type: int
    default: 16
  trim_galore_quality:
    type: int
    default: 20
  trim_galore_clip_r1:
    type: int?
    default: 6
  trim_galore_three_prime_clip_r1:
    type: int?
    default: 6
  bismark_ignore:
    type: int?
  bismark_ignore_3prime:
    type: int?

steps:
  qc_pretrim:
    scatter: [reads]
    scatterMethod: 'dotproduct'
    run: "../tools/fastqc_scPBAT.cwl"
    in:
      reads: reads
    out:
      - fastqc_zip
      - fastqc_html

  trim:
    scatter: [reads]
    scatterMethod: 'dotproduct'
    run: "../tools/trim_galore_scPBAT.cwl"
    in:
      reads: reads
      quality: trim_galore_quality
      clip_r1: trim_galore_clip_r1
      three_prime_clip_r1: trim_galore_three_prime_clip_r1
      threads: threads
    out:
      - log
      - reads_trimmed
  
  qc_posttrim:
    scatter: [reads]
    scatterMethod: 'dotproduct'
    run: "../tools/fastqc_scPBAT.cwl"
    in:
      reads: trim/reads_trimmed
    out:
      - fastqc_zip
      - fastqc_html
      
  align:
    scatter: [reads]
    scatterMethod: 'dotproduct'
    run: "../tools/bismark_align_scPBAT_se.cwl"
    in:
      reads: trim/reads_trimmed
      threads: threads
      genome: genome
    out:
      - aligned_reads
      - log

  merge_and_sort:
    run: "../tools/samtools_merge_and_sort.cwl"
    in:
      bams:
        source: align/aligned_reads
      name_sort:
        valueFrom: $(true)
      threads: threads
    out:
       - bam_merged

  remove_duplicates:
    run: "../tools/bismark_deduplicate.cwl"
    in:
      aligned_reads: merge_and_sort/bam_merged
      paired_end:
        valueFrom: ${return(false)}
    out:
      - dedup_reads
      - log

      
  extract_methylation:
    run: "../tools/bismark_methylation_extractor.cwl"
    in:
      aligned_reads: remove_duplicates/dedup_reads
      no_overlap: 
        valueFrom: ${return(false)}
      ignore: bismark_ignore
      ignore_3prime: bismark_ignore_3prime
      threads: threads
      genome: genome
      paired_end:
        valueFrom: ${return(false)}
    out:
      - methylation_calls_bedgraph
      - methylation_calls_bismark
      - mbias_report
      - splitting_report
      - genome_wide_methylation_report
      - context_specific_methylation_reports
  
  bismark_report:
    scatter: [alignment_report]
    run: "../tools/bismark_report.cwl"
    in:
      alignment_report: align/log
      dedup_report: remove_duplicates/log
      splitting_report: extract_methylation/splitting_report
      mbias_report: extract_methylation/mbias_report
    out:
      - report

  qc_post_mapping:
    doc: |
      samtools flagstat
    run: "../tools/samtools_flagstat.cwl"
    in:
      bam: remove_duplicates/dedup_reads
    out:
       - flagstat_output

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
  reads_trimmed:
    type: File[]
    outputSource: trim/reads_trimmed
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
    type: File
    outputSource: remove_duplicates/dedup_reads
  dedup_log:
    type: File
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
    type: File
    outputSource: qc_post_mapping/flagstat_output