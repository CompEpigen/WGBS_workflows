cwlVersion: v1.0
class: Workflow

requirements:
  ScatterFeatureRequirement: {}
  MultipleInputFeatureRequirement: {}
  SubworkflowFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}

inputs:
  bams:
    type:
      type: array
      items: File
  genome: 
    type: Directory
  threads:
    type: int
    default: 16
  bismark_pbat:
    type: boolean
    default: false
  bismark_ignore:
    type: int?
  bismark_ignore_r2:
    type: int?
  bismark_ignore_3prime:
    type: int?
  bismark_ignore_3prime_r2:
    type: int?
  bismark_no_overlap:
    type: boolean
    default: true

steps:
  merge_and_sort:
    run: "../tools/samtools_merge_and_sort.cwl"
    in:
      bams:
        source: bams
      name_sort:
        valueFrom: $(true)
      threads: threads
    out:
       - bam_merged

  remove_duplicates:
    run: "../tools/bismark_deduplicate.cwl"
    in:
      aligned_reads: merge_and_sort/bam_merged
    out:
      - dedup_reads
      - log

      
  extract_methylation:
    run: "../tools/bismark_methylation_extractor.cwl"
    in:
      aligned_reads: remove_duplicates/dedup_reads
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
  

  qc_post_mapping:
    doc: |
      samtools flagstat
    run: "../tools/samtools_flagstat.cwl"
    in:
      bam: remove_duplicates/dedup_reads
    out:
       - flagstat_output

outputs:
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
  flagstats_post_mapping:
    type: File
    outputSource: qc_post_mapping/flagstat_output