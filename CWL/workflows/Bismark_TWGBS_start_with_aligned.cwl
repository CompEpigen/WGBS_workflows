cwlVersion: v1.0
class: Workflow

requirements:
  ScatterFeatureRequirement: {}
  MultipleInputFeatureRequirement: {}
  SubworkflowFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}

inputs:
  bams_doublex1:
    type: File[]
  bams_doublex2:
    type: File[]
  genome: 
    type: Directory
  threads:
    type: int
    default: 16
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

steps:
  merge_and_sort_doublexes:
    run: "../tools/samtools_merge_and_sort_doublexes.cwl"
    scatter: [bams_doublex1, bams_doublex2]
    scatterMethod: 'dotproduct'
    in:
      bams_doublex1:
        source: bams_doublex1
      bams_doublex2:
        source: bams_doublex2
      name_sort:
        valueFrom: $(true)
      threads: threads
    out:
       - bam_merged

  remove_duplicates:
    run: "../tools/bismark_deduplicate.cwl"
    scatter: [aligned_reads]
    scatterMethod: 'dotproduct'
    in:
      aligned_reads: merge_and_sort_doublexes/bam_merged
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

outputs:
  merged_reads:
    type: File
    outputSource: merge_and_sort/bam_merged
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
  flagstats_post_mapping:
    type: File[]
    outputSource: qc_post_mapping/flagstat_output