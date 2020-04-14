class: Workflow
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: fastq1
    type: 'File[]'
    'sbg:x': -82.05426025390625
    'sbg:y': -708.0413818359375
  - id: fastq2
    type: 'File[]'
    'sbg:x': -86.61282348632812
    'sbg:y': -541.4498291015625
  - id: is_non_directional
    type: boolean
    default: false
    'sbg:x': -86.61282348632812
    'sbg:y': -365.7083740234375
  - id: max_threads
    type: int
    default: 10
    'sbg:x': -78.44940948486328
    'sbg:y': 89.95835876464844
  - id: methyldackel_ctob
    type: string?
    'sbg:x': 0
    'sbg:y': 2280.761474609375
  - id: methyldackel_ctot
    type: string?
    'sbg:x': -18.234277725219727
    'sbg:y': 2155.620849609375
  - id: methyldackel_min_depth
    type: int
    default: 1
    'sbg:x': 7.168551405811741e-7
    'sbg:y': 2016.8048095703125
  - id: methyldackel_min_mapq
    type: int
    default: 0
    'sbg:x': -5.7605383574355074e-8
    'sbg:y': 1864.312744140625
  - id: methyldackel_min_phred
    type: int
    default: 0
    'sbg:x': 9.11713981628418
    'sbg:y': 1693.5865478515625
  - id: methyldackel_nctob
    type: string?
    default: '2,2,2,2'
    'sbg:x': 0.0000010176772775594145
    'sbg:y': 1554.770263671875
  - id: methyldackel_nctot
    type: string?
    default: '2,2,2,2'
    'sbg:x': -4.558569431304932
    'sbg:y': 1415.9539794921875
  - id: methyldackel_noCG
    type: boolean
    default: false
    'sbg:x': 0
    'sbg:y': 1124.6458740234375
  - id: methyldackel_nob
    type: string?
    default: '2,2,2,2'
    'sbg:x': 4.558570384979248
    'sbg:y': 1254.344970703125
  - id: methyldackel_not
    type: string?
    default: '2,2,2,2'
    'sbg:x': -9.11713981628418
    'sbg:y': 985.8296508789062
  - id: methyldackel_ob
    type: string?
    'sbg:x': 0
    'sbg:y': 847.013427734375
  - id: methyldackel_ot
    type: string?
    'sbg:x': -9.681075096130371
    'sbg:y': 702.4104614257812
  - id: reference
    type: File
    secondaryFiles:
      - .fai
      - .bwameth.c2t
      - .bwameth.c2t.amb
      - .bwameth.c2t.ann
      - .bwameth.c2t.bwt
      - .bwameth.c2t.pac
      - .bwameth.c2t.sa
    'sbg:x': -94.5317153930664
    'sbg:y': -168.46800231933594
  - id: sample_id
    type: string
    'sbg:x': -73.00712585449219
    'sbg:y': 248.1683807373047
  - id: spikein_chr_name
    type: string
    doc: >-
      chromosome name in the reference genome which corresponds to the
      unmethylated spikein DNA
    default: Lamda
    'sbg:x': -136.75709533691406
    'sbg:y': -1331.1024169921875
outputs:
  - id: bam
    outputSource:
      - index_bam/bam_sorted_indexed
    type: File
    'sbg:x': 4508
    'sbg:y': 658
  - id: bam_spike_in
    outputSource:
      - index_spike_in_bam/bam_sorted_indexed
    type: File
    'sbg:x': 4755.5654296875
    'sbg:y': 1112.589111328125
  - id: bisulfite_conversion_file
    outputSource:
      - conversion_estimation_spike_in/bisulfite_conversion_file
    type: File
    'sbg:x': 4739.7578125
    'sbg:y': 936.0545654296875
  - id: mbias_output
    outputSource:
      - mbias_calculation/mbias_output
    type: 'File[]'
    'sbg:x': 4634.173828125
    'sbg:y': 1724.246826171875
  - id: mbias_output_wo_trimming
    outputSource:
      - mbias_calculation_wo_trimming/mbias_output
    type: 'File[]'
    'sbg:x': 4529.32666015625
    'sbg:y': -46.497440338134766
  - id: mcall_bedgraph
    outputSource:
      - methylation_calling/mcall_bedgraph
    type: File
    'sbg:x': 4697.99365234375
    'sbg:y': 2162.25732421875
  - id: mcall_bedgraph_spike_in
    outputSource:
      - methylation_calling_spike_in/mcall_bedgraph
    type: File
    'sbg:x': 4748.85400390625
    'sbg:y': 1281.314208984375
  - id: multiqc_html
    outputSource:
      - create_summary_qc_report/multiqc_html
    type: File
    'sbg:x': 4540.59521484375
    'sbg:y': 492.2149353027344
  - id: multiqc_zip
    outputSource:
      - create_summary_qc_report/multiqc_zip
    type: File
    'sbg:x': 4572.505859375
    'sbg:y': 285.0201416015625
  - id: picard_markdup_log
    outputSource:
      - remove_duplicates/picard_markdup_log
    type: File
    'sbg:x': 2336.66064453125
    'sbg:y': -639.9028930664062
  - id: post_mapping_fastqc_html
    outputSource:
      - qc_post_mapping/fastqc_html
    type: 'File[]'
    'sbg:x': 4597.70458984375
    'sbg:y': -652.3868408203125
  - id: post_mapping_fastqc_zip
    outputSource:
      - qc_post_mapping/fastqc_zip
    type: 'File[]'
    'sbg:x': 4538.44384765625
    'sbg:y': -750.1759643554688
  - id: post_mapping_flagstats
    outputSource:
      - flagstats_post_mapping/flagstat_output
    type: File
    'sbg:x': 4442.7138671875
    'sbg:y': -392.1081237792969
steps:
  - id: conversion_estimation_spike_in
    in:
      - id: output_basename
        source: sample_id
      - id: spike_in_mcall_bedgraph
        source: methylation_calling_spike_in/mcall_bedgraph
    out:
      - id: bisulfite_conversion_file
    run: ../tools/bisulfite_conversion_spike_in.cwl
    'sbg:x': 4128.88916015625
    'sbg:y': 757.0206298828125
  - id: create_summary_qc_report
    in:
      - id: qc_files_array
        linkMerge: merge_flattened
        source:
          - remove_duplicates/picard_markdup_log
          - qc_post_mapping/fastqc_zip
          - qc_post_mapping/fastqc_html
          - flagstats_post_mapping/flagstat_output
      - id: report_name
        source: sample_id
    out:
      - id: multiqc_html
      - id: multiqc_zip
    run: ../tools/multiqc_hack.cwl
    doc: |
      multiqc summarizes the qc results from fastqc 
      and other tools
    'sbg:x': 4182.88134765625
    'sbg:y': 385.27703857421875
  - id: flagstats_post_mapping
    in:
      - id: bam
        source: index_bam/bam_sorted_indexed
    out:
      - id: flagstat_output
    run: ../tools/samtools_flagstat.cwl
    doc: |
      samtools flagstat
    'sbg:x': 3396.394775390625
    'sbg:y': -404.9775695800781
  - id: get_spike_in_reads
    in:
      - id: bam
        source: index_bam/bam_sorted_indexed
      - id: region
        source: spikein_chr_name
    out:
      - id: bam_spike_in
    run: ../tools/samtools_view_extract_spike_in.cwl
    doc: extracts reads mapping to the unmethylated spike in
    'sbg:x': 2843
    'sbg:y': 886
  - id: index_bam
    in:
      - id: bam_sorted
        source: remove_duplicates/bam_duprem
    out:
      - id: bam_sorted_indexed
    run: ../tools/samtools_index_hack.cwl
    doc: |
      samtools index - indexes sorted bam
    'sbg:x': 1901
    'sbg:y': 817
  - id: index_spike_in_bam
    in:
      - id: bam_sorted
        source: get_spike_in_reads/bam_spike_in
    out:
      - id: bam_sorted_indexed
    run: ../tools/samtools_index_hack.cwl
    doc: |
      samtools index - indexes sorted bam
    'sbg:x': 3262.044189453125
    'sbg:y': 905.8668212890625
  - id: lane_replicate_merging
    in:
      - id: output_name
        source: sample_id
        valueFrom: $(self + ".bam")
      - id: bams
        source:
          - map/sam
    out:
      - id: bam_merged
    run: ../tools/samtools_merge.cwl
    doc: samtools merge - merging bam files of lane replicates
    'sbg:x': 1070.7197265625
    'sbg:y': -374.6963806152344
  - id: mbias_calculation
    in:
      - id: bam
        source: index_bam/bam_sorted_indexed
      - id: nCTOB
        source: methyldackel_nctob
      - id: nCTOT
        source: methyldackel_nctot
      - id: nOB
        source: methyldackel_nob
      - id: nOT
        source: methyldackel_not
      - id: noCG
        source: methyldackel_noCG
      - id: output_basename
        source: sample_id
      - id: reference
        source: reference
      - id: threads
        source: max_threads
    out:
      - id: mbias_output
    run: ../tools/methyldackel_mbias.cwl
    'sbg:x': 2867.6005859375
    'sbg:y': 1690.6920166015625
  - id: mbias_calculation_wo_trimming
    in:
      - id: bam
        source: index_bam/bam_sorted_indexed
      - id: output_basename
        source: sample_id
        valueFrom: $(self + "_wo_trimming")
      - id: reference
        source: reference
      - id: threads
        source: max_threads
    out:
      - id: mbias_output
    run: ../tools/methyldackel_mbias.cwl
    'sbg:x': 3018.033447265625
    'sbg:y': 48.106319427490234
  - id: methylation_calling
    in:
      - id: CTOB
        source: methyldackel_ctob
      - id: CTOT
        source: methyldackel_ctot
      - id: OB
        source: methyldackel_ob
      - id: OT
        source: methyldackel_ot
      - id: bam
        source: index_bam/bam_sorted_indexed
      - id: min_depth
        source: methyldackel_min_depth
      - id: min_mapq
        source: methyldackel_min_mapq
      - id: min_phred
        source: methyldackel_min_phred
      - id: nCTOB
        source: methyldackel_nctob
      - id: nCTOT
        source: methyldackel_nctot
      - id: nOB
        source: methyldackel_nob
      - id: nOT
        source: methyldackel_not
      - id: noCG
        source: methyldackel_noCG
      - id: output_basename
        source: sample_id
      - id: reference
        source: reference
      - id: threads
        source: max_threads
    out:
      - id: mcall_bedgraph
    run: ../tools/methyldackel_extract.cwl
    'sbg:x': 2831.132080078125
    'sbg:y': 2103.622802734375
  - id: methylation_calling_spike_in
    in:
      - id: CTOB
        source: methyldackel_ctob
      - id: CTOT
        source: methyldackel_ctot
      - id: OB
        source: methyldackel_ob
      - id: OT
        source: methyldackel_ot
      - id: bam
        source: index_spike_in_bam/bam_sorted_indexed
      - id: min_depth
        source: methyldackel_min_depth
      - id: min_mapq
        source: methyldackel_min_mapq
      - id: min_phred
        source: methyldackel_min_phred
      - id: nCTOB
        source: methyldackel_nctob
      - id: nCTOT
        source: methyldackel_nctot
      - id: nOB
        source: methyldackel_nob
      - id: nOT
        source: methyldackel_not
      - id: noCG
        source: methyldackel_noCG
      - id: output_basename
        source: sample_id
        valueFrom: $(self + "_spike_in")
      - id: reference
        source: reference
      - id: threads
        source: max_threads
    out:
      - id: mcall_bedgraph
    run: ../tools/methyldackel_extract.cwl
    'sbg:x': 3848.409912109375
    'sbg:y': 1283.5234375
  - id: qc_post_mapping
    in:
      - id: bam
        source: index_bam/bam_sorted_indexed
    out:
      - id: fastqc_html
      - id: fastqc_zip
    run: ../tools/fastqc.cwl
    doc: fastqc - quality control for reads after mapping and duplication removal
    'sbg:x': 3405.511962890625
    'sbg:y': -539.9010620117188
  - id: remove_duplicates
    in:
      - id: bam_sorted
        source: sorting_merged_bam/bam_sorted
    out:
      - id: bam_duprem
      - id: picard_markdup_log
    run: ../tools/picard_markdup.cwl
    doc: picard markdup - emoves duplicates from a single sorted bam file.
    'sbg:x': 1912.4864501953125
    'sbg:y': -388.3721008300781
  - id: sorting_merged_bam
    in:
      - id: bam_unsorted
        source: lane_replicate_merging/bam_merged
    out:
      - id: bam_sorted
    run: ../tools/samtools_sort.cwl
    doc: samtools sort - sorting of merged bam
    'sbg:x': 1469.843017578125
    'sbg:y': -376.81353759765625
  - id: map
    in:
      - id: fastq1
        source: fastq1
      - id: fastq2
        source: fastq2
      - id: is_non_directional
        default: false
        source: is_non_directional
      - id: reference
        source: reference
      - id: threads
        source: max_threads
    out:
      - id: sam
    run: ../tools/bwameth.cwl
    scatter:
      - fastq1
      - fastq2
    scatterMethod: dotproduct
    'sbg:x': 670.3844604492188
    'sbg:y': -368.344970703125
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: SubworkflowFeatureRequirement
