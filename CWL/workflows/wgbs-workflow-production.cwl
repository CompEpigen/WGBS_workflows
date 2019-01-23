cwlVersion: v1.0
class: Workflow
requirements:
  ScatterFeatureRequirement: {}
  MultipleInputFeatureRequirement: {}
  SubworkflowFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement:
    expressionLib:
    - var new_ext = function() { var ext=inputs.bai?'.bai':inputs.csi?'.csi':'.bai';
      return inputs.input.split('/').slice(-1)[0]+ext; };
    - var prepend = function(array,prefix) {  var file_array = []; for(var i=0; i<array.length;
      i++){ file_array[i] = prefix + '/' + array[i]; } return file_array; };
inputs:
  analysis_dir: string
  inp_read1:
    items: string
    type: array
  inp_read2:
    items: string
    type: array
  fastq_batch_size: int
  reference_fasta: string
  max_reads: int
  chromosomes:
    items: string
    type: array
  chr_prefix: string
  conversion_chr_name: string
  lib_type: string
  trimmomatic_adapters_file: string
  illuminaclip: string
  trimmomatic_phred: string
  trimmomatic_leading: int
  trimmomatic_trailing: int
  trimmomatic_crop: int
  trimmomatic_headcrop: int
  trimmomatic_tailcrop: int
  trimmomatic_minlen: int
  trimmomatic_avgqual: int
  pileometh_min_phred: int
  pileometh_min_depth: int
  pileometh_min_mapq: int
  pileometh_ot: string
  pileometh_ob: string
  pileometh_ctot: string
  pileometh_ctob: string
  pileometh_not: string
  pileometh_nob: string
  pileometh_nctot: string
  pileometh_nctob: string
  temp_dir: string
  trimmomatic_jar_file: string
  clean_raw_fastqs: boolean
  clean_trimmed_fastqs: boolean
  clean_primary_sams: boolean
  clean_primary_bams: boolean
  clean_chr_bams: boolean
  clean_merged_bams: boolean
  clean_dup_rm_bams: boolean
steps:
  split_read1_files:
    scatter: '#split_read1_files/file'
    run: ../tools_string/split-compressed-files.yml
    out:
    - output
    in:
      file: inp_read1
      file_dir: temp_dir
      size: fastq_batch_size
      suffix:
        source: inp_read1
        valueFrom: _$(inputs.file.substr(inputs.file.lastIndexOf('/')+1,inputs.file.lastIndexOf('.fastq')))_$(inputs.size/1000)k_R1.fastq
      max_reads: max_reads
  split_read2_files:
    scatter: '#split_read2_files/file'
    run: ../tools_string/split-compressed-files.yml
    out:
    - output
    in:
      file: inp_read2
      file_dir: temp_dir
      size: fastq_batch_size
      suffix:
        source: inp_read2
        valueFrom: _$(inputs.file.substr(inputs.file.lastIndexOf('/')+1,inputs.file.lastIndexOf('.fastq')))_$(inputs.size/1000)k_R2.fastq
      max_reads: max_reads
  flatten1:
    run: ../tools_string/flatten_fastq_arrays.yml
    out:
    - flattened_fastq_array
    in:
      fastq_arrays: split_read1_files/output
  flatten2:
    run: ../tools_string/flatten_fastq_arrays.yml
    out:
    - flattened_fastq_array
    in:
      fastq_arrays: split_read2_files/output
  adaptor_trimming:
    run: ../tools_string/trimmomatic.yml
    scatterMethod: dotproduct
    scatter: []
    out:
    - output_read1_trimmed_file
    - output_read2_trimmed_paired_file
    - output_log_file
    in:
      file_dir: temp_dir
      input_read1_fastq_file: flatten1/flattened_fastq_array
      input_read2_fastq_file: flatten2/flattened_fastq_array
      java_opts:
        default: -XX:-UseCompressedClassPointers -Xmx3000M -verbose
      trimmomatic_jar_path: trimmomatic_jar_file
      end_mode:
        default: PE
      input_adapters_file: trimmomatic_adapters_file
      phred: trimmomatic_phred
      illuminaclip: illuminaclip
      leading: trimmomatic_leading
      trailing: trimmomatic_trailing
      crop: trimmomatic_crop
      headcrop: trimmomatic_headcrop
      tailcrop: trimmomatic_tailcrop
      minlen: trimmomatic_minlen
      avgqual: trimmomatic_avgqual
      nthreads:
        default: 10
      log_filename:
        source: flatten1/flattened_fastq_array
        valueFrom: $(inputs.input_read1_fastq_file.substr(0, inputs.input_read1_fastq_file.lastIndexOf('.')))_trimming.log
  clean_raw_fastqs_1:
    run: ../tools_string/cleanup.yml
    out:
    - output
    in:
      files: flatten1/flattened_fastq_array
      data_dir: temp_dir
      do_clean: clean_raw_fastqs
      previous_step: adaptor_trimming/output_log_file
  clean_raw_fastqs_2:
    run: ../tools_string/cleanup.yml
    out:
    - output
    in:
      files: flatten2/flattened_fastq_array
      data_dir: temp_dir
      do_clean: clean_raw_fastqs
      previous_step: adaptor_trimming/output_log_file
  alignment:
    run: ../tools_string/bwameth.yml
    scatterMethod: dotproduct
    scatter: []
    out:
    - ''
    in:
      file_dir: temp_dir
      reference: reference_fasta
      read1: adaptor_trimming/output_read1_trimmed_file
      read2: adaptor_trimming/output_read2_trimmed_paired_file
      filename:
        source: adaptor_trimming/output_read1_trimmed_file
        valueFrom: $(inputs.read1.substr(0, inputs.read1.lastIndexOf('.')))
      cleanup_check1: clean_raw_fastqs_1/output
      cleanup_check2: clean_raw_fastqs_2/output
      threads:
        default: 4
      lib_type: lib_type
  clean_trimmed_fastqs_1:
    run: ../tools_string/cleanup.yml
    out:
    - output
    in:
      files: adaptor_trimming/output_read1_trimmed_file
      data_dir: temp_dir
      do_clean: clean_trimmed_fastqs
      previous_step: alignment/alignment
  clean_trimmed_fastqs_2:
    run: ../tools_string/cleanup.yml
    out:
    - output
    in:
      files: adaptor_trimming/output_read2_trimmed_paired_file
      data_dir: temp_dir
      do_clean: clean_trimmed_fastqs
      previous_step: alignment/alignment
  flag_stat_aligned:
    scatter: '#flag_stat_aligned/input_bam_file'
    run: ../tools_string/samtools-flagstat.yml
    out:
    - output
    in:
      file_dir: temp_dir
      input_bam_file: alignment/alignment
      cleanup_check1: clean_trimmed_fastqs_1/output
      cleanup_check2: clean_trimmed_fastqs_2/output
  sam_to_bam:
    scatter: '#sam_to_bam/input'
    run: ../tools_string/samtools-view.yml
    out:
    - output
    in:
      input: alignment/alignment
      isbam:
        default: 'true'
      output_name:
        source: alignment/alignment
        valueFrom: $(inputs.input.substr(0, inputs.input.lastIndexOf('.'))).bam
  clean_sams:
    run: ../tools_string/cleanup.yml
    out:
    - output
    in:
      files: alignment/alignment
      data_dir: temp_dir
      do_clean: clean_primary_sams
      previous_step1: sam_to_bam/output
      previous_step2: flag_stat_aligned/output
  split_by_chromosome:
    scatter: '#split_by_chromosome/input_bam_file'
    run: ../tools_string/bamtools-split.yml
    out:
    - output_bam_files
    in:
      file_dir: temp_dir
      input_bam_file: sam_to_bam/output
      split_options:
        default: reference
      ref_prefix: chr_prefix
      cleanup_check: clean_sams/output
  clean_bams:
    run: ../tools_string/cleanup.yml
    out:
    - output
    in:
      files: sam_to_bam/output
      data_dir: temp_dir
      do_clean: clean_primary_bams
      previous_step: split_by_chromosome/output_bam_files
  fix_all_bams:
    scatter: '#fix_all_bams/array_of_bams'
    run: ../tools_string/fix-all-bam-files.yml
    out:
    - array_of_fixed_bams
    in:
      file_dir: temp_dir
      array_of_bams: split_by_chromosome/output_bam_files
      cleanup_check: clean_bams/output
  clean_split_bams:
    scatter: '#clean_split_bams/files'
    run: ../tools_string/cleanup.yml
    out:
    - output
    in:
      files: split_by_chromosome/output_bam_files
      data_dir: temp_dir
      do_clean: clean_chr_bams
      previous_step: fix_all_bams/array_of_fixed_bams
  rearrange_bams:
    run: ../tools_string/rearrange_bams.yml
    out:
    - bam_arrays_per_chr
    - chrom_names
    in:
      bam_arrays: fix_all_bams/array_of_fixed_bams
      chromosomes: chromosomes
      cleanup_check: clean_split_bams/output
  bam_merging:
    run: ../tools_string/picard-MergeSamFiles.yml
    scatterMethod: dotproduct
    scatter: []
    out:
    - mergeSam_output
    in:
      file_dir: temp_dir
      inputFileName_mergedSam: rearrange_bams/bam_arrays_per_chr
      outputFileName_mergedSam:
        source: rearrange_bams/chrom_names
        valueFrom: ${ return inputs.file_dir + '/' + self + '_merged.bam'; }
      tmpdir: temp_dir
      createIndex:
        default: 'true'
  clean_fixed_bams:
    scatter: '#clean_fixed_bams/files'
    run: ../tools_string/cleanup.yml
    out:
    - output
    in:
      files: fix_all_bams/array_of_fixed_bams
      data_dir: temp_dir
      do_clean: clean_chr_bams
      previous_step: bam_merging/mergeSam_output
  index_bam_file:
    scatter: '#index_bam_file/input'
    run: ../tools_string/samtools-index.yml
    out:
    - index
    in:
      file_dir: temp_dir
      input: bam_merging/mergeSam_output
      bai:
        default: true
      cleanup_check: clean_fixed_bams/output
  duplicates_removal:
    scatter: '#duplicates_removal/inputFileName_markDups'
    run: ../tools_string/picard-MarkDuplicates.yml
    out:
    - markDups_output
    in:
      inputFileName_markDups:
        source: bam_merging/mergeSam_output
        valueFrom: $([self])
      outputFileName_markDups:
        source: bam_merging/mergeSam_output
        valueFrom: $(inputs.inputFileName_markDups.substr(0, inputs.inputFileName_markDups.lastIndexOf('.')))_dupl_removed.bam
      tmpdir: temp_dir
      removeDuplicates:
        default: 'true'
      metricsFile:
        source: bam_merging/mergeSam_output
        valueFrom: $(inputs.inputFileName_markDups.substr(0, inputs.inputFileName_markDups.lastIndexOf('.')))_duplicate_metrics.txt
      createIndex:
        default: 'true'
  clean_merged_bams:
    run: ../tools_string/cleanup.yml
    out:
    - output
    in:
      files: bam_merging/mergeSam_output
      data_dir: temp_dir
      do_clean: clean_merged_bams
      previous_step1: duplicates_removal/markDups_output
      previous_step2: index_bam_file/index
  flag_stat_dup_removed:
    scatter: '#flag_stat_dup_removed/input_bam_file'
    run: ../tools_string/samtools-flagstat.yml
    out:
    - output
    in:
      file_dir: temp_dir
      input_bam_file: duplicates_removal/markDups_output
  insert_size_dist:
    scatter: '#insert_size_dist/inputFileName_insertSize'
    run: ../tools_string/picard-InsertSizeMetric.yml
    out:
    - insertSize_output
    in:
      file_dir: temp_dir
      inputFileName_insertSize: duplicates_removal/markDups_output
      outputFileName_insertSize:
        source: duplicates_removal/markDups_output
        valueFrom: $( inputs.inputFileName_insertSize.substr(0, inputs.inputFileName_insertSize.lastIndexOf('.'))
          + '_insert_size_metrix.txt')
      histogramFile:
        source: duplicates_removal/markDups_output
        valueFrom: $( inputs.inputFileName_insertSize.substr(0, inputs.inputFileName_insertSize.lastIndexOf('.'))
          + '_insert_size_histogram.pdf')
      tmpdir: temp_dir
      createIndex:
        default: 'true'
  mbias_calculation:
    scatter: '#mbias_calculation/bam_file'
    run: ../tools_string/pileometh-mbias.yml
    out:
    - mbias_file
    in:
      bam_file: duplicates_removal/markDups_output
      reference: reference_fasta
      mbiasfile_name:
        source: duplicates_removal/markDups_output
        valueFrom: $(inputs.bam_file.substr(inputs.bam_file.lastIndexOf('/')+1, inputs.bam_file.lastIndexOf('.')))_mbias
  mbias_calculation_trimmed:
    scatter: '#mbias_calculation_trimmed/bam_file'
    run: ../tools_string/methyldackel-mbias.yml
    out:
    - mbias_file
    in:
      bam_file: duplicates_removal/markDups_output
      reference: reference_fasta
      mbiasfile_name:
        source: duplicates_removal/markDups_output
        valueFrom: $(inputs.bam_file.substr(inputs.bam_file.lastIndexOf('/')+1, inputs.bam_file.lastIndexOf('.')))_mbiasTrimmed
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
  methylation_calling:
    scatter: '#methylation_calling/bam_file'
    run: ../tools_string/methyldackel-extract.yml
    out:
    - methcall_bed
    in:
      file_dir: temp_dir
      bam_file: duplicates_removal/markDups_output
      bedfile_name:
        source: duplicates_removal/markDups_output
        valueFrom: $(inputs.bam_file.substr(inputs.bam_file.lastIndexOf('/')+1, inputs.bam_file.lastIndexOf('.')))
      reference: reference_fasta
      noCG:
        default: false
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
  merge_meth_calls:
    run: ../tools_string/methcall-merger.yml
    out:
    - merged_bed_file
    in:
      file_dir: temp_dir
      input_bed_files:
        source: methylation_calling/methcall_bed
        valueFrom: $(prepend(self, inputs.file_dir))
      output_file_name:
        source: temp_dir
        valueFrom: methylation_calls_CpG.bedGraph
  find_lambda_file:
    run: ../tools_string/find_lambda.yml
    out:
    - lambda_bam
    in:
      split_files: duplicates_removal/markDups_output
      filename: conversion_chr_name
  methylation_calling_lambda:
    run: ../tools_string/methyldackel-extract.yml
    out:
    - methcall_bed
    in:
      file_dir: temp_dir
      reference: reference_fasta
      noCG:
        default: true
      bedfile_name: conversion_chr_name
      OT: pileometh_ot
      OB: pileometh_ob
      CTOT: pileometh_ctot
      CTOB: pileometh_ctob
      nOT: pileometh_not
      nOB: pileometh_nob
      nCTOT: pileometh_nctot
      nCTOB: pileometh_nctob
      bam_file: find_lambda_file/lambda_bam
  clean_rmdup_bams:
    run: ../tools_string/cleanup.yml
    out:
    - output
    in:
      files: duplicates_removal/markDups_output
      data_dir: temp_dir
      do_clean: clean_dup_rm_bams
      previous_step1: methylation_calling/methcall_bed
      previous_step2: methylation_calling_lambda/methcall_bed
  conversion_estimation_lambda:
    run: ../tools_string/bisulfite-conversion-lambda.yml
    out:
    - bisulfite_conversion_file
    in:
      file_dir: temp_dir
      input_bed_file:
        source: methylation_calling_lambda/methcall_bed
        valueFrom: $(inputs.file_dir + '/' + self)
      output_file_name:
        source: temp_dir
        valueFrom: bisulfite_conversion.txt
outputs:
  trimming_report:
    type:
      items:
      - 'null'
      - File
      type: array
    outputSource: adaptor_trimming/output_log_file
  alignment_flagstat:
    type:
      items: string
      type: array
    outputSource: flag_stat_aligned/output
  duplicate_removed_flagstat:
    type:
      items: string
      type: array
    outputSource: flag_stat_dup_removed/output
  insert_size_metrics:
    type:
      items: string
      type: array
    outputSource: insert_size_dist/insertSize_output
  mbias_file:
    type:
      items:
        items: File
        type: array
      type: array
    outputSource: mbias_calculation/mbias_file
  mbias_file_trimmed:
    type:
      items:
        items: File
        type: array
      type: array
    outputSource: mbias_calculation_trimmed/mbias_file
  methylation_calls:
    type: File
    outputSource: merge_meth_calls/merged_bed_file
  bisulfite_conversion_estimation:
    type: File
    outputSource: conversion_estimation_lambda/bisulfite_conversion_file

