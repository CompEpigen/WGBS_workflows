cwlVersion: v1.0
class: Workflow

requirements:
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}

inputs:
  sc_id:
    type: string[]
  fastq1:
    type: File[]
  fastq2: 
    type: File[]
  reference_index:
    type: File
    secondaryFiles:
      - _strands
  fastqs_are_gzipped:
    type: boolean
    default: true
  # adapter1: 
  #   type: string?
  # adapter2:
  #   type: string?
        
steps:
  # qc_raw:
  #   doc: fastqc - quality control for trimmed fastq
  #   run: "../tools/fastqc.cwl"
  #   in:
  #     fastq1:
  #       source: fastq1
  #     fastq2:
  #       source: fastq2
  #   out:
  #     - fastqc_zip
  #     - fastqc_html

  # adaptor_trimming_and_qc_trimmed:
  #   doc: trim galore - adapter trimming using trim_galore
  #   run: "../tools/trim_galore.cwl"
  #   in:
  #     fastq1:
  #       source: fastq1
  #     fastq2:
  #       source: fastq2
  #     adapter1:
  #       source: adapter1
  #     adapter2:
  #       source: adapter2   
  #   out:
  #     - fastq1_trimmed
  #     - fastq2_trimmed
  #     - fastq1_trimmed_unpaired
  #     - fastq2_trimmed_unpaired
  #     - trim_galore_log
  #     - trimmed_fastqc_html
  #     - trimmed_fastqc_zip

  build_meta_tsv:
    in:
      sc_id: sc_id
      fastq1: fastq1
      fastq2: fastq2
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      hints:
        # ResourceRequirement:
          # coresMin: 1
          # ramMin: 15000
          #tmpdirMin: 10000
        DockerRequirement:
          dockerPull: python:3.7.4
      requirements:
        InitialWorkDirRequirement:
          listing: 
            - entryname: build_meta_tsv.py
              entry: | 
                import sys
                args = sys.argv
                sc_ids = args[1].split(",")
                fastq1 = args[2].split(",")
                fastq2 = args[3].split(",")
                for i, sc_id in enumerate(sc_ids):
                  print(sc_id + "\t" + fastq1[i] + "\t" + fastq2[i])
      baseCommand: ["python3", "build_meta_tsv.py"]
      inputs:
        sc_id:
          type: string[]
          inputBinding:
            position: 1
            itemSeparator: ","
        fastq1:
          type: File[]
          inputBinding:
            position: 2
            itemSeparator: ","
        fastq2:
          type: File[]
          inputBinding:
            position: 3
            itemSeparator: ","
      stdout: meta.tsv
      outputs:
        meta_tsv:
          type: stdout
    out:
      - meta_tsv

  align_and_call_meth:
    run: ../tools/fame_sc.cwl
    in:
      load_index: reference_index
      sc_input: build_meta_tsv/meta_tsv
      gzip_reads: fastqs_are_gzipped
    out:
      - sc_methylation_calls
      - bulk_methylation_calls   

outputs:
  # raw_fastqc_zip:
  #   type:
  #     type: array
  #     items: File
  #   outputSource: qc_raw/fastqc_zip
  # raw_fastqc_html:
  #   type:
  #     type: array
  #     items: File
  #   outputSource: qc_raw/fastqc_html
  # fastq1_trimmed:
  #   type: File
  #   outputSource: adaptor_trimming_and_qc_trimmed/fastq1_trimmed
  # fastq2_trimmed:
  #   type: File
  #   outputSource: adaptor_trimming_and_qc_trimmed/fastq2_trimmed
  # trim_galore_log:
  #   type:
  #     type: array
  #     items: File
  #   outputSource: adaptor_trimming_and_qc_trimmed/trim_galore_log
  # trimmed_fastqc_html:
  #   type:
  #     type: array
  #     items: File
  #   outputSource: adaptor_trimming_and_qc_trimmed/trimmed_fastqc_html
  # trimmed_fastqc_zip:
  #   type:
  #     type: array
  #     items: File
  #   outputSource: adaptor_trimming_and_qc_trimmed/trimmed_fastqc_zip
  #   type: File
  #   outputSource: sort_bam/bam_sorted
  meta_tsv:
    type: File
    outputSource: build_meta_tsv/meta_tsv
  sc_methylation_calls:
    type: File
    outputSource: align_and_call_meth/sc_methylation_calls
  bulk_methylation_calls:
    type: File
    outputSource: align_and_call_meth/bulk_methylation_calls


    