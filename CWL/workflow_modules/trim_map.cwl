cwlVersion: v1.0
class: Workflow

inputs:
  fastq1:
    type: File
    # type: array
    # items: File
  fastq2:
    type: File
    # type: array
    # items: File
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

  # # optional parameter that likely need revision:
  is_non_directional:
    type: boolean
    default: false
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

steps:
  qc_raw:
    doc: fastqc - quality control for raw fastqs
    run: "../tools/fastqc.cwl"
    in:
      fastq1: fastq1
      fastq2: fastq2
    out:
      - fastqc_zip
      - fastqc_html

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

  qc_trimmed:
    doc: fastqc - quality control for raw fastqs
    run: "../tools/fastqc.cwl"
    in:
      fastq1: adaptor_trimming/fastq1_trimmed
      fastq2: adaptor_trimming/fastq2_trimmed
    out:
      - fastqc_zip
      - fastqc_html

  mapping:
    run: "../tools/bwameth.cwl"
    in:
      fastq1: adaptor_trimming/fastq1_trimmed
      fastq2: adaptor_trimming/fastq2_trimmed
      reference: reference
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

  bam:
    type: File
    outputSource: sort_bam/bam_sorted

  raw_fastqc_zip:
    type:
      type: array
      items: File
    outputSource: qc_raw/fastqc_zip
  raw_fastqc_html:
    type:
      type: array
      items: File
    outputSource: qc_raw/fastqc_html
  trimmed_fastqc_zip:
    type:
      type: array
      items: File
    outputSource: qc_trimmed/fastqc_zip
  trimmed_fastqc_html:
    type:
      type: array
      items: File
    outputSource: qc_trimmed/fastqc_html