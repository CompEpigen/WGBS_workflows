cwlVersion: v1.0
class: Workflow

requirements:
  ScatterFeatureRequirement: {}
  MultipleInputFeatureRequirement: {}
  SubworkflowFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}

inputs:
  # raw data:
  fastq1:
    type: File
    # type: array
    # items: File
  fastq2:
    type: File
    # type: array
    # items: File
#   reference_fasta: 
#     type: File
#     secondaryFiles:
#       - .fai
#       - .bwameth.c2t
#       - .bwameth.c2t.amb
#       - .bwameth.c2t.ann
#       - .bwameth.c2t.bwt
#       - .bwameth.c2t.pac
#       - .bwameth.c2t.sa
  trimmomatic_adapters_file: File

  # # optional parameter that likely need revision:
  # lib_type:
  #   type: string
  #   default: "directional"
  max_threads:
    type: int
    default: 10
  # chr_prefix:
  #   doc: #!
  #   type: string
  #   default: "chr"
  # chromosomes: 
  #   type: array
  #   items: string
  #   default: [ 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'MT','Lambda' ]


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
  # pileometh_min_phred: 
  #   type: int
  #   default: 0
  # pileometh_min_depth: 
  #   type: int
  #   default: 1
  # pileometh_min_mapq: 
  #   type: int
  #   default: 0
  # pileometh_ot: 
  #   type: int
  #   default: "0,0,0,0"
  # pileometh_ob: 
  #   type: int
  #   default: "0,0,0,0"
  # pileometh_ctot: 
  #   type: int
  #   default: "0,0,0,0"
  # pileometh_ctob: 
  #   type: int
  #   default: "0,0,0,0"
  # pileometh_not: 
  #   type: int
  #   default: "10,10,10,10"
  # pileometh_nob: 
  #   type: int
  #   default: "10,10,10,10"
  # pileometh_nctot: 
  #   type: int
  #   default: "10,10,10,10"
  # pileometh_nctob: 
  #   type: int
  #   default: "10,10,10,10"

# conversion_chr_name: string

steps:
  adaptor_trimming:
    run: ../tools/trimmomatic.cwl
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
      headcrop: trimmomatic_headcrop
      tailcrop: trimmomatic_tailcrop
      minlen: trimmomatic_minlen
      avgqual: trimmomatic_avgqual
      nthreads: max_threads
    out:
    - fastq1_trimmed
    - fastq2_trimmed
    - fastq1_trimmed_unpaired
    - fastq2_trimmed_unpaired
    - trimmomatic_log

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

