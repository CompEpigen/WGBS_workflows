cwlVersion: v1.0
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}
hints:
  ResourceRequirement:
    coresMin: $( inputs.threads )
    ramMin: 28000
    tmpdirMin: 10000
  DockerRequirement:
    dockerPull: dukegcb/trimmomatic:latest

baseCommand: ["java"]
arguments: 
  - valueFrom: |
      ${
        if ( inputs.fastq2 == null){
          return "SE"
        } else {
          return "PE"
        }
      }
    position: 3
  - valueFrom: $(inputs.fastq1.nameroot + ".trimmomatic.log")
    position: 4
    prefix: -trimlog
  - valueFrom: $(inputs.fastq1.nameroot + "_trimmed.fastq")
    position: 7
  - valueFrom: $(inputs.fastq1.nameroot + "_trimmed_unpaired.fastq")
    position: 8
  - valueFrom: |
      ${ 
        if( inputs.fastq2 != null){
          return inputs.fastq2.nameroot + "_trimmed.fastq"
        }
        else{
          return null
        }
      }
    position: 9
  - valueFrom: |
      ${ 
        if( inputs.fastq2 != null){
          return inputs.fastq2.nameroot + "_trimmed_unpaired.fastq"
        }
        else{
          return null
        }
      }
    position: 10
  - valueFrom: $("ILLUMINACLIP:" + inputs.adapters_file.path + ":"+ inputs.illuminaclip)
    position: 11

inputs:
  fastq1:
    doc: FASTQ file for input read (read R1 in Paired End mode)
    type: File
    inputBinding:
      position: 5
  fastq2:
    doc: FASTQ file for read R2 in Paired End mode
    type: File?
    inputBinding:
      position: 6
  adapters_file:
    doc: FASTA file containing adapters, PCR sequences, etc. It is used to search
      for and remove these sequences in the input FASTQ file(s)
    type: File

  path_to_trimmomatic:
    doc: |
      default path matching the applied docker container; 
      if the container is not used please adapt
    type: string
    default: "/usr/share/java/trimmomatic.jar"
    inputBinding:
      position: 2
      prefix: "-jar"
  
  # additional arguments
  java_opts:
    type: string?
    inputBinding:
      shellQuote: false
      position: 1
  threads:
    doc: Number of threads
    type: int
    default: 10
    inputBinding:
      position: 4
      prefix: -threads
  phred:
    type: string
    default: '64'
    inputBinding:
      position: 4
      prefix: -phred
      separate: false
  log_filename:
    type: string?
    inputBinding:
      position: 4
      prefix: -trimlog
  illuminaclip:
    type: string
  slidingwindow:
    type: string?
    inputBinding:
      position: 15
      prefix: 'SLIDINGWINDOW:'
      separate: false
  leading:
    type: int?
    inputBinding:
      position: 14
      prefix: 'LEADING:'
      separate: false
  trailing:
    type: int?
    inputBinding:
      position: 14
      prefix: 'TRAILING:'
      separate: false
  crop:
    type: int?
    inputBinding:
      position: 13
      prefix: 'CROP:'
      separate: false
  # headcrop:
  #   type: int?
  #   inputBinding:
  #     position: 13
  #     prefix: 'HEADCROP:'
  #     separate: false
  # tailcrop:
  #   type: int?
  #   inputBinding:
  #     position: 13
  #     prefix: 'TAILCROP:'
  #     separate: false
  minlen:
    type: int?
    inputBinding:
      position: 100
      prefix: 'MINLEN:'
      separate: false
  avgqual:
    type: int?
    inputBinding:
      position: 101
      prefix: 'AVGQUAL:'
      separate: false

outputs:
  fastq1_trimmed:
    type: File
    outputBinding:
      glob: $(inputs.fastq1.nameroot + "_trimmed.fastq")
  fastq2_trimmed:    
    type: File?
    outputBinding:
      glob: |
        ${ 
          if( inputs.fastq2 != null){
            return inputs.fastq2.nameroot + "_trimmed.fastq"
          }
          else{
            return null
          }
        }
  fastq1_trimmed_unpaired:    
    type: File?
    outputBinding:
      glob: $(inputs.fastq1.nameroot + "_trimmed_unpaired.fastq")
  fastq2_trimmed_unpaired:    
    type: File?
    outputBinding:
      glob: |
        ${ 
          if( inputs.fastq2 != null){
            return inputs.fastq2.nameroot + "_trimmed_unpaired.fastq"
          }
          else{
            return null
          }
        }
  trimmomatic_log: # can be used by multiqc
    type: File
    outputBinding:
      glob:  $(inputs.fastq1.nameroot + ".trimmomatic.log")
