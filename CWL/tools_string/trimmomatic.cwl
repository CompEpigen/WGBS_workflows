cwlVersion: v1.0
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}
hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 28000
inputs:
  file_dir: string
  java_opts:
    doc: JVM arguments should be a quoted, space separated list (e.g. "-Xms128m -Xmx512m")
    type: string?
    inputBinding:
      shellQuote: false
      position: 1
  trimmomatic_jar_path:
    type: File
    inputBinding:
      position: 2
      prefix: -jar
  end_mode:
    doc: "SE|PE\nSingle End (SE) or Paired End (PE) mode\n"
    type: string
    inputBinding:
      position: 3
  nthreads:
    doc: Number of threads
    type: int
    default: 12
    inputBinding:
      position: 4
      prefix: -threads
  phred:
    doc: "\"33\"|\"64\"\n-phred33 (\"33\") or -phred64 (\"64\") specifies the base\
      \ quality encoding. Default: -phred64\n"
    type: string
    default: '64'
    inputBinding:
      position: 4
      prefix: -phred
      separate: false
  log_filename:
    doc: "<ouptut log file name>\nSpecifying a trimlog file creates a log of all read\
      \ trimmings, indicating the following details:\n  the read name\n  the surviving\
      \ sequence length\n  the location of the first surviving base, aka. the amount\
      \ trimmed from the start\n  the location of the last surviving base in the original\
      \ read\n  the amount trimmed from the end\n<ouptut log file name>: filename\
      \ for the generated output log file.\n"
    type: string?
    inputBinding:
      position: 4
      prefix: -trimlog
  input_read1_fastq_file:
    doc: FASTQ file for input read (read R1 in Paired End mode)
    type: string
    inputBinding:
      position: 5
  input_read2_fastq_file:
    doc: FASTQ file for read R2 in Paired End mode
    type: string?
    inputBinding:
      position: 6
  input_adapters_file:
    doc: FASTA file containing adapters, PCR sequences, etc. It is used to search
      for and remove these sequences in the input FASTQ file(s)
    type: string
  illuminaclip:
    doc: "<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple\
      \ clip threshold>:<minAdapterLength>:<keepBothReads>\nFind and remove Illumina\
      \ adapters.\nREQUIRED:\n<fastaWithAdaptersEtc>: specifies the path to a fasta\
      \ file containing all the adapters, PCR sequences etc.\nThe naming of the various\
      \ sequences within this file determines how they are used. See below.\n<seedMismatches>:\
      \ specifies the maximum mismatch count which will still allow a full match to\
      \ be performed\n<palindromeClipThreshold>: specifies how accurate the match\
      \ between the two 'adapter ligated' reads must be\nfor PE palindrome read alignment.\n\
      <simpleClipThreshold>: specifies how accurate the match between any adapter\
      \ etc. sequence must be against a read\nOPTIONAL:\n<minAdapterLength>: In addition\
      \ to the alignment score, palindrome mode can verify\nthat a minimum length\
      \ of adapter has been detected. If unspecified, this defaults to 8 bases,\n\
      for historical reasons. However, since palindrome mode has a very low false\
      \ positive rate, this\ncan be safely reduced, even down to 1, to allow shorter\
      \ adapter fragments to be removed.\n<keepBothReads>: After read-though has been\
      \ detected by palindrome mode, and the\nadapter sequence removed, the reverse\
      \ read contains the same sequence information as the\nforward read, albeit in\
      \ reverse complement. For this reason, the default behaviour is to\nentirely\
      \ drop the reverse read. By specifying „true‟ for this parameter, the reverse\
      \ read will\nalso be retained, which may be useful e.g. if the downstream tools\
      \ cannot handle a\ncombination of paired and unpaired reads.\n"
    type: string
  slidingwindow:
    doc: "<windowSize>:<requiredQuality>\nPerform a sliding window trimming, cutting\
      \ once the average quality within the window falls\nbelow a threshold. By considering\
      \ multiple bases, a single poor quality base will not cause the\nremoval of\
      \ high quality data later in the read.\n<windowSize>: specifies the number of\
      \ bases to average across\n<requiredQuality>: specifies the average quality\
      \ required\n"
    type: string?
    inputBinding:
      position: 15
      prefix: 'SLIDINGWINDOW:'
      separate: false
  maxinfo:
    doc: "<targetLength>:<strictness>\nPerforms an adaptive quality trim, balancing\
      \ the benefits of retaining longer reads against the\ncosts of retaining bases\
      \ with errors.\n<targetLength>: This specifies the read length which is likely\
      \ to allow the\nlocation of the read within the target sequence to be determined.\n\
      <strictness>: This value, which should be set between 0 and 1, specifies the\n\
      balance between preserving as much read length as possible vs. removal of incorrect\n\
      bases. A low value of this parameter (<0.2) favours longer reads, while a high\
      \ value\n(>0.8) favours read correctness.\n"
    type: string?
    inputBinding:
      position: 15
      prefix: 'MAXINFO:'
      separate: false
  leading:
    doc: "<quality>\nRemove low quality bases from the beginning. As long as a base\
      \ has a value below this\nthreshold the base is removed and the next base will\
      \ be investigated.\n<quality>: Specifies the minimum quality required to keep\
      \ a base.\n"
    type: int?
    inputBinding:
      position: 14
      prefix: 'LEADING:'
      separate: false
  trailing:
    doc: "<quality>\nRemove low quality bases from the end. As long as a base has\
      \ a value below this threshold\nthe base is removed and the next base (which\
      \ as trimmomatic is starting from the 3‟ prime end\nwould be base preceding\
      \ the just removed base) will be investigated. This approach can be\nused removing\
      \ the special illumina „low quality segment‟ regions (which are marked with\n\
      quality score of 2), but we recommend Sliding Window or MaxInfo instead\n<quality>:\
      \ Specifies the minimum quality required to keep a base.\n"
    type: int?
    inputBinding:
      position: 14
      prefix: 'TRAILING:'
      separate: false
  crop:
    doc: "<length>\nRemoves bases regardless of quality from the end of the read,\
      \ so that the read has maximally\nthe specified length after this step has been\
      \ performed. Steps performed after CROP might of\ncourse further shorten the\
      \ read.\n<length>: The number of bases to keep, from the start of the read.\n"
    type: int?
    inputBinding:
      position: 13
      prefix: 'CROP:'
      separate: false
  headcrop:
    doc: "<length>\nRemoves the specified number of bases, regardless of quality,\
      \ from the beginning of the read.\n<length>: The number of bases to remove,\
      \ from the start of the read.\n"
    type: int?
    inputBinding:
      position: 13
      prefix: 'HEADCROP:'
      separate: false
  tailcrop:
    doc: "<length>\nRemoves the specified number of bases, regardless of quality,\
      \ from the end of the read.\n<length>: The number of bases to remove, from the\
      \ end of the read.\n"
    type: int?
    inputBinding:
      position: 13
      prefix: 'TAILCROP:'
      separate: false
  minlen:
    doc: "<length>\nThis module removes reads that fall below the specified minimal\
      \ length. If required, it should\nnormally be after all other processing steps.\
      \ Reads removed by this step will be counted and\nincluded in the „dropped reads‟\
      \ count presented in the trimmomatic summary.\n<length>:  Specifies the minimum\
      \ length of reads to be kept\n"
    type: int?
    inputBinding:
      position: 100
      prefix: 'MINLEN:'
      separate: false
  avgqual:
    doc: "<quality>\nDrop the read if the average quality is below the specified level\n\
      <quality>: Specifies the minimum average quality required to keep a read.\n"
    type: int?
    inputBinding:
      position: 101
      prefix: 'AVGQUAL:'
      separate: false
  tophred33:
    doc: This (re)encodes the quality part of the FASTQ file to base 33.
    type: boolean?
    inputBinding:
      position: 12
      prefix: TOPHRED33
      separate: false
  tophred64:
    doc: This (re)encodes the quality part of the FASTQ file to base 64.
    type: boolean?
    inputBinding:
      position: 12
      prefix: TOPHRED64
      separate: false
baseCommand: cd
arguments:
- position: -2
  valueFrom: $(inputs.file_dir)
- position: -1
  valueFrom: ;
- position: 0
  valueFrom: java
- position: 7
  valueFrom: $(inputs.input_read1_fastq_file.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/,
    '') + '_trimmed.fastq')
- position: 8
  valueFrom: "${\n  if (inputs.end_mode == \"PE\" && inputs.input_read2_fastq_file)\n\
    \    return inputs.input_read1_fastq_file.replace(/^.*[\\\\\\/]/, '').replace(/\\\
    .[^/.]+$/, '') + '_trimmed_unpaired.fastq';\n  return null;\n}\n"
- position: 9
  valueFrom: "${\n  if (inputs.end_mode == \"PE\" && inputs.input_read2_fastq_file)\n\
    \    return inputs.input_read2_fastq_file.replace(/^.*[\\\\\\/]/, '').replace(/\\\
    .[^/.]+$/, '') + '_trimmed.fastq';\n  return null;\n}\n"
- position: 10
  valueFrom: "${\n  if (inputs.end_mode == \"PE\" && inputs.input_read2_fastq_file)\n\
    \    return inputs.input_read2_fastq_file.replace(/^.*[\\\\\\/]/, '').replace(/\\\
    .[^/.]+$/, '') + '_trimmed_unpaired.fastq';\n  return null;\n}\n"
- position: 11
  valueFrom: $("ILLUMINACLIP:" + inputs.input_adapters_file + ":"+ inputs.illuminaclip)
- shellQuote: false
  position: 12
  valueFrom: 2>
- position: 13
  valueFrom: $(inputs.file_dir + '/' + inputs.log_filename.replace('.log','_stdout.log'))
outputs:
  output_read1_trimmed_file:
    type: string
    outputBinding:
      outputEval: $(inputs.input_read1_fastq_file.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/,
        '') + '_trimmed.fastq')
  output_read1_trimmed_unpaired_file:
    type: string?
    outputBinding:
      outputEval: "${\n  if (inputs.end_mode == \"PE\")\n    return inputs.file_dir\
        \ + '/' +inputs.input_read1_fastq_file.replace(/^.*[\\\\\\/]/, '').replace(/\\\
        .[^/.]+$/, '') + '_unpaired_trimmed.fastq';\n  return null;\n}\n"
  output_read2_trimmed_paired_file:
    type: string?
    outputBinding:
      outputEval: "${\n  if (inputs.end_mode == \"PE\" && inputs.input_read2_fastq_file)\n\
        \    return inputs.file_dir + '/' + inputs.input_read2_fastq_file.replace(/^.*[\\\
        \\\\/]/, '').replace(/\\.[^/.]+$/, '') + '_trimmed.fastq';\n  return null;\n\
        }\n"
  output_read2_trimmed_unpaired_file:
    type: string?
    outputBinding:
      outputEval: "${\n  if (inputs.end_mode == \"PE\" && inputs.input_read2_fastq_file)\n\
        \    return inputs.file_dir + '/' +inputs.input_read2_fastq_file.replace(/^.*[\\\
        \\\\/]/, '').replace(/\\.[^/.]+$/, '') + '_unpaired_trimmed.fastq';\n  return\
        \ null;\n}\n"
  output_log_file:
    doc: Trimmomatic Log file.
    type: string?
    outputBinding:
      outputEval: ${ return inputs.file_dir + '/' + inputs.log_filename}

