cwlVersion: v1.0
class: CommandLineTool
hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 10000
inputs:
  file_dir: string
  java_arg:
    type: string
    default: -Xmx4g
    inputBinding:
      position: 1
  inputBAMFile:
    doc: One or more input SAM or BAM files to fix
    type: string
  outputFileName: string
baseCommand: java
arguments:
- position: 2
  valueFrom: /ngs_share/tools/htsjdk/build/libs/htsjdk-2.7.0-3-g1c66107-SNAPSHOT-all.jar
  prefix: -cp
- position: 3
  valueFrom: htsjdk.samtools.FixBAMFile
- position: 4
  valueFrom: $(inputs.file_dir + '/' + inputs.inputBAMFile)
- position: 5
  valueFrom: $(inputs.file_dir + '/' + inputs.outputFileName)
outputs:
  fixBam_output:
    type: string
    outputBinding:
      outputEval: $(inputs.file_dir + '/' + inputs.outputFileName)

