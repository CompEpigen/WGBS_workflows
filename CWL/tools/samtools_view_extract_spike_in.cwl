cwlVersion: v1.0
class: CommandLineTool
hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 10000
    tmpdirMin: 10000
  DockerRequirement:
    dockerPull: kerstenbreuer/samtools:1.7

baseCommand: ["samtools", "view"]

inputs:
  bam:
    doc: aligned reads in sam or bam format
    type: File
    secondaryFiles: .bai
    inputBinding:
      position: 2
  region:
    doc: region to extract
    type: string?
    inputBinding:
      position: 3

arguments:
  - valueFrom: -h
    position: 1
    # include the headers
  - valueFrom: -b
    position: 1
    # output in bam format

stdout: $( inputs.bam.nameroot + "_spike_in.bam" )

outputs:
  bam_spike_in:
    type: stdout
  
  
