cwlVersion: v1.0
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 20000
    tmpdirMin: 10000
  DockerRequirement:
    dockerPull: kerstenbreuer/picard_tools:2.17.4
  
baseCommand: ["java", "-jar"]
arguments:
  - valueFrom: "AddOrReplaceReadGroups"
    position: 2
  - valueFrom: $(inputs.bam_withoutRG.nameroot + "_RG.bam")
    prefix: "O="
    position: 3
    separate: false
  - valueFrom: "ID=readGroup_name"
    position: 4
  - valueFrom: "LB=readGroup_name"
    position: 5
  - valueFrom: "PL=illumina"
    position: 6
  - valueFrom: "PU=run"
    position: 7    
  - valueFrom: "SM=sample_name"
    position: 8     

inputs:
  bam_withoutRG:
    type: File
    inputBinding:
      prefix: "I="
      separate: false
      position: 3
  path_to_picards:
    type: string
    default: "/bin/picard.jar"
    inputBinding:
      position: 1
      
outputs:
  bam_withRG:
    type: File
    outputBinding:
      glob: $(inputs.bam_withoutRG.nameroot + "_RG.bam")
    