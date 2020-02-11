cwlVersion: v1.0
class: CommandLineTool
hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 5000
    tmpdirMin: 10000
  DockerRequirement:
    dockerPull: kerstenbreuer/trim_galore:0.6.4_2.6_0.11.8_scPBAT
  
baseCommand: "fastqc"
arguments: 
  - valueFrom: $(runtime.outdir)
    prefix: "-o"
    # specifies output directory
  - valueFrom: "--noextract"
    # reported file will be zipped

inputs:
  reads:
    type: File?
    inputBinding:
      position: 1
  bam:
    type: File?
    inputBinding:
      position: 1

outputs:
  fastqc_zip:
    doc: all data e.g. figures
    type:
      type: array
      items: File
    outputBinding:
      glob: "*_fastqc.zip"
  fastqc_html:
    doc: html report showing results from zip
    type:
      type: array
      items: File
    outputBinding:
      glob: "*_fastqc.html"
    
