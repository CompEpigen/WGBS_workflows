cwlVersion: v1.0
class: Workflow
requirements:
  ScatterFeatureRequirement: {}
  MultipleInputFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement:
    expressionLib:
    - var new_ext = function() { var ext=inputs.bai?'.bai':inputs.csi?'.csi':'.bai';
      return inputs.input.split('/').slice(-1)[0]+ext; };
inputs:
  file_dir: string
  array_of_bams:
    items: string
    type: array
steps:
  fix_bams:
    scatter: '#fix_bams/inputBAMFile'
    run: fix-bam-file.yml
    out:
    - fixBam_output
    in:
      file_dir: file_dir
      inputBAMFile: array_of_bams
      outputFileName:
        source: array_of_bams
        valueFrom: $(inputs.inputBAMFile.substr(0,inputs.inputBAMFile.lastIndexOf('.'))
          + '.fixed').bam
outputs:
  array_of_fixed_bams:
    type:
      items: string
      type: array
    outputSource: fix_bams/fixBam_output

