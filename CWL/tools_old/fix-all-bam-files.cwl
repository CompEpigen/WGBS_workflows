cwlVersion: v1.0
class: Workflow
requirements:
  ScatterFeatureRequirement: {}
  MultipleInputFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement:
    expressionLib:
    - var new_ext = function() { var ext=inputs.bai?'.bai':inputs.csi?'.csi':'.bai';
      return inputs.input.path.split('/').slice(-1)[0]+ext; };
inputs:
  array_of_bams:
    items: File
    type: array
steps:
  fix_bams:
    scatter: '#fix_bams/inputBAMFile'
    run: fix-bam-file.yml
    out:
    - fixBam_output
    in:
      inputBAMFile: array_of_bams
      outputFileName:
        source: array_of_bams
        valueFrom: $(inputs.inputBAMFile.basename.substr(0,inputs.inputBAMFile.basename.lastIndexOf('.'))
          + '.fixed').bam
outputs:
  array_of_fixed_bams:
    type:
      items: File
      type: array
    outputSource: fix_bams/fixBam_output

