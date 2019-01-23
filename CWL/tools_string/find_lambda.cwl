cwlVersion: v1.0
class: ExpressionTool
requirements:
- class: InlineJavascriptRequirement
inputs:
- type:
    items: string
    type: array
  id: split_files
- type: string
  id: filename
expression: ${ var index = inputs.split_files.length - 1; for(var i=0; i<inputs.split_files.length;
  i++){ if(inputs.split_files[i].indexOf(inputs.filename) !== -1){ index = i; } }
  var output = {}; output['lambda_bam'] = inputs.split_files[index]; return output;
  }
outputs:
- type: string
  id: lambda_bam

