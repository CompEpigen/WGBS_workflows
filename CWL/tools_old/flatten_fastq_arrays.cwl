cwlVersion: v1.0
class: ExpressionTool
requirements:
- class: InlineJavascriptRequirement
inputs:
- type:
    items:
      items: File
      type: array
    type: array
  id: fastq_arrays
expression: "${\nvar rearranged = []; for(var i=0; i<inputs.fastq_arrays.length; i++){\
  \ for(var j=0; j<inputs.fastq_arrays[i].length; j++){ rearranged[rearranged.length]\
  \ = inputs.fastq_arrays[i][j]; } }\nvar output = {}; output['flattened_fastq_array']\
  \ = rearranged; return output; }"
outputs:
- type:
    items: File
    type: array
  id: flattened_fastq_array

