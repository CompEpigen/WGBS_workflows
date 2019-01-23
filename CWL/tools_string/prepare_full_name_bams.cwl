cwlVersion: v1.0
class: ExpressionTool
requirements:
- class: InlineJavascriptRequirement
inputs:
- type:
    items: string
    type: array
  id: bam_array
- type: string
  id: dir
expression: "${\nvar dir = inputs.dir; var full_named = [];\nfor(var fi=0; fi<inputs.bam_array.length;\
  \ fi++){ full_named[fi] = dir + '/' + inputs.bam_array[fi]; }\nvar output = {};\
  \ output['bam_array_full_name'] = full_named; return output; }"
outputs:
- type:
    items: string
    type: array
  id: bam_array_full_name

