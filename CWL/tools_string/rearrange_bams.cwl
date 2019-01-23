cwlVersion: v1.0
class: ExpressionTool
requirements:
- class: InlineJavascriptRequirement
inputs:
- type:
    items:
      items: string
      type: array
    type: array
  id: bam_arrays
- type:
    items: string
    type: array
  id: chromosomes
expression: "${\nvar chroms = inputs.chromosomes; var rearranged = []; var found_chroms\
  \ = []; var any_found = false;\nfor(var ci=0; ci<chroms.length; ci++){ var chr_output\
  \ = []; any_found = false; for(var i=0; i<inputs.bam_arrays.length; i++){ for(var\
  \ j=0; j<inputs.bam_arrays[i].length; j++){ var this_file = inputs.bam_arrays[i][j];\
  \ if(this_file.indexOf('.' + chroms[ci] + '.') !== -1){ chr_output[chr_output.length]\
  \ = this_file; any_found = true; } } }\nif (any_found) { rearranged[rearranged.length]\
  \ = chr_output; found_chroms[found_chroms.length] = chroms[ci]; } } var output =\
  \ {}; output['bam_arrays_per_chr'] = rearranged; output['chrom_names'] = found_chroms;\
  \ return output; }"
outputs:
- type:
    items:
      items: string
      type: array
    type: array
  id: bam_arrays_per_chr
- type:
    items: string
    type: array
  id: chrom_names

