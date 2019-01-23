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
  id: bam_arrays
expression: "${\nvar chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',\
  \ 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16',\
  \ 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']; var rearranged\
  \ = [];\nfor(var ci=0; ci<chroms.length; ci++){ var chr_output = []; for(var i=0;\
  \ i<inputs.bam_arrays.length; i++){ for(var j=0; j<inputs.bam_arrays[i].length;\
  \ j++){ var this_file = inputs.bam_arrays[i][j]; if(this_file.basename.indexOf('.'\
  \ + chroms[ci] + '.') !== -1){ chr_output[chr_output.length] = this_file; } } }\
  \ rearranged[ci] = chr_output; } var output = {}; output['bam_arrays_per_chr'] =\
  \ rearranged; output['chrom_names'] = chroms; return output; }"
outputs:
- type:
    items:
      items: File
      type: array
    type: array
  id: bam_arrays_per_chr
- type:
    items: string
    type: array
  id: chrom_names

