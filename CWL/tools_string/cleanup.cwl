cwlVersion: v1.0
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}
inputs:
  files:
    doc: "the input files to remove\n"
    type:
      items: string
      type: array
  data_dir:
    doc: "the working directory with files\n"
    type: string
  do_clean:
    doc: "action switch\n"
    type: boolean
baseCommand: cd
arguments:
- position: 1
  valueFrom: $(inputs.data_dir)
- shellQuote: false
  position: 2
  valueFrom: ;
- position: 3
  valueFrom: "${ \n    if(inputs.do_clean){\n       return \"rm\";\n    }else{\n \
    \      return \"ls\";\n    }\n}\n"
- position: 100000
  valueFrom: $(inputs.files)
outputs:
  output:
    type:
      items: string
      type: array
    outputBinding:
      outputEval: $(inputs.files)

