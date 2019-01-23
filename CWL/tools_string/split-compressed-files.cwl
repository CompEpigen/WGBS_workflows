cwlVersion: v1.0
class: CommandLineTool
requirements:
  ShellCommandRequirement: {}
hints:
  ResourceRequirement:
    coresMin: 1
    diskMin: 2500
    ramMin: 2500
inputs:
  file_dir:
    type: string
    inputBinding:
      position: 1
  file:
    doc: "the input fastq file\n"
    type: string
    inputBinding:
      position: 1002
  size:
    doc: "size of a chunck in lines\n"
    type: int
    inputBinding:
      position: 10004
      prefix: --lines=
      valueFrom: $(self * 4)
      separate: false
  max_reads:
    doc: "take at most max_reads reads\n"
    type: int
  suffix:
    doc: "a suffix for chuncks\n"
    type: string
    inputBinding:
      position: 10003
      prefix: --additional-suffix=
      separate: false
baseCommand: cd
arguments:
- position: 1000
  valueFrom: ;
- position: 1001
  valueFrom: zcat
- shellQuote: false
  position: 2000
  valueFrom: "${\n     if(inputs.max_reads>0){\n           return \"|\";\n     }else{\n\
    \           return \"\";\n     }\n }\n"
- position: 2001
  valueFrom: "${\n     if(inputs.max_reads>0){\n           return \"head\";\n    \
    \ }else{\n           return \"\";\n     }\n }\n"
- shellQuote: false
  position: 2002
  valueFrom: "${\n     if(inputs.max_reads>0){\n           return \"-n\";\n     }else{\n\
    \           return \"\";\n     }\n }\n"
- shellQuote: false
  position: 2003
  valueFrom: "${\n     if(inputs.max_reads>0){\n           return inputs.max_reads*4;\n\
    \     }else{\n           return \"\";\n     }\n }\n"
- shellQuote: false
  position: 10000
  valueFrom: '|'
- position: 10001
  valueFrom: split
- shellQuote: false
  position: 10002
  valueFrom: '-'
- shellQuote: true
  position: 100000
  valueFrom: ;
- shellQuote: false
  position: 100001
  valueFrom: ls -1
- shellQuote: false
  position: 100002
  valueFrom: '|'
- shellQuote: false
  position: 100003
  valueFrom: grep
- shellQuote: false
  position: 100004
  valueFrom: $(inputs.suffix)
stdout: list_of_files
outputs:
  output:
    type:
      items: string
      type: array
    outputBinding:
      glob: list_of_files
      loadContents: true
      outputEval: $(self[0].contents.split("\n").filter(function(n){ return n != ""
        }))

