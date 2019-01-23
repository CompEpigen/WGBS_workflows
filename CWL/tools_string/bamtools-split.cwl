cwlVersion: v1.0
class: CommandLineTool
label: BamToolsSplit
doc: "BamTools Split splits a BAM file based on a user-specified property. It creates\
  \ a new BAM output file for each value found.\n\n**Warning:** Splitting  by tags\
  \ or reference can output a large number of files.\n\n**Common issues:** Splitting\
  \ by tag can produce no output if the selected tag doesn't exist in the BAM file."
requirements:
  ShellCommandRequirement: {}
inputs:
  file_dir: string
  input_bam_file:
    label: Input BAM file
    doc: The input BAM file.
    type: string
    inputBinding:
      position: 2
      prefix: -in
      separate: true
    streamable: false
  ref_prefix:
    label: Reference prefix
    doc: Custom prefix for splitting by references. Currently files end with REF_<refName>.bam.
      This option allows you to replace "REF_" with a prefix of your choosing.
    type: string?
    inputBinding:
      position: 4
      prefix: -refPrefix
    streamable: false
  tag_prefix:
    label: Tag prefix
    doc: Custom prefix for splitting by tags. Current files end with TAG_<tagname>_<tagvalue>.bam.
      This option allows you to replace "TAG_" with a prefix of your choosing.
    type: string?
    inputBinding:
      position: 4
      prefix: ''
    streamable: false
  tag:
    label: Tag split
    doc: Splits alignments based on all values of TAG encountered (i.e. -tag RG creates
      a BAM file for each read group in original BAM file).
    type: string?
    streamable: false
  split_options:
    label: Split Options
    doc: 'Property upon which the BAM splitting is performed. Warning: Splitting  by
      tags or reference can output a large number of files.'
    type:
      type: enum
      name: split_options
      symbols:
      - mapped
      - paired
      - reference
      - tag
    streamable: false
stdin: ''
baseCommand: cd
arguments:
- position: -4
  valueFrom: $(inputs.file_dir)
- position: -3
  valueFrom: ;
- position: -2
  valueFrom: bamtools
- position: -1
  valueFrom: split
- position: 3
  valueFrom: "${\n   var filepath = inputs.input_bam_file;\n   var file_path_sep =\
    \ filepath.split(\"/\");\n   var filename = file_path_sep[file_path_sep.length-1];\n\
    \   var file_dot_sep = filename.split(\".\");\n   var base_name = file_dot_sep.slice(0,-1).join(\"\
    .\");\n   return base_name;\n}\n"
  prefix: -stub
- position: 4
  valueFrom: "${\n  var line = \"\";\n  \n  if (inputs.split_options == 'mapped'){\n\
    \    line = line.concat(\"-mapped\");\n  } else if (inputs.split_options == 'paired'){\n\
    \    line = line.concat(\"-paired\");\n  } else if (inputs.split_options == 'reference'){\n\
    \    line = line.concat(\"-reference\");\n  } else if (inputs.split_options ==\
    \ 'tag'){\n      line = line.concat(\"-tag \");\n      line = line.concat(inputs.tag)\n\
    \      if (inputs.tag_prefix){\n          line = line.concat(\" -tagPrefix \"\
    );\n          line = line.concat(inputs.tag_prefix);\n      }\n  }\n  return line;\n\
    }\n"
  separate: true
- position: 100
  valueFrom: ;
- position: 1000
  valueFrom: ls
- position: 1001
  valueFrom: '-1'
- position: 1002
  valueFrom: '|'
- position: 1003
  valueFrom: grep
- position: 1004
  valueFrom: ${ var filepath = inputs.input_bam_file; var file_path_sep = filepath.split("/");
    var filename = file_path_sep[file_path_sep.length-1]; var file_dot_sep = filename.split(".bam");
    var base_name = file_dot_sep.slice(0,1); return base_name + "." + inputs.ref_prefix
    +".*.bam"; }
stdout: split_chr_files
outputs:
  output_bam_files:
    label: Output BAM files
    doc: Output BAM files.
    type: string[]
    outputBinding:
      glob: split_chr_files
      loadContents: true
      outputEval: $(self[0].contents.split("\n").filter(function(n){ return n != ""
        }))
    streamable: false
successCodes: []
temporaryFailCodes: []

