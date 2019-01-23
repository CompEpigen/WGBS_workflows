cwlVersion: v1.0
class: CommandLineTool
doc: "picard-MergeSamFiles.cwl is developed for CWL consortium\n"
requirements:
  InlineJavascriptRequirement: {}
hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 10000
inputs:
  java_arg:
    type: string
    default: -Xmx4g
    inputBinding:
      position: 1
  outputFileName_mergedSam:
    doc: "SAM or BAM file to write merged result to Required\n"
    type: string
    inputBinding:
      position: 4
      prefix: OUTPUT=
  inputFileName_mergedSam:
    doc: "SAM or BAM input file Default value null. This option must be specified\
      \ at least 1 times\n"
    type:
      items: File
      inputBinding:
        prefix: INPUT=
      type: array
    inputBinding:
      position: 5
  readSorted:
    doc: "If true, assume that the input files are in the same sort order as the requested\
      \ output sort order, even if their headers say otherwise. Default value false.\
      \ This option can be set to 'null' to clear the default value. Possible values\
      \ {true, false}\n"
    type: boolean?
    inputBinding:
      position: 6
      prefix: ASSUME_SORTED=
  mergeSequenceDictionaries:
    doc: "Merge the sequence dictionaries Default value false. This option can be\
      \ set to null to clear the default value. Possible values {true, false}\n"
    type: boolean?
    inputBinding:
      position: 7
      prefix: MERGE_SEQUENCE_DICTIONARIES=
  useThreading:
    doc: "Option to create a background thread to encode, compress and write to disk\
      \ the output file. The threaded version uses about 20% more CPU and decreases\
      \ runtime by ~20% when writing out a compressed BAM file. Default value false.\
      \ This option can be set to 'null' to clear the default value. Possible values\
      \ {true, false}\n"
    type: boolean?
    inputBinding:
      position: 8
      prefix: USE_THREADING=
  comment:
    doc: "Comment(s) to include in the merged output files header. Default value null.\
      \ This option may be specified 0 or more times\n"
    type: string?
    inputBinding:
      position: 9
      prefix: COMMENT=
  tmpdir:
    doc: "Default value null. This option may be specified 0 or more times.\n"
    type: string
    inputBinding:
      position: 10
      prefix: TMP_DIR=
  createIndex:
    doc: Whether to create a BAM index when writing a coordinate-sorted BAM file.
      Default value false. This option can be set to 'null' to clear the default value.
      Possible values {true, false}
    type: string?
    default: 'true'
    inputBinding:
      position: 20
      prefix: CREATE_INDEX=
      separate: false
baseCommand: java
arguments:
- position: 2
  valueFrom: /ngs_share/tools/miniconda3/envs/py27/share/picard-2.5.0-1/picard.jar
  prefix: -jar
- position: 3
  valueFrom: MergeSamFiles
outputs:
  mergeSam_output:
    type: File
    outputBinding:
      glob: $(inputs.outputFileName_mergedSam)
dct:creator:
- foaf:name: THE UNIVERSITY OF MELBOURNE
  class: foaf:Organization
  foaf:member:
  - id: farahk@student.unimelb.edu.au
    foaf:name: Farah Zaib Khan
    class: foaf:Person
    foaf:mbox: mailto:farahk@student.unimelb.edu.au
  - id: skanwal@student.unimelb.edu.au
    foaf:name: Sehrish Kanwal
    class: foaf:Person
    foaf:mbox: mailto:skanwal@student.unimelb.edu.au
adms:includedAsset:
  doap:homepage: http://broadinstitute.github.io/picard/
  doap:programming-language: JAVA
  doap:category: commandline tool
  doap:release:
  - doap:revision: '1.141'
    class: doap:Version
  doap:repository:
  - doap:location: https://github.com/broadinstitute/picard.git
    class: doap:GitRepository
  doap:developer:
  - foaf:name: Broad Institute
    class: foaf:Organization
  doap:description: "A set of Java command line tools for manipulating high-throughput\
    \ sequencing data (HTS) data and formats. Picard is implemented using the HTSJDK\
    \ Java library HTSJDK, supporting accessing of common file formats, such as SAM\
    \ and VCF, used for high-throughput sequencing data. http://broadinstitute.github.io/picard/command-line-overview.html#MergeSamFiles\n"
  doap:name: picard
  doap:license: MIT, Apache2
dcat:downloadURL: https://github.com/common-workflow-language/workflows/blob/master/tools/picard-MergeSamFiles.cwl
$schemas:
- http://dublincore.org/2012/06/14/dcterms.rdf
- http://xmlns.com/foaf/spec/20140114.rdf
- http://usefulinc.com/ns/doap#
- http://www.w3.org/ns/adms#
- http://www.w3.org/ns/dcat.rdf
doap:maintainer:
- foaf:name: Barski Lab, Cincinnati Children's Hospital Medical Center
  class: foaf:Organization
  foaf:member:
  - id: http://orcid.org/0000-0001-9102-5681
    foaf:name: Andrey Kartashov
    class: foaf:Person
    foaf:openid: http://orcid.org/0000-0001-9102-5681
    foaf:mbox: mailto:Andrey.Kartashov@cchmc.org
$namespaces:
  foaf: http://xmlns.com/foaf/0.1/
  doap: http://usefulinc.com/ns/doap#
  dct: http://purl.org/dc/terms/
  dcat: http://www.w3.org/ns/dcat#
  adms: http://www.w3.org/ns/adms#
doap:name: picard-MergeSamFiles.cwl

