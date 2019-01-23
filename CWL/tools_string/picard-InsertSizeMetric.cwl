cwlVersion: v1.0
class: CommandLineTool
doc: "picard-CollectInsertSizeMetrics.cwl is developed for CWL consortium\n Collect\
  \ metrics about the insert size distribution of a paired-end library.\n"
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
  inputFileName_insertSize:
    doc: One or more input SAM or BAM files to analyze. Must be coordinate sorted.
      Default value null. This option may be specified 0 or more times
    type: string
    inputBinding:
      position: 4
      prefix: INPUT=
  outputFileName_insertSize:
    doc: The output file to write marked records to Required
    type: string
    inputBinding:
      position: 5
      prefix: OUTPUT=
  histogramFile:
    doc: File to write duplication metrics to Required
    type: string?
    inputBinding:
      position: 6
      prefix: HISTOGRAM_FILE=
  includeDuplicates:
    doc: If true do not write duplicates to the output file instead of writing them
      with appropriate flags set. Default value false. This option can be set to 'null'
      to clear the default value. Possible values {true, false}
    type: string?
    inputBinding:
      position: 7
      prefix: INCLUDE_DUPLICATES=
      separate: false
  histogramWidth:
    doc: Explicitly sets the Histogram width, overriding automatic truncation of Histogram
      tail. Also, when calculating mean and standard deviation, only bins less or
      equal Histogram_WIDTH will be included. Default value null.
    type: int?
    inputBinding:
      position: 8
      prefix: HISTOGRAM_WIDTH=
  deviations:
    doc: Generate mean, sd and plots by trimming the data down to MEDIAN + DEVIATIONS*MEDIAN_ABSOLUTE_DEVIATION.
      This is done because insert size data typically includes enough anomalous values
      from chimeras and other artifacts to make the mean and sd grossly misleading
      regarding the real distribution. Default value 10.0. This option can be set
      to 'null' to clear the default value.
    type: double?
    inputBinding:
      position: 9
      prefix: DEVIATIONS=
  minPCT:
    doc: When generating the Histogram, discard any data categories (out of FR, TANDEM,
      RF) that have fewer than this percentage of overall reads. (Range 0 to 1). Default
      value 0.05. This option can be set to 'null' to clear the default value.
    type: double?
    inputBinding:
      position: 9
      prefix: MINIMUM_PCT=
  metricAccumulationLevel:
    doc: The level(s) at which to accumulate metrics. Default value ALL_READS. This
      option can be set to 'null' to clear the default value. Possible values ALL_READS,
      SAMPLE, LIBRARY, READ_GROUP This option may be specified 0 or more times. This
      option can be set to 'null' to clear the default list.
    type: string?
    inputBinding:
      position: 10
      prefix: METRIC_ACCUMULATION_LEVEL =
  readSorted:
    doc: If true, assume that the input file is coordinate sorted even if the header
      says otherwise. Default value false. This option can be set to 'null' to clear
      the default value. Possible values {true, false}
    type: boolean?
    inputBinding:
      position: 22
      prefix: ASSUME_SORTED=
  stopAfter:
    doc: Stop after processing N reads, mainly for debugging. Default value 0. This
      option can be set to 'null' to clear the default value.
    type: long?
    inputBinding:
      position: 8
      prefix: STOP_AFTER=
baseCommand: java
arguments:
- position: 2
  valueFrom: /ngs_share/tools/miniconda3/envs/py27/share/picard-2.5.0-1/picard.jar
  prefix: -jar
- position: 3
  valueFrom: CollectInsertSizeMetrics
outputs:
  insertSize_output:
    type: string
    outputBinding:
      outputEval: $(inputs.outputFileName_insertSize)
dct:creator:
- foaf:name: DKFZ
  class: foaf:Organization
  foaf:member:
  - id: p.lutsik@dkfz.de
    class: foaf:Person
    foaf:mbox: mailto:p.lutsik@dkfz.de
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
  doap:description: A set of Java command line tools for manipulating high-throughput
    sequencing data (HTS) data and formats. Picard is implemented using the HTSJDK
    Java library HTSJDK, supporting accessing of common file formats, such as SAM
    and VCF, used for high-throughput sequencing data. http://broadinstitute.github.io/picard/command-line-overview.html#BuildBamIndex
  doap:name: picard
  doap:license: MIT, Apache2
dcat:downloadURL: https://github.com/common-workflow-language/workflows/blob/master/tools/picard-MarkDuplicates.cwl
$schemas:
- http://dublincore.org/2012/06/14/dcterms.rdf
- http://xmlns.com/foaf/spec/20140114.rdf
- http://usefulinc.com/ns/doap#
- http://www.w3.org/ns/adms#
- http://www.w3.org/ns/dcat.rdf
doap:maintainer:
- foaf:name: DKFZ
  class: foaf:Organization
  foaf:member:
  - id: p.lutsik@dkfz.de
    foaf:name: Pavlo Lutsik
    class: foaf:Person
    foaf:mbox: mailto:p.lutsik@dkfz.de
$namespaces:
  foaf: http://xmlns.com/foaf/0.1/
  doap: http://usefulinc.com/ns/doap#
  dct: http://purl.org/dc/terms/
  dcat: http://www.w3.org/ns/dcat#
  adms: http://www.w3.org/ns/adms#
doap:name: picard-InsertSizeMetric.cwl

