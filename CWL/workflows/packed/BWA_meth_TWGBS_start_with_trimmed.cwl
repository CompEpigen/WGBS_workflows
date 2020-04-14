{
    "$graph": [
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "hints": [
                {
                    "dockerPull": "kerstenbreuer/samtools:1.7",
                    "class": "DockerRequirement"
                },
                {
                    "coresMin": 1,
                    "ramMin": 4000,
                    "tmpdirMin": 10000,
                    "class": "ResourceRequirement"
                }
            ],
            "baseCommand": "awk",
            "arguments": [
                {
                    "shellQuote": false,
                    "position": 1,
                    "valueFrom": "${ return \"{SUM1 += $4/100; SUM2 +=1} END {print 1-SUM1/SUM2}\" }"
                }
            ],
            "stdout": "$(inputs.output_basename + \".bs_convers.txt\")",
            "inputs": [
                {
                    "type": "string",
                    "id": "#bisulfite_conversion_spike_in.cwl/output_basename"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "position": 10
                    },
                    "id": "#bisulfite_conversion_spike_in.cwl/spike_in_mcall_bedgraph"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.output_basename + \".bs_convers.txt\")"
                    },
                    "id": "#bisulfite_conversion_spike_in.cwl/bisulfite_conversion_file"
                }
            ],
            "id": "#bisulfite_conversion_spike_in.cwl",
            "$namespaces": {
                "sbg": "https://www.sevenbridges.com/"
            }
        },
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "ShellCommandRequirement"
                }
            ],
            "hints": [
                {
                    "dockerPull": "kerstenbreuer/bwameth:latest",
                    "class": "DockerRequirement"
                },
                {
                    "coresMin": "$( inputs.threads )",
                    "ramMin": 28000,
                    "tmpdirMin": 30000,
                    "class": "ResourceRequirement"
                }
            ],
            "baseCommand": [
                "bwameth.py"
            ],
            "stdout": "$(inputs.fastq1.nameroot).sam",
            "inputs": [
                {
                    "doc": "the input fastq file with the first mate",
                    "type": "File",
                    "inputBinding": {
                        "position": 11
                    },
                    "id": "#bwameth.cwl/fastq1"
                },
                {
                    "doc": "the input fastq file with the second mate",
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 12
                    },
                    "id": "#bwameth.cwl/fastq2"
                },
                {
                    "doc": "Is library type non-directional",
                    "type": "boolean",
                    "default": false,
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--non-directional"
                    },
                    "id": "#bwameth.cwl/is_non_directional"
                },
                {
                    "doc": "the reference fasta file location",
                    "type": "File",
                    "secondaryFiles": [
                        ".fai",
                        ".bwameth.c2t",
                        ".bwameth.c2t.amb",
                        ".bwameth.c2t.ann",
                        ".bwameth.c2t.bwt",
                        ".bwameth.c2t.pac",
                        ".bwameth.c2t.sa"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--reference",
                        "separate": true
                    },
                    "id": "#bwameth.cwl/reference"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "default": 1,
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--threads"
                    },
                    "id": "#bwameth.cwl/threads"
                }
            ],
            "outputs": [
                {
                    "type": "stdout",
                    "id": "#bwameth.cwl/sam"
                }
            ],
            "id": "#bwameth.cwl"
        },
        {
            "class": "CommandLineTool",
            "hints": [
                {
                    "dockerPull": "kerstenbreuer/trim_galore:0.6.4_2.6_0.11.8",
                    "class": "DockerRequirement"
                },
                {
                    "coresMin": 1,
                    "ramMin": 5000,
                    "tmpdirMin": 10000,
                    "class": "ResourceRequirement"
                }
            ],
            "baseCommand": "fastqc",
            "arguments": [
                {
                    "valueFrom": "$(runtime.outdir)",
                    "prefix": "-o"
                },
                {
                    "valueFrom": "--noextract"
                }
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 1
                    },
                    "id": "#fastqc.cwl/bam"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 1
                    },
                    "id": "#fastqc.cwl/read1"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 2
                    },
                    "id": "#fastqc.cwl/read2"
                }
            ],
            "outputs": [
                {
                    "doc": "html report showing results from zip",
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputBinding": {
                        "glob": "*_fastqc.html"
                    },
                    "id": "#fastqc.cwl/fastqc_html"
                },
                {
                    "doc": "all data e.g. figures",
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputBinding": {
                        "glob": "*_fastqc.zip"
                    },
                    "id": "#fastqc.cwl/fastqc_zip"
                }
            ],
            "id": "#fastqc.cwl"
        },
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "ShellCommandRequirement"
                }
            ],
            "hints": [
                {
                    "dockerPull": "nfcore/methylseq",
                    "class": "DockerRequirement"
                },
                {
                    "coresMin": "$( inputs.threads )",
                    "ramMin": 28000,
                    "tmpdirMin": 30000,
                    "class": "ResourceRequirement"
                }
            ],
            "baseCommand": [
                "MethylDackel",
                "extract"
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--CTOB"
                    },
                    "id": "#methyldackel_extract.cwl/CTOB"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--CTOT"
                    },
                    "id": "#methyldackel_extract.cwl/CTOT"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--OB"
                    },
                    "id": "#methyldackel_extract.cwl/OB"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--OT"
                    },
                    "id": "#methyldackel_extract.cwl/OT"
                },
                {
                    "type": "File",
                    "secondaryFiles": ".bai",
                    "inputBinding": {
                        "position": 12
                    },
                    "id": "#methyldackel_extract.cwl/bam"
                },
                {
                    "type": "int",
                    "default": 1,
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--minDepth"
                    },
                    "id": "#methyldackel_extract.cwl/min_depth"
                },
                {
                    "type": "int",
                    "default": 0,
                    "inputBinding": {
                        "position": 1,
                        "prefix": "-q"
                    },
                    "id": "#methyldackel_extract.cwl/min_mapq"
                },
                {
                    "type": "int",
                    "default": 1,
                    "inputBinding": {
                        "position": 1,
                        "prefix": "-p"
                    },
                    "id": "#methyldackel_extract.cwl/min_phred"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--nCTOB"
                    },
                    "id": "#methyldackel_extract.cwl/nCTOB"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--nCTOT"
                    },
                    "id": "#methyldackel_extract.cwl/nCTOT"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--nOB"
                    },
                    "id": "#methyldackel_extract.cwl/nOB"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--nOT"
                    },
                    "id": "#methyldackel_extract.cwl/nOT"
                },
                {
                    "type": "boolean",
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--noCpG"
                    },
                    "id": "#methyldackel_extract.cwl/noCG"
                },
                {
                    "type": "string",
                    "inputBinding": {
                        "valueFrom": "$(self + \".mcall\")",
                        "position": 10,
                        "prefix": "-o"
                    },
                    "id": "#methyldackel_extract.cwl/output_basename"
                },
                {
                    "type": "File",
                    "secondaryFiles": [
                        ".fai",
                        ".bwameth.c2t",
                        ".bwameth.c2t.amb",
                        ".bwameth.c2t.ann",
                        ".bwameth.c2t.bwt",
                        ".bwameth.c2t.pac",
                        ".bwameth.c2t.sa"
                    ],
                    "inputBinding": {
                        "position": 11
                    },
                    "id": "#methyldackel_extract.cwl/reference"
                },
                {
                    "type": "int",
                    "default": 1,
                    "inputBinding": {
                        "position": 1,
                        "prefix": "-@"
                    },
                    "id": "#methyldackel_extract.cwl/threads"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.bedGraph"
                    },
                    "id": "#methyldackel_extract.cwl/mcall_bedgraph"
                }
            ],
            "id": "#methyldackel_extract.cwl"
        },
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "ShellCommandRequirement"
                }
            ],
            "hints": [
                {
                    "dockerPull": "nfcore/methylseq",
                    "class": "DockerRequirement"
                },
                {
                    "coresMin": "$( inputs.threads )",
                    "ramMin": 28000,
                    "tmpdirMin": 30000,
                    "class": "ResourceRequirement"
                }
            ],
            "baseCommand": [
                "MethylDackel",
                "mbias"
            ],
            "arguments": [
                {
                    "valueFrom": "--txt",
                    "position": 1
                }
            ],
            "stdout": "$(inputs.output_basename + \".mbias.txt\")",
            "inputs": [
                {
                    "type": "File",
                    "secondaryFiles": ".bai",
                    "inputBinding": {
                        "position": 12
                    },
                    "id": "#methyldackel_mbias.cwl/bam"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--nCTOB"
                    },
                    "id": "#methyldackel_mbias.cwl/nCTOB"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--nCTOT"
                    },
                    "id": "#methyldackel_mbias.cwl/nCTOT"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--nOB"
                    },
                    "id": "#methyldackel_mbias.cwl/nOB"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--nOT"
                    },
                    "id": "#methyldackel_mbias.cwl/nOT"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--noCpG"
                    },
                    "id": "#methyldackel_mbias.cwl/noCG"
                },
                {
                    "type": "string",
                    "inputBinding": {
                        "valueFrom": "$(self + \".mbias_plot\")",
                        "position": 13
                    },
                    "id": "#methyldackel_mbias.cwl/output_basename"
                },
                {
                    "type": "File",
                    "secondaryFiles": [
                        ".fai",
                        ".bwameth.c2t",
                        ".bwameth.c2t.amb",
                        ".bwameth.c2t.ann",
                        ".bwameth.c2t.bwt",
                        ".bwameth.c2t.pac",
                        ".bwameth.c2t.sa"
                    ],
                    "inputBinding": {
                        "position": 11
                    },
                    "id": "#methyldackel_mbias.cwl/reference"
                },
                {
                    "type": "int",
                    "default": 1,
                    "inputBinding": {
                        "position": 1,
                        "prefix": "-@"
                    },
                    "id": "#methyldackel_mbias.cwl/threads"
                }
            ],
            "outputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputBinding": {
                        "glob": "*$(inputs.output_basename)*"
                    },
                    "id": "#methyldackel_mbias.cwl/mbias_output"
                }
            ],
            "id": "#methyldackel_mbias.cwl"
        },
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                }
            ],
            "hints": [
                {
                    "dockerPull": "kerstenbreuer/multiqc:1.7",
                    "class": "DockerRequirement"
                },
                {
                    "coresMin": 1,
                    "ramMin": 10000,
                    "tmpdirMin": 10000,
                    "class": "ResourceRequirement"
                }
            ],
            "baseCommand": [
                "bash",
                "-c"
            ],
            "arguments": [
                {
                    "valueFrom": "${\n    var qc_files_array = inputs.qc_files_array;\n    var qc_files_array_of_array = inputs.qc_files_array_of_array;\n    var cmdline = \"echo 'copying input file ...'\";\n\n    if ( qc_files_array != null ){\n      for (var i=0; i<qc_files_array.length; i++){\n        if( qc_files_array[i] != null ){\n          cmdline += \"; cp \" + qc_files_array[i].path + \" .\";\n        }\n      }\n    }\n\n    if ( qc_files_array_of_array != null ){\n      for (var i=0; i<qc_files_array_of_array.length; i++){ \n        for (var ii=0; ii<qc_files_array_of_array[i].length; ii++){\n          if( qc_files_array_of_array[i][ii] != null ){\n            cmdline += \"; cp \" + qc_files_array_of_array[i][ii].path + \" .\";\n          }\n        }\n      }\n    }\n    \n    cmdline += \"; echo \\'copying done\\'\" +\n        \"; multiqc --zip-data-dir --cl_config \\'log_filesize_limit: 100000000\\' \" +\n        \"--outdir \" + runtime.outdir +\n        \" --filename \" + inputs.report_name + \"_report .\";\n\n    return cmdline\n  }\n"
                }
            ],
            "inputs": [
                {
                    "doc": "qc files which shall be part of the multiqc summary;\noptional, only one of qc_files_array or qc_files_array_of_array \nmust be provided\n",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": [
                                "File",
                                "null"
                            ]
                        }
                    ],
                    "id": "#multiqc_hack.cwl/qc_files_array"
                },
                {
                    "doc": "qc files which shall be part of the multiqc summary;\noptional, only one of qc_files_array or qc_files_array_of_array \nmust be provided\n",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": {
                                "type": "array",
                                "items": [
                                    "File",
                                    "null"
                                ]
                            }
                        }
                    ],
                    "id": "#multiqc_hack.cwl/qc_files_array_of_array"
                },
                {
                    "doc": "name used for the html report and the corresponding zip file",
                    "type": "string",
                    "default": "multiqc",
                    "id": "#multiqc_hack.cwl/report_name"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.report_name + \"_report.html\")"
                    },
                    "id": "#multiqc_hack.cwl/multiqc_html"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.report_name + \"_report_data.zip\")"
                    },
                    "id": "#multiqc_hack.cwl/multiqc_zip"
                }
            ],
            "id": "#multiqc_hack.cwl"
        },
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                }
            ],
            "hints": [
                {
                    "dockerPull": "kerstenbreuer/picard_tools:2.17.4",
                    "class": "DockerRequirement"
                },
                {
                    "coresMin": 1,
                    "ramMin": 20000,
                    "tmpdirMin": 10000,
                    "class": "ResourceRequirement"
                }
            ],
            "baseCommand": [
                "java",
                "-jar"
            ],
            "arguments": [
                {
                    "valueFrom": "MarkDuplicates",
                    "position": 2
                },
                {
                    "valueFrom": "$(inputs.bam_sorted.nameroot + \"_duprem.bam\")",
                    "prefix": "OUTPUT=",
                    "separate": false,
                    "position": 13
                },
                {
                    "valueFrom": "$(inputs.bam_sorted.nameroot + \"_duprem.log\")",
                    "prefix": "METRICS_FILE=",
                    "separate": false,
                    "position": 13
                },
                {
                    "valueFrom": "REMOVE_DUPLICATES=TRUE",
                    "position": 14
                },
                {
                    "valueFrom": "ASSUME_SORTED=TRUE",
                    "position": 15
                },
                {
                    "valueFrom": "VALIDATION_STRINGENCY=SILENT",
                    "position": 16
                }
            ],
            "stdout": "$(inputs.bam_sorted.nameroot + \".picard_markdup.stdout\")",
            "inputs": [
                {
                    "doc": "sorted bam input file",
                    "type": "File",
                    "inputBinding": {
                        "prefix": "INPUT=",
                        "separate": false,
                        "position": 11
                    },
                    "id": "#picard_markdup.cwl/bam_sorted"
                },
                {
                    "type": "string",
                    "default": "/bin/picard.jar",
                    "inputBinding": {
                        "position": 1
                    },
                    "id": "#picard_markdup.cwl/path_to_picards"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.bam_sorted.nameroot + \"_duprem.bam\")"
                    },
                    "id": "#picard_markdup.cwl/bam_duprem"
                },
                {
                    "type": "stdout",
                    "id": "#picard_markdup.cwl/picard_markdup_log"
                }
            ],
            "id": "#picard_markdup.cwl"
        },
        {
            "class": "CommandLineTool",
            "hints": [
                {
                    "dockerPull": "kerstenbreuer/samtools:1.7",
                    "class": "DockerRequirement"
                },
                {
                    "coresMin": 4,
                    "ramMin": 15000,
                    "tmpdirMin": 10000,
                    "class": "ResourceRequirement"
                }
            ],
            "baseCommand": [
                "samtools",
                "flagstat"
            ],
            "stdout": "$(inputs.bam.nameroot + inputs.output_suffix)",
            "inputs": [
                {
                    "type": "File",
                    "inputBinding": {
                        "position": 2
                    },
                    "id": "#samtools_flagstat.cwl/bam"
                },
                {
                    "type": "string",
                    "default": ".flagStat",
                    "id": "#samtools_flagstat.cwl/output_suffix"
                }
            ],
            "outputs": [
                {
                    "type": "stdout",
                    "id": "#samtools_flagstat.cwl/flagstat_output"
                }
            ],
            "id": "#samtools_flagstat.cwl"
        },
        {
            "class": "CommandLineTool",
            "hints": [
                {
                    "dockerPull": "kerstenbreuer/samtools:1.7",
                    "class": "DockerRequirement"
                },
                {
                    "coresMin": 1,
                    "ramMin": 20000,
                    "tmpdirMin": 10000,
                    "class": "ResourceRequirement"
                }
            ],
            "baseCommand": [
                "bash",
                "-c"
            ],
            "arguments": [
                {
                    "valueFrom": "$(\"cp \" + inputs.bam_sorted.path + \" . && samtools index -b \" + inputs.bam_sorted.basename )"
                }
            ],
            "inputs": [
                {
                    "doc": "sorted bam input file",
                    "type": "File",
                    "id": "#samtools_index_hack.cwl/bam_sorted"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "secondaryFiles": ".bai",
                    "outputBinding": {
                        "glob": "$(inputs.bam_sorted.basename)"
                    },
                    "id": "#samtools_index_hack.cwl/bam_sorted_indexed"
                }
            ],
            "id": "#samtools_index_hack.cwl"
        },
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                }
            ],
            "hints": [
                {
                    "dockerPull": "kerstenbreuer/samtools:1.7",
                    "class": "DockerRequirement"
                },
                {
                    "coresMin": 1,
                    "ramMin": 20000,
                    "tmpdirMin": 30000,
                    "class": "ResourceRequirement"
                }
            ],
            "baseCommand": [
                "samtools",
                "merge"
            ],
            "inputs": [
                {
                    "id": "#samtools_merge.cwl/output_name",
                    "doc": "name of merged bam file",
                    "type": "string",
                    "inputBinding": {
                        "position": 1
                    }
                },
                {
                    "id": "#samtools_merge.cwl/bams",
                    "doc": "bam files to be merged",
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "inputBinding": {
                        "position": 2
                    }
                }
            ],
            "outputs": [
                {
                    "id": "#samtools_merge.cwl/bam_merged",
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.output_name)"
                    }
                }
            ],
            "id": "#samtools_merge.cwl"
        },
        {
            "doc": "Sort a bam file by read names.",
            "class": "CommandLineTool",
            "hints": [
                {
                    "dockerPull": "kerstenbreuer/samtools:1.7",
                    "class": "DockerRequirement"
                },
                {
                    "coresMin": 4,
                    "ramMin": 15000,
                    "tmpdirMin": 30000,
                    "class": "ResourceRequirement"
                }
            ],
            "baseCommand": [
                "samtools",
                "sort"
            ],
            "arguments": [
                {
                    "valueFrom": "$(runtime.cores)",
                    "prefix": "-@"
                }
            ],
            "inputs": [
                {
                    "doc": "aligned reads to be checked in sam or bam format",
                    "type": "File",
                    "inputBinding": {
                        "position": 2
                    },
                    "id": "#samtools_sort.cwl/bam_unsorted"
                }
            ],
            "stdout": "$(inputs.bam_unsorted.basename)",
            "outputs": [
                {
                    "type": "stdout",
                    "id": "#samtools_sort.cwl/bam_sorted"
                }
            ],
            "id": "#samtools_sort.cwl"
        },
        {
            "class": "CommandLineTool",
            "hints": [
                {
                    "dockerPull": "kerstenbreuer/samtools:1.7",
                    "class": "DockerRequirement"
                },
                {
                    "coresMin": 1,
                    "ramMin": 10000,
                    "tmpdirMin": 10000,
                    "class": "ResourceRequirement"
                }
            ],
            "baseCommand": [
                "samtools",
                "view"
            ],
            "inputs": [
                {
                    "doc": "aligned reads in sam or bam format",
                    "type": "File",
                    "secondaryFiles": ".bai",
                    "inputBinding": {
                        "position": 2
                    },
                    "id": "#samtools_view_extract_spike_in.cwl/bam"
                },
                {
                    "doc": "region to extract",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 3
                    },
                    "id": "#samtools_view_extract_spike_in.cwl/region"
                }
            ],
            "arguments": [
                {
                    "valueFrom": "-h",
                    "position": 1
                },
                {
                    "valueFrom": "-b",
                    "position": 1
                }
            ],
            "stdout": "$( inputs.bam.nameroot + \"_spike_in.bam\" )",
            "outputs": [
                {
                    "type": "stdout",
                    "id": "#samtools_view_extract_spike_in.cwl/bam_spike_in"
                }
            ],
            "id": "#samtools_view_extract_spike_in.cwl"
        },
        {
            "class": "Workflow",
            "inputs": [
                {
                    "id": "#main/fastq1",
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "https://www.sevenbridges.com/x": -82.05426025390625,
                    "https://www.sevenbridges.com/y": -708.0413818359375
                },
                {
                    "id": "#main/fastq2",
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "https://www.sevenbridges.com/x": -86.61282348632812,
                    "https://www.sevenbridges.com/y": -541.4498291015625
                },
                {
                    "id": "#main/is_non_directional",
                    "type": "boolean",
                    "default": false,
                    "https://www.sevenbridges.com/x": -86.61282348632812,
                    "https://www.sevenbridges.com/y": -365.7083740234375
                },
                {
                    "id": "#main/max_threads",
                    "type": "int",
                    "default": 10,
                    "https://www.sevenbridges.com/x": -78.44940948486328,
                    "https://www.sevenbridges.com/y": 89.95835876464844
                },
                {
                    "id": "#main/methyldackel_ctob",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2280.761474609375
                },
                {
                    "id": "#main/methyldackel_ctot",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": -18.234277725219727,
                    "https://www.sevenbridges.com/y": 2155.620849609375
                },
                {
                    "id": "#main/methyldackel_min_depth",
                    "type": "int",
                    "default": 1,
                    "https://www.sevenbridges.com/x": 7.168551405811741e-07,
                    "https://www.sevenbridges.com/y": 2016.8048095703125
                },
                {
                    "id": "#main/methyldackel_min_mapq",
                    "type": "int",
                    "default": 0,
                    "https://www.sevenbridges.com/x": -5.7605383574355074e-08,
                    "https://www.sevenbridges.com/y": 1864.312744140625
                },
                {
                    "id": "#main/methyldackel_min_phred",
                    "type": "int",
                    "default": 0,
                    "https://www.sevenbridges.com/x": 9.11713981628418,
                    "https://www.sevenbridges.com/y": 1693.5865478515625
                },
                {
                    "id": "#main/methyldackel_nctob",
                    "type": [
                        "null",
                        "string"
                    ],
                    "default": "0,0,0,0",
                    "https://www.sevenbridges.com/x": 1.0176772775594145e-06,
                    "https://www.sevenbridges.com/y": 1554.770263671875
                },
                {
                    "id": "#main/methyldackel_nctot",
                    "type": [
                        "null",
                        "string"
                    ],
                    "default": "0,0,0,0",
                    "https://www.sevenbridges.com/x": -4.558569431304932,
                    "https://www.sevenbridges.com/y": 1415.9539794921875
                },
                {
                    "id": "#main/methyldackel_noCG",
                    "type": "boolean",
                    "default": false,
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1124.6458740234375
                },
                {
                    "id": "#main/methyldackel_nob",
                    "type": [
                        "null",
                        "string"
                    ],
                    "default": "2,9,2,12",
                    "https://www.sevenbridges.com/x": 4.558570384979248,
                    "https://www.sevenbridges.com/y": 1254.344970703125
                },
                {
                    "id": "#main/methyldackel_not",
                    "type": [
                        "null",
                        "string"
                    ],
                    "default": "2,9,2,12",
                    "https://www.sevenbridges.com/x": -9.11713981628418,
                    "https://www.sevenbridges.com/y": 985.8296508789062
                },
                {
                    "id": "#main/methyldackel_ob",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 847.013427734375
                },
                {
                    "id": "#main/methyldackel_ot",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": -9.681075096130371,
                    "https://www.sevenbridges.com/y": 702.4104614257812
                },
                {
                    "id": "#main/reference",
                    "type": "File",
                    "secondaryFiles": [
                        ".fai",
                        ".bwameth.c2t",
                        ".bwameth.c2t.amb",
                        ".bwameth.c2t.ann",
                        ".bwameth.c2t.bwt",
                        ".bwameth.c2t.pac",
                        ".bwameth.c2t.sa"
                    ],
                    "https://www.sevenbridges.com/x": -94.5317153930664,
                    "https://www.sevenbridges.com/y": -168.46800231933594
                },
                {
                    "id": "#main/sample_id",
                    "type": "string",
                    "https://www.sevenbridges.com/x": -73.00712585449219,
                    "https://www.sevenbridges.com/y": 248.1683807373047
                },
                {
                    "id": "#main/spikein_chr_name",
                    "type": "string",
                    "doc": "chromosome name in the reference genome which corresponds to the unmethylated spikein DNA",
                    "default": "Lamda",
                    "https://www.sevenbridges.com/x": -136.75709533691406,
                    "https://www.sevenbridges.com/y": -1331.1024169921875
                }
            ],
            "outputs": [
                {
                    "id": "#main/bam",
                    "outputSource": [
                        "#main/index_bam/bam_sorted_indexed"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 4508,
                    "https://www.sevenbridges.com/y": 658
                },
                {
                    "id": "#main/bam_spike_in",
                    "outputSource": [
                        "#main/index_spike_in_bam/bam_sorted_indexed"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 4755.5654296875,
                    "https://www.sevenbridges.com/y": 1112.589111328125
                },
                {
                    "id": "#main/bisulfite_conversion_file",
                    "outputSource": [
                        "#main/conversion_estimation_spike_in/bisulfite_conversion_file"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 4739.7578125,
                    "https://www.sevenbridges.com/y": 936.0545654296875
                },
                {
                    "id": "#main/mbias_output",
                    "outputSource": [
                        "#main/mbias_calculation/mbias_output"
                    ],
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "https://www.sevenbridges.com/x": 4634.173828125,
                    "https://www.sevenbridges.com/y": 1724.246826171875
                },
                {
                    "id": "#main/mbias_output_wo_trimming",
                    "outputSource": [
                        "#main/mbias_calculation_wo_trimming/mbias_output"
                    ],
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "https://www.sevenbridges.com/x": 4529.32666015625,
                    "https://www.sevenbridges.com/y": -46.497440338134766
                },
                {
                    "id": "#main/mcall_bedgraph",
                    "outputSource": [
                        "#main/methylation_calling/mcall_bedgraph"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 4697.99365234375,
                    "https://www.sevenbridges.com/y": 2162.25732421875
                },
                {
                    "id": "#main/mcall_bedgraph_spike_in",
                    "outputSource": [
                        "#main/methylation_calling_spike_in/mcall_bedgraph"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 4748.85400390625,
                    "https://www.sevenbridges.com/y": 1281.314208984375
                },
                {
                    "id": "#main/multiqc_html",
                    "outputSource": [
                        "#main/create_summary_qc_report/multiqc_html"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 4540.59521484375,
                    "https://www.sevenbridges.com/y": 492.2149353027344
                },
                {
                    "id": "#main/multiqc_zip",
                    "outputSource": [
                        "#main/create_summary_qc_report/multiqc_zip"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 4572.505859375,
                    "https://www.sevenbridges.com/y": 285.0201416015625
                },
                {
                    "id": "#main/picard_markdup_log",
                    "outputSource": [
                        "#main/remove_duplicates/picard_markdup_log"
                    ],
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "https://www.sevenbridges.com/x": 2336.66064453125,
                    "https://www.sevenbridges.com/y": -639.9028930664062
                },
                {
                    "id": "#main/post_mapping_fastqc_html",
                    "outputSource": [
                        "#main/qc_post_mapping/fastqc_html"
                    ],
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "https://www.sevenbridges.com/x": 4597.70458984375,
                    "https://www.sevenbridges.com/y": -652.3868408203125
                },
                {
                    "id": "#main/post_mapping_fastqc_zip",
                    "outputSource": [
                        "#main/qc_post_mapping/fastqc_zip"
                    ],
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "https://www.sevenbridges.com/x": 4538.44384765625,
                    "https://www.sevenbridges.com/y": -750.1759643554688
                },
                {
                    "id": "#main/post_mapping_flagstats",
                    "outputSource": [
                        "#main/flagstats_post_mapping/flagstat_output"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 4442.7138671875,
                    "https://www.sevenbridges.com/y": -392.1081237792969
                }
            ],
            "steps": [
                {
                    "id": "#main/conversion_estimation_spike_in",
                    "in": [
                        {
                            "id": "#main/conversion_estimation_spike_in/output_basename",
                            "source": "#main/sample_id"
                        },
                        {
                            "id": "#main/conversion_estimation_spike_in/spike_in_mcall_bedgraph",
                            "source": "#main/methylation_calling_spike_in/mcall_bedgraph"
                        }
                    ],
                    "out": [
                        {
                            "id": "#main/conversion_estimation_spike_in/bisulfite_conversion_file"
                        }
                    ],
                    "run": "#bisulfite_conversion_spike_in.cwl",
                    "https://www.sevenbridges.com/x": 4128.88916015625,
                    "https://www.sevenbridges.com/y": 757.0206298828125
                },
                {
                    "id": "#main/create_summary_qc_report",
                    "in": [
                        {
                            "id": "#main/create_summary_qc_report/qc_files_array",
                            "linkMerge": "merge_flattened",
                            "source": [
                                "#main/remove_duplicates/picard_markdup_log",
                                "#main/qc_post_mapping/fastqc_zip",
                                "#main/qc_post_mapping/fastqc_html",
                                "#main/flagstats_post_mapping/flagstat_output"
                            ]
                        },
                        {
                            "id": "#main/create_summary_qc_report/report_name",
                            "source": "#main/sample_id"
                        }
                    ],
                    "out": [
                        {
                            "id": "#main/create_summary_qc_report/multiqc_html"
                        },
                        {
                            "id": "#main/create_summary_qc_report/multiqc_zip"
                        }
                    ],
                    "run": "#multiqc_hack.cwl",
                    "doc": "multiqc summarizes the qc results from fastqc \nand other tools\n",
                    "https://www.sevenbridges.com/x": 4182.88134765625,
                    "https://www.sevenbridges.com/y": 385.27703857421875
                },
                {
                    "id": "#main/flagstats_post_mapping",
                    "in": [
                        {
                            "id": "#main/flagstats_post_mapping/bam",
                            "source": "#main/index_bam/bam_sorted_indexed"
                        }
                    ],
                    "out": [
                        {
                            "id": "#main/flagstats_post_mapping/flagstat_output"
                        }
                    ],
                    "run": "#samtools_flagstat.cwl",
                    "doc": "samtools flagstat\n",
                    "https://www.sevenbridges.com/x": 3396.394775390625,
                    "https://www.sevenbridges.com/y": -404.9775695800781
                },
                {
                    "id": "#main/get_spike_in_reads",
                    "in": [
                        {
                            "id": "#main/get_spike_in_reads/bam",
                            "source": "#main/index_bam/bam_sorted_indexed"
                        },
                        {
                            "id": "#main/get_spike_in_reads/region",
                            "source": "#main/spikein_chr_name"
                        }
                    ],
                    "out": [
                        {
                            "id": "#main/get_spike_in_reads/bam_spike_in"
                        }
                    ],
                    "run": "#samtools_view_extract_spike_in.cwl",
                    "doc": "extracts reads mapping to the unmethylated spike in",
                    "https://www.sevenbridges.com/x": 2843,
                    "https://www.sevenbridges.com/y": 886
                },
                {
                    "id": "#main/index_bam",
                    "in": [
                        {
                            "id": "#main/index_bam/bam_sorted",
                            "source": "#main/sorting_merged_bam/bam_sorted"
                        }
                    ],
                    "out": [
                        {
                            "id": "#main/index_bam/bam_sorted_indexed"
                        }
                    ],
                    "run": "#samtools_index_hack.cwl",
                    "doc": "samtools index - indexes sorted bam\n",
                    "https://www.sevenbridges.com/x": 1901,
                    "https://www.sevenbridges.com/y": 817
                },
                {
                    "id": "#main/index_spike_in_bam",
                    "in": [
                        {
                            "id": "#main/index_spike_in_bam/bam_sorted",
                            "source": "#main/get_spike_in_reads/bam_spike_in"
                        }
                    ],
                    "out": [
                        {
                            "id": "#main/index_spike_in_bam/bam_sorted_indexed"
                        }
                    ],
                    "run": "#samtools_index_hack.cwl",
                    "doc": "samtools index - indexes sorted bam\n",
                    "https://www.sevenbridges.com/x": 3262.044189453125,
                    "https://www.sevenbridges.com/y": 905.8668212890625
                },
                {
                    "id": "#main/lane_replicate_merging",
                    "in": [
                        {
                            "id": "#main/lane_replicate_merging/output_name",
                            "source": "#main/sample_id",
                            "valueFrom": "$(self + \".bam\")"
                        },
                        {
                            "id": "#main/lane_replicate_merging/bams",
                            "source": [
                                "#main/remove_duplicates/bam_duprem"
                            ]
                        }
                    ],
                    "out": [
                        {
                            "id": "#main/lane_replicate_merging/bam_merged"
                        }
                    ],
                    "run": "#samtools_merge.cwl",
                    "doc": "samtools merge - merging bam files of lane replicates",
                    "https://www.sevenbridges.com/x": 1076.6666259765625,
                    "https://www.sevenbridges.com/y": -376.22222900390625
                },
                {
                    "id": "#main/mbias_calculation",
                    "in": [
                        {
                            "id": "#main/mbias_calculation/bam",
                            "source": "#main/index_bam/bam_sorted_indexed"
                        },
                        {
                            "id": "#main/mbias_calculation/nCTOB",
                            "source": "#main/methyldackel_nctob"
                        },
                        {
                            "id": "#main/mbias_calculation/nCTOT",
                            "source": "#main/methyldackel_nctot"
                        },
                        {
                            "id": "#main/mbias_calculation/nOB",
                            "source": "#main/methyldackel_nob"
                        },
                        {
                            "id": "#main/mbias_calculation/nOT",
                            "source": "#main/methyldackel_not"
                        },
                        {
                            "id": "#main/mbias_calculation/noCG",
                            "source": "#main/methyldackel_noCG"
                        },
                        {
                            "id": "#main/mbias_calculation/output_basename",
                            "source": "#main/sample_id"
                        },
                        {
                            "id": "#main/mbias_calculation/reference",
                            "source": "#main/reference"
                        },
                        {
                            "id": "#main/mbias_calculation/threads",
                            "source": "#main/max_threads"
                        }
                    ],
                    "out": [
                        {
                            "id": "#main/mbias_calculation/mbias_output"
                        }
                    ],
                    "run": "#methyldackel_mbias.cwl",
                    "https://www.sevenbridges.com/x": 2867.6005859375,
                    "https://www.sevenbridges.com/y": 1690.6920166015625
                },
                {
                    "id": "#main/mbias_calculation_wo_trimming",
                    "in": [
                        {
                            "id": "#main/mbias_calculation_wo_trimming/bam",
                            "source": "#main/index_bam/bam_sorted_indexed"
                        },
                        {
                            "id": "#main/mbias_calculation_wo_trimming/output_basename",
                            "source": "#main/sample_id",
                            "valueFrom": "$(self + \"_wo_trimming\")"
                        },
                        {
                            "id": "#main/mbias_calculation_wo_trimming/reference",
                            "source": "#main/reference"
                        },
                        {
                            "id": "#main/mbias_calculation_wo_trimming/threads",
                            "source": "#main/max_threads"
                        }
                    ],
                    "out": [
                        {
                            "id": "#main/mbias_calculation_wo_trimming/mbias_output"
                        }
                    ],
                    "run": "#methyldackel_mbias.cwl",
                    "https://www.sevenbridges.com/x": 3018.033447265625,
                    "https://www.sevenbridges.com/y": 48.106319427490234
                },
                {
                    "id": "#main/methylation_calling",
                    "in": [
                        {
                            "id": "#main/methylation_calling/CTOB",
                            "source": "#main/methyldackel_ctob"
                        },
                        {
                            "id": "#main/methylation_calling/CTOT",
                            "source": "#main/methyldackel_ctot"
                        },
                        {
                            "id": "#main/methylation_calling/OB",
                            "source": "#main/methyldackel_ob"
                        },
                        {
                            "id": "#main/methylation_calling/OT",
                            "source": "#main/methyldackel_ot"
                        },
                        {
                            "id": "#main/methylation_calling/bam",
                            "source": "#main/index_bam/bam_sorted_indexed"
                        },
                        {
                            "id": "#main/methylation_calling/min_depth",
                            "source": "#main/methyldackel_min_depth"
                        },
                        {
                            "id": "#main/methylation_calling/min_mapq",
                            "source": "#main/methyldackel_min_mapq"
                        },
                        {
                            "id": "#main/methylation_calling/min_phred",
                            "source": "#main/methyldackel_min_phred"
                        },
                        {
                            "id": "#main/methylation_calling/nCTOB",
                            "source": "#main/methyldackel_nctob"
                        },
                        {
                            "id": "#main/methylation_calling/nCTOT",
                            "source": "#main/methyldackel_nctot"
                        },
                        {
                            "id": "#main/methylation_calling/nOB",
                            "source": "#main/methyldackel_nob"
                        },
                        {
                            "id": "#main/methylation_calling/nOT",
                            "source": "#main/methyldackel_not"
                        },
                        {
                            "id": "#main/methylation_calling/noCG",
                            "source": "#main/methyldackel_noCG"
                        },
                        {
                            "id": "#main/methylation_calling/output_basename",
                            "source": "#main/sample_id"
                        },
                        {
                            "id": "#main/methylation_calling/reference",
                            "source": "#main/reference"
                        },
                        {
                            "id": "#main/methylation_calling/threads",
                            "source": "#main/max_threads"
                        }
                    ],
                    "out": [
                        {
                            "id": "#main/methylation_calling/mcall_bedgraph"
                        }
                    ],
                    "run": "#methyldackel_extract.cwl",
                    "https://www.sevenbridges.com/x": 2831.132080078125,
                    "https://www.sevenbridges.com/y": 2103.622802734375
                },
                {
                    "id": "#main/methylation_calling_spike_in",
                    "in": [
                        {
                            "id": "#main/methylation_calling_spike_in/CTOB",
                            "source": "#main/methyldackel_ctob"
                        },
                        {
                            "id": "#main/methylation_calling_spike_in/CTOT",
                            "source": "#main/methyldackel_ctot"
                        },
                        {
                            "id": "#main/methylation_calling_spike_in/OB",
                            "source": "#main/methyldackel_ob"
                        },
                        {
                            "id": "#main/methylation_calling_spike_in/OT",
                            "source": "#main/methyldackel_ot"
                        },
                        {
                            "id": "#main/methylation_calling_spike_in/bam",
                            "source": "#main/index_spike_in_bam/bam_sorted_indexed"
                        },
                        {
                            "id": "#main/methylation_calling_spike_in/min_depth",
                            "source": "#main/methyldackel_min_depth"
                        },
                        {
                            "id": "#main/methylation_calling_spike_in/min_mapq",
                            "source": "#main/methyldackel_min_mapq"
                        },
                        {
                            "id": "#main/methylation_calling_spike_in/min_phred",
                            "source": "#main/methyldackel_min_phred"
                        },
                        {
                            "id": "#main/methylation_calling_spike_in/nCTOB",
                            "source": "#main/methyldackel_nctob"
                        },
                        {
                            "id": "#main/methylation_calling_spike_in/nCTOT",
                            "source": "#main/methyldackel_nctot"
                        },
                        {
                            "id": "#main/methylation_calling_spike_in/nOB",
                            "source": "#main/methyldackel_nob"
                        },
                        {
                            "id": "#main/methylation_calling_spike_in/nOT",
                            "source": "#main/methyldackel_not"
                        },
                        {
                            "id": "#main/methylation_calling_spike_in/noCG",
                            "source": "#main/methyldackel_noCG"
                        },
                        {
                            "id": "#main/methylation_calling_spike_in/output_basename",
                            "source": "#main/sample_id",
                            "valueFrom": "$(self + \"_spike_in\")"
                        },
                        {
                            "id": "#main/methylation_calling_spike_in/reference",
                            "source": "#main/reference"
                        },
                        {
                            "id": "#main/methylation_calling_spike_in/threads",
                            "source": "#main/max_threads"
                        }
                    ],
                    "out": [
                        {
                            "id": "#main/methylation_calling_spike_in/mcall_bedgraph"
                        }
                    ],
                    "run": "#methyldackel_extract.cwl",
                    "https://www.sevenbridges.com/x": 3848.409912109375,
                    "https://www.sevenbridges.com/y": 1283.5234375
                },
                {
                    "id": "#main/qc_post_mapping",
                    "in": [
                        {
                            "id": "#main/qc_post_mapping/bam",
                            "source": "#main/index_bam/bam_sorted_indexed"
                        }
                    ],
                    "out": [
                        {
                            "id": "#main/qc_post_mapping/fastqc_html"
                        },
                        {
                            "id": "#main/qc_post_mapping/fastqc_zip"
                        }
                    ],
                    "run": "#fastqc.cwl",
                    "doc": "fastqc - quality control for reads after mapping and duplication removal",
                    "https://www.sevenbridges.com/x": 3405.511962890625,
                    "https://www.sevenbridges.com/y": -539.9010620117188
                },
                {
                    "id": "#main/remove_duplicates",
                    "in": [
                        {
                            "id": "#main/remove_duplicates/bam_sorted",
                            "source": "#main/samtools_sort/bam_sorted"
                        }
                    ],
                    "out": [
                        {
                            "id": "#main/remove_duplicates/bam_duprem"
                        },
                        {
                            "id": "#main/remove_duplicates/picard_markdup_log"
                        }
                    ],
                    "run": "#picard_markdup.cwl",
                    "doc": "picard markdup - emoves duplicates from a single sorted bam file.",
                    "scatter": [
                        "#main/remove_duplicates/bam_sorted"
                    ],
                    "scatterMethod": "dotproduct",
                    "https://www.sevenbridges.com/x": 776,
                    "https://www.sevenbridges.com/y": -538.3333129882812
                },
                {
                    "id": "#main/sorting_merged_bam",
                    "in": [
                        {
                            "id": "#main/sorting_merged_bam/bam_unsorted",
                            "source": "#main/lane_replicate_merging/bam_merged"
                        }
                    ],
                    "out": [
                        {
                            "id": "#main/sorting_merged_bam/bam_sorted"
                        }
                    ],
                    "run": "#samtools_sort.cwl",
                    "doc": "samtools sort - sorting of merged bam",
                    "https://www.sevenbridges.com/x": 1469.843017578125,
                    "https://www.sevenbridges.com/y": -376.81353759765625
                },
                {
                    "id": "#main/map",
                    "in": [
                        {
                            "id": "#main/map/fastq1",
                            "source": "#main/fastq1"
                        },
                        {
                            "id": "#main/map/fastq2",
                            "source": "#main/fastq2"
                        },
                        {
                            "id": "#main/map/is_non_directional",
                            "default": false,
                            "source": "#main/is_non_directional"
                        },
                        {
                            "id": "#main/map/reference",
                            "source": "#main/reference"
                        },
                        {
                            "id": "#main/map/threads",
                            "source": "#main/max_threads"
                        }
                    ],
                    "out": [
                        {
                            "id": "#main/map/sam"
                        }
                    ],
                    "run": "#bwameth.cwl",
                    "scatter": [
                        "#main/map/fastq1",
                        "#main/map/fastq2"
                    ],
                    "scatterMethod": "dotproduct",
                    "https://www.sevenbridges.com/x": 303.3333435058594,
                    "https://www.sevenbridges.com/y": -534.6666870117188
                },
                {
                    "id": "#main/samtools_sort",
                    "in": [
                        {
                            "id": "#main/samtools_sort/bam_unsorted",
                            "source": "#main/map/sam"
                        }
                    ],
                    "out": [
                        {
                            "id": "#main/samtools_sort/bam_sorted"
                        }
                    ],
                    "run": "#samtools_sort.cwl",
                    "scatter": [
                        "#main/samtools_sort/bam_unsorted"
                    ],
                    "scatterMethod": "dotproduct",
                    "https://www.sevenbridges.com/x": 518.7777709960938,
                    "https://www.sevenbridges.com/y": -536.4444580078125
                }
            ],
            "requirements": [
                {
                    "class": "ScatterFeatureRequirement"
                },
                {
                    "class": "MultipleInputFeatureRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "SubworkflowFeatureRequirement"
                }
            ],
            "id": "#main"
        }
    ],
    "cwlVersion": "v1.0"
}