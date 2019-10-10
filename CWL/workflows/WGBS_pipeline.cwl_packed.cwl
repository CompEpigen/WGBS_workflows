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
            "id": "#bisulfite_conversion_spike_in.cwl"
        },
        {
            "class": "CommandLineTool",
            "requirements": [
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
                    "class": "ResourceRequirement"
                }
            ],
            "baseCommand": [
                "bwameth.py"
            ],
            "stdout": "$(inputs.fastq1.nameroot + \".sam\")",
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
                    "dockerPull": "kerstenbreuer/trim_galore:0.4.4_1.14_0.11.7",
                    "class": "DockerRequirement"
                },
                {
                    "coresMin": 1,
                    "ramMin": 5000,
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
                    "id": "#fastqc.cwl/fastq1"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 2
                    },
                    "id": "#fastqc.cwl/fastq2"
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
                    "class": "ResourceRequirement"
                }
            ],
            "baseCommand": [
                "MethylDackel",
                "extract"
            ],
            "inputs": [
                {
                    "type": "string",
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--CTOB"
                    },
                    "id": "#methyldackel_extract.cwl/CTOB"
                },
                {
                    "type": "string",
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--CTOT"
                    },
                    "id": "#methyldackel_extract.cwl/CTOT"
                },
                {
                    "type": "string",
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--OB"
                    },
                    "id": "#methyldackel_extract.cwl/OB"
                },
                {
                    "type": "string",
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
                    "type": "string",
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--nCTOB"
                    },
                    "id": "#methyldackel_extract.cwl/nCTOB"
                },
                {
                    "type": "string",
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--nCTOT"
                    },
                    "id": "#methyldackel_extract.cwl/nCTOT"
                },
                {
                    "type": "string",
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--nOB"
                    },
                    "id": "#methyldackel_extract.cwl/nOB"
                },
                {
                    "type": "string",
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
            "class": "CommandLineTool",
            "hints": [
                {
                    "dockerPull": "kerstenbreuer/samtools:1.7",
                    "class": "DockerRequirement"
                },
                {
                    "coresMin": 1,
                    "ramMin": 10000,
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
                    "inputBinding": {
                        "position": 2
                    },
                    "id": "#samtools_view_sam2bam.cwl/sam"
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
            "stdout": "$( inputs.sam.nameroot + \".bam\" )",
            "outputs": [
                {
                    "type": "stdout",
                    "id": "#samtools_view_sam2bam.cwl/bam_unsorted"
                }
            ],
            "id": "#samtools_view_sam2bam.cwl"
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
                    "dockerPull": "dukegcb/trimmomatic:latest",
                    "class": "DockerRequirement"
                },
                {
                    "coresMin": "$( inputs.threads )",
                    "ramMin": 28000,
                    "class": "ResourceRequirement"
                }
            ],
            "baseCommand": [
                "java"
            ],
            "arguments": [
                {
                    "valueFrom": "${\n  if ( inputs.fastq2 == null){\n    return \"SE\"\n  } else {\n    return \"PE\"\n  }\n}\n",
                    "position": 3
                },
                {
                    "valueFrom": "$(inputs.fastq1.nameroot + \".trimmomatic.log\")",
                    "position": 4,
                    "prefix": "-trimlog"
                },
                {
                    "valueFrom": "$(inputs.fastq1.nameroot + \"_trimmed.fastq\")",
                    "position": 7
                },
                {
                    "valueFrom": "$(inputs.fastq1.nameroot + \"_trimmed_unpaired.fastq\")",
                    "position": 8
                },
                {
                    "valueFrom": "${ \n  if( inputs.fastq2 != null){\n    return inputs.fastq2.nameroot + \"_trimmed.fastq\"\n  }\n  else{\n    return null\n  }\n}\n",
                    "position": 9
                },
                {
                    "valueFrom": "${ \n  if( inputs.fastq2 != null){\n    return inputs.fastq2.nameroot + \"_trimmed_unpaired.fastq\"\n  }\n  else{\n    return null\n  }\n}\n",
                    "position": 10
                },
                {
                    "valueFrom": "$(\"ILLUMINACLIP:\" + inputs.adapters_file.path + \":\"+ inputs.illuminaclip)",
                    "position": 11
                }
            ],
            "inputs": [
                {
                    "doc": "FASTA file containing adapters, PCR sequences, etc. It is used to search for and remove these sequences in the input FASTQ file(s)",
                    "type": "File",
                    "id": "#trimmomatic.cwl/adapters_file"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 101,
                        "prefix": "AVGQUAL:",
                        "separate": false
                    },
                    "id": "#trimmomatic.cwl/avgqual"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 13,
                        "prefix": "CROP:",
                        "separate": false
                    },
                    "id": "#trimmomatic.cwl/crop"
                },
                {
                    "doc": "FASTQ file for input read (read R1 in Paired End mode)",
                    "type": "File",
                    "inputBinding": {
                        "position": 5
                    },
                    "id": "#trimmomatic.cwl/fastq1"
                },
                {
                    "doc": "FASTQ file for read R2 in Paired End mode",
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 6
                    },
                    "id": "#trimmomatic.cwl/fastq2"
                },
                {
                    "type": "string",
                    "id": "#trimmomatic.cwl/illuminaclip"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "shellQuote": false,
                        "position": 1
                    },
                    "id": "#trimmomatic.cwl/java_opts"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 14,
                        "prefix": "LEADING:",
                        "separate": false
                    },
                    "id": "#trimmomatic.cwl/leading"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 4,
                        "prefix": "-trimlog"
                    },
                    "id": "#trimmomatic.cwl/log_filename"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 100,
                        "prefix": "MINLEN:",
                        "separate": false
                    },
                    "id": "#trimmomatic.cwl/minlen"
                },
                {
                    "doc": "default path matching the applied docker container; \nif the container is not used please adapt\n",
                    "type": "string",
                    "default": "/usr/share/java/trimmomatic.jar",
                    "inputBinding": {
                        "position": 2,
                        "prefix": "-jar"
                    },
                    "id": "#trimmomatic.cwl/path_to_trimmomatic"
                },
                {
                    "type": "string",
                    "default": "64",
                    "inputBinding": {
                        "position": 4,
                        "prefix": "-phred",
                        "separate": false
                    },
                    "id": "#trimmomatic.cwl/phred"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 15,
                        "prefix": "SLIDINGWINDOW:",
                        "separate": false
                    },
                    "id": "#trimmomatic.cwl/slidingwindow"
                },
                {
                    "doc": "Number of threads",
                    "type": "int",
                    "default": 10,
                    "inputBinding": {
                        "position": 4,
                        "prefix": "-threads"
                    },
                    "id": "#trimmomatic.cwl/threads"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 14,
                        "prefix": "TRAILING:",
                        "separate": false
                    },
                    "id": "#trimmomatic.cwl/trailing"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.fastq1.nameroot + \"_trimmed.fastq\")"
                    },
                    "id": "#trimmomatic.cwl/fastq1_trimmed"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "$(inputs.fastq1.nameroot + \"_trimmed_unpaired.fastq\")"
                    },
                    "id": "#trimmomatic.cwl/fastq1_trimmed_unpaired"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "${ \n  if( inputs.fastq2 != null){\n    return inputs.fastq2.nameroot + \"_trimmed.fastq\"\n  }\n  else{\n    return null\n  }\n}\n"
                    },
                    "id": "#trimmomatic.cwl/fastq2_trimmed"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "${ \n  if( inputs.fastq2 != null){\n    return inputs.fastq2.nameroot + \"_trimmed_unpaired.fastq\"\n  }\n  else{\n    return null\n  }\n}\n"
                    },
                    "id": "#trimmomatic.cwl/fastq2_trimmed_unpaired"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.fastq1.nameroot + \".trimmomatic.log\")"
                    },
                    "id": "#trimmomatic.cwl/trimmomatic_log"
                }
            ],
            "id": "#trimmomatic.cwl"
        },
        {
            "class": "Workflow",
            "inputs": [
                {
                    "type": "File",
                    "id": "#trim_map.cwl/fastq1"
                },
                {
                    "type": "File",
                    "id": "#trim_map.cwl/fastq2"
                },
                {
                    "type": "boolean",
                    "default": false,
                    "id": "#trim_map.cwl/is_non_directional"
                },
                {
                    "type": "int",
                    "default": 10,
                    "id": "#trim_map.cwl/max_threads"
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
                    "id": "#trim_map.cwl/reference"
                },
                {
                    "type": "File",
                    "id": "#trim_map.cwl/trimmomatic_adapters_file"
                },
                {
                    "type": "int",
                    "default": 1,
                    "id": "#trim_map.cwl/trimmomatic_avgqual"
                },
                {
                    "type": "int",
                    "default": 1000,
                    "id": "#trim_map.cwl/trimmomatic_crop"
                },
                {
                    "type": "int",
                    "default": 10,
                    "id": "#trim_map.cwl/trimmomatic_headcrop"
                },
                {
                    "type": "string",
                    "default": "2:30:10:8:true",
                    "id": "#trim_map.cwl/trimmomatic_illuminaclip"
                },
                {
                    "type": "int",
                    "default": 0,
                    "id": "#trim_map.cwl/trimmomatic_leading"
                },
                {
                    "type": "int",
                    "default": 0,
                    "id": "#trim_map.cwl/trimmomatic_minlen"
                },
                {
                    "type": "string",
                    "default": "64",
                    "id": "#trim_map.cwl/trimmomatic_phred"
                },
                {
                    "type": "int",
                    "default": 10,
                    "id": "#trim_map.cwl/trimmomatic_tailcrop"
                },
                {
                    "type": "int",
                    "default": 0,
                    "id": "#trim_map.cwl/trimmomatic_trailing"
                }
            ],
            "steps": [
                {
                    "run": "#trimmomatic.cwl",
                    "in": [
                        {
                            "source": "#trim_map.cwl/trimmomatic_adapters_file",
                            "id": "#trim_map.cwl/adaptor_trimming/adapters_file"
                        },
                        {
                            "source": "#trim_map.cwl/trimmomatic_avgqual",
                            "id": "#trim_map.cwl/adaptor_trimming/avgqual"
                        },
                        {
                            "source": "#trim_map.cwl/trimmomatic_crop",
                            "id": "#trim_map.cwl/adaptor_trimming/crop"
                        },
                        {
                            "source": "#trim_map.cwl/fastq1",
                            "id": "#trim_map.cwl/adaptor_trimming/fastq1"
                        },
                        {
                            "source": "#trim_map.cwl/fastq2",
                            "id": "#trim_map.cwl/adaptor_trimming/fastq2"
                        },
                        {
                            "source": "#trim_map.cwl/trimmomatic_illuminaclip",
                            "id": "#trim_map.cwl/adaptor_trimming/illuminaclip"
                        },
                        {
                            "source": "#trim_map.cwl/trimmomatic_leading",
                            "id": "#trim_map.cwl/adaptor_trimming/leading"
                        },
                        {
                            "source": "#trim_map.cwl/trimmomatic_minlen",
                            "id": "#trim_map.cwl/adaptor_trimming/minlen"
                        },
                        {
                            "source": "#trim_map.cwl/trimmomatic_phred",
                            "id": "#trim_map.cwl/adaptor_trimming/phred"
                        },
                        {
                            "source": "#trim_map.cwl/max_threads",
                            "id": "#trim_map.cwl/adaptor_trimming/threads"
                        },
                        {
                            "source": "#trim_map.cwl/trimmomatic_trailing",
                            "id": "#trim_map.cwl/adaptor_trimming/trailing"
                        }
                    ],
                    "out": [
                        "#trim_map.cwl/adaptor_trimming/fastq1_trimmed",
                        "#trim_map.cwl/adaptor_trimming/fastq2_trimmed",
                        "#trim_map.cwl/adaptor_trimming/fastq1_trimmed_unpaired",
                        "#trim_map.cwl/adaptor_trimming/fastq2_trimmed_unpaired",
                        "#trim_map.cwl/adaptor_trimming/trimmomatic_log"
                    ],
                    "id": "#trim_map.cwl/adaptor_trimming"
                },
                {
                    "run": "#bwameth.cwl",
                    "in": [
                        {
                            "source": "#trim_map.cwl/adaptor_trimming/fastq1_trimmed",
                            "id": "#trim_map.cwl/mapping/fastq1"
                        },
                        {
                            "source": "#trim_map.cwl/adaptor_trimming/fastq2_trimmed",
                            "id": "#trim_map.cwl/mapping/fastq2"
                        },
                        {
                            "source": "#trim_map.cwl/is_non_directional",
                            "id": "#trim_map.cwl/mapping/is_non_directional"
                        },
                        {
                            "source": "#trim_map.cwl/reference",
                            "id": "#trim_map.cwl/mapping/reference"
                        },
                        {
                            "source": "#trim_map.cwl/max_threads",
                            "id": "#trim_map.cwl/mapping/threads"
                        }
                    ],
                    "out": [
                        "#trim_map.cwl/mapping/sam"
                    ],
                    "id": "#trim_map.cwl/mapping"
                },
                {
                    "doc": "fastqc - quality control for raw fastqs",
                    "run": "#fastqc.cwl",
                    "in": [
                        {
                            "source": "#trim_map.cwl/fastq1",
                            "id": "#trim_map.cwl/qc_raw/fastq1"
                        },
                        {
                            "source": "#trim_map.cwl/fastq2",
                            "id": "#trim_map.cwl/qc_raw/fastq2"
                        }
                    ],
                    "out": [
                        "#trim_map.cwl/qc_raw/fastqc_zip",
                        "#trim_map.cwl/qc_raw/fastqc_html"
                    ],
                    "id": "#trim_map.cwl/qc_raw"
                },
                {
                    "doc": "fastqc - quality control for raw fastqs",
                    "run": "#fastqc.cwl",
                    "in": [
                        {
                            "source": "#trim_map.cwl/adaptor_trimming/fastq1_trimmed",
                            "id": "#trim_map.cwl/qc_trimmed/fastq1"
                        },
                        {
                            "source": "#trim_map.cwl/adaptor_trimming/fastq2_trimmed",
                            "id": "#trim_map.cwl/qc_trimmed/fastq2"
                        }
                    ],
                    "out": [
                        "#trim_map.cwl/qc_trimmed/fastqc_zip",
                        "#trim_map.cwl/qc_trimmed/fastqc_html"
                    ],
                    "id": "#trim_map.cwl/qc_trimmed"
                },
                {
                    "doc": "samtools view - convert sam to bam",
                    "run": "#samtools_view_sam2bam.cwl",
                    "in": [
                        {
                            "source": "#trim_map.cwl/mapping/sam",
                            "id": "#trim_map.cwl/sam2bam/sam"
                        }
                    ],
                    "out": [
                        "#trim_map.cwl/sam2bam/bam_unsorted"
                    ],
                    "id": "#trim_map.cwl/sam2bam"
                },
                {
                    "doc": "samtools sort - sorts unsorted bam file by coordinates.",
                    "run": "#samtools_sort.cwl",
                    "in": [
                        {
                            "source": "#trim_map.cwl/sam2bam/bam_unsorted",
                            "id": "#trim_map.cwl/sort_bam/bam_unsorted"
                        }
                    ],
                    "out": [
                        "#trim_map.cwl/sort_bam/bam_sorted"
                    ],
                    "id": "#trim_map.cwl/sort_bam"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputSource": "#trim_map.cwl/sort_bam/bam_sorted",
                    "id": "#trim_map.cwl/bam"
                },
                {
                    "type": "File",
                    "outputSource": "#trim_map.cwl/adaptor_trimming/fastq1_trimmed",
                    "id": "#trim_map.cwl/fastq1_trimmed"
                },
                {
                    "type": "File",
                    "outputSource": "#trim_map.cwl/adaptor_trimming/fastq1_trimmed_unpaired",
                    "id": "#trim_map.cwl/fastq1_trimmed_unpaired"
                },
                {
                    "type": "File",
                    "outputSource": "#trim_map.cwl/adaptor_trimming/fastq2_trimmed",
                    "id": "#trim_map.cwl/fastq2_trimmed"
                },
                {
                    "type": "File",
                    "outputSource": "#trim_map.cwl/adaptor_trimming/fastq2_trimmed_unpaired",
                    "id": "#trim_map.cwl/fastq2_trimmed_unpaired"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputSource": "#trim_map.cwl/qc_raw/fastqc_html",
                    "id": "#trim_map.cwl/raw_fastqc_html"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputSource": "#trim_map.cwl/qc_raw/fastqc_zip",
                    "id": "#trim_map.cwl/raw_fastqc_zip"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputSource": "#trim_map.cwl/qc_trimmed/fastqc_html",
                    "id": "#trim_map.cwl/trimmed_fastqc_html"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputSource": "#trim_map.cwl/qc_trimmed/fastqc_zip",
                    "id": "#trim_map.cwl/trimmed_fastqc_zip"
                },
                {
                    "type": "File",
                    "outputSource": "#trim_map.cwl/adaptor_trimming/trimmomatic_log",
                    "id": "#trim_map.cwl/trimmomatic_log"
                }
            ],
            "id": "#trim_map.cwl"
        },
        {
            "class": "Workflow",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "MultipleInputFeatureRequirement"
                },
                {
                    "class": "ScatterFeatureRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                },
                {
                    "class": "SubworkflowFeatureRequirement"
                }
            ],
            "inputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "id": "#main/fastq1"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "id": "#main/fastq2"
                },
                {
                    "type": "boolean",
                    "default": false,
                    "id": "#main/is_non_directional"
                },
                {
                    "type": "int",
                    "default": 10,
                    "id": "#main/max_threads"
                },
                {
                    "type": "string",
                    "default": "0,0,0,0",
                    "id": "#main/methyldackel_ctob"
                },
                {
                    "type": "string",
                    "default": "0,0,0,0",
                    "id": "#main/methyldackel_ctot"
                },
                {
                    "type": "int",
                    "default": 1,
                    "id": "#main/methyldackel_min_depth"
                },
                {
                    "type": "int",
                    "default": 0,
                    "id": "#main/methyldackel_min_mapq"
                },
                {
                    "type": "int",
                    "default": 0,
                    "id": "#main/methyldackel_min_phred"
                },
                {
                    "type": "string",
                    "default": "10,10,10,10",
                    "id": "#main/methyldackel_nctob"
                },
                {
                    "type": "string",
                    "default": "10,10,10,10",
                    "id": "#main/methyldackel_nctot"
                },
                {
                    "type": "boolean",
                    "default": false,
                    "id": "#main/methyldackel_noCG"
                },
                {
                    "type": "string",
                    "default": "10,10,10,10",
                    "id": "#main/methyldackel_nob"
                },
                {
                    "type": "string",
                    "default": "10,10,10,10",
                    "id": "#main/methyldackel_not"
                },
                {
                    "type": "string",
                    "default": "0,0,0,0",
                    "id": "#main/methyldackel_ob"
                },
                {
                    "type": "string",
                    "default": "0,0,0,0",
                    "id": "#main/methyldackel_ot"
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
                    "id": "#main/reference"
                },
                {
                    "type": "string",
                    "id": "#main/sample_id"
                },
                {
                    "doc": "chromosome name in the reference genome which corresponds to the unmethylated spikein DNA",
                    "type": "string",
                    "default": "Lamda",
                    "id": "#main/spikein_chr_name"
                },
                {
                    "type": "File",
                    "id": "#main/trimmomatic_adapters_file"
                },
                {
                    "type": "int",
                    "default": 1,
                    "id": "#main/trimmomatic_avgqual"
                },
                {
                    "type": "int",
                    "default": 1000,
                    "id": "#main/trimmomatic_crop"
                },
                {
                    "type": "string",
                    "default": "2:30:10:8:true",
                    "id": "#main/trimmomatic_illuminaclip"
                },
                {
                    "type": "int",
                    "default": 0,
                    "id": "#main/trimmomatic_leading"
                },
                {
                    "type": "int",
                    "default": 0,
                    "id": "#main/trimmomatic_minlen"
                },
                {
                    "type": "string",
                    "default": "64",
                    "id": "#main/trimmomatic_phred"
                },
                {
                    "type": "int",
                    "default": 0,
                    "id": "#main/trimmomatic_trailing"
                }
            ],
            "steps": [
                {
                    "run": "#bisulfite_conversion_spike_in.cwl",
                    "in": [
                        {
                            "source": "#main/sample_id",
                            "id": "#main/conversion_estimation_spike_in/output_basename"
                        },
                        {
                            "source": "#main/methylation_calling_spike_in/mcall_bedgraph",
                            "id": "#main/conversion_estimation_spike_in/spike_in_mcall_bedgraph"
                        }
                    ],
                    "out": [
                        "#main/conversion_estimation_spike_in/bisulfite_conversion_file"
                    ],
                    "id": "#main/conversion_estimation_spike_in"
                },
                {
                    "doc": "multiqc summarizes the qc results from fastqc \nand other tools\n",
                    "run": "#multiqc_hack.cwl",
                    "in": [
                        {
                            "source": [
                                "#main/remove_duplicates/picard_markdup_log",
                                "#main/trim_map/trimmomatic_log",
                                "#main/qc_post_mapping/fastqc_zip",
                                "#main/qc_post_mapping/fastqc_html",
                                "#main/flagstats_post_mapping/flagstat_output"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#main/create_summary_qc_report/qc_files_array"
                        },
                        {
                            "source": [
                                "#main/trim_map/raw_fastqc_zip",
                                "#main/trim_map/raw_fastqc_html",
                                "#main/trim_map/trimmed_fastqc_html",
                                "#main/trim_map/trimmed_fastqc_zip"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#main/create_summary_qc_report/qc_files_array_of_array"
                        },
                        {
                            "source": "#main/sample_id",
                            "id": "#main/create_summary_qc_report/report_name"
                        }
                    ],
                    "out": [
                        "#main/create_summary_qc_report/multiqc_zip",
                        "#main/create_summary_qc_report/multiqc_html"
                    ],
                    "id": "#main/create_summary_qc_report"
                },
                {
                    "doc": "samtools flagstat\n",
                    "run": "#samtools_flagstat.cwl",
                    "in": [
                        {
                            "source": "#main/index_bam/bam_sorted_indexed",
                            "id": "#main/flagstats_post_mapping/bam"
                        }
                    ],
                    "out": [
                        "#main/flagstats_post_mapping/flagstat_output"
                    ],
                    "id": "#main/flagstats_post_mapping"
                },
                {
                    "doc": "extracts reads mapping to the unmethylated spike in",
                    "run": "#samtools_view_extract_spike_in.cwl",
                    "in": [
                        {
                            "source": "#main/index_bam/bam_sorted_indexed",
                            "id": "#main/get_spike_in_reads/bam"
                        },
                        {
                            "source": "#main/spikein_chr_name",
                            "id": "#main/get_spike_in_reads/region"
                        }
                    ],
                    "out": [
                        "#main/get_spike_in_reads/bam_spike_in"
                    ],
                    "id": "#main/get_spike_in_reads"
                },
                {
                    "doc": "samtools index - indexes sorted bam\n",
                    "run": "#samtools_index_hack.cwl",
                    "in": [
                        {
                            "source": "#main/remove_duplicates/bam_duprem",
                            "id": "#main/index_bam/bam_sorted"
                        }
                    ],
                    "out": [
                        "#main/index_bam/bam_sorted_indexed"
                    ],
                    "id": "#main/index_bam"
                },
                {
                    "doc": "samtools index - indexes sorted bam\n",
                    "run": "#samtools_index_hack.cwl",
                    "in": [
                        {
                            "source": "#main/get_spike_in_reads/bam_spike_in",
                            "id": "#main/index_spike_in_bam/bam_sorted"
                        }
                    ],
                    "out": [
                        "#main/index_spike_in_bam/bam_sorted_indexed"
                    ],
                    "id": "#main/index_spike_in_bam"
                },
                {
                    "doc": "samtools merge - merging bam files of lane replicates",
                    "run": "#samtools_merge.cwl",
                    "in": [
                        {
                            "source": "#main/trim_map/bam",
                            "id": "#main/lane_replicate_merging/bams"
                        },
                        {
                            "source": "#main/sample_id",
                            "valueFrom": "$(self + \".bam\")",
                            "id": "#main/lane_replicate_merging/output_name"
                        }
                    ],
                    "out": [
                        "#main/lane_replicate_merging/bam_merged"
                    ],
                    "id": "#main/lane_replicate_merging"
                },
                {
                    "run": "#methyldackel_mbias.cwl",
                    "in": [
                        {
                            "source": "#main/index_bam/bam_sorted_indexed",
                            "id": "#main/mbias_calculation/bam"
                        },
                        {
                            "source": "#main/methyldackel_nctob",
                            "id": "#main/mbias_calculation/nCTOB"
                        },
                        {
                            "source": "#main/methyldackel_nctot",
                            "id": "#main/mbias_calculation/nCTOT"
                        },
                        {
                            "source": "#main/methyldackel_nob",
                            "id": "#main/mbias_calculation/nOB"
                        },
                        {
                            "source": "#main/methyldackel_not",
                            "id": "#main/mbias_calculation/nOT"
                        },
                        {
                            "source": "#main/methyldackel_noCG",
                            "id": "#main/mbias_calculation/noCG"
                        },
                        {
                            "source": "#main/sample_id",
                            "id": "#main/mbias_calculation/output_basename"
                        },
                        {
                            "source": "#main/reference",
                            "id": "#main/mbias_calculation/reference"
                        },
                        {
                            "source": "#main/max_threads",
                            "id": "#main/mbias_calculation/threads"
                        }
                    ],
                    "out": [
                        "#main/mbias_calculation/mbias_output"
                    ],
                    "id": "#main/mbias_calculation"
                },
                {
                    "run": "#methyldackel_mbias.cwl",
                    "in": [
                        {
                            "source": "#main/index_bam/bam_sorted_indexed",
                            "id": "#main/mbias_calculation_wo_trimming/bam"
                        },
                        {
                            "source": "#main/sample_id",
                            "valueFrom": "$(self + \"_wo_trimming\")",
                            "id": "#main/mbias_calculation_wo_trimming/output_basename"
                        },
                        {
                            "source": "#main/reference",
                            "id": "#main/mbias_calculation_wo_trimming/reference"
                        },
                        {
                            "source": "#main/max_threads",
                            "id": "#main/mbias_calculation_wo_trimming/threads"
                        }
                    ],
                    "out": [
                        "#main/mbias_calculation_wo_trimming/mbias_output"
                    ],
                    "id": "#main/mbias_calculation_wo_trimming"
                },
                {
                    "run": "#methyldackel_extract.cwl",
                    "in": [
                        {
                            "source": "#main/methyldackel_ctob",
                            "id": "#main/methylation_calling/CTOB"
                        },
                        {
                            "source": "#main/methyldackel_ctot",
                            "id": "#main/methylation_calling/CTOT"
                        },
                        {
                            "source": "#main/methyldackel_ob",
                            "id": "#main/methylation_calling/OB"
                        },
                        {
                            "source": "#main/methyldackel_ot",
                            "id": "#main/methylation_calling/OT"
                        },
                        {
                            "source": "#main/index_bam/bam_sorted_indexed",
                            "id": "#main/methylation_calling/bam"
                        },
                        {
                            "source": "#main/methyldackel_min_depth",
                            "id": "#main/methylation_calling/min_depth"
                        },
                        {
                            "source": "#main/methyldackel_min_mapq",
                            "id": "#main/methylation_calling/min_mapq"
                        },
                        {
                            "source": "#main/methyldackel_min_phred",
                            "id": "#main/methylation_calling/min_phred"
                        },
                        {
                            "source": "#main/methyldackel_nctob",
                            "id": "#main/methylation_calling/nCTOB"
                        },
                        {
                            "source": "#main/methyldackel_nctot",
                            "id": "#main/methylation_calling/nCTOT"
                        },
                        {
                            "source": "#main/methyldackel_nob",
                            "id": "#main/methylation_calling/nOB"
                        },
                        {
                            "source": "#main/methyldackel_not",
                            "id": "#main/methylation_calling/nOT"
                        },
                        {
                            "source": "#main/methyldackel_noCG",
                            "id": "#main/methylation_calling/noCG"
                        },
                        {
                            "source": "#main/sample_id",
                            "id": "#main/methylation_calling/output_basename"
                        },
                        {
                            "source": "#main/reference",
                            "id": "#main/methylation_calling/reference"
                        },
                        {
                            "source": "#main/max_threads",
                            "id": "#main/methylation_calling/threads"
                        }
                    ],
                    "out": [
                        "#main/methylation_calling/mcall_bedgraph"
                    ],
                    "id": "#main/methylation_calling"
                },
                {
                    "run": "#methyldackel_extract.cwl",
                    "in": [
                        {
                            "source": "#main/methyldackel_ctob",
                            "id": "#main/methylation_calling_spike_in/CTOB"
                        },
                        {
                            "source": "#main/methyldackel_ctot",
                            "id": "#main/methylation_calling_spike_in/CTOT"
                        },
                        {
                            "source": "#main/methyldackel_ob",
                            "id": "#main/methylation_calling_spike_in/OB"
                        },
                        {
                            "source": "#main/methyldackel_ot",
                            "id": "#main/methylation_calling_spike_in/OT"
                        },
                        {
                            "source": "#main/index_spike_in_bam/bam_sorted_indexed",
                            "id": "#main/methylation_calling_spike_in/bam"
                        },
                        {
                            "source": "#main/methyldackel_min_depth",
                            "id": "#main/methylation_calling_spike_in/min_depth"
                        },
                        {
                            "source": "#main/methyldackel_min_mapq",
                            "id": "#main/methylation_calling_spike_in/min_mapq"
                        },
                        {
                            "source": "#main/methyldackel_min_phred",
                            "id": "#main/methylation_calling_spike_in/min_phred"
                        },
                        {
                            "source": "#main/methyldackel_nctob",
                            "id": "#main/methylation_calling_spike_in/nCTOB"
                        },
                        {
                            "source": "#main/methyldackel_nctot",
                            "id": "#main/methylation_calling_spike_in/nCTOT"
                        },
                        {
                            "source": "#main/methyldackel_nob",
                            "id": "#main/methylation_calling_spike_in/nOB"
                        },
                        {
                            "source": "#main/methyldackel_not",
                            "id": "#main/methylation_calling_spike_in/nOT"
                        },
                        {
                            "source": "#main/methyldackel_noCG",
                            "id": "#main/methylation_calling_spike_in/noCG"
                        },
                        {
                            "source": "#main/sample_id",
                            "valueFrom": "$(self + \"_spike_in\")",
                            "id": "#main/methylation_calling_spike_in/output_basename"
                        },
                        {
                            "source": "#main/reference",
                            "id": "#main/methylation_calling_spike_in/reference"
                        },
                        {
                            "source": "#main/max_threads",
                            "id": "#main/methylation_calling_spike_in/threads"
                        }
                    ],
                    "out": [
                        "#main/methylation_calling_spike_in/mcall_bedgraph"
                    ],
                    "id": "#main/methylation_calling_spike_in"
                },
                {
                    "doc": "fastqc - quality control for reads after mapping and duplication removal",
                    "run": "#fastqc.cwl",
                    "in": [
                        {
                            "source": "#main/index_bam/bam_sorted_indexed",
                            "id": "#main/qc_post_mapping/bam"
                        }
                    ],
                    "out": [
                        "#main/qc_post_mapping/fastqc_zip",
                        "#main/qc_post_mapping/fastqc_html"
                    ],
                    "id": "#main/qc_post_mapping"
                },
                {
                    "doc": "picard markdup - emoves duplicates from a single sorted bam file.",
                    "run": "#picard_markdup.cwl",
                    "in": [
                        {
                            "source": "#main/sorting_merged_bam/bam_sorted",
                            "id": "#main/remove_duplicates/bam_sorted"
                        }
                    ],
                    "out": [
                        "#main/remove_duplicates/bam_duprem",
                        "#main/remove_duplicates/picard_markdup_log"
                    ],
                    "id": "#main/remove_duplicates"
                },
                {
                    "doc": "samtools sort - sorting of merged bam",
                    "run": "#samtools_sort.cwl",
                    "in": [
                        {
                            "source": "#main/lane_replicate_merging/bam_merged",
                            "id": "#main/sorting_merged_bam/bam_unsorted"
                        }
                    ],
                    "out": [
                        "#main/sorting_merged_bam/bam_sorted"
                    ],
                    "id": "#main/sorting_merged_bam"
                },
                {
                    "scatter": [
                        "#main/trim_map/fastq1",
                        "#main/trim_map/fastq2"
                    ],
                    "scatterMethod": "dotproduct",
                    "run": "#trim_map.cwl",
                    "in": [
                        {
                            "source": "#main/fastq1",
                            "id": "#main/trim_map/fastq1"
                        },
                        {
                            "source": "#main/fastq2",
                            "id": "#main/trim_map/fastq2"
                        },
                        {
                            "source": "#main/is_non_directional",
                            "id": "#main/trim_map/is_non_directional"
                        },
                        {
                            "source": "#main/max_threads",
                            "id": "#main/trim_map/max_threads"
                        },
                        {
                            "source": "#main/reference",
                            "id": "#main/trim_map/reference"
                        },
                        {
                            "source": "#main/trimmomatic_adapters_file",
                            "id": "#main/trim_map/trimmomatic_adapters_file"
                        },
                        {
                            "source": "#main/trimmomatic_avgqual",
                            "id": "#main/trim_map/trimmomatic_avgqual"
                        },
                        {
                            "source": "#main/trimmomatic_crop",
                            "id": "#main/trim_map/trimmomatic_crop"
                        },
                        {
                            "source": "#main/trimmomatic_illuminaclip",
                            "id": "#main/trim_map/trimmomatic_illuminaclip"
                        },
                        {
                            "source": "#main/trimmomatic_leading",
                            "id": "#main/trim_map/trimmomatic_leading"
                        },
                        {
                            "source": "#main/trimmomatic_minlen",
                            "id": "#main/trim_map/trimmomatic_minlen"
                        },
                        {
                            "source": "#main/trimmomatic_phred",
                            "id": "#main/trim_map/trimmomatic_phred"
                        },
                        {
                            "source": "#main/trimmomatic_trailing",
                            "id": "#main/trim_map/trimmomatic_trailing"
                        }
                    ],
                    "out": [
                        "#main/trim_map/trimmomatic_log",
                        "#main/trim_map/bam",
                        "#main/trim_map/raw_fastqc_zip",
                        "#main/trim_map/raw_fastqc_html",
                        "#main/trim_map/trimmed_fastqc_zip",
                        "#main/trim_map/trimmed_fastqc_html"
                    ],
                    "id": "#main/trim_map"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputSource": "#main/index_bam/bam_sorted_indexed",
                    "id": "#main/bam"
                },
                {
                    "type": "File",
                    "outputSource": "#main/index_spike_in_bam/bam_sorted_indexed",
                    "id": "#main/bam_spike_in"
                },
                {
                    "type": "File",
                    "outputSource": "#main/conversion_estimation_spike_in/bisulfite_conversion_file",
                    "id": "#main/bisulfite_conversion_file"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputSource": "#main/mbias_calculation/mbias_output",
                    "id": "#main/mbias_output"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputSource": "#main/mbias_calculation_wo_trimming/mbias_output",
                    "id": "#main/mbias_output_wo_trimming"
                },
                {
                    "type": "File",
                    "outputSource": "#main/methylation_calling/mcall_bedgraph",
                    "id": "#main/mcall_bedgraph"
                },
                {
                    "type": "File",
                    "outputSource": "#main/methylation_calling_spike_in/mcall_bedgraph",
                    "id": "#main/mcall_bedgraph_spike_in"
                },
                {
                    "type": "File",
                    "outputSource": "#main/create_summary_qc_report/multiqc_html",
                    "id": "#main/multiqc_html"
                },
                {
                    "type": "File",
                    "outputSource": "#main/create_summary_qc_report/multiqc_zip",
                    "id": "#main/multiqc_zip"
                },
                {
                    "type": "File",
                    "outputSource": "#main/remove_duplicates/picard_markdup_log",
                    "id": "#main/picard_markdup_log"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputSource": "#main/qc_post_mapping/fastqc_html",
                    "id": "#main/post_mapping_fastqc_html"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputSource": "#main/qc_post_mapping/fastqc_zip",
                    "id": "#main/post_mapping_fastqc_zip"
                },
                {
                    "type": "File",
                    "outputSource": "#main/flagstats_post_mapping/flagstat_output",
                    "id": "#main/post_mapping_flagstats"
                },
                {
                    "type": {
                        "type": "array",
                        "items": {
                            "type": "array",
                            "items": "File"
                        }
                    },
                    "outputSource": "#main/trim_map/raw_fastqc_html",
                    "id": "#main/raw_fastqc_html"
                },
                {
                    "type": {
                        "type": "array",
                        "items": {
                            "type": "array",
                            "items": "File"
                        }
                    },
                    "outputSource": "#main/trim_map/raw_fastqc_zip",
                    "id": "#main/raw_fastqc_zip"
                },
                {
                    "type": {
                        "type": "array",
                        "items": {
                            "type": "array",
                            "items": "File"
                        }
                    },
                    "outputSource": "#main/trim_map/trimmed_fastqc_html",
                    "id": "#main/trimmed_fastqc_html"
                },
                {
                    "type": {
                        "type": "array",
                        "items": {
                            "type": "array",
                            "items": "File"
                        }
                    },
                    "outputSource": "#main/trim_map/trimmed_fastqc_zip",
                    "id": "#main/trimmed_fastqc_zip"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputSource": "#main/trim_map/trimmomatic_log",
                    "id": "#main/trimmomatic_log"
                }
            ],
            "id": "#main"
        }
    ],
    "cwlVersion": "v1.0"
}