{
    "$graph": [
        {
            "class": "CommandLineTool",
            "hints": [
                {
                    "dockerPull": "kerstenbreuer/bismark:0.22.3",
                    "class": "DockerRequirement"
                },
                {
                    "coresMin": 4,
                    "ramMin": 20000,
                    "tmpdirMin": 30000,
                    "class": "ResourceRequirement"
                }
            ],
            "baseCommand": "deduplicate_bismark",
            "arguments": [
                {
                    "valueFrom": "${\n  if (inputs.paired_end){\n    return(\"-p\")\n  }\n  else {\n    return(\"-s\")\n  }\n}\n"
                }
            ],
            "inputs": [
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--bam",
                        "position": 10
                    },
                    "id": "#bismark_deduplicate.cwl/aligned_reads"
                },
                {
                    "type": "boolean",
                    "default": true,
                    "id": "#bismark_deduplicate.cwl/paired_end"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.deduplicated.bam"
                    },
                    "id": "#bismark_deduplicate.cwl/dedup_reads"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.deduplication_report.txt"
                    },
                    "id": "#bismark_deduplicate.cwl/log"
                }
            ],
            "id": "#bismark_deduplicate.cwl"
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
                    "dockerPull": "kerstenbreuer/bismark:0.22.3",
                    "class": "DockerRequirement"
                },
                {
                    "coresMin": "$(inputs.threads)",
                    "ramMin": "${return(inputs.threads*30000)}",
                    "tmpdirMin": 30000,
                    "class": "ResourceRequirement"
                }
            ],
            "baseCommand": "bismark_methylation_extractor",
            "arguments": [
                {
                    "valueFrom": "--bedGraph",
                    "position": 1
                },
                {
                    "valueFrom": "--gzip",
                    "position": 1
                },
                {
                    "valueFrom": "--counts",
                    "position": 1
                },
                {
                    "valueFrom": "${\n  if (inputs.paired_end){\n    return(\"-p\")\n  }\n  else {\n    return(\"-s\")\n  }\n}\n",
                    "position": 1
                },
                {
                    "valueFrom": "--report",
                    "position": 1
                },
                {
                    "valueFrom": "--cytosine_report",
                    "position": 1
                }
            ],
            "inputs": [
                {
                    "type": "File",
                    "inputBinding": {
                        "position": 10
                    },
                    "id": "#bismark_methylation_extractor.cwl/aligned_reads"
                },
                {
                    "type": "Directory",
                    "inputBinding": {
                        "prefix": "--genome",
                        "position": 10
                    },
                    "id": "#bismark_methylation_extractor.cwl/genome"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "prefix": "--ignore",
                        "position": 3
                    },
                    "id": "#bismark_methylation_extractor.cwl/ignore"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "prefix": "--ignore_3prime",
                        "position": 3
                    },
                    "id": "#bismark_methylation_extractor.cwl/ignore_3prime"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "prefix": "--ignore_3prime_r2",
                        "position": 3
                    },
                    "id": "#bismark_methylation_extractor.cwl/ignore_3prime_r2"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "prefix": "--ignore_r2",
                        "position": 3
                    },
                    "id": "#bismark_methylation_extractor.cwl/ignore_r2"
                },
                {
                    "type": "boolean",
                    "default": true,
                    "inputBinding": {
                        "prefix": "--no_overlap",
                        "position": 2
                    },
                    "id": "#bismark_methylation_extractor.cwl/no_overlap"
                },
                {
                    "type": "boolean",
                    "default": true,
                    "id": "#bismark_methylation_extractor.cwl/paired_end"
                },
                {
                    "type": "int",
                    "default": 1,
                    "inputBinding": {
                        "prefix": "--multicore",
                        "valueFrom": "${return(Math.ceil(self/5))}",
                        "position": 1
                    },
                    "id": "#bismark_methylation_extractor.cwl/threads"
                }
            ],
            "outputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputBinding": {
                        "glob": "C*_O*.txt*"
                    },
                    "id": "#bismark_methylation_extractor.cwl/context_specific_methylation_reports"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*CpG_report.txt*"
                    },
                    "id": "#bismark_methylation_extractor.cwl/genome_wide_methylation_report"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.M-bias.txt*"
                    },
                    "id": "#bismark_methylation_extractor.cwl/mbias_report"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*bedGraph.gz"
                    },
                    "id": "#bismark_methylation_extractor.cwl/methylation_calls_bedgraph"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*bismark.cov.gz"
                    },
                    "id": "#bismark_methylation_extractor.cwl/methylation_calls_bismark"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*splitting_report.txt*"
                    },
                    "id": "#bismark_methylation_extractor.cwl/splitting_report"
                }
            ],
            "id": "#bismark_methylation_extractor.cwl"
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
            "requirements": [
                {
                    "class": "ShellCommandRequirement"
                }
            ],
            "hints": [
                {
                    "dockerPull": "kerstenbreuer/samtools:1.7",
                    "class": "DockerRequirement"
                },
                {
                    "coresMin": "$(inputs.threads)",
                    "ramMin": 20000,
                    "tmpdirMin": 30000,
                    "class": "ResourceRequirement"
                }
            ],
            "baseCommand": [
                "samtools",
                "merge"
            ],
            "arguments": [
                {
                    "valueFrom": "-",
                    "position": 1
                },
                {
                    "valueFrom": "|",
                    "position": 3,
                    "shellQuote": false
                },
                {
                    "valueFrom": "samtools",
                    "position": 4
                },
                {
                    "valueFrom": "sort",
                    "position": 5
                },
                {
                    "prefix": "-@",
                    "valueFrom": "$(inputs.threads)",
                    "position": 6
                },
                {
                    "prefix": "-o",
                    "valueFrom": "$(inputs.output_name)",
                    "position": 7
                },
                {
                    "valueFrom": "-",
                    "position": 8
                }
            ],
            "inputs": [
                {
                    "doc": "bam files to be merged",
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "inputBinding": {
                        "position": 2
                    },
                    "id": "#samtools_merge_and_sort.cwl/bams"
                },
                {
                    "type": "boolean",
                    "default": false,
                    "inputBinding": {
                        "prefix": "-n",
                        "position": 7
                    },
                    "id": "#samtools_merge_and_sort.cwl/name_sort"
                },
                {
                    "type": "string",
                    "default": "merged_reads.bam",
                    "id": "#samtools_merge_and_sort.cwl/output_name"
                },
                {
                    "type": "int",
                    "default": 1,
                    "id": "#samtools_merge_and_sort.cwl/threads"
                }
            ],
            "outputs": [
                {
                    "id": "#samtools_merge_and_sort.cwl/bam_merged",
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.output_name)"
                    }
                }
            ],
            "id": "#samtools_merge_and_sort.cwl"
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
                    "id": "#main/bams"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#main/bismark_ignore"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#main/bismark_ignore_3prime"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#main/bismark_ignore_3prime_r2"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#main/bismark_ignore_r2"
                },
                {
                    "type": "boolean",
                    "default": true,
                    "id": "#main/bismark_no_overlap"
                },
                {
                    "type": "boolean",
                    "default": false,
                    "id": "#main/bismark_pbat"
                },
                {
                    "type": "Directory",
                    "id": "#main/genome"
                },
                {
                    "type": "int",
                    "default": 16,
                    "id": "#main/threads"
                }
            ],
            "steps": [
                {
                    "run": "#bismark_methylation_extractor.cwl",
                    "in": [
                        {
                            "source": "#main/remove_duplicates/dedup_reads",
                            "id": "#main/extract_methylation/aligned_reads"
                        },
                        {
                            "source": "#main/genome",
                            "id": "#main/extract_methylation/genome"
                        },
                        {
                            "source": "#main/bismark_ignore",
                            "id": "#main/extract_methylation/ignore"
                        },
                        {
                            "source": "#main/bismark_ignore_3prime",
                            "id": "#main/extract_methylation/ignore_3prime"
                        },
                        {
                            "source": "#main/bismark_ignore_3prime_r2",
                            "id": "#main/extract_methylation/ignore_3prime_r2"
                        },
                        {
                            "source": "#main/bismark_ignore_r2",
                            "id": "#main/extract_methylation/ignore_r2"
                        },
                        {
                            "source": "#main/bismark_no_overlap",
                            "id": "#main/extract_methylation/no_overlap"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/extract_methylation/threads"
                        }
                    ],
                    "out": [
                        "#main/extract_methylation/methylation_calls_bedgraph",
                        "#main/extract_methylation/methylation_calls_bismark",
                        "#main/extract_methylation/mbias_report",
                        "#main/extract_methylation/splitting_report",
                        "#main/extract_methylation/genome_wide_methylation_report",
                        "#main/extract_methylation/context_specific_methylation_reports"
                    ],
                    "id": "#main/extract_methylation"
                },
                {
                    "run": "#samtools_merge_and_sort.cwl",
                    "in": [
                        {
                            "source": "#main/bams",
                            "id": "#main/merge_and_sort/bams"
                        },
                        {
                            "valueFrom": "$(true)",
                            "id": "#main/merge_and_sort/name_sort"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/merge_and_sort/threads"
                        }
                    ],
                    "out": [
                        "#main/merge_and_sort/bam_merged"
                    ],
                    "id": "#main/merge_and_sort"
                },
                {
                    "doc": "samtools flagstat\n",
                    "run": "#samtools_flagstat.cwl",
                    "in": [
                        {
                            "source": "#main/remove_duplicates/dedup_reads",
                            "id": "#main/qc_post_mapping/bam"
                        }
                    ],
                    "out": [
                        "#main/qc_post_mapping/flagstat_output"
                    ],
                    "id": "#main/qc_post_mapping"
                },
                {
                    "run": "#bismark_deduplicate.cwl",
                    "in": [
                        {
                            "source": "#main/merge_and_sort/bam_merged",
                            "id": "#main/remove_duplicates/aligned_reads"
                        }
                    ],
                    "out": [
                        "#main/remove_duplicates/dedup_reads",
                        "#main/remove_duplicates/log"
                    ],
                    "id": "#main/remove_duplicates"
                }
            ],
            "outputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputSource": "#main/extract_methylation/context_specific_methylation_reports",
                    "id": "#main/context_specific_methylation_reports"
                },
                {
                    "type": "File",
                    "outputSource": "#main/remove_duplicates/log",
                    "id": "#main/dedup_log"
                },
                {
                    "type": "File",
                    "outputSource": "#main/remove_duplicates/dedup_reads",
                    "id": "#main/dedup_reads"
                },
                {
                    "type": "File",
                    "outputSource": "#main/qc_post_mapping/flagstat_output",
                    "id": "#main/flagstats_post_mapping"
                },
                {
                    "type": "File",
                    "outputSource": "#main/extract_methylation/genome_wide_methylation_report",
                    "id": "#main/genome_wide_methylation_report"
                },
                {
                    "type": "File",
                    "outputSource": "#main/extract_methylation/mbias_report",
                    "id": "#main/mbias_report"
                },
                {
                    "type": "File",
                    "outputSource": "#main/extract_methylation/methylation_calls_bedgraph",
                    "id": "#main/methylation_calls_bedgraph"
                },
                {
                    "type": "File",
                    "outputSource": "#main/extract_methylation/methylation_calls_bismark",
                    "id": "#main/methylation_calls_bismark"
                },
                {
                    "type": "File",
                    "outputSource": "#main/extract_methylation/splitting_report",
                    "id": "#main/splitting_report"
                }
            ],
            "id": "#main"
        }
    ],
    "cwlVersion": "v1.0"
}