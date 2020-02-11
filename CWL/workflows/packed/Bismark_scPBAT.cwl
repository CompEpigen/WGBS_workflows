{
    "$graph": [
        {
            "class": "CommandLineTool",
            "requirements": {
                "InlineJavascriptRequirement": {},
                "StepInputExpressionRequirement": {}
            },
            "hints": {
                "ResourceRequirement": {
                    "coresMin": "$(inputs.threads)",
                    "ramMin": "${return(Math.ceil(inputs.threads/5)*14000)}"
                },
                "DockerRequirement": {
                    "dockerPull": "quay.io/biocontainers/bismark:0.22.3--0"
                }
            },
            "baseCommand": "bismark",
            "arguments": [
                {
                    "valueFrom": "--bam",
                    "position": 1
                },
                {
                    "valueFrom": "--non_directional",
                    "position": 1
                }
            ],
            "inputs": {
                "genome": {
                    "type": "Directory",
                    "inputBinding": {
                        "prefix": "--genome",
                        "position": 10
                    }
                },
                "reads": {
                    "type": "File",
                    "inputBinding": {
                        "position": 11
                    }
                },
                "threads": {
                    "type": "int",
                    "default": 1,
                    "inputBinding": {
                        "prefix": "--multicore",
                        "valueFrom": "${return(Math.ceil(self/5))}",
                        "position": 1
                    }
                }
            },
            "outputs": {
                "aligned_reads": {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.bam"
                    }
                },
                "log": {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.txt"
                    }
                }
            },
            "id": "#bismark_align_scPBAT.cwl"
        },
        {
            "class": "CommandLineTool",
            "hints": {
                "ResourceRequirement": {
                    "coresMin": 4,
                    "ramMin": 20000,
                    "tmpdirMin": 30000
                },
                "DockerRequirement": {
                    "dockerPull": "kerstenbreuer/bismark:0.22.2"
                }
            },
            "baseCommand": "deduplicate_bismark",
            "arguments": [
                {
                    "valueFrom": "${\n  if (inputs.paired_end){\n    return(\"-p\")\n  }\n  else {\n    return(\"-s\")\n  }\n}\n"
                }
            ],
            "inputs": {
                "aligned_reads": {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--bam",
                        "position": 10
                    }
                },
                "paired_end": {
                    "type": "boolean",
                    "default": true
                }
            },
            "outputs": {
                "dedup_reads": {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.deduplicated.bam"
                    }
                },
                "log": {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.deduplication_report.txt"
                    }
                }
            },
            "id": "#bismark_deduplicate.cwl"
        },
        {
            "class": "CommandLineTool",
            "requirements": {
                "InlineJavascriptRequirement": {},
                "StepInputExpressionRequirement": {}
            },
            "hints": {
                "ResourceRequirement": {
                    "coresMin": "$(inputs.threads)",
                    "ramMin": "${return(Math.ceil(inputs.threads/5)*14)}",
                    "tmpdirMin": 30000
                },
                "DockerRequirement": {
                    "dockerPull": "quay.io/biocontainers/bismark:0.22.3--0"
                }
            },
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
            "inputs": {
                "aligned_reads": {
                    "type": "File",
                    "inputBinding": {
                        "position": 10
                    }
                },
                "genome": {
                    "type": "Directory",
                    "inputBinding": {
                        "prefix": "--genome",
                        "position": 10
                    }
                },
                "no_overlap": {
                    "type": "boolean",
                    "default": true,
                    "inputBinding": {
                        "prefix": "--no_overlap",
                        "position": 2
                    }
                },
                "paired_end": {
                    "type": "boolean",
                    "default": true
                },
                "ignore": {
                    "type": "int?",
                    "inputBinding": {
                        "prefix": "--ignore",
                        "position": 3
                    }
                },
                "ignore_r2": {
                    "type": "int?",
                    "inputBinding": {
                        "prefix": "--ignore_r2",
                        "position": 3
                    }
                },
                "ignore_3prime": {
                    "type": "int?",
                    "inputBinding": {
                        "prefix": "--ignore_3prime",
                        "position": 3
                    }
                },
                "ignore_3prime_r2": {
                    "type": "int?",
                    "inputBinding": {
                        "prefix": "--ignore_3prime_r2",
                        "position": 3
                    }
                },
                "threads": {
                    "type": "int",
                    "default": 1,
                    "inputBinding": {
                        "prefix": "--multicore",
                        "valueFrom": "${return(Math.ceil(self/5))}",
                        "position": 1
                    }
                }
            },
            "outputs": {
                "methylation_calls_bedgraph": {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*bedGraph.gz"
                    }
                },
                "methylation_calls_bismark": {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*bismark.cov.gz"
                    }
                },
                "mbias_report": {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.M-bias.txt*"
                    }
                },
                "splitting_report": {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*splitting_report.txt*"
                    }
                },
                "genome_wide_methylation_report": {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*CpG_report.txt*"
                    }
                },
                "context_specific_methylation_reports": {
                    "type": "File[]",
                    "outputBinding": {
                        "glob": "C*_O*.txt*"
                    }
                }
            },
            "id": "#bismark_methylation_extractor.cwl"
        },
        {
            "class": "CommandLineTool",
            "hints": {
                "ResourceRequirement": {
                    "coresMin": 4,
                    "ramMin": 20000,
                    "tmpdirMin": 10000
                },
                "DockerRequirement": {
                    "dockerPull": "quay.io/biocontainers/bismark:0.22.3--0"
                }
            },
            "baseCommand": "bismark2report",
            "inputs": {
                "alignment_report": {
                    "type": "File?",
                    "inputBinding": {
                        "prefix": "--alignment_report",
                        "position": 1
                    }
                },
                "dedup_report": {
                    "type": "File?",
                    "inputBinding": {
                        "prefix": "--dedup_report",
                        "position": 1
                    }
                },
                "splitting_report": {
                    "type": "File?",
                    "inputBinding": {
                        "prefix": "--splitting_report",
                        "position": 1
                    }
                },
                "mbias_report": {
                    "type": "File?",
                    "inputBinding": {
                        "prefix": "--mbias_report",
                        "position": 1
                    }
                }
            },
            "outputs": {
                "report": {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.html"
                    }
                }
            },
            "id": "#bismark_report.cwl"
        },
        {
            "class": "CommandLineTool",
            "hints": {
                "ResourceRequirement": {
                    "coresMin": 1,
                    "ramMin": 5000,
                    "tmpdirMin": 10000
                },
                "DockerRequirement": {
                    "dockerPull": "kerstenbreuer/trim_galore:0.6.4_2.6_0.11.8_scPBAT"
                }
            },
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
            "inputs": {
                "reads": {
                    "type": "File?",
                    "inputBinding": {
                        "position": 1
                    }
                },
                "bam": {
                    "type": "File?",
                    "inputBinding": {
                        "position": 1
                    }
                }
            },
            "outputs": {
                "fastqc_zip": {
                    "doc": "all data e.g. figures",
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputBinding": {
                        "glob": "*_fastqc.zip"
                    }
                },
                "fastqc_html": {
                    "doc": "html report showing results from zip",
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputBinding": {
                        "glob": "*_fastqc.html"
                    }
                }
            },
            "id": "#fastqc_scPBAT.cwl"
        },
        {
            "class": "CommandLineTool",
            "hints": {
                "ResourceRequirement": {
                    "coresMin": 4,
                    "ramMin": 15000,
                    "tmpdirMin": 10000
                },
                "DockerRequirement": {
                    "dockerPull": "kerstenbreuer/samtools:1.7"
                }
            },
            "baseCommand": [
                "samtools",
                "flagstat"
            ],
            "stdout": "$(inputs.bam.nameroot + inputs.output_suffix)",
            "inputs": {
                "bam": {
                    "type": "File",
                    "inputBinding": {
                        "position": 2
                    }
                },
                "output_suffix": {
                    "type": "string",
                    "default": ".flagStat"
                }
            },
            "outputs": {
                "flagstat_output": {
                    "type": "stdout"
                }
            },
            "id": "#samtools_flagstat.cwl"
        },
        {
            "class": "CommandLineTool",
            "requirements": {
                "ShellCommandRequirement": {}
            },
            "hints": {
                "ResourceRequirement": {
                    "coresMin": "$(inputs.threads)",
                    "ramMin": 20000,
                    "tmpdirMin": 30000
                },
                "DockerRequirement": {
                    "dockerPull": "kerstenbreuer/samtools:1.7"
                }
            },
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
            "inputs": {
                "output_name": {
                    "type": "string",
                    "default": "merged_reads.bam"
                },
                "bams": {
                    "doc": "bam files to be merged",
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "inputBinding": {
                        "position": 2
                    }
                },
                "name_sort": {
                    "type": "boolean",
                    "default": false,
                    "inputBinding": {
                        "prefix": "-n",
                        "position": 7
                    }
                },
                "threads": {
                    "type": "int",
                    "default": 1
                }
            },
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
            "class": "CommandLineTool",
            "requirements": {
                "InlineJavascriptRequirement": {}
            },
            "hints": {
                "ResourceRequirement": {
                    "coresMin": "$(Math.min(inputs.threads, 8))",
                    "ramMin": 7000,
                    "tmpdirMin": 10000
                },
                "DockerRequirement": {
                    "dockerPull": "kerstenbreuer/trim_galore:0.6.4_2.6_0.11.8_scPBAT"
                }
            },
            "baseCommand": "trim_galore",
            "arguments": [
                {
                    "valueFrom": "--gzip",
                    "position": 1
                }
            ],
            "inputs": {
                "reads": {
                    "type": "File",
                    "inputBinding": {
                        "position": 10
                    }
                },
                "quality": {
                    "type": "int",
                    "default": 20,
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--quality"
                    }
                },
                "clip_r1": {
                    "type": "int?",
                    "default": 6,
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--clip_r1"
                    }
                },
                "three_prime_clip_r1": {
                    "type": "int?",
                    "default": 6,
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--three_prime_clip_r1"
                    }
                },
                "threads": {
                    "type": "int?",
                    "inputBinding": {
                        "valueFrom": "$(Math.min(self, 8))",
                        "prefix": "--cores",
                        "position": 1
                    }
                }
            },
            "outputs": {
                "reads_trimmed": {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*trimmed.fq.gz"
                    }
                },
                "log": {
                    "type": "File[]",
                    "outputBinding": {
                        "glob": "*trimming_report.txt"
                    }
                }
            },
            "id": "#trim_galore_scPBAT.cwl"
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
                    "type": "Directory",
                    "id": "#main/genome"
                },
                {
                    "doc": "Even though seqeuncing has been done in PE mode, trimming and alignment are performed \nin SE-mode. Please specify a list of all fastq files belonging to one sample \n(first and second read, multiple lanes).\n",
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "id": "#main/reads"
                },
                {
                    "type": "int",
                    "default": 16,
                    "id": "#main/threads"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "default": 6,
                    "id": "#main/trim_galore_clip_r1"
                },
                {
                    "type": "int",
                    "default": 20,
                    "id": "#main/trim_galore_quality"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "default": 6,
                    "id": "#main/trim_galore_three_prime_clip_r1"
                }
            ],
            "steps": [
                {
                    "scatter": [
                        "#main/align/reads"
                    ],
                    "scatterMethod": "dotproduct",
                    "run": "#bismark_align_scPBAT.cwl",
                    "in": [
                        {
                            "source": "#main/genome",
                            "id": "#main/align/genome"
                        },
                        {
                            "source": "#main/trim/reads_trimmed",
                            "id": "#main/align/reads"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/align/threads"
                        }
                    ],
                    "out": [
                        "#main/align/aligned_reads",
                        "#main/align/log"
                    ],
                    "id": "#main/align"
                },
                {
                    "scatter": [
                        "#main/bismark_report/alignment_report"
                    ],
                    "run": "#bismark_report.cwl",
                    "in": [
                        {
                            "source": "#main/align/log",
                            "id": "#main/bismark_report/alignment_report"
                        },
                        {
                            "source": "#main/remove_duplicates/log",
                            "id": "#main/bismark_report/dedup_report"
                        },
                        {
                            "source": "#main/extract_methylation/mbias_report",
                            "id": "#main/bismark_report/mbias_report"
                        },
                        {
                            "source": "#main/extract_methylation/splitting_report",
                            "id": "#main/bismark_report/splitting_report"
                        }
                    ],
                    "out": [
                        "#main/bismark_report/report"
                    ],
                    "id": "#main/bismark_report"
                },
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
                            "valueFrom": "${return(false)}",
                            "id": "#main/extract_methylation/no_overlap"
                        },
                        {
                            "valueFrom": "${return(false)}",
                            "id": "#main/extract_methylation/paired_end"
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
                            "source": "#main/align/aligned_reads",
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
                    "scatter": [
                        "#main/qc_posttrim/reads"
                    ],
                    "scatterMethod": "dotproduct",
                    "run": "#fastqc_scPBAT.cwl",
                    "in": [
                        {
                            "source": "#main/trim/reads_trimmed",
                            "id": "#main/qc_posttrim/reads"
                        }
                    ],
                    "out": [
                        "#main/qc_posttrim/fastqc_zip",
                        "#main/qc_posttrim/fastqc_html"
                    ],
                    "id": "#main/qc_posttrim"
                },
                {
                    "scatter": [
                        "#main/qc_pretrim/reads"
                    ],
                    "scatterMethod": "dotproduct",
                    "run": "#fastqc_scPBAT.cwl",
                    "in": [
                        {
                            "source": "#main/reads",
                            "id": "#main/qc_pretrim/reads"
                        }
                    ],
                    "out": [
                        "#main/qc_pretrim/fastqc_zip",
                        "#main/qc_pretrim/fastqc_html"
                    ],
                    "id": "#main/qc_pretrim"
                },
                {
                    "run": "#bismark_deduplicate.cwl",
                    "in": [
                        {
                            "source": "#main/merge_and_sort/bam_merged",
                            "id": "#main/remove_duplicates/aligned_reads"
                        },
                        {
                            "valueFrom": "${return(false)}",
                            "id": "#main/remove_duplicates/paired_end"
                        }
                    ],
                    "out": [
                        "#main/remove_duplicates/dedup_reads",
                        "#main/remove_duplicates/log"
                    ],
                    "id": "#main/remove_duplicates"
                },
                {
                    "scatter": [
                        "#main/trim/reads"
                    ],
                    "scatterMethod": "dotproduct",
                    "run": "#trim_galore_scPBAT.cwl",
                    "in": [
                        {
                            "source": "#main/trim_galore_clip_r1",
                            "id": "#main/trim/clip_r1"
                        },
                        {
                            "source": "#main/trim_galore_quality",
                            "id": "#main/trim/quality"
                        },
                        {
                            "source": "#main/reads",
                            "id": "#main/trim/reads"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/trim/threads"
                        },
                        {
                            "source": "#main/trim_galore_three_prime_clip_r1",
                            "id": "#main/trim/three_prime_clip_r1"
                        }
                    ],
                    "out": [
                        "#main/trim/log",
                        "#main/trim/reads_trimmed"
                    ],
                    "id": "#main/trim"
                }
            ],
            "outputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputSource": "#main/align/log",
                    "id": "#main/align_log"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputSource": "#main/bismark_report/report",
                    "id": "#main/bismark_report_html"
                },
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
                    "type": {
                        "type": "array",
                        "items": {
                            "type": "array",
                            "items": "File"
                        }
                    },
                    "outputSource": "#main/qc_posttrim/fastqc_html",
                    "id": "#main/qc_posttrim_fastqc_html"
                },
                {
                    "type": {
                        "type": "array",
                        "items": {
                            "type": "array",
                            "items": "File"
                        }
                    },
                    "outputSource": "#main/qc_posttrim/fastqc_zip",
                    "id": "#main/qc_posttrim_fastqc_zip"
                },
                {
                    "type": {
                        "type": "array",
                        "items": {
                            "type": "array",
                            "items": "File"
                        }
                    },
                    "outputSource": "#main/qc_pretrim/fastqc_html",
                    "id": "#main/qc_pretrim_fastqc_html"
                },
                {
                    "type": {
                        "type": "array",
                        "items": {
                            "type": "array",
                            "items": "File"
                        }
                    },
                    "outputSource": "#main/qc_pretrim/fastqc_zip",
                    "id": "#main/qc_pretrim_fastqc_zip"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputSource": "#main/trim/reads_trimmed",
                    "id": "#main/reads_trimmed"
                },
                {
                    "type": "File",
                    "outputSource": "#main/extract_methylation/splitting_report",
                    "id": "#main/splitting_report"
                },
                {
                    "type": {
                        "type": "array",
                        "items": {
                            "type": "array",
                            "items": "File"
                        }
                    },
                    "outputSource": "#main/trim/log",
                    "id": "#main/trim_log"
                }
            ],
            "id": "#main"
        }
    ],
    "cwlVersion": "v1.0"
}