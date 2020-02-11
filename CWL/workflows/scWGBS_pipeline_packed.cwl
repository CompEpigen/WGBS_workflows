{
    "$graph": [
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "ShellCommandRequirement"
                }
            ],
            "hints": [
                {
                    "dockerPull": "kerstenbreuer/fame:0.2_8cores_150bp_reads",
                    "class": "DockerRequirement"
                },
                {
                    "coresMin": 8,
                    "ramMin": 28000,
                    "class": "ResourceRequirement"
                }
            ],
            "baseCommand": [
                "fame"
            ],
            "arguments": [
                {
                    "valueFrom": "--paired"
                }
            ],
            "inputs": [
                {
                    "type": "boolean",
                    "default": true,
                    "inputBinding": {
                        "prefix": "--gzip_reads"
                    },
                    "id": "#fame_sc.cwl/gzip_reads"
                },
                {
                    "type": "File",
                    "secondaryFiles": [
                        "_strands"
                    ],
                    "inputBinding": {
                        "prefix": "--load_index"
                    },
                    "id": "#fame_sc.cwl/load_index"
                },
                {
                    "type": "string",
                    "default": "bulk_methylation_calls",
                    "inputBinding": {
                        "prefix": "--out_basename"
                    },
                    "id": "#fame_sc.cwl/out_basename"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--sc_input"
                    },
                    "id": "#fame_sc.cwl/sc_input"
                },
                {
                    "type": "string",
                    "default": "sc_methylation_calls",
                    "inputBinding": {
                        "prefix": "--sc_output"
                    },
                    "id": "#fame_sc.cwl/sc_output"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.out_basename)*"
                    },
                    "id": "#fame_sc.cwl/bulk_methylation_calls"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.sc_output)*"
                    },
                    "id": "#fame_sc.cwl/sc_methylation_calls"
                }
            ],
            "id": "#fame_sc.cwl"
        },
        {
            "class": "Workflow",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
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
                    "default": true,
                    "id": "#main/fastqs_are_gzipped"
                },
                {
                    "type": "File",
                    "secondaryFiles": [
                        "_strands"
                    ],
                    "id": "#main/reference_index"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "string"
                    },
                    "id": "#main/sc_id"
                }
            ],
            "steps": [
                {
                    "run": "#fame_sc.cwl",
                    "in": [
                        {
                            "source": "#main/fastqs_are_gzipped",
                            "id": "#main/align_and_call_meth/gzip_reads"
                        },
                        {
                            "source": "#main/reference_index",
                            "id": "#main/align_and_call_meth/load_index"
                        },
                        {
                            "source": "#main/build_meta_tsv/meta_tsv",
                            "id": "#main/align_and_call_meth/sc_input"
                        }
                    ],
                    "out": [
                        "#main/align_and_call_meth/sc_methylation_calls",
                        "#main/align_and_call_meth/bulk_methylation_calls"
                    ],
                    "id": "#main/align_and_call_meth"
                },
                {
                    "in": [
                        {
                            "source": "#main/fastq1",
                            "id": "#main/build_meta_tsv/fastq1"
                        },
                        {
                            "source": "#main/fastq2",
                            "id": "#main/build_meta_tsv/fastq2"
                        },
                        {
                            "source": "#main/sc_id",
                            "id": "#main/build_meta_tsv/sc_id"
                        }
                    ],
                    "run": {
                        "cwlVersion": "v1.0",
                        "class": "CommandLineTool",
                        "hints": [
                            {
                                "dockerPull": "python:3.7.4",
                                "class": "DockerRequirement"
                            }
                        ],
                        "requirements": [
                            {
                                "listing": [
                                    {
                                        "entryname": "build_meta_tsv.py",
                                        "entry": "import sys\nargs = sys.argv\nsc_ids = args[1].split(\",\")\nfastq1 = args[2].split(\",\")\nfastq2 = args[3].split(\",\")\nfor i, sc_id in enumerate(sc_ids):\n  print(sc_id + \"\\t\" + fastq1[i] + \"\\t\" + fastq2[i])\n"
                                    }
                                ],
                                "class": "InitialWorkDirRequirement"
                            }
                        ],
                        "baseCommand": [
                            "python3",
                            "build_meta_tsv.py"
                        ],
                        "inputs": [
                            {
                                "type": {
                                    "type": "array",
                                    "items": "File"
                                },
                                "inputBinding": {
                                    "position": 2,
                                    "itemSeparator": ","
                                },
                                "id": "#main/build_meta_tsv/b7aac25d-2377-46ec-9723-37682863ae2f/fastq1"
                            },
                            {
                                "type": {
                                    "type": "array",
                                    "items": "File"
                                },
                                "inputBinding": {
                                    "position": 3,
                                    "itemSeparator": ","
                                },
                                "id": "#main/build_meta_tsv/b7aac25d-2377-46ec-9723-37682863ae2f/fastq2"
                            },
                            {
                                "type": {
                                    "type": "array",
                                    "items": "string"
                                },
                                "inputBinding": {
                                    "position": 1,
                                    "itemSeparator": ","
                                },
                                "id": "#main/build_meta_tsv/b7aac25d-2377-46ec-9723-37682863ae2f/sc_id"
                            }
                        ],
                        "stdout": "meta.tsv",
                        "outputs": [
                            {
                                "type": "File",
                                "id": "#main/build_meta_tsv/b7aac25d-2377-46ec-9723-37682863ae2f/meta_tsv",
                                "outputBinding": {
                                    "glob": "meta.tsv"
                                }
                            }
                        ],
                        "id": "#main/build_meta_tsv/b7aac25d-2377-46ec-9723-37682863ae2f"
                    },
                    "out": [
                        "#main/build_meta_tsv/meta_tsv"
                    ],
                    "id": "#main/build_meta_tsv"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputSource": "#main/align_and_call_meth/bulk_methylation_calls",
                    "id": "#main/bulk_methylation_calls"
                },
                {
                    "type": "File",
                    "outputSource": "#main/build_meta_tsv/meta_tsv",
                    "id": "#main/meta_tsv"
                },
                {
                    "type": "File",
                    "outputSource": "#main/align_and_call_meth/sc_methylation_calls",
                    "id": "#main/sc_methylation_calls"
                }
            ],
            "id": "#main"
        }
    ],
    "cwlVersion": "v1.0"
}