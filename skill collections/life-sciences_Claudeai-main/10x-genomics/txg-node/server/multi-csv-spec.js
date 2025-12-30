export function getMultiCsvConfigSpec() {
    return {
        format: "Cell Ranger multi config CSV for 10x Cloud",
        structure: "INI-style sections with CSV data mixed in",
        notes: "Use txg:// prefix for files already uploaded to cloud, or local paths for automatic upload. Supported prefixes: txg://fastqs/<filename> or txg://fastqs/<uuid> for FASTQ files, txg://references/<custom_reference_name> for custom references, txg://files/<file_id> for general files, txg://analyses/<analysis_id> for aggregation analyses.",
        sections: {
            "[gene-expression]": {
                description: "Gene expression analysis parameters",
                fields: {
                    "reference": {
                        description: "Reference transcriptome to use (e.g., 'refdata-cellranger-GRCh38-2024-A', or custom reference ID). Use the 'list_custom_references' or 'list_prebuilt_references' tools to find available references.",
                        required: true,
                        example: "refdata-cellranger-GRCh38-2024-A"
                    },
                    "create-bam": {
                        description: "Generate BAM file. See https://10xgen.com/create-bam for additional guidance.",
                        required: true,
                        values: ["true", "false"],
                        default: "true"
                    },
                    "cell-annotation-model": {
                        description: "Cell annotation model to use. If auto, uses the default model for the species. If not given, does not run cell annotation.",
                        required: false,
                        values: ["auto", "human_pca_v1_beta", "mouse_pca_v1_beta"],
                        example: "human_pca_v1_beta"
                    },
                    "r1-length": {
                        description: "Hard trim R1 reads to this length",
                        required: false,
                        type: "integer"
                    },
                    "r2-length": {
                        description: "Hard trim R2 reads to this length",
                        required: false,
                        type: "integer"
                    },
                    "chemistry": {
                        description: "Assay configuration. Auto-detected by default (recommended). Only specify manually if auto-detection fails or for specific requirements.",
                        required: false,
                        default: "auto",
                        values: [
                            // General
                            "auto",
                            "threeprime",
                            "fiveprime",
                            // 3' chemistries
                            "SC3Pv1",
                            "SC3Pv2",
                            "SC3Pv3",
                            "SC3Pv3-polyA",
                            "SC3Pv4",
                            "SC3Pv4-polyA",
                            "SC3Pv3HT",
                            "SC3Pv3HT-polyA",
                            "SC-FB",
                            // 5' chemistries
                            "SC5P-PE",
                            "SC5P-PE-v3",
                            "SC5P-R2",
                            "SC5P-R2-v3",
                            "SC5PHT",
                            // Flex chemistries
                            "SFRP",
                            "MFRP",
                            "MFRP-R1",
                            "MFRP-RNA",
                            "MFRP-Ab",
                            "MFRP-Ab-R2pos50",
                            "MFRP-RNA-R1",
                            "MFRP-Ab-R1",
                            // Multiome
                            "ARC-v1"
                        ],
                        notes: [
                            "The -polyA suffix is only needed when manually specifying chemistry for antibody libraries with poly-A capture (e.g., TotalSeq-A)",
                            "For most cases, specifying the base chemistry (e.g., SC3Pv3) automatically maps library types to appropriate chemistries",
                            "ARC-v1 is for Gene Expression portion of Multiome data; auto-detection of ARC-v1 triggers an error"
                        ]
                    },
                    "expect-cells": {
                        description: "Expected number of recovered cells. If specified, override the pipeline’s auto-estimation of cells",
                        required: false,
                        type: "integer",
                        example: "10000"
                    },
                    "force-cells": {
                        description: "Force cell calling to this number",
                        required: false,
                        type: "integer"
                    },
                    "include-introns": {
                        description: "Include intronic reads in count",
                        required: false,
                        values: ["true", "false"],
                        default: "true"
                    },
                    "no-secondary": {
                        description: "Disable secondary analysis",
                        required: false,
                        values: ["true", "false"],
                        default: "false"
                    }
                },
                additional_fields: {
                    threePrime: {
                        description: "These [gene-expression] options only apply to 3' Cell Multiplexing data analysis.",
                        fields: {
                            "min-assignment-confidence": {
                                description: "Optional. The minimum estimated likelihood to call a sample as tagged with a Cell Multiplexing Oligo (CMO) instead of 'Unassigned'",
                                required: false,
                                type: "float",
                                default: "0.9"
                            },
                            "cmo-set": {
                                description: "Optional. Path to a custom CMO set CSV file declaring CMO constructs and associated barcodes. Default CMO reference IDs are built into Cell Ranger and don't need to be specified.",
                                required: false,
                                type: "string"
                            },
                            "barcode-sample-assignment": {
                                description: "Optional. Absolute path to a barcode-sample assignment CSV file that specifies the barcodes belonging to each sample. Also applicable to cell hashing with Antibody Capture.",
                                required: false,
                                type: "string"
                            }
                        }
                    },
                    flex: {
                        description: "These [gene-expression] options only apply to Flex data analysis.",
                        fields: {
                            "probe-set": {
                                description: "Probe set CSV file (Fixed RNA Profiling/Flex only)",
                                required: "For Flex assays",
                            },
                            "filter-probes": {
                                description: "Filter probes not in reference (Flex only)",
                                required: false,
                                values: ["true", "false"],
                                default: "true"
                            },
                        }
                    }
                }
            },
            "[feature]": {
                description: "Feature barcode analysis parameters",
                required: "When using Antibody/CRISPR libraries",
                fields: {
                    "reference": {
                        description: "Feature reference CSV file",
                        required: "Required only for Antibody Capture, Antigen Capture, or CRISPR Guide Capture libraries.",
                    },
                    "r1-length": {
                        description: "Optional. Limit the length of input Read 1 sequence of Feature Barcode libraries to the first N bases. Note: length includes the 10x Barcode and UMI sequences so do not set below 26. Useful for determining optimal read length for sequencing. Default: do not trim Read 1.",
                        required: false,
                        type: "integer",
                        minimum: 26
                    },
                    "r2-length": {
                        description: "Optional. Limit the length of input Read 2 sequence of Feature Barcode libraries to the first N bases. Trimming occurs before sequencing metrics are computed, so limiting Read 2 length may affect Q30 scores. Default: do not trim Read 2.",
                        required: false,
                        type: "integer"
                    },
                    "min-crispr-umi": {
                        description: "Optional. Set the minimum number of CRISPR guide RNA UMIs required for protospacer detection. Can be customized according to specific experimental needs. Only applicable to datasets with a CRISPR Guide Capture library.",
                        required: false,
                        type: "integer",
                        default: "3"
                    }
                }
            },
            "[libraries]": {
                description: "Input library definitions (CSV table)",
                required: true,
                format: "CSV with headers",
                columns: {
                    "fastq_id": {
                        description: "Name prefix of FASTQ files",
                        required: true,
                        example: "sample1"
                    },
                    "fastqs": {
                        description: "Path to FASTQ directory or file",
                        required: true,
                        example: "/path/to/fastqs"
                    },
                    "feature_types": {
                        description: "Library type",
                        required: true,
                        values: [
                            "Gene Expression",
                            "Antibody Capture",
                            "CRISPR Guide Capture",
                            "Multiplexing Capture",
                            "VDJ",
                            "VDJ-B",
                            "VDJ-T",
                            "VDJ-T-GD",
                            "Antigen Capture"
                        ]
                    },
                    "lanes": {
                        description: "Lanes to use",
                        required: false,
                        default: "all",
                        example: "1|2|3"
                    },
                    "physical_library_id": {
                        description: "Physical library ID for grouping",
                        required: false
                    },
                    "subsample_rate": {
                        description: "Fraction of reads to use",
                        required: false,
                        type: "float",
                        range: "0.0-1.0"
                    },
                    "chemistry": {
                        description: "Optional (only applicable to Flex). Library-specific assay configuration. Auto-detected by default (recommended). Typically users will not need to specify a chemistry. has same options as [gene-expression] chemistry field.",
                        required: false,
                        default: "auto",
                        type: "string"
                    }
                }
            },
            "[vdj]": {
                description: "V(D)J analysis parameters",
                required: "When using VDJ libraries",
                fields: {
                    "reference": {
                        description: "Path to Cell Ranger V(D)J reference",
                        required: true,
                        example: "/path/to/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0"
                    },
                    "inner-enrichment-primers": {
                        description: "Inner enrichment primers file",
                        required: false
                    },
                    "r1-length": {
                        description: "Hard trim R1 to this length",
                        required: false,
                        type: "integer"
                    },
                    "r2-length": {
                        description: "Hard trim R2 to this length",
                        required: false,
                        type: "integer"
                    }
                }
            },
            "[samples]": {
                description: "Sample definitions for multiplexing (CSV table)",
                required: "For multiplexed experiments",
                format: "CSV with headers",
                columns: {
                    "sample_id": {
                        description: "Required. A name to identify a multiplexed sample. Must be alphanumeric with hyphens and/or underscores, and less than 64 characters.",
                        required: true,
                        type: "string",
                        pattern: "^[a-zA-Z0-9_-]{1,63}$"
                    },
                    "expect_cells": {
                        description: "Optional. Override the pipeline's auto-estimation of cells with the expected number of recovered cells. For Flex, only valid for multiplex configuration. Note: column name uses underscore (_).",
                        required: false,
                        type: "integer"
                    },
                    "force_cells": {
                        description: "Optional. Force pipeline to use this number of cells, bypassing cell detection. Default: detect cells using EmptyDrops. For Flex, only valid for multiplex configuration. Note: column name uses underscore (_).",
                        required: false,
                        type: "integer"
                    },
                    "description": {
                        description: "Optional. A description for the sample.",
                        required: false,
                        type: "string"
                    },
                    "ocm_barcode_ids": {
                        description: "Required for 3' and 5' GEM-X On-chip multiplexing (OCM) only. The OCM barcode IDs used to multiplex this sample (e.g., OB1, OB2, OB3, OB4). Separate multiple IDs with pipe (e.g., OB1|OB2).",
                        required: false,
                        type: "string",
                        pattern: "^OB[1-4](\\|OB[1-4])*$"
                    },
                    "hashtag_ids": {
                        description: "Required for cell or sample hashing with Antibody Capture only. The hashtag IDs used to multiplex this sample. Separate multiple IDs with pipe (e.g., ABHT-1|ABHT-2).",
                        required: false,
                        type: "string"
                    },
                    "cmo_ids": {
                        description: "Required for 3' Cell Multiplexing only. The Cell Multiplexing oligo IDs used to multiplex this sample. Only input CMOs used in the experiment. Separate multiple IDs with pipe (e.g., CMO301|CMO302).",
                        required: false,
                        type: "string"
                    },
                    "probe_barcode_ids": {
                        description: "Required for Flex only. Fixed RNA Probe Barcode IDs, Antibody Multiplexing Barcode IDs, and CRISPR Multiplexing Barcode IDs. Recommended format when Antibody Capture present: BC001+AB001 (BC+AB, no spaces). For CRISPR: BC001+CR001. Separate multiple IDs with pipe (e.g., BC001|BC002).",
                        required: false,
                        type: "string"
                    },
                    "emptydrops_minimum_umis": {
                        description: "Optional for multiplex Flex only. Adjust UMI cutoff during second cell calling step. Only evaluates barcodes above this threshold. Recommended starting value: 100. Cannot be used with force-cells. Cell Ranger v7.1.0+.",
                        required: false,
                        type: "integer"
                    }
                }
            }
        },
        rules: [
            "Sections start with [section-name]",
            "Parameters use key,value format",
            "[libraries] and [samples] are CSV tables with headers",
            "Multiple values separated by | (pipe)",
            "Paths can be absolute or relative",
        ],
        examples: {
            "flex_gex_crispr": {
                description: "Flex with Gene Expression and CRISPR (cloud files)",
                csv: `[gene-expression]
ref,refdata-cellranger-GRCh38-2024-A
probe-set,txg://references/probe-set_reference.csv
create-bam,false

[feature]
ref,txg://references/feature_reference.csv

[libraries]
fastq_id,fastqs,lanes,feature_types
crispr_fastq,txg://fastqs/crispr_S1_L001,any,CRISPR Guide Capture
gex_fastq,txg://fastqs/gex_S1_L001,any,Gene Expression

[samples]
sample_id,probe_barcode_ids
sample1,BC001+CR001|BC002+CR002|BC003+CR003|BC004+CR004
sample2,BC005+CR005|BC006+CR006|BC007+CR007|BC008+CR008
sample3,BC009+CR009|BC010+CR010|BC011+CR011|BC012+CR012
sample4,BC013+CR013|BC014+CR014|BC015+CR015|BC016+CR016`
            },
            "flex_with_prebuilt": {
                description: "Flex using prebuilt references (auto-upload local FASTQs)",
                csv: `[gene-expression]
reference,refdata-cellranger-GRCh38-2024-A
create-bam,false
probe-set,Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2024-A.csv

[libraries]
fastq_id,fastqs,feature_types
GEX_fastq,~/path_to/fastqs/GEX_fastqs/,Gene Expression
CRISPR_fastq,~/path_to/CRISPR_fastqs/,CRISPR Guide Capture

[feature]
reference,~/path_to/references/RTL_lnRNA.csv`
            },
            "5prime_gex_crispr": {
                description: "5' Gene Expression with CRISPR",
                csv: `[gene-expression]
ref,refdata-cellranger-GRCh38-2024-A
create-bam,false

[feature]
ref,txg://references/feature_reference_5p_crispr.csv

[libraries]
fastq_id,fastqs,lanes,feature_types
gex_fastq,txg://fastqs/gex_S1_L001,1,Gene Expression
cr_fastq,txg://fastqs/crispr_S1_L003,3,CRISPR Guide Capture`
            },
            "3prime_cmo_multiplex": {
                description: "3' Cell Multiplexing with CMO",
                csv: `[gene-expression]
reference,refdata-cellranger-GRCh38-2024-A
expect-cells,20000

[feature]
reference,cmo_reference.csv

[libraries]
fastq_id,fastqs,feature_types
GEX,txg://fastqs/gex_S1_L001,Gene Expression
CMO,txg://fastqs/cmo_S2_L001,Multiplexing Capture

[samples]
sample_id,cmo_ids
Sample1,CMO301
Sample2,CMO302
Sample3,CMO303
Sample4,CMO304`
            },
            "aggr_config": {
                description: "Aggregation CSV for combining analyses",
                csv: `sample_id,molecule_h5
5p_multi,txg://analyses/<analysisID1>
5p_multi_2,txg://analyses/<analysisID2>
5p_count,txg://analyses/<analysisID3>`
            }
        }
    };
}
