const PIPELINE_SELECTION_GUIDANCE = "Pipeline selection: Use 'multi' for Flex assays, multimodal experiments, feature barcoding, or multiplexed samples. Use 'count' for simple single-sample gene expression only. If unsure, ask user about assay type.";
const ANALYSIS_CONFIRMATION_INSTRUCTION = "CRITICAL: Always explicitly confirm parameters before creating analysis (can be expensive, long-running compute). " +
    "Ask about experimental setup, explain parameter choices, show final parameters in table, and ASK 'Please review and confirm these parameters before I create the analysis.'";
const ACCEPT_PAYMENT_PARAMETER = {
    anyOf: [{ type: "string" }, { type: "null" }],
    default: null,
    description: 'CRITICAL: Do NOT set to "true" unless user has EXPLICITLY confirmed they accept the charges. When set to "false", omitted, or any other value, the command will fail with a cost estimate if payment is required. Only set to "true" after user explicitly agrees to the charges shown in the error message.',
};
export const TOOL_DEFINITIONS = [
    // Core Tools
    {
        name: "get_tool_version",
        description: "Get the version of the TXG CLI tool.",
        inputSchema: {
            type: "object",
            properties: {},
        },
        annotations: {
            readOnlyHint: true,
        },
    },
    {
        name: "verify_auth",
        description: "Verify authentication for the TXG CLI tool. Returns email of the authenticated user. If authentication fails, instruct the user to update the token in Claude Desktop settings > Extensions > 10x Genomics Cloud.",
        inputSchema: {
            type: "object",
            properties: {},
        },
        annotations: {
            readOnlyHint: true,
        },
    },
    // Analysis Tools
    {
        name: "get_multi_csv_config_spec",
        description: "Get the CSV configuration specification for Cell Ranger multi analysis. IMPORTANT: Always use this tool BEFORE creating a csv file to use with multi analysis to understand the correct CSV format for 10x Cloud (which differs from Cell Ranger CLI format). Returns detailed format requirements, field descriptions, and examples",
        inputSchema: {
            type: "object",
            properties: {},
        },
        annotations: {
            readOnlyHint: true,
        },
    },
    {
        name: "create_cellranger_multi_analysis",
        description: "Creates a new Cell Ranger 'multi' analysis. REQUIRED for Flex assays, multimodal experiments, feature barcoding, multiplexed samples. IMPORTANT: if csv file is not already provided, use 'get_multi_csv_config_spec' tool FIRST to get the correct CSV format for 10x Cloud. " +
            PIPELINE_SELECTION_GUIDANCE +
            " " +
            ANALYSIS_CONFIRMATION_INSTRUCTION,
        inputSchema: {
            type: "object",
            properties: {
                analysis_name: {
                    type: "string",
                    description: "Name of the analysis to create",
                },
                csv_path: {
                    type: "string",
                    description: "Path to the multi config CSV file. CRITICAL: if you need to generate the file, you must use 'get_multi_csv_config_spec' tool to get the correct format specification and examples.",
                },
                project_id: {
                    anyOf: [{ type: "string" }, { type: "null" }],
                    default: null,
                    description: "ID of existing project to create the analysis in. Either project_id or project_name must be specified",
                },
                project_name: {
                    anyOf: [{ type: "string" }, { type: "null" }],
                    default: null,
                    description: "Name of a new project to create for this analysis (alternative to project_id)",
                },
                description: {
                    anyOf: [{ type: "string" }, { type: "null" }],
                    default: null,
                    description: "Analysis description",
                },
                accept_payment: ACCEPT_PAYMENT_PARAMETER,
            },
            required: ["analysis_name", "csv_path"],
        },
    },
    {
        name: "create_cellranger_count_analysis",
        description: "Creates a new Cell Ranger 'count' analysis. ONLY for simple single-sample gene expression without additional modalities. NOT for Flex assays. " +
            PIPELINE_SELECTION_GUIDANCE +
            " " +
            ANALYSIS_CONFIRMATION_INSTRUCTION,
        inputSchema: {
            type: "object",
            properties: {
                analysis_name: {
                    type: "string",
                    description: "Name of the analysis to create",
                },
                transcriptome: {
                    type: "string",
                    description: "Reference transcriptome to use (e.g., 'refdata-cellranger-GRCh38-2024-A', or custom reference ID). Use the 'list_custom_references' or 'list_prebuilt_references' tools to find available references.",
                },
                fastqs: {
                    type: "array",
                    items: { type: "string" },
                    description: "List of FASTQ file paths or FASTQ URIs to analyze. Can be local paths or previously uploaded FASTQs in the format txg://fastqs/<filename> or txg://fastqs/<fastq_uuid>",
                },
                project_id: {
                    anyOf: [{ type: "string" }, { type: "null" }],
                    default: null,
                    description: "ID of existing project to create the analysis in",
                },
                project_name: {
                    anyOf: [{ type: "string" }, { type: "null" }],
                    default: null,
                    description: "Name of a new project to create for this analysis (alternative to project_id)",
                },
                expect_cells: {
                    anyOf: [{ type: "integer" }, { type: "null" }],
                    default: null,
                    description: "Expected number of recovered cells, used as input to cell calling algorithm. Mutually exclusive with force_cells.",
                },
                force_cells: {
                    anyOf: [{ type: "integer" }, { type: "null" }],
                    default: null,
                    description: "Force pipeline to use this number of cells, bypassing cell calling algorithm. Minimum: 10. Mutually exclusive with expect_cells.",
                },
                chemistry: {
                    anyOf: [{ type: "string" }, { type: "null" }],
                    default: "auto",
                    description: "Assay chemistry version (default: 'auto' for automatic detection)",
                },
                cell_annotation_model: {
                    anyOf: [{ type: "string" }, { type: "null" }],
                    default: null,
                    description: "Cell annotation model to use. If not provided, does not run cell annotation. Use 'auto' for default model for the species. or specify a particular model ID from 'list_annotation_models' tool.",
                },
                include_introns: {
                    type: "boolean",
                    default: true,
                    description: "Whether to include intronic reads in count (default: true)",
                },
                create_bam: {
                    type: "boolean",
                    default: true,
                    description: "Enable or disable BAM file generation (default: true). Setting to false reduces computation time and output size.",
                },
                nosecondary: {
                    type: "boolean",
                    default: false,
                    description: "Disable secondary analysis (e.g., clustering) (default: false)",
                },
                r1_length: {
                    anyOf: [{ type: "integer" }, { type: "null" }],
                    default: null,
                    description: "Hard trim the input Read 1 to this length before analysis",
                },
                r2_length: {
                    anyOf: [{ type: "integer" }, { type: "null" }],
                    default: null,
                    description: "Hard trim the input Read 2 to this length before analysis",
                },
                accept_payment: ACCEPT_PAYMENT_PARAMETER,
            },
            required: ["analysis_name", "transcriptome", "fastqs"],
        },
    },
    {
        name: "get_aggr_csv_config_spec",
        description: "Get the CSV configuration specification for Cell Ranger aggr (aggregation) analysis. IMPORTANT: Always use this tool BEFORE creating a csv file to use with aggr analysis to understand the correct CSV format for 10x Cloud. Returns detailed format requirements, field descriptions, and examples",
        inputSchema: {
            type: "object",
            properties: {},
        },
        annotations: {
            readOnlyHint: true,
        },
    },
    {
        name: "create_cellranger_aggr_analysis",
        description: "Creates a new Cell Ranger 'aggr' analysis to aggregate multiple runs. " +
            ANALYSIS_CONFIRMATION_INSTRUCTION,
        inputSchema: {
            type: "object",
            properties: {
                analysis_name: {
                    type: "string",
                    description: "Name of the aggregation analysis",
                },
                csv_path: {
                    type: "string",
                    description: "Path to aggr config CSV file. CRITICAL: if you need to generate the file, you must use 'get_aggr_csv_config_spec' tool to get the correct format specification and examples.",
                },
                project_id: {
                    anyOf: [{ type: "string" }, { type: "null" }],
                    default: null,
                    description: "ID of existing project to create the analysis in",
                },
                project_name: {
                    anyOf: [{ type: "string" }, { type: "null" }],
                    default: null,
                    description: "Name of a new project to create for this analysis (alternative to project_id)",
                },
                normalize: {
                    anyOf: [{ type: "string" }, { type: "null" }],
                    default: "mapped",
                    description: "Normalization method: 'mapped' (default) or 'none'",
                },
                description: {
                    anyOf: [{ type: "string" }, { type: "null" }],
                    default: null,
                    description: "Description for the aggregation analysis",
                },
                accept_payment: ACCEPT_PAYMENT_PARAMETER,
            },
            required: ["analysis_name", "csv_path"],
        },
    },
    {
        name: "list_analyses",
        description: "Lists all analyses in a specific project.",
        inputSchema: {
            type: "object",
            properties: {
                project_id: {
                    type: "string",
                    description: "Project ID to list analyses from",
                },
            },
            required: ["project_id"],
        },
        annotations: {
            readOnlyHint: true,
        },
    },
    {
        name: "get_analysis_details",
        description: "Shows details about a single analysis.",
        inputSchema: {
            type: "object",
            properties: {
                analysis_id: {
                    type: "string",
                    description: "Analysis ID to get details for",
                },
            },
            required: ["analysis_id"],
        },
        annotations: {
            readOnlyHint: true,
        },
    },
    {
        name: "list_analysis_files",
        description: "Lists all files within a single analysis.",
        inputSchema: {
            type: "object",
            properties: {
                analysis_id: {
                    type: "string",
                    description: "Analysis ID to list files from",
                },
            },
            required: ["analysis_id"],
        },
        annotations: {
            readOnlyHint: true,
        },
    },
    {
        name: "download_analysis_files",
        description: "Downloads all files from an analysis to a specified local path.",
        inputSchema: {
            type: "object",
            properties: {
                analysis_id: {
                    type: "string",
                    description: "Analysis ID to download files from",
                },
                output_path: {
                    type: "string",
                    description: "Local directory path where files will be saved",
                },
            },
            required: ["analysis_id", "output_path"],
        },
    },
    // Annotation Tools
    {
        name: "list_annotation_models",
        description: "Lists available cell annotation models.",
        inputSchema: {
            type: "object",
            properties: {},
        },
        annotations: {
            readOnlyHint: true,
        },
    },
    // FASTQ Tools
    {
        name: "list_fastqs",
        description: "Lists FASTQ files for a given project.",
        inputSchema: {
            type: "object",
            properties: {
                project_id: {
                    type: "string",
                    description: "Project ID to list FASTQ files from",
                },
            },
            required: ["project_id"],
        },
        annotations: {
            readOnlyHint: true,
        },
    },
    {
        name: "upload_fastqs",
        description: "Uploads FASTQ files from a local path to a project. This can take a while depending on file size and network speed.",
        inputSchema: {
            type: "object",
            properties: {
                project_id: {
                    type: "string",
                    description: "Project ID to upload FASTQ files to",
                },
                file_path: {
                    type: "string",
                    description: "Path to FASTQ file(s) or directory containing FASTQ files. Files must follow Illumina naming convention (e.g., sample_S1_L001_R1_001.fastq.gz)",
                },
            },
            required: ["project_id", "file_path"],
        },
    },
    {
        name: "list_libraries",
        description: "List all available FASTQ library types. Returns both user-friendly names (TEXT) and internal identifiers (NAME) for each library type.",
        inputSchema: {
            type: "object",
            properties: {},
        },
        annotations: {
            readOnlyHint: true,
        },
    },
    {
        name: "set_library",
        description: "Set the library type for FASTQ sets. Use the internal NAME identifier (e.g., 'singleThreeExpression', 'antibodyCapture') not the user-friendly TEXT description. When discussing with users, reference the friendly names like 'Next GEM 3' Gene Expression v3' or 'Antibody Capture', but always pass the corresponding NAME value to this tool.",
        inputSchema: {
            type: "object",
            properties: {
                project_id: {
                    type: "string",
                    description: "Project ID containing the FASTQ set",
                },
                fastq_id: { type: "string", description: "FASTQ set ID to update" },
                library_type: {
                    type: "string",
                    description: "Internal library type identifier (NAME field from list_libraries). Examples: 'singleThreeExpression' for 'Next GEM 3' Gene Expression v3', 'antibodyCapture' for 'Antibody Capture', 'gemXSingleThreeExpression' for 'GEM-X 3' Gene Expression v4'. Always use the NAME value, not the TEXT description.",
                },
            },
            required: ["project_id", "fastq_id", "library_type"],
        },
    },
    // File Management Tools
    {
        name: "list_project_files",
        description: "Lists files for a given project.",
        inputSchema: {
            type: "object",
            properties: {
                project_id: {
                    type: "string",
                    description: "Project ID to list files from",
                },
            },
            required: ["project_id"],
        },
        annotations: {
            readOnlyHint: true,
        },
    },
    {
        name: "upload_project_file",
        description: "Uploads a local file to a project.",
        inputSchema: {
            type: "object",
            properties: {
                project_id: {
                    type: "string",
                    description: "Project ID to upload file to",
                },
                file_path: {
                    type: "string",
                    description: "Path to the local file to upload",
                },
            },
            required: ["project_id", "file_path"],
        },
    },
    {
        name: "download_project_file",
        description: "Downloads a specific file from a project to a local path.",
        inputSchema: {
            type: "object",
            properties: {
                project_id: {
                    type: "string",
                    description: "Project ID containing the file",
                },
                file_id: {
                    type: "string",
                    description: "IDs of the files to download. Use the 'list_project_files' tool to find file IDs.",
                },
                output_path: {
                    type: "string",
                    description: "Local directory path where file will be saved",
                },
            },
            required: ["project_id", "file_id", "output_path"],
        },
    },
    // Project Management Tools
    {
        name: "list_projects",
        description: "Lists all available projects.",
        inputSchema: {
            type: "object",
            properties: {},
        },
        annotations: {
            readOnlyHint: true,
        },
    },
    {
        name: "create_project",
        description: "Creates a new project with a given name and optional description.",
        inputSchema: {
            type: "object",
            properties: {
                name: { type: "string", description: "Name for the new project" },
                description: {
                    anyOf: [{ type: "string" }, { type: "null" }],
                    default: null,
                    description: "Optional description for the project",
                },
            },
            required: ["name"],
        },
    },
    {
        name: "update_project",
        description: "Updates the name and/or description of an existing project.",
        inputSchema: {
            type: "object",
            properties: {
                project_id: { type: "string", description: "Project ID to update" },
                name: {
                    anyOf: [{ type: "string" }, { type: "null" }],
                    default: null,
                    description: "New name for the project",
                },
                description: {
                    anyOf: [{ type: "string" }, { type: "null" }],
                    default: null,
                    description: "New description for the project",
                },
            },
            required: ["project_id"],
        },
    },
    // Reference Management Tools
    {
        name: "list_custom_references",
        description: "Lists all custom references.",
        inputSchema: {
            type: "object",
            properties: {},
        },
        annotations: {
            readOnlyHint: true,
        },
    },
    {
        name: "list_prebuilt_references",
        description: "Lists all prebuilt references.",
        inputSchema: {
            type: "object",
            properties: {},
        },
        annotations: {
            readOnlyHint: true,
        },
    },
    {
        name: "get_reference",
        description: "Shows details for a specific custom reference.",
        inputSchema: {
            type: "object",
            properties: {
                reference_id: {
                    type: "string",
                    description: "Reference ID to get details for",
                },
            },
            required: ["reference_id"],
        },
        annotations: {
            readOnlyHint: true,
        },
    },
    {
        name: "update_reference",
        description: "Updates the name and/or description for a custom reference.",
        inputSchema: {
            type: "object",
            properties: {
                reference_id: { type: "string", description: "Reference ID to update" },
                name: {
                    anyOf: [{ type: "string" }, { type: "null" }],
                    default: null,
                    description: "New name for the reference",
                },
                description: {
                    anyOf: [{ type: "string" }, { type: "null" }],
                    default: null,
                    description: "New description for the reference",
                },
                organism: {
                    anyOf: [{ type: "string" }, { type: "null" }],
                    default: null,
                    description: "Scientific name of the organism followed by the common name in parentheses (optional)",
                },
            },
            required: ["reference_id"],
        },
    },
    {
        name: "upload_reference",
        description: "Uploads a custom reference from a local file path. Depending on size this can take a while.",
        inputSchema: {
            type: "object",
            properties: {
                file_path: {
                    type: "string",
                    description: "Path to reference file: can be a cellranger mkref/mkvdjref folder, .tar.gz file, or feature/probeset .CSV file",
                },
                name: { type: "string", description: "Name for the custom reference" },
            },
            required: ["file_path", "name"],
        },
    },
];
