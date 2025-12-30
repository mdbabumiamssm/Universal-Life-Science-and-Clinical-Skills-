import { CallToolRequestSchema } from "@modelcontextprotocol/sdk/types.js";
import { z } from "zod";
import { txgCli, LONG_OPERATION_TIMEOUT, LONG_OPERATION_SOFT_TIMEOUT, } from "./txg-cli-manager.js";
import { toResponse, toAnalysisResponse, wrapResponse } from "./middleware.js";
import { getMultiCsvConfigSpec } from "./multi-csv-spec.js";
import { getAggrCsvConfigSpec } from "./aggr-csv-spec.js";
/**
 * Add payment acceptance flag to command.
 * Only passes --accept-payment=true if explicitly set to "true" (string) or true (boolean).
 * Otherwise, passes --accept-payment=false for safety.
 */
function addAcceptPaymentFlag(command, acceptPayment) {
    if (acceptPayment === "true" || acceptPayment === true) {
        command.push("--accept-payment=true");
    }
    else {
        command.push("--accept-payment=false");
    }
}
export function registerTools(server) {
    // Register the tool call handler
    server.setRequestHandler(CallToolRequestSchema, async (request) => {
        const { name, arguments: args } = request.params;
        switch (name) {
            case "get_tool_version":
                return wrapResponse(toResponse(await txgCli.runCommand(["--version"])));
            case "verify_auth":
                return wrapResponse(toResponse(await txgCli.runCommand(["auth", "verify"])));
            // --- Analysis Tools ---
            case "get_multi_csv_config_spec": {
                const spec = getMultiCsvConfigSpec();
                return {
                    content: [
                        {
                            type: "text",
                            text: JSON.stringify(spec, null, 2),
                        },
                    ],
                };
            }
            case "get_aggr_csv_config_spec": {
                const spec = getAggrCsvConfigSpec();
                return {
                    content: [
                        {
                            type: "text",
                            text: JSON.stringify(spec, null, 2),
                        },
                    ],
                };
            }
            case "create_cellranger_multi_analysis": {
                const schema = z.object({
                    analysis_name: z.string().describe("Name of the analysis to create"),
                    csv_path: z.string().describe("Path to the multi config CSV file"),
                    project_id: z.string().optional().describe("ID of existing project"),
                    project_name: z.string().optional().describe("Name of a new project"),
                    description: z.string().optional(),
                    accept_payment: z.union([z.string(), z.boolean()]).optional(),
                });
                const params = schema.parse(args);
                const command = [
                    "analyses",
                    "create",
                    "cellranger",
                    "multi",
                    "--product-version",
                    "latest",
                    "--analysis-name",
                    params.analysis_name,
                    "--csv",
                    params.csv_path,
                    "--wait-completion=false",
                ];
                if (params.project_id) {
                    command.push("--project-id", params.project_id);
                }
                else if (params.project_name) {
                    command.push("--project-name", params.project_name);
                }
                if (params.description) {
                    command.push("--description", params.description);
                }
                addAcceptPaymentFlag(command, params.accept_payment);
                return wrapResponse(toAnalysisResponse(await txgCli.runCommand(command)));
            }
            case "create_cellranger_count_analysis": {
                const schema = z.object({
                    analysis_name: z.string().describe("Name of the analysis to create"),
                    transcriptome: z.string().describe("Reference transcriptome to use"),
                    fastqs: z
                        .array(z.string())
                        .describe("List of FASTQ file paths or IDs"),
                    project_id: z.string().optional(),
                    project_name: z.string().optional(),
                    expect_cells: z
                        .union([z.number(), z.string().transform((v) => parseInt(v, 10))])
                        .optional(),
                    force_cells: z
                        .union([z.number(), z.string().transform((v) => parseInt(v, 10))])
                        .optional(),
                    chemistry: z.string().optional().default("auto"),
                    include_introns: z.union([z.string(), z.boolean()]).optional(),
                    cell_annotation_model: z.string().optional(),
                    create_bam: z.union([z.string(), z.boolean()]).optional(),
                    nosecondary: z.union([z.string(), z.boolean()]).optional(),
                    r1_length: z
                        .union([z.number(), z.string().transform((v) => parseInt(v, 10))])
                        .optional(),
                    r2_length: z
                        .union([z.number(), z.string().transform((v) => parseInt(v, 10))])
                        .optional(),
                    accept_payment: z.union([z.string(), z.boolean()]).optional(),
                });
                const params = schema.parse(args);
                const command = [
                    "analyses",
                    "create",
                    "cellranger",
                    "count",
                    "--product-version",
                    "latest",
                    "--analysis-name",
                    params.analysis_name,
                    "--transcriptome",
                    params.transcriptome,
                    "--wait-completion=false",
                ];
                for (const fq of params.fastqs) {
                    command.push("--fastqs", fq);
                }
                if (params.project_id) {
                    command.push("--project-id", params.project_id);
                }
                else if (params.project_name) {
                    command.push("--project-name", params.project_name);
                }
                if (params.expect_cells) {
                    command.push("--expect-cells", params.expect_cells.toString());
                }
                if (params.force_cells) {
                    command.push("--force-cells", params.force_cells.toString());
                }
                if (params.chemistry) {
                    command.push("--chemistry", params.chemistry);
                }
                if (params.include_introns !== undefined) {
                    if (params.include_introns === true ||
                        params.include_introns === "true" ||
                        params.include_introns === "yes") {
                        command.push("--include-introns=true");
                    }
                    else if (params.include_introns === false ||
                        params.include_introns === "false" ||
                        params.include_introns === "no") {
                        command.push("--include-introns=false");
                    }
                }
                if (params.create_bam !== undefined) {
                    if (params.create_bam === true ||
                        params.create_bam === "true" ||
                        params.create_bam === "yes") {
                        command.push("--create-bam=true");
                    }
                    else if (params.create_bam === false ||
                        params.create_bam === "false" ||
                        params.create_bam === "no") {
                        command.push("--create-bam=false");
                    }
                }
                if (params.nosecondary === true ||
                    params.nosecondary === "true" ||
                    params.nosecondary === "yes") {
                    command.push("--nosecondary");
                }
                if (params.r1_length) {
                    command.push("--r1-length", params.r1_length.toString());
                }
                if (params.r2_length) {
                    command.push("--r2-length", params.r2_length.toString());
                }
                if (params.cell_annotation_model) {
                    command.push("--cell-annotation-model", params.cell_annotation_model);
                }
                addAcceptPaymentFlag(command, params.accept_payment);
                return wrapResponse(toAnalysisResponse(await txgCli.runCommand(command)));
            }
            case "create_cellranger_aggr_analysis": {
                const schema = z.object({
                    analysis_name: z.string(),
                    csv_path: z.string(),
                    project_id: z.string().optional(),
                    project_name: z.string().optional(),
                    normalize: z.string().optional().default("mapped"),
                    description: z.string().optional(),
                    accept_payment: z.union([z.string(), z.boolean()]).optional(),
                });
                const params = schema.parse(args);
                const command = [
                    "analyses",
                    "create",
                    "cellranger",
                    "aggr",
                    "--product-version",
                    "latest",
                    "--analysis-name",
                    params.analysis_name,
                    "--csv",
                    params.csv_path,
                    "--wait-completion=false",
                ];
                if (params.project_id) {
                    command.push("--project-id", params.project_id);
                }
                else if (params.project_name) {
                    command.push("--project-name", params.project_name);
                }
                if (params.normalize) {
                    command.push("--normalize", params.normalize);
                }
                if (params.description) {
                    command.push("--description", params.description);
                }
                addAcceptPaymentFlag(command, params.accept_payment);
                return wrapResponse(toAnalysisResponse(await txgCli.runCommand(command)));
            }
            case "list_analyses": {
                const schema = z.object({
                    project_id: z.string().describe("Project ID to list analyses from"),
                });
                const params = schema.parse(args);
                return wrapResponse(toResponse(await txgCli.runCommand(["analyses", "list", params.project_id])));
            }
            case "get_analysis_details": {
                const schema = z.object({
                    analysis_id: z.string().describe("Analysis ID to get details for"),
                });
                const params = schema.parse(args);
                return wrapResponse(toResponse(await txgCli.runCommand(["analyses", "get", params.analysis_id])));
            }
            case "list_analysis_files": {
                const schema = z.object({
                    analysis_id: z.string(),
                });
                const params = schema.parse(args);
                return wrapResponse(toResponse(await txgCli.runCommand(["analyses", "files", params.analysis_id])));
            }
            case "download_analysis_files": {
                const schema = z.object({
                    analysis_id: z.string(),
                    output_path: z.string(),
                });
                const params = schema.parse(args);
                return wrapResponse(toResponse(await txgCli.runCommand([
                    "analyses",
                    "download",
                    params.analysis_id,
                    "--target-dir",
                    params.output_path,
                ], LONG_OPERATION_TIMEOUT, LONG_OPERATION_SOFT_TIMEOUT)));
            }
            // --- Annotation Tools ---
            case "list_annotation_models":
                return wrapResponse(toResponse(await txgCli.runCommand(["annotation", "models", "list"])));
            // --- FASTQ Tools ---
            case "list_fastqs": {
                const schema = z.object({
                    project_id: z.string(),
                });
                const params = schema.parse(args);
                return wrapResponse(toResponse(await txgCli.runCommand(["fastqs", "list", params.project_id])));
            }
            case "upload_fastqs": {
                const schema = z.object({
                    project_id: z.string(),
                    file_path: z.string(),
                });
                const params = schema.parse(args);
                return wrapResponse(toResponse(await txgCli.runCommand([
                    "fastqs",
                    "upload",
                    "--project-id",
                    params.project_id,
                    params.file_path,
                ], LONG_OPERATION_TIMEOUT, LONG_OPERATION_SOFT_TIMEOUT)));
            }
            case "list_libraries": {
                return wrapResponse(toResponse(await txgCli.runCommand(["fastqs", "library"])));
            }
            case "set_library": {
                const schema = z.object({
                    project_id: z.string(),
                    fastq_id: z.string(),
                    library_type: z.string(),
                });
                const params = schema.parse(args);
                return wrapResponse(toResponse(await txgCli.runCommand([
                    "fastqs",
                    "library",
                    "--project-id",
                    params.project_id,
                    params.fastq_id,
                    params.library_type,
                ])));
            }
            // --- File Management Tools ---
            case "list_project_files": {
                const schema = z.object({
                    project_id: z.string(),
                });
                const params = schema.parse(args);
                return wrapResponse(toResponse(await txgCli.runCommand(["files", "list", params.project_id])));
            }
            case "upload_project_file": {
                const schema = z.object({
                    project_id: z.string(),
                    file_path: z.string(),
                });
                const params = schema.parse(args);
                return wrapResponse(toResponse(await txgCli.runCommand([
                    "files",
                    "upload",
                    "--project-id",
                    params.project_id,
                    params.file_path,
                ], LONG_OPERATION_TIMEOUT, LONG_OPERATION_SOFT_TIMEOUT)));
            }
            case "download_project_file": {
                const schema = z.object({
                    project_id: z.string(),
                    file_id: z.string(),
                    output_path: z.string(),
                });
                const params = schema.parse(args);
                return wrapResponse(toResponse(await txgCli.runCommand([
                    "files",
                    "download",
                    params.project_id,
                    "--file-id",
                    params.file_id,
                    "--target-dir",
                    params.output_path,
                ], LONG_OPERATION_TIMEOUT, LONG_OPERATION_SOFT_TIMEOUT)));
            }
            // --- Project Management Tools ---
            case "list_projects":
                return wrapResponse(toResponse(await txgCli.runCommand(["projects", "list"])));
            case "create_project": {
                const schema = z.object({
                    name: z.string(),
                    description: z.string().optional(),
                });
                const params = schema.parse(args);
                const command = ["projects", "create", "--name", params.name];
                if (params.description) {
                    command.push("--description", params.description);
                }
                return wrapResponse(toResponse(await txgCli.runCommand(command)));
            }
            case "update_project": {
                const schema = z.object({
                    project_id: z.string(),
                    name: z.string().optional(),
                    description: z.string().optional(),
                });
                const params = schema.parse(args);
                const command = ["projects", "update", params.project_id];
                if (params.name) {
                    command.push("--name", params.name);
                }
                if (params.description) {
                    command.push("--description", params.description);
                }
                return wrapResponse(toResponse(await txgCli.runCommand(command)));
            }
            // --- Reference Management Tools ---
            case "list_custom_references":
                return wrapResponse(toResponse(await txgCli.runCommand(["references", "list"])));
            case "list_prebuilt_references":
                return wrapResponse(toResponse(await txgCli.runCommand(["prebuilt-references"])));
            case "get_reference": {
                const schema = z.object({
                    reference_id: z.string(),
                });
                const params = schema.parse(args);
                return wrapResponse(toResponse(await txgCli.runCommand(["references", "get", params.reference_id])));
            }
            case "update_reference": {
                const schema = z.object({
                    reference_id: z.string(),
                    name: z.string().optional(),
                    description: z.string().optional(),
                    organism: z.string().optional(),
                });
                const params = schema.parse(args);
                const command = ["references", "update", params.reference_id];
                if (params.name) {
                    command.push("--name", params.name);
                }
                if (params.description) {
                    command.push("--description", params.description);
                }
                if (params.organism) {
                    command.push("--organism", params.organism);
                }
                return wrapResponse(toResponse(await txgCli.runCommand(command)));
            }
            case "upload_reference": {
                const schema = z.object({
                    file_path: z.string(),
                    name: z.string(),
                });
                const params = schema.parse(args);
                return wrapResponse(toResponse(await txgCli.runCommand(["references", "upload", params.file_path, "--name", params.name], LONG_OPERATION_TIMEOUT, LONG_OPERATION_SOFT_TIMEOUT)));
            }
            default:
                throw new Error(`Unknown tool: ${name}`);
        }
    });
}
