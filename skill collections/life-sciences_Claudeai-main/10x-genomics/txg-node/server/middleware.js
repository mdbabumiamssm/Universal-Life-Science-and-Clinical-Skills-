// Wrap response in MCP protocol structure
// This is separate from toResponse so that:
// 1. Tests can test toResponse logic without MCP wrapper
// 2. If MCP format changes or we switch to MCPServer, only this function changes
export function wrapResponse(response) {
    return {
        content: [
            {
                type: "text",
                text: JSON.stringify(response),
            },
        ],
    };
}
/**
 * Redacts the access token from the command string for security.
 */
function redactAccessToken(command) {
    return command.replace(/--access-token\s+\S+/, "--access-token ***");
}
export function toResponse(result) {
    // Handle in-progress result
    if (result.inProgress) {
        return {
            content: (result.stdout || result.stderr).trim(),
            error: "",
            status: "in_progress",
            fullCommand: redactAccessToken(result.fullCommand),
            message: "Operation in progress. Depending on file size, this may take a while. The operation will continue in the background.",
        };
    }
    // Handle completed result
    // When command succeeds, treat stderr as additional output, not an error
    // Some CLI tools output success messages to stderr
    const isSuccess = result.exitCode === 0;
    return {
        content: isSuccess
            ? (result.stdout || result.stderr).trim()
            : result.stdout,
        error: isSuccess ? "" : result.stderr,
        status: isSuccess ? "success" : "error",
        returncode: result.exitCode,
        fullCommand: redactAccessToken(result.fullCommand),
    };
}
export function toAnalysisResponse(result) {
    const response = toResponse(result);
    if (response.status === "success") {
        return {
            ...response,
            message: "Analysis started successfully. Expected completion time can be several hours depending on data size. The job will continue running even if this conversation ends.",
            next_steps: "Use 'get_analysis_details' tool to check progress. You will receive an email notification when the analysis is complete.",
        };
    }
    // Return error response as-is
    return response;
}
