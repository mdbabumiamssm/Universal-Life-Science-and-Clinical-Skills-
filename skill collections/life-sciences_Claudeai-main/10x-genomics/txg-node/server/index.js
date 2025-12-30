#!/usr/bin/env node
import { Server } from "@modelcontextprotocol/sdk/server/index.js";
import { StdioServerTransport } from "@modelcontextprotocol/sdk/server/stdio.js";
import { ListToolsRequestSchema, ListPromptsRequestSchema, ListResourcesRequestSchema, } from "@modelcontextprotocol/sdk/types.js";
import { registerTools } from "./tools.js";
import { registerPrompts } from "./prompts.js";
import { TOOL_DEFINITIONS } from "./tool-definitions.js";
// Factory function for creating server (useful for testing)
export function createServer() {
    const server = new Server({
        name: "10x-genomics",
        version: "1.0.1",
    }, {
        capabilities: {
            tools: {},
            prompts: {},
            resources: {},
        },
    });
    // Register tools list handler
    server.setRequestHandler(ListToolsRequestSchema, async () => {
        return {
            tools: TOOL_DEFINITIONS,
        };
    });
    // Register prompts list handler
    server.setRequestHandler(ListPromptsRequestSchema, async () => {
        return {
            prompts: [], //PROMPT_DEFINITIONS // Hiding prompts for now untill we have a better use case (e.g. analysis)
        };
    });
    // Register resources list handler (returns empty array as we don't have resources)
    server.setRequestHandler(ListResourcesRequestSchema, async () => {
        return {
            resources: [],
        };
    });
    // Register tool and prompt handlers
    registerTools(server);
    registerPrompts(server);
    // Error handling
    server.onerror = (error) => {
        console.error("[MCP Error]", error);
    };
    return server;
}
// Only run main when not imported (not in test mode)
if (process.env.NODE_ENV !== "test") {
    // Validate environment
    if (!process.env.TXG_CLI_ACCESS_TOKEN) {
        console.error("Warning: TXG_CLI_ACCESS_TOKEN environment variable not set");
    }
    // Create server instance
    const server = createServer();
    process.on("SIGINT", async () => {
        await server.close();
        process.exit(0);
    });
    // Start server
    const main = async () => {
        const transport = new StdioServerTransport();
        await server.connect(transport);
        console.error("10x Genomics MCP Server running...");
    };
    main().catch((error) => {
        console.error("Failed to start server:", error);
        process.exit(1);
    });
}
