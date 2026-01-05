# BioMCP Reference Server

## Overview
This is a reference implementation of a **Model Context Protocol (MCP)** server tailored for biomedical tasks. It allows AI clients (like Claude Desktop or MCP-enabled IDEs) to connect to specialized biomedical "tools" locally.

## What is MCP?
The Model Context Protocol (MCP) is an open standard that enables AI models to interact with external data and tools. This server exposes biomedical capabilities as MCP "tools".

## Available Tools
1.  **search_pubmed(query)**: Simulates searching PubMed for literature.
2.  **get_gene_info(gene_symbol)**: Retrieves summaries for key genes (TP53, BRCA1, etc.).

## Installation & Usage

### 1. Prerequisites
- Python 3.10+
- An MCP Client (e.g., Claude Desktop App)

### 2. Configuration (Claude Desktop)
Add the following to your `claude_desktop_config.json` (usually in `~/Library/Application Support/Claude/` on Mac or `%APPDATA%\Claude\` on Windows, or `~/.config/Claude/` on Linux):

```json
{
  "mcpServers": {
    "biomcp": {
      "command": "python3",
      "args": [
        "/absolute/path/to/Skills/MCP_Servers/BioMCP_Reference/server.py"
      ]
    }
  }
}
```

### 3. Testing
Once configured, restart Claude Desktop. You should see a "socket" icon indicating the server is connected. You can then ask Claude:
> "Check PubMed for recent papers on CRISPR."
> "What is the function of the TP53 gene?"

The AI will call the tools provided by this server to answer.

## Logs
Server logs are written to `biomcp.log` in the same directory as the script.
