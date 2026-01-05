import sys
import json
import logging

# Configure logging to file since stdout is used for MCP protocol
logging.basicConfig(filename='biomcp.log', level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def log(msg):
    logging.info(msg)

class BioMCPServer:
    def __init__(self):
        self.tools = {
            "search_pubmed": {
                "description": "Search for biomedical literature on PubMed.",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "query": {"type": "string", "description": "Search terms"}
                    },
                    "required": ["query"]
                },
                "handler": self.handle_pubmed_search
            },
            "get_gene_info": {
                "description": "Retrieve summary information for a specific gene.",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "gene_symbol": {"type": "string", "description": "Gene symbol (e.g., TP53, BRCA1)"}
                    },
                    "required": ["gene_symbol"]
                },
                "handler": self.handle_gene_info
            }
        }

    def handle_pubmed_search(self, args):
        query = args.get('query')
        log(f"Executing PubMed search for: {query}")
        # Mock response
        return {
            "content": [
                {
                    "type": "text",
                    "text": f"Found 3 relevant papers for '{query}':\n1. Recent advances in {query} therapy (2025).\n2. {query} expression analysis in oncology (2024).\n3. Clinical guidelines for {query} (2024)."
                }
            ]
        }

    def handle_gene_info(self, args):
        gene = args.get('gene_symbol', '').upper()
        log(f"Fetching info for gene: {gene}")
        
        # Mock database
        db = {
            "TP53": "Tumor protein p53. Builds a protein that suppresses tumor growth. Mutations are found in ~50% of cancers.",
            "BRCA1": "Breast cancer type 1 susceptibility protein. Involved in DNA repair. Mutations increase risk of breast/ovarian cancer.",
            "EGFR": "Epidermal Growth Factor Receptor. A target for cancer therapy (e.g., in NSCLC)."
        }
        
        info = db.get(gene, f"No detailed information found for {gene} in local cache.")
        
        return {
            "content": [
                {
                    "type": "text",
                    "text": info
                }
            ]
        }

    def run(self):
        log("BioMCP Server Started. Listening on Stdio...")
        
        while True:
            try:
                line = sys.stdin.readline()
                if not line:
                    break
                
                request = json.loads(line)
                log(f"Received request: {request.get('method')}")
                
                response = None
                
                if request.get("method") == "tools/list":
                    # Return list of tools
                    tool_list = []
                    for name, details in self.tools.items():
                        tool_list.append({
                            "name": name,
                            "description": details["description"],
                            "inputSchema": details["parameters"]
                        })
                    
                    response = {
                        "jsonrpc": "2.0",
                        "id": request["id"],
                        "result": {
                            "tools": tool_list
                        }
                    }

                elif request.get("method") == "tools/call":
                    # Execute tool
                    params = request.get("params", {})
                    tool_name = params.get("name")
                    args = params.get("arguments", {})
                    
                    if tool_name in self.tools:
                        result = self.tools[tool_name]["handler"](args)
                        response = {
                            "jsonrpc": "2.0",
                            "id": request["id"],
                            "result": result
                        }
                    else:
                         response = {
                            "jsonrpc": "2.0",
                            "id": request["id"],
                            "error": {"code": -32601, "message": "Tool not found"}
                        }
                
                elif request.get("method") == "initialize":
                     response = {
                        "jsonrpc": "2.0",
                        "id": request["id"],
                        "result": {
                            "protocolVersion": "2024-11-05",
                            "capabilities": {
                                "tools": {}
                            },
                            "serverInfo": {
                                "name": "biomcp-reference",
                                "version": "0.1.0"
                            }
                        }
                    }

                elif request.get("method") == "notifications/initialized":
                    # No response needed for notifications
                    continue

                else:
                    # Ignore other methods or ping
                     pass

                if response:
                    sys.stdout.write(json.dumps(response) + "\n")
                    sys.stdout.flush()

            except Exception as e:
                log(f"Error: {e}")
                # Try to send error if possible, or just continue
                continue

if __name__ == "__main__":
    server = BioMCPServer()
    server.run()
