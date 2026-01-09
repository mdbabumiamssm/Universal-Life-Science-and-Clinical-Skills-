import argparse
import sys
import json
from typing import Any, Dict

# Mocking MCP components for demonstration if library is missing
try:
    from mcp.server.fastmcp import FastMCP
except ImportError:
    class FastMCP:
        def __init__(self, name, version="0.1.0"):
            self.name = name
            self.version = version
            self.tools = {}
        
        def tool(self):
            def decorator(func):
                self.tools[func.__name__] = func
                return func
            return decorator
        
        def run(self):
            print(f"[{self.name}] Server running (Mock Mode). Available tools: {list(self.tools.keys())}")
            # In a real scenario, this would start a stdio/SSE server
            # For this script, we'll just demonstrate one tool call
            print("Simulating tool call 'sequence_length('ATCG')'...")
            res = self.tools['sequence_length']("ATCG")
            print(f"Result: {res}")

# Initialize Server
mcp = FastMCP("BioMCP", version="1.0.0")

@mcp.tool()
def sequence_length(sequence: str) -> int:
    """Calculates the length of a DNA/RNA/Protein sequence."""
    return len(sequence)

@mcp.tool()
def reverse_complement(sequence: str) -> str:
    """Returns the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(sequence.upper()))

@mcp.tool()
def calculate_gc_content(sequence: str) -> float:
    """Calculates GC content percentage."""
    if not sequence:
        return 0.0
    g = sequence.upper().count('G')
    c = sequence.upper().count('C')
    return (g + c) / len(sequence) * 100.0

if __name__ == "__main__":
    # If run directly, start the server
    mcp.run()
