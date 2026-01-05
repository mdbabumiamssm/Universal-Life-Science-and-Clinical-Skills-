import sys
import random

# Mock RDKit if not available
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

class ChemCrowAgent:
    """
    A simplified agent inspired by ChemCrow that provides tools for 
    chemical property analysis, safety checks, and basic synthesis planning.
    """
    
    def __init__(self):
        self.tools = {
            "MolWeight": self.get_molecular_weight,
            "Safety": self.check_safety,
            "SMILES2Name": self.smiles_to_name,
            "Validity": self.check_validity
        }
        
        if not HAS_RDKIT:
            print("Warning: 'rdkit' not found. Agent running in MOCK mode.")

    def run_tool(self, tool_name, input_data):
        if tool_name not in self.tools:
            return f"Error: Tool '{tool_name}' not found."
        return self.tools[tool_name](input_data)

    def check_validity(self, smiles):
        """Checks if a SMILES string is valid."""
        if HAS_RDKIT:
            mol = Chem.MolFromSmiles(smiles)
            return mol is not None
        else:
            # Mock validation: just check if it looks vaguely like SMILES
            return len(smiles) > 0 and not " " in smiles

    def get_molecular_weight(self, smiles):
        """Calculates Molecular Weight."""
        if not self.check_validity(smiles):
            return "Invalid SMILES"

        if HAS_RDKIT:
            mol = Chem.MolFromSmiles(smiles)
            return Descriptors.MolWt(mol)
        else:
            # Mock: Random weight based on length to be deterministic-ish
            return round(len(smiles) * 12.5 + 2.0, 2)

    def check_safety(self, smiles):
        """
        Checks for basic safety flags (Explosive, Toxic).
        Real implementation would check substructures or query a database.
        """
        if not self.check_validity(smiles):
            return "Invalid SMILES"
            
        # Basic heuristic flags
        flags = []
        if "N(=O)=O" in smiles or "nitro" in smiles.lower():
            flags.append("High Energy/Explosive Potential (Nitro group)")
        if "Hg" in smiles:
            flags.append("Toxic Heavy Metal (Mercury)")
        if "P" in smiles and "F" in smiles: # G-series nerve agents often have P-F
            flags.append("Potential Toxicity (Organophosphate-like)")
            
        if not flags:
            return "No obvious structural alerts found (Basic Check)."
        return "WARNING: " + "; ".join(flags)

    def smiles_to_name(self, smiles):
        """
        Mocks converting SMILES to a common name.
        Real implementation would query PubChem or similar.
        """
        # Hardcoded dictionary for demo
        lookup = {
            "CC(=O)Oc1ccccc1C(=O)O": "Aspirin",
            "CC(=O)Nc1ccc(O)cc1": "Acetaminophen (Tylenol)",
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C": "Caffeine"
        }
        return lookup.get(smiles, "Unknown Molecule (External Search Required)")

def main():
    agent = ChemCrowAgent()
    
    print("--- ChemCrow Lite Agent ---")
    print("Tools: " + ", ".join(agent.tools.keys()))
    
    # Interactive Loop
    if len(sys.argv) > 1:
        # CLI Mode
        cmd = sys.argv[1]
        arg = sys.argv[2] if len(sys.argv) > 2 else ""
        print(agent.run_tool(cmd, arg))
    else:
        # Demo Mode
        test_smiles = [
            ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
            ("TNT", "Cc1c(N(=O)=O)cc(N(=O)=O)cc1N(=O)=O"), # Trinitrotoluene
            ("Water", "O")
        ]
        
        for name, smiles in test_smiles:
            print(f"\nAnalyzing {name} ({smiles}):")
            print(f"  Valid: {agent.run_tool('Validity', smiles)}")
            print(f"  MW:    {agent.run_tool('MolWeight', smiles)}")
            print(f"  Safe:  {agent.run_tool('Safety', smiles)}")

if __name__ == "__main__":
    main()
