from typing import List, Dict, Any
from dataclasses import dataclass
import random

@dataclass
class PlanNode:
    id: str
    instruction: str
    status: str = "pending" # pending, in_progress, completed, failed
    output: str = ""

class PlanAndSolveAgent:
    """
    Implements the 'Plan-and-Solve' prompting strategy.
    Unlike ReAct (which interleaves thought/action), this creates a full plan upfront.
    
    Paper: 'Plan-and-Solve Prompting: Improving Zero-Shot Chain-of-Thought Reasoning by Large Language Models'
    """
    
    def __init__(self, name: str):
        self.name = name
        self.current_plan: List[PlanNode] = []

    def generate_plan(self, complex_query: str) -> List[PlanNode]:
        """
        Uses a 'Planner' prompt to break down the query.
        """
        print(f"[{self.name}] Generating plan for: {complex_query}")
        
        # Mocking the LLM decomposition
        # Query: "Find proteins associated with Diabetes and check if Aspirin binds to them."
        
        steps = [
            PlanNode("1", "Search knowledge graph for proteins associated with Diabetes."),
            PlanNode("2", "Filter list for druggable targets."),
            PlanNode("3", "For each target, check binding affinity of Aspirin."),
            PlanNode("4", "Summarize findings.")
        ]
        self.current_plan = steps
        return steps

    def execute_plan(self):
        """
        Executes the plan sequentially.
        """
        print(f"[{self.name}] Executing Plan ({len(self.current_plan)} steps)...")
        
        for step in self.current_plan:
            step.status = "in_progress"
            print(f"  > Step {step.id}: {step.instruction}")
            
            # Mock Execution
            # In prod, this would call tools or other agents
            step.output = f"Result for '{step.instruction}'"
            step.status = "completed"
            
    def get_summary(self) -> str:
        return "\n".join([f"{n.id}. {n.instruction} -> {n.output}" for n in self.current_plan])

# --- Example Usage ---
if __name__ == "__main__":
    agent = PlanAndSolveAgent("PlannerBot")
    query = "Analyze the effect of inhibitor X on Pathway Y considering mutation Z."
    
    agent.generate_plan(query)
    agent.execute_plan()
    
    print("\n--- Final Report ---")
    print(agent.get_summary())
