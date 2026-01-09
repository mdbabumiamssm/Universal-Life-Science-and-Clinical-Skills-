from typing import List, Dict, Any, Optional
from dataclasses import dataclass
import json
import re

@dataclass
class AgentMessage:
    role: str
    content: str
    sender: str

class BaseAgent:
    """Mock base agent for demonstration."""
    def __init__(self, name: str, expertise: str):
        self.name = name
        self.expertise = expertise

    def execute(self, task: str) -> str:
        # In a real system, this calls an LLM
        return f"[{self.name}]: Analyzed '{task}' using {self.expertise}."

class SupervisorAgent:
    """
    An Agent Orchestrator that plans, delegates, and synthesizes.
    Implements a 'Router' pattern using structured JSON outputs.
    """
    
    def __init__(self, name: str, sub_agents: List[BaseAgent]):
        self.name = name
        self.sub_agents = {agent.name: agent for agent in sub_agents}
        self.memory: List[AgentMessage] = []

    def plan(self, user_query: str) -> List[Dict[str, str]]:
        """
        Decomposes a query into a structured plan.
        Returns a list of steps: [{"agent": "Coder", "instruction": "Write python script..."}]
        """
        print(f"[{self.name}] Planning for: {user_query}")
        
        # MOCK LLM OUTPUT (Simulation of a Router Chain)
        # Real implementation would prompt LLM to return JSON
        
        plan = []
        # Robust keyword matching for the demo
        q = user_query.lower()
        if any(k in q for k in ["drug", "molecule", "inhibitor", "chemistry"]):
            plan.append({"agent": "Chemist", "instruction": "Analyze molecular properties."})
            plan.append({"agent": "SafetyOfficer", "instruction": "Check for toxicity risks."})
        elif any(k in q for k in ["code", "implement", "python", "script"]):
             plan.append({"agent": "Coder", "instruction": "Generate the requested code."})
             plan.append({"agent": "Reviewer", "instruction": "Review code for bugs."})
        else:
            # Fallback if no specific domain detected, try to use a default or just log
            plan.append({"agent": "Chemist", "instruction": "General scientific analysis."})
            
        return plan

    def execute_workflow(self, user_query: str) -> Dict[str, Any]:
        """
        Executes the planned workflow.
        """
        self.memory.append(AgentMessage(role="user", content=user_query, sender="User"))
        
        # 1. Plan
        plan = self.plan(user_query)
        results = {}
        
        # 2. Delegate
        for step in plan:
            agent_name = step["agent"]
            instruction = step["instruction"]
            
            if agent_name in self.sub_agents:
                print(f"[{self.name}] Delegating to -> {agent_name}: {instruction}")
                output = self.sub_agents[agent_name].execute(instruction)
                results[agent_name] = output
                self.memory.append(AgentMessage(role="function", content=output, sender=agent_name))
            else:
                print(f"[{self.name}] Error: Agent {agent_name} not found.")

        # 3. Synthesize (Mock LLM synthesis)
        final_response = "Synthesized Answer: Based on the team's input, here is the result...\n"
        for agent, res in results.items():
            final_response += f"- {res}\n"
            
        self.memory.append(AgentMessage(role="assistant", content=final_response, sender=self.name))
        return {
            "query": user_query,
            "plan": plan,
            "trace": results,
            "final_answer": final_response
        }

# --- Example Usage ---
if __name__ == "__main__":
    # Define specialized agents
    chemist = BaseAgent("Chemist", "Cheminformatics & RDKit")
    safety = BaseAgent("SafetyOfficer", "Toxicology & FDA Regulations")
    coder = BaseAgent("Coder", "Python & PyTorch")
    reviewer = BaseAgent("Reviewer", "Security & Performance")
    
    # Initialize Supervisor
    supervisor = SupervisorAgent("ChiefScientist", [chemist, safety, coder, reviewer])
    
    # Run a Bio-task
    print("--- Workflow 1: Drug Discovery ---")
    result_bio = supervisor.execute_workflow("I need to design a new inhibitor for EGFR.")
    print(result_bio["final_answer"])
    
    print("\n--- Workflow 2: Coding Task ---")
    result_code = supervisor.execute_workflow("Implement a PyTorch training loop.")
    print(result_code["final_answer"])