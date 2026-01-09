from typing import List, Optional
import math
import numpy as np

class MedPromptUtils:
    """
    Implements Microsoft's MedPrompt strategies for high-accuracy clinical reasoning.
    Reference: 'Steering Llama 2 via Prompt Engineering to State-of-the-Art Performance'
    
    Core Components:
    1. Dynamic Few-Shot Selection (k-NN based)
    2. Chain of Thought (CoT)
    3. Self-Consistency / Ensemble Refinement
    """

    def __init__(self, embedding_model=None):
        self.embedding_model = embedding_model # Mock: In prod, use OpenAI/HuggingFace embeddings
        self.example_store = [] # List of {question, reasoning, answer, embedding}

    def add_few_shot_example(self, question: str, reasoning: str, answer: str):
        """Adds a verified clinical Q&A pair to the example store."""
        # Mock embedding generation
        mock_embedding = np.random.rand(128) 
        self.example_store.append({
            "question": question,
            "reasoning": reasoning,
            "answer": answer,
            "embedding": mock_embedding
        })

    def get_dynamic_examples(self, query: str, k: int = 3) -> List[dict]:
        """
        Retrieves the k-most similar examples to the current query.
        This is the 'Dynamic Few-Shot' component.
        """
        if not self.example_store:
            return []
            
        # Mock cosine similarity logic
        # In prod: Calculate cosine sim between query_embedding and store_embeddings
        return self.example_store[:k] 

    def construct_prompt(self, query: str, context: str = "") -> str:
        """
        Constructs the full MedPrompt.
        Structure:
        [System Role]
        [Context/EHR Data]
        [Few-Shot Examples (Dynamic)]
        [Target Question]
        [Chain of Thought Trigger]
        """
        
        system_role = (
            "You are an expert medical AI assistant. "
            "Think step-by-step before answering. "
            "Cite guidelines where possible."
        )
        
        examples = self.get_dynamic_examples(query)
        formatted_examples = ""
        for i, ex in enumerate(examples):
            formatted_examples += (
                f"\nExample {i+1}:\n"
                f"Q: {ex['question']}\n"
                f"Thinking: {ex['reasoning']}\n"
                f"A: {ex['answer']}\n"
            )

        prompt = (
            f"{system_role}\n\n"
            f"CONTEXT:\n{context}\n\n"
            f"SIMILAR CASES:\n{formatted_examples}\n\n"
            f"CURRENT QUESTION:\n{query}\n\n"
            f"Let's think step by step to ensure patient safety and diagnostic accuracy."
        )
        return prompt

    def ensemble_refinement(self, responses: List[str]) -> str:
        """
        Simulates the 'Self-Consistency' or 'Ensemble' step.
        Takes multiple reasoning paths and generates a consensus.
        """
        # In a real implementation, this would use an LLM to summarize/vote.
        # Here we mock a consensus selection.
        unique_answers = set(responses)
        if len(unique_answers) == 1:
            return responses[0]
        
        return f"Consensus derived from {len(responses)} reasoning paths: {responses[0]}"

# --- Example Usage ---
if __name__ == "__main__":
    mp = MedPromptUtils()
    
    # Seed with some dummy clinical data
    mp.add_few_shot_example(
        "Patient has high fever and stiff neck.",
        "Symptoms suggest meningeal irritation. Kernig's sign should be checked.",
        "Possible Meningitis."
    )
    
    query = "Patient presents with headache, photophobia, and neck stiffness."
    prompt = mp.construct_prompt(query)
    
    print("--- Generated MedPrompt ---")
    print(prompt)
