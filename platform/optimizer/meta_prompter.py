"""
Meta-Prompter: AI-Driven Prompt Optimizer for Biomedical Skills

Optimizes USDL skill prompts for each target platform (Claude, OpenAI, Gemini)
using a "teacher" LLM to rewrite prompts according to platform-specific best practices.
"""

import yaml
import json
import os
from pathlib import Path
from typing import Dict, Any, List, Optional
from dataclasses import dataclass, field
from abc import ABC, abstractmethod


# Platform-specific style guides for the optimizer
PLATFORM_STYLE_GUIDES = {
    "claude": """
You are optimizing a biomedical prompt for Claude (Anthropic).

Claude's strengths and preferences:
1. XML-style tags for structure: <input>, <instructions>, <context>, <examples>
2. Long context window (200K tokens) - can include extensive context
3. Assistant prefilling - start the assistant response to steer output
4. Explicit reasoning steps - Claude excels with step-by-step thinking
5. Safety-conscious - include appropriate medical disclaimers

Optimization rules:
- Wrap user inputs in <input> tags
- Use <instructions> for clear directives
- Add <thinking> prompts for complex reasoning
- Include medical safety disclaimers where appropriate
- Use markdown formatting for structure
- Prefill assistant response when steering is needed
""",

    "openai": """
You are optimizing a biomedical prompt for GPT-4 (OpenAI).

GPT-4's strengths and preferences:
1. Concise system messages - be direct and clear
2. Strong function calling - leverage structured outputs
3. JSON mode - request specific output formats
4. Role-play adherence - persona instructions work well
5. Temperature sensitivity - specify when determinism matters

Optimization rules:
- Keep system prompts focused and concise
- Use numbered lists for multi-step instructions
- Specify output format explicitly (JSON, markdown, etc.)
- Include few-shot examples inline
- Add "You must..." for critical constraints
- Use markdown headers sparingly
""",

    "gemini": """
You are optimizing a biomedical prompt for Gemini (Google).

Gemini's strengths and preferences:
1. Multimodal native - reference images/documents naturally
2. Long context (1M+ tokens) - massive context window
3. Grounding with Google Search - can verify facts
4. Structured output mode - JSON schemas work well
5. Code execution - can run Python directly

Optimization rules:
- Structure prompts for potential search grounding
- Use clear section headers
- Specify when factual verification is needed
- Leverage code execution for calculations
- Include explicit output schemas
- Reference Google's knowledge when appropriate
"""
}


@dataclass
class OptimizationResult:
    """Container for optimization results"""
    original_prompt: str
    optimized_prompt: str
    platform: str
    changes_made: List[str]
    confidence_score: float
    metadata: Dict[str, Any] = field(default_factory=dict)


class LLMBackend(ABC):
    """Abstract base class for LLM backends"""

    @abstractmethod
    def generate(self, prompt: str, system: str = None) -> str:
        pass


class MockLLMBackend(LLMBackend):
    """Mock backend for testing without API calls"""

    def generate(self, prompt: str, system: str = None) -> str:
        # Return a mock optimization (for testing)
        return f"[OPTIMIZED]\n{prompt}\n[/OPTIMIZED]"


class AnthropicBackend(LLMBackend):
    """Anthropic Claude backend"""

    def __init__(self, api_key: str = None):
        self.api_key = api_key or os.environ.get("ANTHROPIC_API_KEY")

    def generate(self, prompt: str, system: str = None) -> str:
        try:
            import anthropic
            client = anthropic.Anthropic(api_key=self.api_key)
            message = client.messages.create(
                model="claude-sonnet-4-20250514",
                max_tokens=4096,
                system=system or "",
                messages=[{"role": "user", "content": prompt}]
            )
            return message.content[0].text
        except ImportError:
            raise ImportError("anthropic package not installed. Run: pip install anthropic")


class OpenAIBackend(LLMBackend):
    """OpenAI GPT backend"""

    def __init__(self, api_key: str = None):
        self.api_key = api_key or os.environ.get("OPENAI_API_KEY")

    def generate(self, prompt: str, system: str = None) -> str:
        try:
            import openai
            client = openai.OpenAI(api_key=self.api_key)
            messages = []
            if system:
                messages.append({"role": "system", "content": system})
            messages.append({"role": "user", "content": prompt})

            response = client.chat.completions.create(
                model="gpt-4-turbo-preview",
                messages=messages,
                max_tokens=4096
            )
            return response.choices[0].message.content
        except ImportError:
            raise ImportError("openai package not installed. Run: pip install openai")


class MetaPrompter:
    """
    AI-driven prompt optimizer for biomedical skills.

    Uses a "teacher" LLM to rewrite prompts according to
    platform-specific best practices.
    """

    def __init__(self, backend: LLMBackend = None):
        """
        Initialize the Meta-Prompter.

        Args:
            backend: LLM backend for optimization. Defaults to MockLLMBackend.
        """
        self.backend = backend or MockLLMBackend()
        self.style_guides = PLATFORM_STYLE_GUIDES

    def optimize_skill(self, usdl_path: str, target_platform: str) -> OptimizationResult:
        """
        Optimize a USDL skill for a specific platform.

        Args:
            usdl_path: Path to USDL YAML file
            target_platform: Target platform (claude, openai, gemini)

        Returns:
            OptimizationResult with optimized prompts
        """
        with open(usdl_path, 'r') as f:
            usdl = yaml.safe_load(f)

        skill = usdl['skill']
        original_prompt = skill.get('prompts', {}).get('system', '')

        # Get platform style guide
        style_guide = self.style_guides.get(target_platform, "")

        # Build optimization prompt
        optimization_prompt = self._build_optimization_prompt(
            original_prompt=original_prompt,
            skill_context=skill,
            target_platform=target_platform
        )

        # Call LLM to optimize
        optimized_prompt = self.backend.generate(
            prompt=optimization_prompt,
            system=style_guide
        )

        # Extract the optimized prompt from response
        optimized_prompt = self._extract_optimized(optimized_prompt)

        # Analyze changes
        changes = self._analyze_changes(original_prompt, optimized_prompt, target_platform)

        return OptimizationResult(
            original_prompt=original_prompt,
            optimized_prompt=optimized_prompt,
            platform=target_platform,
            changes_made=changes,
            confidence_score=0.85,  # Would be calculated from eval results
            metadata={
                "skill_id": skill['id'],
                "skill_version": skill['version']
            }
        )

    def optimize_for_all_platforms(self, usdl_path: str) -> Dict[str, OptimizationResult]:
        """
        Optimize a skill for all supported platforms.

        Args:
            usdl_path: Path to USDL YAML file

        Returns:
            Dictionary mapping platform to OptimizationResult
        """
        results = {}
        for platform in ['claude', 'openai', 'gemini']:
            results[platform] = self.optimize_skill(usdl_path, platform)
        return results

    def _build_optimization_prompt(
        self,
        original_prompt: str,
        skill_context: Dict[str, Any],
        target_platform: str
    ) -> str:
        """Build the prompt for the optimizer LLM."""

        capabilities_text = "\n".join([
            f"- {cap['name']}: {cap['description']}"
            for cap in skill_context.get('capabilities', [])
        ])

        return f"""
Please optimize the following biomedical AI prompt for the {target_platform.upper()} platform.

## Skill Information
- **Name**: {skill_context.get('name', 'Unknown')}
- **Category**: {skill_context.get('category', 'Unknown')}
- **Capabilities**:
{capabilities_text}

## Original System Prompt
```
{original_prompt}
```

## Your Task
1. Rewrite the prompt following {target_platform}-specific best practices
2. Maintain all biomedical accuracy and safety guidelines
3. Preserve the core functionality and intent
4. Optimize for the platform's strengths

## Output Format
Return ONLY the optimized prompt wrapped in <optimized> tags:

<optimized>
[Your optimized prompt here]
</optimized>
"""

    def _extract_optimized(self, response: str) -> str:
        """Extract optimized prompt from LLM response."""
        if "<optimized>" in response and "</optimized>" in response:
            start = response.find("<optimized>") + len("<optimized>")
            end = response.find("</optimized>")
            return response[start:end].strip()
        return response.strip()

    def _analyze_changes(
        self,
        original: str,
        optimized: str,
        platform: str
    ) -> List[str]:
        """Analyze what changes were made during optimization."""
        changes = []

        # Check for platform-specific patterns
        if platform == "claude":
            if "<" in optimized and ">" in optimized:
                changes.append("Added XML-style tags for structure")
            if "step" in optimized.lower() or "think" in optimized.lower():
                changes.append("Added reasoning steps")

        elif platform == "openai":
            if "You must" in optimized:
                changes.append("Added explicit constraints")
            if optimized.count("\n") < original.count("\n"):
                changes.append("Condensed for conciseness")

        elif platform == "gemini":
            if "verify" in optimized.lower() or "search" in optimized.lower():
                changes.append("Added grounding references")
            if "```" in optimized:
                changes.append("Added code formatting")

        # General changes
        if len(optimized) > len(original) * 1.2:
            changes.append("Expanded with additional context")
        elif len(optimized) < len(original) * 0.8:
            changes.append("Condensed for efficiency")

        return changes if changes else ["Minor formatting adjustments"]

    def generate_variants(
        self,
        usdl_path: str,
        platform: str,
        n_variants: int = 3
    ) -> List[str]:
        """
        Generate multiple prompt variants for A/B testing.

        Args:
            usdl_path: Path to USDL file
            platform: Target platform
            n_variants: Number of variants to generate

        Returns:
            List of prompt variants
        """
        variants = []
        for i in range(n_variants):
            result = self.optimize_skill(usdl_path, platform)
            variants.append(result.optimized_prompt)
        return variants

    def save_optimized(
        self,
        result: OptimizationResult,
        output_dir: str
    ) -> str:
        """
        Save optimized prompt to file.

        Args:
            result: OptimizationResult to save
            output_dir: Output directory

        Returns:
            Path to saved file
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        output_file = output_path / f"optimized_{result.platform}.yaml"

        output_data = {
            "platform": result.platform,
            "optimized_prompt": result.optimized_prompt,
            "changes_made": result.changes_made,
            "confidence_score": result.confidence_score,
            "metadata": result.metadata
        }

        with open(output_file, 'w') as f:
            yaml.dump(output_data, f, default_flow_style=False)

        return str(output_file)


# CLI interface
if __name__ == "__main__":
    import sys

    if len(sys.argv) < 3:
        print("Usage: python meta_prompter.py <usdl_file.yaml> <platform> [output_dir]")
        print("Platforms: claude, openai, gemini, all")
        sys.exit(1)

    usdl_file = sys.argv[1]
    platform = sys.argv[2].lower()
    output_dir = sys.argv[3] if len(sys.argv) > 3 else "./optimized"

    # Use mock backend for demo (replace with real backend in production)
    prompter = MetaPrompter(backend=MockLLMBackend())

    if platform == "all":
        results = prompter.optimize_for_all_platforms(usdl_file)
        for plat, result in results.items():
            path = prompter.save_optimized(result, output_dir)
            print(f"Optimized for {plat}: {path}")
    else:
        result = prompter.optimize_skill(usdl_file, platform)
        path = prompter.save_optimized(result, output_dir)
        print(f"Optimized for {platform}: {path}")
        print(f"Changes made: {', '.join(result.changes_made)}")
