"""
BioKernel Runtime: Local Execution Server for Biomedical AI Skills

A unified API server that:
1. Routes requests to optimal LLM (Claude, GPT, Gemini)
2. Executes biomedical code (Python, R, shell)
3. Manages data files and caching
4. Provides MCP and OpenAI-compatible APIs
"""

import os
import sys
import json
import time
import uuid
import asyncio
import logging
from pathlib import Path
from typing import Dict, Any, List, Optional, Union
from dataclasses import dataclass, field, asdict
from datetime import datetime
from enum import Enum
import subprocess
import tempfile


# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("biokernel")


class ModelProvider(Enum):
    """Supported model providers"""
    CLAUDE = "claude"
    OPENAI = "openai"
    GEMINI = "gemini"
    LOCAL = "local"


class TaskComplexity(Enum):
    """Task complexity levels for routing"""
    SIMPLE = "simple"      # Quick summaries, lookups
    MODERATE = "moderate"  # Standard analysis
    COMPLEX = "complex"    # Multi-step reasoning
    EXPERT = "expert"      # Domain expertise required


@dataclass
class RoutingStrategy:
    """Configuration for intelligent model routing"""
    default_provider: ModelProvider = ModelProvider.CLAUDE
    cost_sensitive: bool = False
    latency_sensitive: bool = False
    accuracy_priority: bool = True

    # Provider preferences by task type
    complex_tasks: ModelProvider = ModelProvider.CLAUDE
    simple_tasks: ModelProvider = ModelProvider.GEMINI
    creative_tasks: ModelProvider = ModelProvider.OPENAI
    biomedical_tasks: ModelProvider = ModelProvider.CLAUDE


@dataclass
class ExecutionResult:
    """Result of code execution"""
    success: bool
    output: str
    error: Optional[str] = None
    execution_time_ms: float = 0
    files_created: List[str] = field(default_factory=list)


@dataclass
class JobStatus:
    """Status of an async job"""
    job_id: str
    status: str  # pending, running, completed, failed
    progress: float = 0.0
    result: Optional[Any] = None
    error: Optional[str] = None
    created_at: str = ""
    updated_at: str = ""


class CodeExecutor:
    """
    Secure code execution environment for biomedical computations.

    Supports Python, R, and shell commands with sandboxing.
    """

    def __init__(self, work_dir: str = None, timeout: int = 300):
        self.work_dir = Path(work_dir or tempfile.mkdtemp(prefix="biokernel_"))
        self.timeout = timeout
        self.work_dir.mkdir(parents=True, exist_ok=True)

    def execute_python(self, code: str, inputs: Dict[str, Any] = None) -> ExecutionResult:
        """Execute Python code in isolated environment."""
        start_time = time.time()

        # Create a temporary script
        script_path = self.work_dir / f"script_{uuid.uuid4().hex[:8]}.py"

        # Prepare the script with inputs
        setup_code = ""
        if inputs:
            setup_code = f"inputs = {json.dumps(inputs)}\n"

        full_code = f"""
import sys
import json

# Redirect prints to capture output
_output_buffer = []
_original_print = print
def print(*args, **kwargs):
    import io
    buffer = io.StringIO()
    _original_print(*args, file=buffer, **kwargs)
    _output_buffer.append(buffer.getvalue())
    _original_print(*args, **kwargs)

{setup_code}

# User code
{code}

# Output captured prints
if _output_buffer:
    with open("{script_path}.output", "w") as f:
        f.write("".join(_output_buffer))
"""

        script_path.write_text(full_code)

        try:
            result = subprocess.run(
                [sys.executable, str(script_path)],
                capture_output=True,
                text=True,
                timeout=self.timeout,
                cwd=str(self.work_dir)
            )

            output = result.stdout
            if Path(f"{script_path}.output").exists():
                output = Path(f"{script_path}.output").read_text()

            execution_time = (time.time() - start_time) * 1000

            if result.returncode != 0:
                return ExecutionResult(
                    success=False,
                    output=output,
                    error=result.stderr,
                    execution_time_ms=execution_time
                )

            return ExecutionResult(
                success=True,
                output=output,
                execution_time_ms=execution_time,
                files_created=self._list_created_files()
            )

        except subprocess.TimeoutExpired:
            return ExecutionResult(
                success=False,
                output="",
                error=f"Execution timed out after {self.timeout}s",
                execution_time_ms=self.timeout * 1000
            )
        except Exception as e:
            return ExecutionResult(
                success=False,
                output="",
                error=str(e),
                execution_time_ms=(time.time() - start_time) * 1000
            )

    def execute_shell(self, command: str) -> ExecutionResult:
        """Execute shell command with restrictions."""
        start_time = time.time()

        # Basic command sanitization (in production, use proper sandboxing)
        dangerous = ['rm -rf', 'sudo', 'chmod', '>', '>>', '|']
        for d in dangerous:
            if d in command:
                return ExecutionResult(
                    success=False,
                    output="",
                    error=f"Dangerous command pattern detected: {d}"
                )

        try:
            result = subprocess.run(
                command,
                shell=True,
                capture_output=True,
                text=True,
                timeout=self.timeout,
                cwd=str(self.work_dir)
            )

            return ExecutionResult(
                success=result.returncode == 0,
                output=result.stdout,
                error=result.stderr if result.returncode != 0 else None,
                execution_time_ms=(time.time() - start_time) * 1000
            )
        except Exception as e:
            return ExecutionResult(
                success=False,
                output="",
                error=str(e),
                execution_time_ms=(time.time() - start_time) * 1000
            )

    def _list_created_files(self) -> List[str]:
        """List files created in work directory."""
        return [str(f.relative_to(self.work_dir)) for f in self.work_dir.glob("*")]


class ModelRouter:
    """
    Intelligent routing of requests to optimal LLM provider.

    Considers task complexity, cost, latency, and accuracy requirements.
    """

    def __init__(self, strategy: RoutingStrategy = None):
        self.strategy = strategy or RoutingStrategy()
        self.provider_stats: Dict[str, Dict[str, float]] = {
            "claude": {"latency": 1500, "cost": 0.015, "accuracy": 0.92},
            "openai": {"latency": 1200, "cost": 0.010, "accuracy": 0.90},
            "gemini": {"latency": 800, "cost": 0.005, "accuracy": 0.85},
        }

    def classify_task(self, prompt: str) -> TaskComplexity:
        """Classify task complexity based on prompt analysis."""
        prompt_lower = prompt.lower()

        # Expert-level biomedical tasks
        expert_keywords = ['crispr', 'gene editing', 'clinical trial', 'drug interaction',
                          'variant pathogenicity', 'protein structure']
        if any(kw in prompt_lower for kw in expert_keywords):
            return TaskComplexity.EXPERT

        # Complex reasoning tasks
        complex_keywords = ['analyze', 'compare', 'evaluate', 'design', 'optimize',
                          'differential expression', 'pathway analysis']
        if any(kw in prompt_lower for kw in complex_keywords):
            return TaskComplexity.COMPLEX

        # Moderate tasks
        moderate_keywords = ['summarize', 'explain', 'describe', 'list', 'find']
        if any(kw in prompt_lower for kw in moderate_keywords):
            return TaskComplexity.MODERATE

        return TaskComplexity.SIMPLE

    def select_provider(self, prompt: str, metadata: Dict[str, Any] = None) -> ModelProvider:
        """
        Select optimal provider for a given prompt.

        Args:
            prompt: The user's prompt
            metadata: Optional metadata (e.g., required capabilities)

        Returns:
            Selected ModelProvider
        """
        complexity = self.classify_task(prompt)
        metadata = metadata or {}

        # Check for explicit provider request
        if 'provider' in metadata:
            return ModelProvider(metadata['provider'])

        # Route based on complexity and strategy
        if complexity == TaskComplexity.EXPERT:
            return self.strategy.complex_tasks

        if complexity == TaskComplexity.COMPLEX:
            if self.strategy.accuracy_priority:
                return self.strategy.complex_tasks
            return self.strategy.default_provider

        if complexity == TaskComplexity.SIMPLE:
            if self.strategy.cost_sensitive or self.strategy.latency_sensitive:
                return self.strategy.simple_tasks
            return self.strategy.default_provider

        return self.strategy.default_provider

    def get_provider_config(self, provider: ModelProvider) -> Dict[str, Any]:
        """Get configuration for a provider."""
        configs = {
            ModelProvider.CLAUDE: {
                "model": "claude-sonnet-4-20250514",
                "max_tokens": 4096,
                "api_base": "https://api.anthropic.com"
            },
            ModelProvider.OPENAI: {
                "model": "gpt-4-turbo-preview",
                "max_tokens": 4096,
                "api_base": "https://api.openai.com"
            },
            ModelProvider.GEMINI: {
                "model": "gemini-pro",
                "max_tokens": 4096,
                "api_base": "https://generativelanguage.googleapis.com"
            },
            ModelProvider.LOCAL: {
                "model": "local",
                "max_tokens": 4096,
                "api_base": "http://localhost:11434"  # Ollama default
            }
        }
        return configs.get(provider, configs[ModelProvider.CLAUDE])


class BioKernel:
    """
    Main BioKernel runtime server.

    Provides unified API for:
    - LLM routing and execution
    - Code execution (Python, R, shell)
    - File management
    - Job queue management
    """

    def __init__(
        self,
        routing_strategy: RoutingStrategy = None,
        work_dir: str = None
    ):
        self.router = ModelRouter(strategy=routing_strategy)
        self.executor = CodeExecutor(work_dir=work_dir)
        self.jobs: Dict[str, JobStatus] = {}
        self.skills: Dict[str, Dict[str, Any]] = {}
        self.start_time = datetime.now()

        logger.info("BioKernel initialized")

    def load_skill(self, usdl_path: str) -> str:
        """Load a USDL skill into the kernel."""
        import yaml
        with open(usdl_path, 'r') as f:
            usdl = yaml.safe_load(f)

        skill = usdl['skill']
        skill_id = skill['id']
        self.skills[skill_id] = skill
        logger.info(f"Loaded skill: {skill_id}")
        return skill_id

    def chat(
        self,
        message: str,
        skill_id: str = None,
        provider: str = None,
        **kwargs
    ) -> Dict[str, Any]:
        """
        Process a chat message using optimal routing.

        Args:
            message: User message
            skill_id: Optional skill to use
            provider: Optional provider override

        Returns:
            Response dict with content and metadata
        """
        start_time = time.time()

        # Select provider
        metadata = {"provider": provider} if provider else {}
        selected_provider = self.router.select_provider(message, metadata)
        provider_config = self.router.get_provider_config(selected_provider)

        # Get skill context if specified
        system_prompt = None
        if skill_id and skill_id in self.skills:
            skill = self.skills[skill_id]
            system_prompt = skill.get('prompts', {}).get('system', '')

        # Execute (mock implementation - integrate real APIs in production)
        response_content = self._mock_llm_call(
            message, system_prompt, selected_provider, provider_config
        )

        latency_ms = (time.time() - start_time) * 1000

        return {
            "content": response_content,
            "provider": selected_provider.value,
            "model": provider_config['model'],
            "latency_ms": latency_ms,
            "skill_id": skill_id
        }

    def execute(
        self,
        kernel: str,
        code: str,
        inputs: Dict[str, Any] = None
    ) -> Dict[str, Any]:
        """
        Execute code in specified kernel.

        Args:
            kernel: Kernel type (python, r, shell)
            code: Code to execute
            inputs: Optional input variables

        Returns:
            Execution result dict
        """
        if kernel == "python":
            result = self.executor.execute_python(code, inputs)
        elif kernel == "shell":
            result = self.executor.execute_shell(code)
        else:
            return {"error": f"Unsupported kernel: {kernel}"}

        return asdict(result)

    def submit_job(self, task: Dict[str, Any]) -> str:
        """
        Submit an async job.

        Args:
            task: Job specification

        Returns:
            Job ID
        """
        job_id = str(uuid.uuid4())
        now = datetime.now().isoformat()

        self.jobs[job_id] = JobStatus(
            job_id=job_id,
            status="pending",
            created_at=now,
            updated_at=now
        )

        # In production, this would queue the job
        logger.info(f"Job submitted: {job_id}")
        return job_id

    def get_job_status(self, job_id: str) -> Optional[Dict[str, Any]]:
        """Get status of a job."""
        job = self.jobs.get(job_id)
        if job:
            return asdict(job)
        return None

    def health_check(self) -> Dict[str, Any]:
        """Return server health status."""
        uptime = (datetime.now() - self.start_time).total_seconds()
        return {
            "status": "healthy",
            "uptime_seconds": uptime,
            "skills_loaded": len(self.skills),
            "jobs_pending": sum(1 for j in self.jobs.values() if j.status == "pending"),
            "providers": {
                "claude": self._check_provider("claude"),
                "openai": self._check_provider("openai"),
                "gemini": self._check_provider("gemini")
            }
        }

    def _mock_llm_call(
        self,
        message: str,
        system_prompt: str,
        provider: ModelProvider,
        config: Dict[str, Any]
    ) -> str:
        """Mock LLM call for demo purposes."""
        # In production, integrate with actual provider APIs
        return f"""[BioKernel Response]
Provider: {provider.value}
Model: {config['model']}

Based on your query about: "{message[:100]}..."

This is a mock response from BioKernel. In production, this would be routed to the actual {provider.value} API.

System context applied: {'Yes' if system_prompt else 'No'}
"""

    def _check_provider(self, provider: str) -> str:
        """Check if provider API key is configured."""
        env_vars = {
            "claude": "ANTHROPIC_API_KEY",
            "openai": "OPENAI_API_KEY",
            "gemini": "GOOGLE_API_KEY"
        }
        return "configured" if os.environ.get(env_vars.get(provider, "")) else "not_configured"


# FastAPI Server (optional, requires fastapi and uvicorn)
def create_api_server(kernel: BioKernel):
    """Create FastAPI server for BioKernel."""
    try:
        from fastapi import FastAPI, HTTPException
        from pydantic import BaseModel

        app = FastAPI(
            title="BioKernel API",
            description="Unified API for biomedical AI skills",
            version="1.0.0"
        )

        class ChatRequest(BaseModel):
            message: str
            skill_id: Optional[str] = None
            provider: Optional[str] = None

        class ExecuteRequest(BaseModel):
            kernel: str
            code: str
            inputs: Optional[Dict[str, Any]] = None

        @app.get("/health")
        def health():
            return kernel.health_check()

        @app.post("/chat")
        def chat(request: ChatRequest):
            return kernel.chat(
                message=request.message,
                skill_id=request.skill_id,
                provider=request.provider
            )

        @app.post("/execute")
        def execute(request: ExecuteRequest):
            return kernel.execute(
                kernel=request.kernel,
                code=request.code,
                inputs=request.inputs
            )

        @app.post("/skills/load")
        def load_skill(path: str):
            try:
                skill_id = kernel.load_skill(path)
                return {"skill_id": skill_id}
            except Exception as e:
                raise HTTPException(status_code=400, detail=str(e))

        @app.get("/skills")
        def list_skills():
            return {"skills": list(kernel.skills.keys())}

        return app

    except ImportError:
        logger.warning("FastAPI not installed. API server not available.")
        return None


# CLI interface
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="BioKernel Runtime Server")
    parser.add_argument("--host", default="127.0.0.1", help="Host to bind")
    parser.add_argument("--port", type=int, default=8000, help="Port to bind")
    parser.add_argument("--work-dir", help="Working directory for execution")
    parser.add_argument("--load-skills", nargs="*", help="USDL files to preload")
    parser.add_argument("--demo", action="store_true", help="Run in demo mode")

    args = parser.parse_args()

    # Initialize kernel
    kernel = BioKernel(work_dir=args.work_dir)

    # Load skills if provided
    if args.load_skills:
        for skill_path in args.load_skills:
            try:
                kernel.load_skill(skill_path)
            except Exception as e:
                logger.error(f"Failed to load {skill_path}: {e}")

    if args.demo:
        # Demo mode - run some test operations
        print("\n=== BioKernel Demo ===\n")

        # Health check
        health = kernel.health_check()
        print(f"Health: {health['status']}")
        print(f"Skills loaded: {health['skills_loaded']}")

        # Chat demo
        response = kernel.chat("What genes are involved in breast cancer?")
        print(f"\nChat Response:\n{response['content']}")

        # Execution demo
        result = kernel.execute("python", "print('Hello from BioKernel!')")
        print(f"\nExecution Result: {result['output']}")

    else:
        # Start API server
        app = create_api_server(kernel)
        if app:
            import uvicorn
            print(f"\nStarting BioKernel server on {args.host}:{args.port}")
            uvicorn.run(app, host=args.host, port=args.port)
        else:
            print("Install fastapi and uvicorn to run the API server:")
            print("  pip install fastapi uvicorn")
