import { spawn } from "child_process";
import * as path from "path";
import * as fs from "fs";
import * as os from "os";
import { fileURLToPath } from "url";
// Timeout constants for different operation types
// Quick operations: API calls that should complete in seconds
export const QUICK_OPERATION_TIMEOUT = 120000; // 2 minutes
// Long operations: File uploads/downloads that may take hours
export const LONG_OPERATION_TIMEOUT = 7200000; // 2 hours
// Soft timeout: Return early with progress message for long operations
export const LONG_OPERATION_SOFT_TIMEOUT = 60000; // 1 minute
class TxgCliManager {
    txgPath;
    platform;
    packageDir;
    constructor() {
        this.platform = this.getPlatform();
        // Get directory of current module and go up to package root
        const __filename = fileURLToPath(import.meta.url);
        const __dirname = path.dirname(__filename);
        // Go up from build/server/txg-cli-manager.js to package root
        this.packageDir = path.join(__dirname, "..");
        this.txgPath = this.findBundledBinary();
    }
    getPlatform() {
        const system = os.platform();
        // Map Node.js platform names to our directory names
        let platform;
        if (system === "darwin") {
            platform = "darwin";
        }
        else if (system === "linux") {
            platform = "linux";
        }
        else if (system === "win32") {
            platform = "windows";
        }
        else {
            throw new Error(`Unsupported operating system: ${system}`);
        }
        return platform;
    }
    findBundledBinary() {
        // Skip binary check in test environment
        if (process.env.NODE_ENV === "test") {
            return "mock-txg-binary";
        }
        // Construct executable name based on platform
        const exeName = this.platform === "windows" ? "txg.exe" : "txg";
        // Build path to binary
        const binaryPath = path.join(this.packageDir, "bin", this.platform, exeName);
        // Verify the binary exists
        if (!fs.existsSync(binaryPath)) {
            throw new Error(`txg binary not found at ${binaryPath}`);
        }
        // Ensure it's executable on Unix systems
        if (this.platform !== "windows") {
            try {
                fs.chmodSync(binaryPath, 0o755);
            }
            catch (_) {
                // Ignore chmod errors, binary might already be executable
            }
        }
        return binaryPath;
    }
    async runCommand(args, timeout = QUICK_OPERATION_TIMEOUT, softTimeout) {
        // In test mode, return mock result
        if (process.env.NODE_ENV === "test") {
            return Promise.resolve({
                stdout: "",
                stderr: "",
                exitCode: 0,
                fullCommand: `mock-txg ${args.join(" ")}`,
                inProgress: false,
            });
        }
        return new Promise((resolve, reject) => {
            // Build full command array
            const fullArgs = [...args];
            fullArgs.push("--assumeyes");
            fullArgs.push("--tags", "mcpb");
            const fullCommand = `${this.txgPath} ${fullArgs.join(" ")}`;
            // Log command execution to stderr (MCP best practice for stdio transport)
            console.error(`[TXG CLI] Executing: ${fullCommand.replace(/--access-token\s+\S+/, "--access-token ***")}`);
            const child = spawn(this.txgPath, fullArgs, {
                timeout: timeout,
                env: process.env,
            });
            let stdout = "";
            let stderr = "";
            let resolved = false;
            let timedOut = false;
            // Set hard timeout handler
            const hardTimeoutId = setTimeout(() => {
                timedOut = true;
                child.kill();
                reject(new Error(`Command timed out after ${timeout}ms: ${fullCommand}`));
            }, timeout);
            // Set soft timeout handler if specified
            let softTimeoutId;
            if (softTimeout) {
                softTimeoutId = setTimeout(() => {
                    if (!resolved) {
                        resolved = true;
                        console.error(`[TXG CLI] Soft timeout reached, returning partial result while command continues...`);
                        // Return partial result - command keeps running
                        resolve({
                            stdout: stdout.trim(),
                            stderr: stderr.trim(),
                            exitCode: 0, // Not yet completed
                            inProgress: true,
                            fullCommand,
                        });
                    }
                }, softTimeout);
            }
            child.stdout.on("data", (data) => {
                stdout += data.toString();
            });
            child.stderr.on("data", (data) => {
                stderr += data.toString();
            });
            child.on("error", (error) => {
                clearTimeout(hardTimeoutId);
                if (softTimeoutId) {
                    clearTimeout(softTimeoutId);
                }
                if (!resolved) {
                    resolved = true;
                    console.error(`[TXG CLI] Command error: ${error.message}`);
                    reject(new Error(`Failed to run command: ${error.message}`));
                }
                else {
                    // Already resolved with partial result, just log the error
                    console.error(`[TXG CLI] Background command error: ${error.message}`);
                }
            });
            child.on("close", (code) => {
                clearTimeout(hardTimeoutId);
                if (softTimeoutId) {
                    clearTimeout(softTimeoutId);
                }
                if (!resolved && !timedOut) {
                    resolved = true;
                    const exitCode = code || 0;
                    // Log command result to stderr
                    if (exitCode === 0) {
                        console.error(`[TXG CLI] Command succeeded (exit code: 0)`);
                    }
                    else {
                        console.error(`[TXG CLI] Command failed (exit code: ${exitCode})`);
                    }
                    resolve({
                        stdout: stdout.trim(),
                        stderr: stderr.trim(),
                        exitCode: exitCode,
                        fullCommand,
                        inProgress: false,
                    });
                }
                else if (resolved) {
                    // Command completed after soft timeout - just log it
                    console.error(`[TXG CLI] Background command completed with exit code: ${code || 0}`);
                }
            });
        });
    }
}
export const txgCli = new TxgCliManager();
