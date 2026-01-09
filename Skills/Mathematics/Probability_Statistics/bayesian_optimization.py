import numpy as np
from typing import List, Tuple, Callable, Optional
from dataclasses import dataclass

@dataclass
class OptimizationResult:
    best_params: np.ndarray
    best_value: float
    all_history: List[Tuple[np.ndarray, float]]

class BayesianOptimizer:
    """
    A lightweight implementation of Bayesian Optimization using Gaussian Processes (GP).
    Crucial for 'Self-Driving Labs' to select the next experiment that maximizes information gain
    (Upper Confidence Bound - UCB).
    """
    
    def __init__(self, bounds: List[Tuple[float, float]], kappa: float = 2.5):
        """
        :param bounds: List of (min, max) for each dimension.
        :param kappa: Exploration-exploitation balance (higher = more exploration).
        """
        self.bounds = np.array(bounds)
        self.dim = len(bounds)
        self.kappa = kappa
        self.X_sample: List[np.ndarray] = []
        self.Y_sample: List[float] = []
        
        # Hyperparameters for the RBF Kernel
        self.sigma_f = 1.0
        self.length_scale = 1.0
        self.noise = 1e-5

    def _kernel(self, X1: np.ndarray, X2: np.ndarray) -> np.ndarray:
        """
        Radial Basis Function (RBF) / Squared Exponential Kernel.
        Computes covariance between points.
        """
        sqdist = np.sum(X1**2, 1).reshape(-1, 1) + np.sum(X2**2, 1) - 2 * np.dot(X1, X2.T)
        return self.sigma_f**2 * np.exp(-0.5 / self.length_scale**2 * sqdist)

    def fit(self, X: List[np.ndarray], Y: List[float]):
        """Update the internal Gaussian Process with new data."""
        self.X_sample = np.array(X)
        self.Y_sample = np.array(Y).reshape(-1, 1)

    def predict(self, X: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Predicts mean and variance for a new set of points X using the GP.
        """
        if len(self.X_sample) == 0:
            return np.zeros((len(X), 1)), np.ones((len(X), 1))

        K = self._kernel(self.X_sample, self.X_sample) + self.noise * np.eye(len(self.X_sample))
        K_s = self._kernel(self.X_sample, X)
        K_ss = self._kernel(X, X) + 1e-8 * np.eye(len(X))
        
        K_inv = np.linalg.inv(K)
        
        # Mean prediction
        mu = K_s.T.dot(K_inv).dot(self.Y_sample)
        
        # Variance prediction
        cov = K_ss - K_s.T.dot(K_inv).dot(K_s)
        return mu, np.diag(cov).reshape(-1, 1)

    def acquisition_function(self, X: np.ndarray) -> np.ndarray:
        """
        Upper Confidence Bound (UCB).
        Selects points with high mean (exploitation) or high variance (exploration).
        """
        mu, sigma = self.predict(X)
        return mu + self.kappa * np.sqrt(sigma)

    def propose_next_location(self, num_candidates: int = 100) -> np.ndarray:
        """
        Randomly samples the search space and returns the point with the highest Acquisition Score.
        In production, use L-BFGS-B to maximize the acquisition function properly.
        """
        # Uniform random sampling within bounds
        candidates = np.random.uniform(
            self.bounds[:, 0], self.bounds[:, 1], size=(num_candidates, self.dim)
        )
        scores = self.acquisition_function(candidates)
        best_idx = np.argmax(scores)
        return candidates[best_idx]

# --- Example Usage ---
if __name__ == "__main__":
    # Define a "Black Box" function (e.g., experimental yield)
    # Goal: Maximize this function
    def synthetic_experiment(params):
        # Simple 1D function: f(x) = -(x-2)^2 + 10. Max is at x=2, value=10.
        x = params[0]
        return -(x - 2.0)**2 + 10.0

    # Initialize Optimizer for 1D space between [-5, 5]
    opt = BayesianOptimizer(bounds=[(-5.0, 5.0)])
    
    print("--- Starting Self-Driving Lab Simulation ---")
    
    # Initial random observations
    X_init = [np.array([0.0]), np.array([4.0])]
    Y_init = [synthetic_experiment(x) for x in X_init]
    opt.fit(X_init, Y_init)
    
    # Run 5 optimization steps
    for step in range(5):
        next_param = opt.propose_next_location()
        actual_value = synthetic_experiment(next_param)
        
        print(f"Step {step+1}: Suggested Param {next_param}, Result {actual_value:.4f}")
        
        # Add to known data
        X_init.append(next_param)
        Y_init.append(actual_value)
        opt.fit(X_init, Y_init)

    best_y = max(Y_init)
    print(f"Optimization Complete. Best Value Found: {best_y:.4f}")