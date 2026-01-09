import numpy as np
from typing import Optional, Tuple

class TensorOps:
    """
    Fundamental linear algebra operations for understanding Transformers and attention mechanisms.
    While frameworks like PyTorch handle this, understanding the math is key for 'custom attention'
    in biological sequences (e.g., protein folding).
    """

    @staticmethod
    def scaled_dot_product_attention(
        query: np.ndarray, 
        key: np.ndarray, 
        value: np.ndarray, 
        mask: Optional[np.ndarray] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Computes the core 'Attention' mechanism: softmax(QK^T / sqrt(d_k))
        
        Args:
            query: (Batch, Seq_Len, d_k)
            key:   (Batch, Seq_Len, d_k)
            value: (Batch, Seq_Len, d_v)
            mask:  Optional mask (Batch, Seq_Len, Seq_Len) for causal/padding masking.
            
        Returns:
            output: (Batch, Seq_Len, d_v)
            weights: (Batch, Seq_Len, Seq_Len)
        """
        d_k = query.shape[-1]
        
        # 1. MatMul Q * K^T
        # Swap last two dimensions of Key for transpose
        scores = np.matmul(query, key.swapaxes(-2, -1)) / np.sqrt(d_k)
        
        # 2. Apply Mask (optional)
        if mask is not None:
            # -1e9 represents negative infinity for softmax
            scores = np.where(mask == 0, -1e9, scores)
            
        # 3. Softmax
        # Stabilize by subtracting max
        exp_scores = np.exp(scores - np.max(scores, axis=-1, keepdims=True))
        weights = exp_scores / np.sum(exp_scores, axis=-1, keepdims=True)
        
        # 4. MatMul Weights * V
        output = np.matmul(weights, value)
        
        return output, weights

    @staticmethod
    def create_causal_mask(seq_len: int) -> np.ndarray:
        """
        Creates a triangular mask (lower left = 1) to prevent attending to future tokens.
        Crucial for generative models (GPT).
        """
        # np.tril returns lower triangle
        return np.tril(np.ones((seq_len, seq_len)))

# --- Example Usage ---
if __name__ == "__main__":
    ops = TensorOps()
    
    # Simulate a batch of 1 sequence, length 4, dimension 8
    # "The cat sat on" -> Q, K, V
    Q = np.random.randn(1, 4, 8)
    K = np.random.randn(1, 4, 8)
    V = np.random.randn(1, 4, 8)
    
    # 1. Self-Attention (No Mask)
    out, attn = ops.scaled_dot_product_attention(Q, K, V)
    print("--- Self-Attention Output Shape ---")
    print(out.shape) # (1, 4, 8)
    
    # 2. Causal Attention (Masked)
    mask = ops.create_causal_mask(4)
    out_masked, attn_masked = ops.scaled_dot_product_attention(Q, K, V, mask=mask)
    
    print("\n--- Causal Attention Weights (First Head) ---")
    print(np.round(attn_masked[0], 2))
    # You should see upper triangle as 0.0
