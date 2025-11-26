import torch
import numpy as np
from typing import Callable

def de_mcmc_single_chain_sparse(
    d_target: Callable[[torch.Tensor], torch.Tensor],
    theta_0: torch.Tensor,
    epsilon: float,
    n_it: int,
    p_sparse: float = 0.05,
    scale_sparse: float = 0.5,
    verbose: bool = True
):
    torch.set_default_dtype(torch.float64)

    d = theta_0.shape[0]
    gamma = 2.38 / torch.sqrt(torch.tensor(2.0 * d, dtype=torch.float64))

    chain = torch.zeros((n_it, d + 1), dtype=torch.float64)
    theta_local = theta_0.clone().to(torch.float64)
    target_local = d_target(theta_local).item()

    chain[0, :] = torch.cat((torch.tensor([target_local], dtype=torch.float64), theta_local))

    accepted = 0

    for k in range(1, n_it):
        idx_range = torch.arange(max(k // 2, 1))
        idx_range_np = idx_range.numpy()

        if len(idx_range_np) < 2:
            idx1 = idx2 = int(idx_range_np[0])
        else:
            idx1, idx2 = np.random.choice(idx_range_np, size=2, replace=False)

        delta = chain[idx1, 1:] - chain[idx2, 1:]

        # Standard noise (epsilon)
        epsilon_noise = epsilon * (2 * torch.rand(d, dtype=torch.float64) - 1)

        # Sparse noise (ksi)
        sparse_mask = torch.bernoulli(torch.full((d,), p_sparse, dtype=torch.float64))
        sparse_noise = scale_sparse * (2 * torch.rand(d, dtype=torch.float64) - 1)
        ksi = sparse_mask * sparse_noise

        # Proposal
        theta_new = theta_local + gamma * delta + epsilon_noise + ksi
        target_new = d_target(theta_new).item()

        # Metropolis-Hastings acceptance
        r = np.exp(target_new - target_local)
        if np.random.rand() < r:
            accepted += 1
            theta_local = theta_new
            target_local = target_new

        chain[k, :] = torch.cat((torch.tensor([target_local], dtype=torch.float64), theta_local))

        if verbose and k % (n_it // 10) == 0:
            print(f"{k}/{n_it} | {target_local:.6f} | {theta_local[0]:.6f} | {theta_local[1]:.6f}")

    print(f"Acceptance rate: {accepted / n_it:.4f}")

    return {
        "dTarget": torch.tensor([target_local], dtype=torch.float64)[0],
        "Theta_0": theta_local.to(torch.float64),
        "gamma": gamma,  # keep as tensor
        "epsilon": torch.tensor(epsilon, dtype=torch.float64),
        "p_sparse": torch.tensor(p_sparse, dtype=torch.float64),
        "scale_sparse": torch.tensor(scale_sparse, dtype=torch.float64),
        "chain": chain.to(torch.float64)
    }

def de_mcmc_single_chain_sparse_optimiser(
    d_target: Callable[[torch.Tensor], torch.Tensor],
    theta_0: torch.Tensor,
    epsilon: float,
    n_it: int,
    p_sparse: float = 0.05,
    scale_sparse: float = 0.5,
    memory: int = 100,
    verbose: bool = True
):
    torch.set_default_dtype(torch.float64)

    d = theta_0.shape[0]
    gamma = 2.38 / torch.sqrt(torch.tensor(2.0 * d, dtype=torch.float64))

    chain = torch.zeros((n_it, d + 1), dtype=torch.float64)
    theta_local = theta_0.clone().to(torch.float64)
    target_local = d_target(theta_local).item()

    chain[0, :] = torch.cat((torch.tensor([target_local], dtype=torch.float64), theta_local))

    accepted = 0

    for k in range(1, n_it):        
        idx_range = torch.arange(max(k - memory, 0), max(k, 1))
        idx_range_np = idx_range.numpy()

        if len(idx_range_np) < 2:
            idx1 = idx2 = int(idx_range_np[0])
        else:
            ## pick random elements
            idx1, idx2 = sorted(np.random.choice(idx_range_np, size=2, replace=False))
            # idx1, idx2 = np.random.choice(idx_range_np, size=2, replace=False)

        delta = chain[idx2, 1:] - chain[idx1, 1:] # newer - older
        # delta = chain[idx1, 1:] - chain[idx2, 1:]

        # Standard noise (epsilon)
        epsilon_noise = epsilon * (2 * torch.rand(d, dtype=torch.float64) - 1)

        # Sparse noise (ksi)
        sparse_mask = torch.bernoulli(torch.full((d,), p_sparse, dtype=torch.float64))
        sparse_noise = scale_sparse * (2 * torch.rand(d, dtype=torch.float64) - 1)
        ksi = sparse_mask * sparse_noise

        # Proposal
        theta_new = theta_local + gamma * delta + epsilon_noise + ksi
        target_new = d_target(theta_new).item()

        # Metropolis-Hastings acceptance
        # r = np.exp(target_new - target_local)
        if target_new > target_local:
            accepted += 1
            theta_local = theta_new
            target_local = target_new

        chain[k, :] = torch.cat((torch.tensor([target_local], dtype=torch.float64), theta_local))

        if verbose and k % (n_it // 10) == 0:
            print(f"{k}/{n_it} | {target_local:.6f} | {theta_local[0]:.6f} | {theta_local[1]:.6f}")

    print(f"Acceptance rate: {accepted / n_it:.4f}")

    return {
        "dTarget": torch.tensor([target_local], dtype=torch.float64)[0],
        "Theta_0": theta_local.to(torch.float64),
        "gamma": gamma,  # keep as tensor
        "epsilon": torch.tensor(epsilon, dtype=torch.float64),
        "p_sparse": torch.tensor(p_sparse, dtype=torch.float64),
        "scale_sparse": torch.tensor(scale_sparse, dtype=torch.float64),
        "chain": chain.to(torch.float64)
    }