import torch
import numpy as np
from typing import Callable

def de_mcmc_single_chain(
    d_target: Callable[[torch.Tensor], torch.Tensor],
    theta_0: torch.Tensor,
    epsilon: float,
    n_it: int,
    verbose: bool = True
):
    # Ensure float64 precision
    torch.set_default_dtype(torch.float64)

    d = theta_0.shape[0]
    gamma = 2.38 / torch.sqrt(torch.tensor(2.0 * d, dtype=torch.float64))

    # Chain storage in float64
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

        perturbation = epsilon * (2 * torch.rand(d, dtype=torch.float64) - 1)
        theta_new = theta_local + gamma * delta + perturbation
        target_new = d_target(theta_new).item()

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
        "dTarget": target_local,
        "Theta_0": theta_local,
        "gamma": gamma.item(),
        "epsilon": epsilon,
        "chain": chain
    }
