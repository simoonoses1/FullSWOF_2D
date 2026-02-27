from __future__ import annotations

import time

import numpy as np

from fullswof_oil.infiltration import InfiltrationModel
from fullswof_oil.solver import SolverConfig, run_simulation


def make_case(n: int) -> tuple[np.ndarray, np.ndarray]:
    h0 = np.zeros((n, n), dtype=float)
    c0 = n // 2 - max(2, n // 20)
    c1 = n // 2 + max(2, n // 20)
    h0[c0:c1, c0:c1] = 0.03
    z = np.zeros_like(h0)
    return h0, z


def walltime_for_mode(mode: str, n: int, repeats: int = 1) -> float:
    h0, z = make_case(n)
    elapsed: list[float] = []
    for _ in range(repeats):
        cfg = SolverConfig(
            nx=n,
            ny=n,
            dx=1.0,
            dy=1.0,
            t_end=20.0,
            cfl=0.45,
            execution_mode=mode,
        )
        t0 = time.perf_counter()
        run_simulation(cfg, h0, z=z, infiltration=InfiltrationModel(0.0, 0.0, 0.0))
        elapsed.append(time.perf_counter() - t0)
    return float(np.mean(elapsed))


def main() -> None:
    sizes = [64, 128, 192]
    print("size,numpy_s,numba_s,speedup_numba_vs_numpy")
    for n in sizes:
        t_numpy = walltime_for_mode("numpy", n)
        t_numba = walltime_for_mode("numba", n)
        speedup = t_numpy / t_numba if t_numba > 0 else float("inf")
        print(f"{n}x{n},{t_numpy:.4f},{t_numba:.4f},{speedup:.2f}")


if __name__ == "__main__":
    main()
