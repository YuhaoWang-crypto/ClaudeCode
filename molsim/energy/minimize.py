"""Energy minimization using L-BFGS-B via scipy."""
from __future__ import annotations
import numpy as np
from scipy.optimize import minimize as scipy_minimize
from typing import Optional, Callable, Tuple
from ..core import Molecule
from .calculator import EnergyCalculator, nonbonded_energy


def minimize_energy(
    mol: Molecule,
    receptor: Optional[Molecule] = None,
    max_iter: int = 500,
    tol: float = 1e-4,
    verbose: bool = False,
) -> dict:
    """
    Minimize energy of mol (optionally in field of receptor).

    If receptor is given, minimizes receptor-ligand interaction energy.
    Otherwise minimizes intramolecular non-bonded energy of mol.

    Returns dict: {positions, energy, n_iter, success}.
    """
    x0 = mol.positions.flatten()

    if receptor is not None:
        calc = EnergyCalculator(receptor=receptor, ligand=mol)

        def func_grad(x: np.ndarray) -> Tuple[float, np.ndarray]:
            coords = x.reshape(-1, 3)
            e, g = calc.score_and_grad(coords)
            return e, g.flatten()

    else:
        def func_grad(x: np.ndarray) -> Tuple[float, np.ndarray]:
            coords = x.reshape(-1, 3)
            mol.positions = coords
            e_vdw, e_elec = nonbonded_energy(mol)
            e = e_vdw + e_elec
            # Numerical gradient (fallback for intramolecular)
            eps_fd = 1e-4
            g = np.zeros_like(x)
            for k in range(len(x)):
                x_plus = x.copy(); x_plus[k] += eps_fd
                mol.positions = x_plus.reshape(-1, 3)
                ep, _ = nonbonded_energy(mol)
                g[k] = (ep - e) / eps_fd
            mol.positions = coords
            return e, g

    result = scipy_minimize(
        fun=func_grad,
        x0=x0,
        method='L-BFGS-B',
        jac=True,
        options={'maxiter': max_iter, 'ftol': tol, 'gtol': tol * 0.1,
                 'disp': verbose},
    )

    mol.positions = result.x.reshape(-1, 3)
    return {
        'positions': mol.positions.copy(),
        'energy': float(result.fun),
        'n_iter': result.nit,
        'success': result.success,
        'message': result.message,
    }
