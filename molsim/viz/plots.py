"""Matplotlib-based visualizations for molecular analysis results."""
from __future__ import annotations
import numpy as np
from typing import List, Dict, Optional, Tuple

try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import matplotlib.colors as mcolors
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

_SS_COLORS = {'H': '#FF4444', 'E': '#4444FF', 'T': '#FFAA00', 'C': '#AAAAAA'}
_SS_LABELS = {'H': 'Helix', 'E': 'Sheet', 'T': 'Turn', 'C': 'Coil'}


def _check_mpl():
    if not HAS_MPL:
        raise ImportError("matplotlib is required for visualization. "
                          "Install with: pip install matplotlib")


def plot_ramachandran(
    angles: Dict[int, Dict],
    title: str = "Ramachandran Plot",
    ax=None,
) -> 'plt.Axes':
    """
    Ramachandran plot (phi vs psi for all residues).

    Parameters
    ----------
    angles : dict from ramachandran_angles()
    """
    _check_mpl()
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))

    phi_vals, psi_vals, labels = [], [], []
    for rid, d in angles.items():
        if 'phi' in d and 'psi' in d:
            phi_vals.append(d['phi'])
            psi_vals.append(d['psi'])
            labels.append(d.get('resname', '?'))

    # Background regions (approx favored/allowed)
    # Alpha helix: phi ~-60, psi ~-50  |  Beta sheet: phi ~-120, psi ~120
    from matplotlib.patches import Ellipse
    ax.add_patch(Ellipse((-65, -40), 60, 60, alpha=0.15, color='red', label='_α-helix region'))
    ax.add_patch(Ellipse((-115, 130), 80, 70, alpha=0.15, color='blue', label='_β-sheet region'))

    sc = ax.scatter(phi_vals, psi_vals, s=30, c='#333333', alpha=0.7, zorder=3)

    ax.axhline(0, color='gray', linewidth=0.5, linestyle='--')
    ax.axvline(0, color='gray', linewidth=0.5, linestyle='--')
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.set_xlabel('φ (degrees)', fontsize=12)
    ax.set_ylabel('ψ (degrees)', fontsize=12)
    ax.set_title(title, fontsize=13)
    ax.set_aspect('equal')

    # Custom legend
    alpha_p = mpatches.Patch(color='red', alpha=0.3, label='α-helix')
    beta_p = mpatches.Patch(color='blue', alpha=0.3, label='β-sheet')
    ax.legend(handles=[alpha_p, beta_p], loc='upper right')

    n_total = len(phi_vals)
    n_alpha = sum(1 for ph, ps in zip(phi_vals, psi_vals) if -100 < ph < -30 and -70 < ps < 10)
    n_beta = sum(1 for ph, ps in zip(phi_vals, psi_vals) if -160 < ph < -60 and 80 < ps < 180)
    ax.set_title(f"{title}\n({n_total} residues, {n_alpha} α-helical, {n_beta} β-strand)", fontsize=11)

    return ax


def plot_secondary_structure(
    ss: Dict[int, str],
    title: str = "Secondary Structure",
    ax=None,
) -> 'plt.Axes':
    """
    Linear secondary structure diagram.

    Parameters
    ----------
    ss : dict from secondary_structure()
    """
    _check_mpl()
    if ax is None:
        fig, ax = plt.subplots(figsize=(14, 2))

    res_ids = sorted(ss.keys())
    ss_codes = [ss[r] for r in res_ids]

    for i, (rid, code) in enumerate(zip(res_ids, ss_codes)):
        color = _SS_COLORS.get(code, '#CCCCCC')
        ax.barh(0, 1, left=i, color=color, height=0.8, linewidth=0)

    patches = [mpatches.Patch(color=_SS_COLORS[k], label=_SS_LABELS[k])
               for k in ('H', 'E', 'T', 'C')]
    ax.legend(handles=patches, loc='upper right', fontsize=9)
    ax.set_xlim(0, len(res_ids))
    ax.set_ylim(-0.5, 0.5)
    ax.set_yticks([])
    ax.set_xlabel('Residue index', fontsize=11)
    ax.set_title(title)

    # Label every 10th residue
    ticks = range(0, len(res_ids), max(1, len(res_ids) // 20))
    ax.set_xticks(list(ticks))
    ax.set_xticklabels([str(res_ids[t]) for t in ticks], fontsize=8)

    return ax


def plot_energy_decomposition(
    energy_dict: Dict[str, float],
    title: str = "Energy Decomposition",
    ax=None,
) -> 'plt.Axes':
    """Bar chart of energy components."""
    _check_mpl()
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 4))

    labels = list(energy_dict.keys())
    values = list(energy_dict.values())
    colors = ['#E74C3C' if v > 0 else '#2ECC71' for v in values]

    bars = ax.bar(labels, values, color=colors, edgecolor='black', linewidth=0.5)
    ax.axhline(0, color='black', linewidth=0.8)

    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2,
                val + (0.5 if val >= 0 else -1.5),
                f'{val:.1f}', ha='center', va='bottom', fontsize=9)

    ax.set_ylabel('Energy (kcal/mol)', fontsize=11)
    ax.set_title(title, fontsize=12)
    ax.grid(axis='y', alpha=0.3)
    return ax


def plot_rmsd(
    rmsd_values: np.ndarray,
    timestep: float = 1.0,
    unit: str = 'step',
    title: str = 'RMSD',
    ax=None,
) -> 'plt.Axes':
    """Time series RMSD plot."""
    _check_mpl()
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 3))

    t = np.arange(len(rmsd_values)) * timestep
    ax.plot(t, rmsd_values, color='#2980B9', linewidth=1.5)
    ax.fill_between(t, rmsd_values, alpha=0.2, color='#2980B9')
    ax.set_xlabel(f'Time ({unit})', fontsize=11)
    ax.set_ylabel('RMSD (Å)', fontsize=11)
    ax.set_title(title, fontsize=12)
    ax.grid(alpha=0.3)
    return ax


def plot_rmsf(
    rmsf_values: np.ndarray,
    residue_ids: Optional[List[int]] = None,
    title: str = 'Per-residue RMSF',
    ax=None,
) -> 'plt.Axes':
    """Per-residue RMSF plot."""
    _check_mpl()
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 3))

    x = residue_ids if residue_ids is not None else np.arange(len(rmsf_values))
    ax.bar(x, rmsf_values, color='#E67E22', alpha=0.8, linewidth=0)
    ax.set_xlabel('Residue ID', fontsize=11)
    ax.set_ylabel('RMSF (Å)', fontsize=11)
    ax.set_title(title, fontsize=12)
    ax.grid(axis='y', alpha=0.3)
    return ax


def plot_docking_results(
    results: list,
    top_n: int = 10,
    ax=None,
) -> 'plt.Axes':
    """
    Bar chart of docking scores for top poses.

    Parameters
    ----------
    results : list of DockingResult
    """
    _check_mpl()
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 4))

    top = results[:top_n]
    ranks = [r.rank for r in top]
    scores = [r.score for r in top]
    vdw = [r.vdw for r in top]
    elec = [r.elec for r in top]
    hb = [r.hbond for r in top]

    x = np.arange(len(top))
    width = 0.6
    bars_vdw = ax.bar(x, vdw, width, label='VdW', color='#3498DB')
    bars_elec = ax.bar(x, elec, width, bottom=vdw, label='Electrostatic', color='#E74C3C')
    bars_hb = ax.bar(x, hb, width,
                     bottom=np.array(vdw) + np.array(elec),
                     label='H-bond', color='#2ECC71')

    ax.plot(x, scores, 'ko-', markersize=4, linewidth=1.5, label='Total', zorder=5)
    ax.axhline(0, color='black', linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels([f'Pose {r}' for r in ranks], rotation=30, ha='right')
    ax.set_ylabel('Score (kcal/mol)', fontsize=11)
    ax.set_title(f'Top {len(top)} Docking Poses', fontsize=12)
    ax.legend(loc='best', fontsize=9)
    ax.grid(axis='y', alpha=0.3)
    return ax


def plot_contact_map(
    contacts: list,
    title: str = 'Residue Contact Map',
    ax=None,
) -> 'plt.Axes':
    """
    Contact map from residue_contacts() output.

    Parameters
    ----------
    contacts : list of (rid1, rid2, dist) from residue_contacts()
    """
    _check_mpl()
    if not contacts:
        return

    r1s = [c[0] for c in contacts]
    r2s = [c[1] for c in contacts]
    dists = [c[2] for c in contacts]

    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 5))

    sc = ax.scatter(r1s, r2s, c=dists, cmap='viridis_r', s=15, vmin=0, vmax=8)
    ax.scatter(r2s, r1s, c=dists, cmap='viridis_r', s=15, vmin=0, vmax=8)
    plt.colorbar(sc, ax=ax, label='Distance (Å)')
    ax.set_xlabel('Residue (mol1)', fontsize=11)
    ax.set_ylabel('Residue (mol2)', fontsize=11)
    ax.set_title(title, fontsize=12)
    return ax
