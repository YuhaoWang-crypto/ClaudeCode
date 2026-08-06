r"""
Real signalling-pathway modules: what pattern does a shared pool make?
======================================================================

The two-protein machinery in :mod:`thermo` and :mod:`phases` answers "given
these Kd's and concentrations, what forms?".  This module asks the question the
other way round -- **given the protein combinations that actually sit in a
cell's signalling pathways, what should we expect, and what has to be true for
anything interesting to happen?** -- and answers it with four closed-form
conditions plus an atlas of literature-parameterised systems.

The four conditions
-------------------
Each one is a single inequality, and each is checked against the full Wertheim
solver in ``test_pathways.py``.

**(1) Valence gate (Flory--Stockmayer).**  A network percolates only if

.. math::  (n_A-1)(n_B-1)\,p_A p_B \;\ge\; 1

with :math:`p_a = 1-X_a` the bonded fraction of *a*'s sites.  Since
:math:`p\le1`, a pair of *bivalent* partners can never percolate: it makes
linear polymers, full stop.  **Valence >= 3 on at least one side is not a
tendency, it is a gate.**

**(2) Threshold: sites, not molecules.**  At stoichiometric match the gate above
is crossed when the *site* concentration reaches

.. math::  S^\ast \;=\; K_d\,\frac{p^\ast}{(1-p^\ast)^2},
           \qquad p^\ast=\frac1{\sqrt{(n_A-1)(n_B-1)}}

so the protein threshold is :math:`c_A^\ast = S^\ast/n_A`.  Measured against the
numerical solver this is accurate to ~1--2 %.  The content of the formula is
that **what must reach the Kd is the concentration of binding sites**, which is
:math:`n` times the protein concentration -- which is how 350 uM SH3--PRM bonds
build a network out of low-micromolar proteins.

**(3) Stoichiometry window (re-entrance).**  Even at infinite affinity, bond
conservation :math:`n_A c_A p_A = n_B c_B p_B` caps how lopsided a mixture can
be and still percolate:

.. math::  \frac1{(n_A-1)(n_B-1)} \;\le\; \frac{n_A c_A}{n_B c_B}
           \;\le\; (n_A-1)(n_B-1).

Excess of *either* partner dissolves the network by capping the other's sites
with dead ends.  This is the re-entrant phase boundary, derived rather than
observed, and it is why over-expressing a component can *destroy* the body it
belongs to.

**(4) Miscibility (one body or two).**  Two self-associating modules share a
condensate when cross-binding beats the average of the self-bindings,

.. math::  \frac{2}{K_d^{AB}} \;>\; \frac1{K_d^{AA}}+\frac1{K_d^{BB}}
           \iff K_d^{AB} < \operatorname{harm}(K_d^{AA},K_d^{BB}),

and split into two immiscible bodies otherwise.  Verified against the
common-tangent construction: at :math:`K_d^{AB}` a tenth of the self-Kd the two
proteins share one dense phase; at ten times, the hull returns two essentially
pure dense phases.

And the one that is not a condition on Kd at all
------------------------------------------------
Everything above is equilibrium, and equilibrium forces are **reciprocal**:
:math:`K_d^{AB}=K_d^{BA}` is not a modelling choice, it is detailed balance.  By
the reciprocity theorem proved for Particle Life in ``THEORY.md``, a reciprocal
system has a Lyapunov function and *must come to rest*.  So no combination of
binding constants, however clever, produces the chasing, the snakes or the
rotors.  A pathway that shows persistent directed motion is reporting an
**energy-consuming cycle** -- kinase/phosphatase, GTPase, unfoldase -- and the
condition for that cycle to actually break reciprocity is that it must outrun
the bond it modifies,

.. math::  k_{\rm cat}[E] \;\gtrsim\; k_{\rm off},

otherwise the modification simply equilibrates and the system is reciprocal
again with renormalised constants.

Provenance of the numbers
-------------------------
.. warning::
   The atlas entries are **order-of-magnitude, literature-typical** values, and
   they are here to show which side of each inequality a real system sits on --
   not to be quoted as measurements.  Valence depends on phosphorylation state,
   splice isoform and construct; published Kd's for the same domain pair scatter
   over a decade.  Every entry carries ``source`` and ``confidence`` fields and
   a concentration *range* rather than a point, and :func:`assess` reports the
   verdict across that range so you can see when the call is robust and when it
   is not.  Substitute your own measurements: that is what the dataclass is for.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .sites import Site, SiteMixture, percolates
from .thermo import mode_coherence, spinodal_modes
from .units import conc_to_density, radius_from_mw, radius_idp

__all__ = ["Component", "System", "ATLAS", "get", "to_mixture",
           "percolation_p", "gel_site_conc", "gel_conc", "stoichiometry_window",
           "chi_eff", "miscible", "membrane_boost", "percolates",
           "assess", "boost_to_percolate", "pool_analysis"]


# --------------------------------------------------------------------------- #
#  the four conditions, in closed form
# --------------------------------------------------------------------------- #
def percolation_p(nA, nB=None):
    r"""Bond probability at the gel point, :math:`p^\ast=1/\sqrt{(n_A-1)(n_B-1)}`.

    ``nB=None`` means self-association, for which the classical result is
    :math:`p^\ast = 1/(n-1)`; that is the same formula with ``nB = nA``.
    Returns ``inf`` when the pair can never percolate (``p* >= 1``).
    """
    nA = np.asarray(nA, float)
    nB = nA if nB is None else np.asarray(nB, float)
    g = np.maximum((nA - 1.0) * (nB - 1.0), 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        p = np.where(g > 0, 1.0 / np.sqrt(np.maximum(g, 1e-300)), np.inf)
    return np.where(p >= 1.0, np.inf, p)


def gel_site_conc(kd_uM, nA, nB=None):
    r"""Site concentration (uM) at the gel point: :math:`S^\ast=K_d p^\ast/(1-p^\ast)^2`.

    ``inf`` when the valence gate is shut.  Checked against the Wertheim solver
    over Kd from 1 to 350 uM and valences 2-8: typically ~1 %, worst case 7.4 %
    (at the weakest Kd, where crowding raises ``g^HS`` above the dilute value
    the closed form assumes).
    """
    p = percolation_p(nA, nB)
    kd = np.asarray(kd_uM, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(np.isfinite(p), kd * p / np.maximum((1.0 - p) ** 2, 1e-300), np.inf)


def gel_conc(kd_uM, nA, nB=None):
    """Protein thresholds ``(cA*, cB*)`` in uM at stoichiometric match."""
    S = gel_site_conc(kd_uM, nA, nB)
    nB_ = nA if nB is None else nB
    return S / np.asarray(nA, float), S / np.asarray(nB_, float)


def stoichiometry_window(nA, nB=None):
    r"""Widest site-ratio range :math:`G=(n_A-1)(n_B-1)` that can still percolate.

    Percolation needs :math:`1/G \le n_Ac_A/(n_Bc_B) \le G`; ``G <= 1`` means the
    window is empty.  Measured window edges reach ``G`` to within 4 %.
    """
    nA = np.asarray(nA, float)
    nB = nA if nB is None else np.asarray(nB, float)
    return np.maximum((nA - 1.0) * (nB - 1.0), 0.0)


def chi_eff(kd_aa_uM, kd_bb_uM, kd_ab_uM):
    r"""Demixing parameter :math:`\chi = \tfrac12(1/K_d^{AA}+1/K_d^{BB}) - 1/K_d^{AB}`.

    Positive => self-association wins => two immiscible bodies.
    Negative => cross-association wins => one shared body.
    Units are 1/uM; only the sign and the ratio to ``1/Kd_self`` are meaningful.
    """
    inv = lambda x: np.where(np.isfinite(x) & (np.asarray(x, float) > 0),
                             1.0 / np.maximum(np.asarray(x, float), 1e-300), 0.0)
    return 0.5 * (inv(kd_aa_uM) + inv(kd_bb_uM)) - inv(kd_ab_uM)


def miscible(kd_aa_uM, kd_bb_uM, kd_ab_uM):
    """True when the two modules share one condensate (``chi < 0``)."""
    return chi_eff(kd_aa_uM, kd_bb_uM, kd_ab_uM) < 0


def membrane_boost(cell_volume_um3=2000.0, membrane_area_um2=1500.0, shell_nm=10.0):
    r"""How much recruiting to a surface multiplies the local concentration.

    A protein spread through the cytosol at 1 uM, once captured within
    ``shell_nm`` of the plasma membrane, occupies a shell of volume
    ``area x shell``; the concentration ratio is ``V_cell / V_shell``.  For a
    mammalian cell that is ~130x, and for a small patch of clustered receptor it
    is far more.

    This is the single most important number in the whole module, because the
    threshold condition (2) is a condition on *local* site concentration.  It is
    why LAT/Grb2/SOS phase separates on a membrane at nanomolar bulk protein and
    never in free solution, and why so many pathways begin with "recruit to a
    surface".
    """
    v_shell = np.asarray(membrane_area_um2, float) * np.asarray(shell_nm, float) * 1e-3
    return np.asarray(cell_volume_um3, float) / np.maximum(v_shell, 1e-300)


# --------------------------------------------------------------------------- #
#  system description
# --------------------------------------------------------------------------- #
@dataclass
class Component:
    """One protein (or nucleic acid) in a pathway module.

    ``domains`` is a ``{name: multiplicity}`` map of *distinct* binding-site
    classes -- ``{"SH2": 1, "SH3": 2}`` for Grb2.  Collapsing those into a single
    valence is what makes a model invent networks that do not exist, so the
    architecture is carried explicitly (see :mod:`plife.bio.sites`).
    """

    name: str
    mw_kda: float
    domains: dict                        # {site-class name: multiplicity}
    conc_uM: tuple                       # (lo, hi) physiological range
    idp: bool = False                    # expanded (disordered) rather than globular
    n_res: int = 0                       # residues, for the IDP radius
    note: str = ""

    @property
    def radius_nm(self):
        if self.idp and self.n_res:
            return float(radius_idp(self.n_res))
        return float(radius_from_mw(self.mw_kda))

    @property
    def valence(self):
        """Total binding sites, for reporting."""
        return float(sum(self.domains.values()))

    @property
    def c_mid(self):
        return float(np.sqrt(self.conc_uM[0] * self.conc_uM[1]))   # geometric centre


@dataclass
class System:
    """A pathway module: components, the site-pair Kd's, and provenance.

    ``bonds`` maps ``("Protein.Domain", "Protein.Domain") -> Kd in uM``.  Only
    the pairs listed bind; everything else is infinite.  Writing the bonds out
    domain by domain rather than protein by protein is the whole point: it is
    what stops a one-SH2 adaptor from being modelled as a two-SH2 crosslinker.
    """

    key: str
    name: str
    name_zh: str
    family: str
    components: list
    bonds: dict
    membrane: bool = False               # does the pathway recruit it to a surface?
    source: str = ""
    confidence: str = "order-of-magnitude"
    note_zh: str = ""

    @property
    def names(self):
        return [c.name for c in self.components]

    @property
    def site_list(self):
        return [Site(c.name, d, m) for c in self.components for d, m in c.domains.items()]

    @property
    def valence(self):
        return np.array([c.valence for c in self.components])

    @property
    def kd_uM(self):
        """Site-pair Kd matrix (A,A) in uM built from ``bonds``."""
        labels = [s.label for s in self.site_list]
        idx = {l: i for i, l in enumerate(labels)}
        A = len(labels)
        kd = np.full((A, A), np.inf)
        for (u, v), val in self.bonds.items():
            if u not in idx or v not in idx:
                raise KeyError(f"{self.key}: bond {(u, v)} names an unknown site; "
                               f"have {labels}")
            kd[idx[u], idx[v]] = kd[idx[v], idx[u]] = val
        return kd

    def conc(self, where="mid"):
        """Physiological concentrations (uM): ``lo``, ``mid`` or ``hi``."""
        i = {"lo": 0, "hi": 1}.get(where)
        if i is None:
            return np.array([c.c_mid for c in self.components])
        return np.array([c.conc_uM[i] for c in self.components])


def to_mixture(sys: System) -> SiteMixture:
    """Convert a :class:`System` into a site-resolved mixture."""
    return SiteMixture(names=list(sys.names),
                       diameters=np.array([2 * c.radius_nm for c in sys.components]),
                       sites=sys.site_list, kd=sys.kd_uM * 1e-6)


# --------------------------------------------------------------------------- #
#  the atlas
# --------------------------------------------------------------------------- #
ATLAS = [
    System(
        key="nephrin", name="Nephrin - NCK - N-WASP", name_zh="肾病蛋白 - NCK - N-WASP",
        family="actin nucleation / slit diaphragm",
        components=[
            Component("pNephrin", 25, {"pY": 3}, (1.0, 50.0), idp=True, n_res=180,
                      note="3 phosphotyrosines; the kinase sets the valence"),
            Component("NCK1", 47, {"SH2": 1, "SH3": 3}, (1.0, 10.0),
                      note="1 SH2 for pY, 3 SH3 for PRM"),
            Component("N-WASP", 65, {"PRM": 6}, (0.3, 5.0), note="6 proline-rich motifs"),
        ],
        bonds={("pNephrin.pY", "NCK1.SH2"): 1.0,
               ("NCK1.SH3", "N-WASP.PRM"): 350.0},
        membrane=True,
        source="Li et al. Nature 483:336 (2012); Banjade & Rosen eLife 3:e04123 (2014)",
        note_zh="教科书级例子：SH3-PRM 单键弱到 350 uM，靠 3x6 的价数把网络撑起来。"
                "磷酸化直接改价数（0→3），是最干净的开关。",
    ),
    System(
        key="lat", name="pLAT - Grb2 - SOS1", name_zh="pLAT - Grb2 - SOS1",
        family="T-cell receptor",
        components=[
            Component("pLAT", 25, {"pY": 4}, (0.5, 5.0), idp=True, n_res=210,
                      note="Y132/Y171/Y191/Y226; membrane-anchored"),
            Component("Grb2", 25, {"SH2": 1, "SH3": 2}, (1.0, 5.0),
                      note="SH3-SH2-SH3: ONE SH2, so it cannot bridge two LATs"),
            Component("SOS1", 150, {"PRM": 4}, (0.05, 0.5), note="4 proline-rich motifs"),
        ],
        bonds={("pLAT.pY", "Grb2.SH2"): 1.0,
               ("Grb2.SH3", "SOS1.PRM"): 10.0},
        membrane=True,
        source="Su et al. Science 352:595 (2016); Houtman et al. NSMB 13:798 (2006)",
        note_zh="Grb2 只有一个 SH2——如果把它当成三个等价位点，模型会凭空造出"
                "一个不存在的网络。这条是本模块要做位点分辨的直接原因。",
    ),
    System(
        key="psd", name="SynGAP - PSD-95", name_zh="SynGAP - PSD-95",
        family="postsynaptic density",
        components=[
            Component("SynGAP", 150, {"PBM": 3}, (30.0, 150.0),
                      note="coiled-coil trimer -> 3 PDZ-binding motifs"),
            Component("PSD-95", 95, {"PDZ": 3}, (30.0, 200.0), note="PDZ1-3"),
        ],
        bonds={("SynGAP.PBM", "PSD-95.PDZ"): 1.0},
        source="Zeng et al. Cell 166:1163 (2016)",
        note_zh="3x3 刚好过价数门槛，且突触后致密区局部浓度 ~100 uM，稳稳在阈值之上。",
    ),
    System(
        key="spop", name="SPOP - DAXX", name_zh="SPOP - DAXX",
        family="ubiquitin ligase / nuclear body",
        components=[
            Component("SPOP", 42, {"MATH": 6}, (0.5, 5.0),
                      note="BTB+BACK oligomer presents ~6 substrate sites"),
            Component("DAXX", 81, {"SBC": 5}, (0.5, 5.0), idp=True, n_res=740,
                      note="~5 SPOP-binding consensus motifs"),
        ],
        bonds={("SPOP.MATH", "DAXX.SBC"): 10.0},
        source="Bouchard et al. Mol Cell 72:19 (2018)",
        note_zh="底物自己提供价数：SPOP 的寡聚化 x DAXX 的多个降解子 = 6x5，门槛很低。",
    ),
    System(
        key="wnt", name="Dishevelled DIX - Axin DAX", name_zh="Dishevelled DIX - Axin DAX",
        family="Wnt signalosome",
        components=[
            Component("Dvl2", 78, {"head": 1, "tail": 1}, (0.1, 1.0),
                      note="DIX polymerises head-to-tail: one of each, never more"),
            Component("Axin1", 96, {"head": 1, "tail": 1}, (0.01, 0.1),
                      note="DAX, same head-to-tail geometry"),
        ],
        bonds={("Dvl2.head", "Dvl2.tail"): 20.0,
               ("Axin1.head", "Axin1.tail"): 20.0,
               ("Dvl2.head", "Axin1.tail"): 20.0,
               ("Axin1.head", "Dvl2.tail"): 20.0},
        source="Schwarz-Romond et al. NSMB 14:484 (2007); Bienz TiBS 39:487 (2014)",
        note_zh="价数门槛的反例：头尾聚合每个分子只有一头一尾，走到任何一个分子上"
                "只剩一条出路，分支因子恒 <= 1，永远长不出三维网络。所以 Wnt "
                "signalosome 是线性聚合物/点状 puncta，不是凝聚相。要成网必须另加"
                "交联（DEP 二聚、受体成簇）。",
    ),
    System(
        key="p62", name="p62/SQSTM1 - K63 polyubiquitin", name_zh="p62/SQSTM1 - K63 多聚泛素",
        family="selective autophagy",
        components=[
            Component("p62", 62, {"UBA": 6}, (1.0, 10.0),
                      note="PB1 helical polymer presents many UBA per filament unit"),
            Component("K63-Ub4", 34, {"Ile44": 4}, (0.5, 10.0),
                      note="valence = ubiquitin chain length"),
        ],
        bonds={("p62.UBA", "K63-Ub4.Ile44"): 100.0},
        source="Sun et al. Cell Res 28:405 (2018); Zaffagnini et al. EMBO J 37:e98308 (2018)",
        note_zh="UBA-Ub 单键弱到 100 uM，但 p62 自聚成丝把价数堆得很高。"
                "泛素链越长越容易成体——链长就是价数。",
    ),
    System(
        key="sg", name="G3BP1 - mRNA", name_zh="G3BP1 - mRNA",
        family="stress granule",
        components=[
            Component("G3BP1", 104, {"RRM": 4}, (1.0, 5.0),
                      note="dimer, 2 RNA-binding surfaces per protomer"),
            Component("mRNA-2kb", 660, {"site": 15}, (0.05, 1.0), idp=True, n_res=2000,
                      note="valence = accessible sites, i.e. unfolded length"),
        ],
        bonds={("G3BP1.RRM", "mRNA-2kb.site"): 1.0},
        source="Yang et al. Cell 181:325 (2020); Guillen-Boixet et al. Cell 181:346 (2020)",
        note_zh="RNA 提供巨大价数，这是应激颗粒的核心。翻译停滞、核糖体脱落"
                "释放出裸 mRNA = 价数瞬间上升 = 越过门槛。",
    ),
    System(
        key="nucleolus", name="NPM1 - R-motif partner", name_zh="NPM1 - R-motif 伙伴",
        family="nucleolar granular component",
        components=[
            Component("NPM1", 165, {"acidic": 5}, (5.0, 30.0),
                      note="pentamer -> 5 acidic tracts"),
            Component("SURF6-like", 40, {"Rmotif": 4}, (1.0, 20.0), idp=True, n_res=360,
                      note="multiple arginine-rich motifs"),
        ],
        bonds={("NPM1.acidic", "SURF6-like.Rmotif"): 5.0},
        source="Mitrea et al. eLife 5:e13571 (2016); Feric et al. Cell 165:1686 (2016)",
        note_zh="五聚体 x 多 R-motif，核仁局部浓度又高，属于稳稳成相的一类。",
    ),
    System(
        key="cgas", name="cGAS - dsDNA", name_zh="cGAS - 双链 DNA",
        family="innate immunity / STING",
        components=[
            Component("cGAS", 60, {"DNA": 2}, (0.05, 0.5), note="two DNA-binding surfaces"),
            Component("dsDNA-1kb", 620, {"turn": 50}, (0.001, 0.1),
                      note="valence = length / ~20 bp: THIS is the ligand-length sensor"),
        ],
        bonds={("cGAS.DNA", "dsDNA-1kb.turn"): 0.5},
        source="Du & Chen Science 361:704 (2018); Andreeva et al. Nature 549:394 (2017)",
        note_zh="二价的 cGAS 配上高价数的长 DNA：(2-1)(50-1)=49，门槛极低。"
                "短 DNA 价数低就点不着——DNA 长度依赖的免疫激活，本质就是价数门槛。",
    ),
    System(
        key="pka", name="PKA-RII - AKAP  [negative control]", name_zh="PKA-RII - AKAP【阴性对照】",
        family="cAMP / scaffolding",
        components=[
            Component("PKA-RII2", 90, {"DD": 1}, (0.2, 2.0),
                      note="one AKAP-helix docking surface per R-subunit dimer"),
            Component("AKAP79", 79, {"helix": 1}, (0.05, 0.5), note="one amphipathic helix"),
        ],
        bonds={("PKA-RII2.DD", "AKAP79.helix"): 0.001},
        source="Gold et al. Mol Cell 24:383 (2006)",
        note_zh="亲和力 1 nM，比表里任何一条都强三个数量级——但价数是 1，"
                "所以只有 1:1 复合物，永远没有凝聚体。最能说明：亲和力 != 相分离。",
    ),
    System(
        key="mapk", name="Raf - MEK - ERK core  [negative control]",
        name_zh="Raf - MEK - ERK 核心【阴性对照】",
        family="MAPK cascade",
        components=[
            Component("BRAF", 84, {"cat": 1}, (0.05, 0.5)),
            Component("MEK1", 43, {"nterm": 1, "cterm": 1}, (0.5, 3.0),
                      note="one surface up to Raf, one down to ERK"),
            Component("ERK2", 42, {"dock": 1}, (0.5, 5.0)),
        ],
        bonds={("BRAF.cat", "MEK1.nterm"): 0.3,
               ("MEK1.cterm", "ERK2.dock"): 2.0},
        source="standard MAPK biochemistry; Kd's are typical kinase-substrate values",
        note_zh="经典级联是一条链：每个分子最多一进一出，分支因子 <= 1，"
                "任何浓度都不成相。要成体必须加 KSR 之类的支架把价数抬上去。",
    ),
    System(
        key="mavs", name="MAVS CARD filament", name_zh="MAVS CARD 纤维",
        family="innate immunity / RIG-I",
        components=[
            Component("MAVS-CARD", 10, {"head": 1, "tail": 1}, (0.1, 1.0),
                      note="head-to-tail helical polymer"),
        ],
        bonds={("MAVS-CARD.head", "MAVS-CARD.tail"): 0.01},
        source="Hou et al. Cell 146:448 (2011); Cai et al. Cell 156:1207 (2014)",
        note_zh="一头一尾 + 极强 + 成核依赖 = 一维纤维（prion 样），不是液滴。"
                "价数门槛说它不能成三维网络，实验看到的确实是丝状体和全或无的数字开关。",
    ),
]

_BY_KEY = {s.key: s for s in ATLAS}


def get(key) -> System:
    """Look up an atlas entry by key."""
    if key not in _BY_KEY:
        raise KeyError(f"unknown system {key!r}; have {sorted(_BY_KEY)}")
    return _BY_KEY[key]


# --------------------------------------------------------------------------- #
#  assessment
# --------------------------------------------------------------------------- #
def assess(sys: System, where="mid", local_boost=1.0):
    r"""Run all the conditions on one atlas entry and return a verdict dict.

    The closed-form conditions (1)-(3) are **two-body** statements, so they are
    reported for the *load-bearing bond*: among the site pairs that can branch at
    all (both multivalent, so :math:`(n_a-1)(n_b-1)>1`), the one sitting highest
    relative to its own gel threshold.  Pairs involving a monovalent site class
    are excluded from that contest on purpose -- they extend a chain but cannot
    branch, so their threshold is infinite and they would otherwise always
    "win" while explaining nothing.

    For a module of three or more components the two-body number is only a
    bound, and usually a *loose* one: a scaffold that ties several copies of a
    partner together raises the effective valence of the pair, so the real
    threshold is lower.  Nephrin is the clean case -- the SH3/PRM bond alone
    would need ~30x more material than the assembled three-component network
    does.  The authoritative answer is therefore ``branching``, the spectral
    radius of the full site-resolved matrix, with ``boost_to_percolate`` giving
    the local enrichment that would switch the module on.
    """
    sm = to_mixture(sys)
    c = sys.conc(where) * local_boost
    labels = sm.labels
    site_uM = c[sm.owner] * sm.n            # site concentration, uM

    ok, lam, bonded = percolates(c, sm)

    # load-bearing bond: best margin among pairs that can actually branch
    best, any_bond = None, None
    kd = sys.kd_uM
    for i in range(sm.A):
        for j in range(i, sm.A):
            if not np.isfinite(kd[i, j]):
                continue
            nA, nB = float(sm.n[i]), float(sm.n[j])
            S_star = float(gel_site_conc(kd[i, j], nA, nB))
            avail = float(min(site_uM[i], site_uM[j]))
            G = float(stoichiometry_window(nA, nB))
            ratio = site_uM[i] / site_uM[j] if site_uM[j] > 0 else np.inf
            rec = dict(pair=(labels[i], labels[j]), kd_uM=float(kd[i, j]),
                       valence_pair=(nA, nB), p_star=float(percolation_p(nA, nB)),
                       S_gel_uM=S_star, site_conc_uM=avail,
                       margin=(avail / S_star if np.isfinite(S_star) and S_star > 0 else 0.0),
                       G=G, site_ratio=float(ratio),
                       in_window=bool(G > 1 and 1.0 / G <= ratio <= G))
            if any_bond is None or rec["kd_uM"] < any_bond["kd_uM"]:
                any_bond = rec
            if np.isfinite(S_star) and (best is None or rec["margin"] > best["margin"]):
                best = rec
    pairwise_gate = best is not None
    if best is None:                       # no pair can branch on its own
        best = any_bond or dict(
            pair=("-", "-"), kd_uM=np.inf, valence_pair=(0.0, 0.0), p_star=np.inf,
            S_gel_uM=np.inf, site_conc_uM=0.0, margin=0.0, G=0.0,
            site_ratio=np.inf, in_window=False)

    return dict(
        key=sys.key, name=sys.name, name_zh=sys.name_zh, family=sys.family,
        where=where, local_boost=local_boost,
        conc_uM=c, valence=sys.valence, site_labels=labels, site_uM=site_uM,
        gate_open=bool(pairwise_gate),
        c_gel_uM=(best["S_gel_uM"] / max(best["valence_pair"][0], 1e-9)
                  if np.isfinite(best["S_gel_uM"]) else np.inf),
        percolates=bool(ok), branching=float(lam), bonded=bonded,
        membrane=sys.membrane, source=sys.source, note_zh=sys.note_zh,
        tightest_pair=best["pair"], **{k: v for k, v in best.items() if k != "pair"},
    )


def boost_to_percolate(sys: System, where="mid", lo=1.0, hi=1e6, tol=1e-3):
    r"""**How much local concentration would this module need to switch on?**

    Returns the smallest uniform concentration multiplier at which the branching
    eigenvalue reaches 1, or ``inf`` if the valence gate is shut so that no
    concentration suffices.  Because the multiplier is a *local* enrichment, the
    answer is directly comparable to :func:`membrane_boost` and to measured
    receptor-microcluster densities -- which turns "does it phase separate?"
    into a falsifiable number.
    """
    sm = to_mixture(sys)
    c = sys.conc(where)
    if percolates(c * hi, sm)[1] < 1.0:
        return np.inf
    if percolates(c * lo, sm)[0]:
        return lo
    while hi / lo > 1 + tol:
        mid = np.sqrt(lo * hi)
        if percolates(c * mid, sm)[0]:
            hi = mid
        else:
            lo = mid
    return float(np.sqrt(lo * hi))


def pool_analysis(systems=None, where="mid", local_boost=1.0, cross=None, modes=True):
    r"""**Put every module in one pool: how many bodies, and who is in each?**

    Two complementary readings, because neither alone is enough:

    *Pairwise* -- for each pair of modules, :func:`chi_eff` says whether they
    would share a body.  Modules are then grouped by the connected components of
    the "would share" graph, which predicts the number of distinct bodies.

    *Global* -- the eigen-decomposition of the full multicomponent stability
    matrix at the pool composition.  Every negative eigenvalue is a growing
    composition fluctuation; :func:`~plife.bio.thermo.mode_coherence` on its
    eigenvector says whether that mode condenses everything together or splits
    the pool.  This is the same order parameter used on the Particle Life
    dispersion relation, and it is what makes the two halves of this project one
    project.

    ``cross`` is an optional ``{(key_a, key_b): kd_uM}`` map of cross-pathway
    binding.  Its **default of nothing** is not a shortcut but the hypothesis
    being tested: signalling modules are built from orthogonal domain-motif
    pairs, and if that orthogonality holds then :func:`chi_eff` is positive for
    every pair and the pool must sort into one body per module rather than one
    common blob.  Supplying a cross Kd tight enough to beat the harmonic mean of
    the two self-Kd's merges those modules -- which is the failure mode that
    promiscuous over-expression actually produces.

    ``modes=False`` skips the global eigen-decomposition, which is essentially
    the whole cost of the call (a numerical Hessian over a mixture with dozens
    of species).  Use it when only the body count is wanted -- e.g. when
    sweeping ``cross`` over many values.
    """
    systems = ATLAS if systems is None else systems
    cross = dict(cross or {})
    per = [assess(s, where=where, local_boost=local_boost) for s in systems]

    # ---- build one big mixture: every component of every module ------------
    names, diam, conc, all_sites = [], [], [], []
    for s in systems:
        cc = s.conc(where) * local_boost
        for k, comp in enumerate(s.components):
            names.append(f"{s.key}:{comp.name}")
            diam.append(2 * comp.radius_nm)
            conc.append(cc[k])
        for st in s.site_list:
            all_sites.append(Site(f"{s.key}:{st.species}", st.name, st.n))
    S = len(names)
    A = len(all_sites)
    kd = np.full((A, A), np.inf)
    blocks, o = {}, 0
    for s in systems:
        m = len(s.site_list)
        kd[o:o + m, o:o + m] = s.kd_uM
        blocks[s.key] = (o, o + m)
        o += m
    for (ka, kb), v in cross.items():
        (i0, i1), (j0, j1) = blocks[ka], blocks[kb]
        kd[i0:i1, j0:j1] = v
        kd[j0:j1, i0:i1] = v
    big = SiteMixture(names=names, diameters=np.array(diam), sites=all_sites,
                      kd=kd * 1e-6)
    conc = np.array(conc)
    rho = conc_to_density(conc[None, :])

    mode_list, vals = [], np.array([])
    if modes:
        vals, vecs = spinodal_modes(rho, big)
        vals, vecs = vals[0], vecs[0]
        for i in np.where(vals < 0)[0]:
            v = vecs[:, i]
            coh = float(mode_coherence(v))
            top = np.argsort(-np.abs(v))[:4]
            mode_list.append(dict(eig=float(vals[i]), coherence=coh,
                                  leaders=[(names[t], float(v[t])) for t in top],
                                  kind="condensation" if coh > 0.7 else
                                       ("demixing" if coh < 0.3 else "mixed")))

    # ---- pairwise miscibility over the modules that actually percolate -----
    live = [i for i, p in enumerate(per) if p["percolates"]]
    tight = {s.key: min(s.bonds.values()) if s.bonds else np.inf for s in systems}
    edges = []
    for ii, a in enumerate(live):
        for b in live[ii + 1:]:
            ka, kb = systems[a].key, systems[b].key
            kab = cross.get((ka, kb), cross.get((kb, ka), np.inf))
            chi = float(chi_eff(tight[ka], tight[kb], kab))
            edges.append(dict(a=ka, b=kb, kd_cross_uM=float(kab), chi=chi,
                              shares=bool(chi < 0)))

    # connected components of the "would share a body" graph
    parent = {systems[i].key: systems[i].key for i in live}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for e in edges:
        if e["shares"]:
            ra, rb = find(e["a"]), find(e["b"])
            if ra != rb:
                parent[ra] = rb
    groups = {}
    for i in live:
        groups.setdefault(find(systems[i].key), []).append(systems[i].key)
    groups = sorted(groups.values(), key=lambda g: (-len(g), g[0]))

    return dict(per_system=per, n_species=S, names=names, conc_uM=conc,
                eigenvalues=vals, n_unstable=int((vals < 0).sum()), modes=mode_list,
                live=[systems[i].key for i in live], groups=groups, edges=edges,
                n_bodies=len(groups), mixture=big)
