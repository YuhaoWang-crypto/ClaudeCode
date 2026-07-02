"""AlphaGenome-backed oracle.

Wraps Google DeepMind's ``alphagenome`` client so the design loop can use real
sequence-to-function predictions. The heavy dependency is imported lazily, so the
rest of the package (and the mock loop) works without it installed.

Usage::

    from alphagenome_evo2.oracle.alphagenome import AlphaGenomeOracle
    oracle = AlphaGenomeOracle(api_key=os.environ["ALPHAGENOME_API_KEY"])

Install with ``pip install alphagenome`` and request an API key from
https://deepmind.google.com/science/alphagenome . The adapter maps AlphaGenome
output tracks onto the pipeline's normalised :class:`OraclePrediction` fields:

* ``expression``    ← aggregate RNA_SEQ / CAGE signal over the sequence
* ``accessibility`` ← aggregate ATAC / DNASE signal
* ``tf_binding``    ← aggregate CHIP_TF signal
* ``splice_anomaly``← divergence of SPLICE_SITES from the reference expectation

The exact track names / call signatures track the public ``alphagenome`` API; if
a future version renames them, only this file needs to change.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence as Seq

from ..types import CellContext, OraclePrediction, Sequence
from .base import Oracle

# Sequence lengths AlphaGenome accepts; we pad/centre to the nearest supported one.
_SUPPORTED_LENGTHS = (2048, 16384, 131072, 524288, 1048576)


def _pad_to_supported(seq: str) -> str:
    """Centre ``seq`` in an N-padded window of the smallest supported length."""
    target = next((L for L in _SUPPORTED_LENGTHS if L >= len(seq)), _SUPPORTED_LENGTHS[-1])
    if len(seq) >= target:
        # Centre-crop overly long inputs.
        start = (len(seq) - target) // 2
        return seq[start:start + target]
    pad = target - len(seq)
    left = pad // 2
    return "N" * left + seq + "N" * (pad - left)


def _mean(values) -> float:
    values = list(values)
    return sum(values) / len(values) if values else 0.0


def _normalise(x: float, scale: float) -> float:
    """Squash an unbounded non-negative track mean into [0, 1)."""
    if scale <= 0:
        return 0.0
    return 1.0 - pow(2.718281828, -max(0.0, x) / scale)


class AlphaGenomeOracle(Oracle):
    def __init__(
        self,
        api_key: Optional[str] = None,
        client=None,
        expression_scale: float = 5.0,
        accessibility_scale: float = 3.0,
        tf_scale: float = 3.0,
    ):
        """Create the oracle.

        Pass ``api_key`` to build a client, or inject a pre-built ``client``
        (useful for testing). The ``*_scale`` values calibrate the squashing of
        raw track means into [0, 1]; tune per assay if needed.
        """
        if client is not None:
            self._client = client
        else:
            self._client = self._make_client(api_key)
        self.expression_scale = expression_scale
        self.accessibility_scale = accessibility_scale
        self.tf_scale = tf_scale

    @staticmethod
    def _make_client(api_key: Optional[str]):
        try:
            from alphagenome.models import dna_client  # type: ignore
        except ImportError as exc:  # pragma: no cover - exercised only w/o dep
            raise ImportError(
                "AlphaGenomeOracle requires the 'alphagenome' package. "
                "Install it with `pip install alphagenome` and obtain an API key, "
                "or use MockOracle for offline development."
            ) from exc
        if not api_key:
            raise ValueError("api_key is required to construct an AlphaGenome client")
        return dna_client.create(api_key)

    def predict(
        self,
        sequence: Sequence,
        contexts: Seq[CellContext],
    ) -> Dict[str, OraclePrediction]:
        requested = self._requested_outputs()
        padded = _pad_to_supported(sequence.seq)
        out: Dict[str, OraclePrediction] = {}
        for ctx in contexts:
            output = self._client.predict_sequence(
                sequence=padded,
                requested_outputs=requested,
                ontology_terms=[ctx.ontology_term] if ctx.ontology_term else None,
            )
            out[ctx.name] = self._to_prediction(ctx, output)
        return out

    @staticmethod
    def _requested_outputs():
        """The AlphaGenome output tracks we consume.

        Returns the ``dna_client.OutputType`` enums when the package is present,
        otherwise ``None`` (which injected/fake clients ignore) so tests can run
        without the heavy dependency installed.
        """
        try:
            from alphagenome.models import dna_client  # type: ignore
        except ImportError:
            return None
        return [
            dna_client.OutputType.RNA_SEQ,
            dna_client.OutputType.ATAC,
            dna_client.OutputType.CHIP_TF,
            dna_client.OutputType.SPLICE_SITES,
        ]

    def _to_prediction(self, ctx: CellContext, output) -> OraclePrediction:
        """Reduce AlphaGenome track tensors to normalised scalar read-outs."""
        expr_raw = _track_mean(output, "rna_seq")
        acc_raw = _track_mean(output, "atac")
        tf_raw = _track_mean(output, "chip_tf")
        splice_raw = _track_mean(output, "splice_sites")

        return OraclePrediction(
            context=ctx,
            expression=_normalise(expr_raw, self.expression_scale),
            accessibility=_normalise(acc_raw, self.accessibility_scale),
            tf_binding=_normalise(tf_raw, self.tf_scale),
            # A crude anomaly proxy; refine against a reference distribution.
            splice_anomaly=min(1.0, splice_raw),
            raw={
                "rna_seq_mean": float(expr_raw),
                "atac_mean": float(acc_raw),
                "chip_tf_mean": float(tf_raw),
                "splice_mean": float(splice_raw),
            },
        )


def _track_mean(output, attr: str) -> float:
    """Best-effort extraction of a track's mean signal from an AlphaGenome output.

    Kept defensive because the concrete container shape can vary across
    ``alphagenome`` versions; on any mismatch we fall back to 0.0 rather than
    crashing an entire optimisation run.
    """
    track = getattr(output, attr, None)
    if track is None:
        return 0.0
    values = getattr(track, "values", track)
    try:
        flat = values.reshape(-1)  # numpy / tensor
        return float(flat.mean())
    except AttributeError:
        try:
            return _mean(float(v) for v in values)
        except TypeError:  # pragma: no cover - defensive
            return 0.0
