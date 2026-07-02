"""The real-model adapters must import and construct against injected fakes,
even when the 'alphagenome' / 'evo2' packages are absent."""

from alphagenome_evo2.generator.evo2 import Evo2Generator, _extract_sequences, _sanitise
from alphagenome_evo2.oracle.alphagenome import AlphaGenomeOracle, _pad_to_supported


class _FakeEvo2:
    def generate(self, prompt_seqs, n_tokens, temperature):
        # Echo prompts with a lowercase/junk tail Evo 2 might emit.
        return {"sequences": [p + "nnACGTxx" for p in prompt_seqs]}


class _FakeAGClient:
    def predict_sequence(self, sequence, requested_outputs, ontology_terms=None):
        class _Track:
            def __init__(self, v):
                self.values = v

        class _Out:
            rna_seq = _Track([1.0, 2.0, 3.0])
            atac = _Track([0.5, 0.5])
            chip_tf = _Track([1.0])
            splice_sites = _Track([0.1])

        return _Out()


def test_evo2_adapter_with_injected_model():
    from alphagenome_evo2.types import Sequence

    gen = Evo2Generator(model=_FakeEvo2(), gen_length=8)
    children = gen.generate([Sequence("p", "ACGTACGT")], 3)
    assert len(children) == 3
    assert all(set(c.seq) <= set("ACGT") for c in children)
    assert all(c.metadata["origin"] == "evo2" for c in children)


def test_evo2_helpers():
    assert _sanitise("acGTxx!!") == "ACGT"
    assert _extract_sequences({"sequences": ["AC", "GT"]}) == ["AC", "GT"]
    assert _extract_sequences(["AC", "GT"]) == ["AC", "GT"]


def test_alphagenome_adapter_with_injected_client():
    from alphagenome_evo2.types import CellContext, Sequence

    oracle = AlphaGenomeOracle(client=_FakeAGClient())
    preds = oracle.predict(Sequence("s", "ACGT" * 10), [CellContext("A", "ONT:A")])
    p = preds["A"]
    assert 0.0 <= p.expression <= 1.0
    assert 0.0 <= p.accessibility <= 1.0
    assert p.raw["rna_seq_mean"] == 2.0  # mean of [1,2,3]


def test_pad_to_supported_centers_short_sequences():
    padded = _pad_to_supported("ACGT")
    assert len(padded) == 2048
    assert "ACGT" in padded
