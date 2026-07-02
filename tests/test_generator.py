from alphagenome_evo2 import MockGenerator, Sequence


def test_generate_returns_requested_count_from_parents():
    gen = MockGenerator(seed=1)
    parents = [Sequence("p1", "ACGT" * 25), Sequence("p2", "TTGGCCAA" * 12)]
    children = gen.generate(parents, 20)
    assert len(children) == 20
    assert all(set(c.seq) <= set("ACGT") for c in children)


def test_unconditional_generation_when_no_parents():
    gen = MockGenerator(seed=2, default_length=120)
    children = gen.generate([], 5)
    assert len(children) == 5
    assert all(len(c.seq) == 120 for c in children)


def test_generation_is_reproducible_with_same_seed():
    parents = [Sequence("p1", "ACGTACGT" * 10)]
    a = MockGenerator(seed=7).generate(parents, 10)
    b = MockGenerator(seed=7).generate(parents, 10)
    assert [c.seq for c in a] == [c.seq for c in b]


def test_child_ids_are_unique():
    gen = MockGenerator(seed=3)
    parents = [Sequence("p1", "ACGT" * 30), Sequence("p2", "GGCC" * 30)]
    children = gen.generate(parents, 50)
    assert len({c.id for c in children}) == 50
