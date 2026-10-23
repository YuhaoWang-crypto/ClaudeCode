"""
DCSC — Disentangled Cell State Code: a feasibility prototype.

    synth     synthetic single-cell world with known ground truth
    codec     64-byte factorized barcode (ID | Mod ternary | Mem | Env) + encoder
    index     ternary-Hamming search (numpy / C POPCNT / faiss binary + HNSW)
    delta     Δ-state signatures and residual pursuit (cocktail recipes)
    partial   partial-view queries (gene panels, cytokines) with calibrated error bars
    twostage  two-stage retrieval: compact-code shortlist -> exact rescoring / residual pursuit
    donor_match  donor-level calls: healthy / match disease X / no similar diseased donor in library
    run_demo  end-to-end scoring; writes dcsc/RESULTS.md

See DCSC.md at the repository root for the feasibility assessment.
"""
