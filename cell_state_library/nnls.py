"""K-reference NNLS scoring. Python standard library only.

Same method as v0.1 (and the author's R function): every reference column and the
target are scaled to unit length on the matched genes, then min ||A x - b||^2, x >= 0.
v0.1 handled exactly two references; this version handles K references by checking
every active set, which is exact and fast for the small K a cell type has.
"""
import itertools
import math


def _dot(x, y):
    return sum(u * v for u, v in zip(x, y))


def _solve(m, rhs):
    """Gaussian elimination with partial pivoting. Returns None if singular."""
    n = len(m)
    a = [list(row) + [r] for row, r in zip(m, rhs)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(a[r][c]))
        if abs(a[p][c]) < 1e-12:
            return None
        a[c], a[p] = a[p], a[c]
        for r in range(c + 1, n):
            f = a[r][c] / a[c][c]
            for k in range(c, n + 1):
                a[r][k] -= f * a[c][k]
    x = [0.0] * n
    for r in range(n - 1, -1, -1):
        x[r] = (a[r][n] - sum(a[r][k] * x[k] for k in range(r + 1, n))) / a[r][r]
    return x


def nnls_exact(gram, atb, btb):
    """Exact NNLS from normal equations by enumerating active sets.

    gram = A^T A (K x K), atb = A^T b, btb = b^T b. Returns (x, loss).
    The optimum satisfies the KKT conditions on some active set, so the best
    feasible unconstrained solution over all subsets is the global optimum.
    """
    k = len(atb)
    best_x, best_loss = [0.0] * k, btb
    for size in range(1, k + 1):
        for subset in itertools.combinations(range(k), size):
            g = [[gram[i][j] for j in subset] for i in subset]
            xs = _solve(g, [atb[i] for i in subset])
            if xs is None or min(xs) < -1e-12:
                continue
            x = [0.0] * k
            for i, v in zip(subset, xs):
                x[i] = max(v, 0.0)
            loss = btb - 2 * _dot(atb, x) + sum(x[i] * gram[i][j] * x[j] for i in range(k) for j in range(k))
            if loss < best_loss - 1e-15:
                best_x, best_loss = x, loss
    return best_x, max(best_loss, 0.0)


def _sym_eigvals(m, sweeps=50):
    """Eigenvalues of a small symmetric matrix (cyclic Jacobi)."""
    a = [list(r) for r in m]
    n = len(a)
    for _ in range(sweeps):
        off = sum(a[i][j] ** 2 for i in range(n) for j in range(n) if i != j)
        if off < 1e-24:
            break
        for p in range(n):
            for q in range(p + 1, n):
                if abs(a[p][q]) < 1e-30:
                    continue
                theta = (a[q][q] - a[p][p]) / (2 * a[p][q])
                t = math.copysign(1.0, theta) / (abs(theta) + math.sqrt(theta * theta + 1))
                c = 1 / math.sqrt(t * t + 1)
                s = t * c
                for k in range(n):
                    akp, akq = a[k][p], a[k][q]
                    a[k][p], a[k][q] = c * akp - s * akq, s * akp + c * akq
                for k in range(n):
                    apk, aqk = a[p][k], a[q][k]
                    a[p][k], a[q][k] = c * apk - s * aqk, s * apk + c * aqk
    return [a[i][i] for i in range(n)]


def score_nnls(references, target, min_genes=30, min_coverage=0.5, enforce_gate=True):
    """Score a target log2FC profile against K reference profiles of one cell type.

    references: {name: {gene: log2fc}}; all references must share one gene list.
    target: {gene: log2fc}. Missing genes are dropped, never set to zero.
    """
    names = sorted(references)
    gene_sets = [set(references[n]) for n in names]
    ref_genes = sorted(gene_sets[0])
    if any(s != gene_sets[0] for s in gene_sets):
        raise ValueError("References of one cell type must share one gene list.")
    for v in target.values():
        if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v):
            raise ValueError("All log2FC values must be finite numbers; missing is not zero.")
    genes = [g for g in ref_genes if g in target]
    coverage = len(genes) / len(ref_genes)
    if enforce_gate and (len(genes) < min_genes or coverage < min_coverage):
        raise ValueError(
            f"Insufficient overlap: {len(genes)}/{len(ref_genes)}. Prototype gate: >={min_genes} genes "
            f"and >={min_coverage:.0%} reference coverage; this is an engineering gate, not validated sensitivity.")
    cols = [[references[n][g] for g in genes] for n in names]
    b = [target[g] for g in genes]
    norms = [math.sqrt(_dot(c, c)) for c in cols] + [math.sqrt(_dot(b, b))]
    if min(norms) <= 1e-12:
        raise ValueError("Zero-norm reference or target; no meaningful score.")
    cols = [[t / n for t in c] for c, n in zip(cols, norms)]
    b = [t / norms[-1] for t in b]
    gram = [[_dot(ci, cj) for cj in cols] for ci in cols]
    eig = _sym_eigvals(gram)
    if min(eig) < 1e-8:
        raise ValueError("References are collinear on matched genes; cannot separate scores.")
    x, _ = nnls_exact(gram, [_dot(c, b) for c in cols], 1.0)
    # residual from the vectors themselves (the normal-equation loss loses precision near 0)
    loss = sum((sum(xk * c[i] for xk, c in zip(x, cols)) - bi) ** 2 for i, bi in enumerate(b))
    pairwise = {f"{names[i]}~{names[j]}": gram[i][j] for i in range(len(names)) for j in range(i + 1, len(names))}
    return {
        "scores": dict(zip(names, x)),
        "matched_genes": len(genes),
        "reference_genes": len(ref_genes),
        "coverage": coverage,
        "relative_residual_norm": math.sqrt(loss),
        "reference_cosines": pairwise,
        # sqrt of eigenvalue ratio = condition number of the normalized reference matrix
        "reference_condition_number": math.sqrt(max(eig) / min(eig)),
    }
