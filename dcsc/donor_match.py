"""
Donor-level matching with an explicit "no similar diseased donor" answer.

Input: lab-centred donor profiles (each lab's healthy-control mean already
subtracted; see dcsc/real/donor_search.py). For a query donor the matcher
returns one of three calls:

    "healthy"             profile is inside the normal range of healthy donors
    "match"               abnormal, and its nearest diseased library donors are
                          about as close as same-disease donors from DIFFERENT
                          labs usually are -> disease label = majority of them
    "no_similar_disease"  abnormal, but nothing in the library is that close:
                          the library does not contain this disease state

Both decisions are conformal-style tests calibrated on the library itself,
always across labs (the independence unit, e.g. a Census collection):

  abnormality   Mahalanobis distance to the healthy distribution (Ledoit-Wolf
                covariance of healthy profiles). Reference = scores of healthy
                library donors, each computed with a covariance fitted WITHOUT
                its own lab. p_abnormal = share of reference scores >= query.
  match         s = best cosine similarity to a library donor of the predicted
                disease. Reference = for every diseased library donor whose
                disease occurs in >= 2 labs, its best cosine to a same-disease
                donor of ANOTHER lab. p_match = share of reference <= s.
"""
import numpy as np
from sklearn.covariance import LedoitWolf


class DonorMatcher:
    def __init__(self, X, disease, lab, healthy_label="normal", k=10,
                 alpha_abnormal=0.10, alpha_match=0.20):
        """Defaults: ~10% false alarms on healthy donors; on real Census donors
        this withholds ~4% of correct matches and turns ~40% of would-be wrong
        labels for unseen diseases into "no_similar_disease" (REAL_RESULTS §6.3).
        Raise alpha_match to be stricter about matches."""
        self.X = np.asarray(X, float)
        self.disease = np.asarray(disease).astype(str)
        self.lab = np.asarray(lab).astype(str)
        self.k, self.alpha_abnormal, self.alpha_match = k, alpha_abnormal, alpha_match
        self.healthy = self.disease == healthy_label
        self.Xn = self.X / (np.linalg.norm(self.X, axis=1, keepdims=True) + 1e-9)
        H = self.X[self.healthy]
        self.prec = LedoitWolf().fit(H).precision_
        # abnormality reference: cross-fitted by lab
        ref = []
        for L in np.unique(self.lab[self.healthy]):
            inL = self.healthy & (self.lab == L)
            out = self.healthy & (self.lab != L)
            if out.sum() < 10:
                continue
            P = LedoitWolf().fit(self.X[out]).precision_
            Xi = self.X[inL]
            ref.append(np.einsum("ij,jk,ik->i", Xi, P, Xi))
        self.ref_abnormal = np.sort(np.concatenate(ref))
        # match reference: cross-lab nearest same-disease similarity
        sick = ~self.healthy
        ref = []
        for i in np.where(sick)[0]:
            same = sick & (self.disease == self.disease[i]) & (self.lab != self.lab[i])
            if same.any():
                ref.append((self.Xn[same] @ self.Xn[i]).max())
        self.ref_match = np.sort(np.array(ref))
        self.n_diseases_multi_lab = len({d for d in self.disease[sick]
                                         if len(set(self.lab[sick & (self.disease == d)])) >= 2})

    def p_abnormal(self, x):
        d2 = float(x @ self.prec @ x)
        n = len(self.ref_abnormal)
        return (1 + n - np.searchsorted(self.ref_abnormal, d2, side="left")) / (1 + n), d2

    def query(self, x, lab=None):
        x = np.asarray(x, float)
        p_abn, d2 = self.p_abnormal(x)
        out = dict(p_abnormal=p_abn, mahalanobis2=d2)
        pool = ~self.healthy if lab is None else (~self.healthy & (self.lab != str(lab)))
        idx = np.where(pool)[0]
        xn = x / (np.linalg.norm(x) + 1e-9)
        s = self.Xn[idx] @ xn
        nb = idx[np.argsort(-s)[:self.k]]
        votes = {}
        for j in nb:
            votes[self.disease[j]] = votes.get(self.disease[j], 0) + max(self.Xn[j] @ xn, 0)
        pred = max(votes, key=votes.get)
        s_best = float(s[self.disease[idx] == pred].max())
        n = len(self.ref_match)
        p_match = (1 + np.searchsorted(self.ref_match, s_best, side="right")) / (1 + n)
        out.update(nearest_disease=pred, best_similarity=s_best, p_match=float(p_match),
                   neighbours=[(self.lab[j], self.disease[j], float(self.Xn[j] @ xn)) for j in nb[:5]])
        if p_abn >= self.alpha_abnormal:
            out["call"] = "healthy"
        elif p_match < self.alpha_match:
            out["call"] = "no_similar_disease"
        else:
            out["call"] = "match"
        return out
