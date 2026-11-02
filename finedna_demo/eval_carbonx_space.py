"""Quiz the live Carbon-X model on a Hugging Face Space with FineDNA questions.

    python finedna_demo/eval_carbonx_space.py --per-type 3

Carbon-X has no hosted inference endpoint, so this calls a public Gradio Space
that runs it on ZeroGPU (default: mbcredits/carbon-x-demo, endpoint /ask_carbon).
Set HF_TOKEN to use your own ZeroGPU quota.

For each task type it asks two kinds of question:
  * seen  - the rl_mix prompt as-is. Carbon-X was RL-trained on rl_mix, so these
            are training questions.
  * fresh - the same prompt with its DNA swapped for a random sequence of the same
            length; the gold answer is recomputed by verifiers.py. The model has
            never seen these, so they test computing, not recall.
Replies are cached in results/carbonx_space_eval.jsonl, so a run that hits the
GPU quota can be resumed later.
"""
import argparse
import collections
import json
import os
import random
import re
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from verifiers import CHECKERS, compare, longest  # noqa: E402

SPACE, API = "mbcredits/carbon-x-demo", "/ask_carbon"
RL_CACHE = HERE / "cache" / "rl_mix_5000_0.jsonl"
OUT = HERE / "results" / "carbonx_space_eval.jsonl"
# Task types whose prompt makes no claim about the sequence, so any random DNA is valid.
FRESH_OK = ["gc_content", "reverse_complement", "aligned_complement", "transcription_coding",
            "pyrimidine_count", "purine_count", "cpg_count", "pam_count",
            "hydrogen_bond_count", "wallace_tm", "long_oligo_tm", "complete_codons"]
SEQ_ANSWER = {"reverse_complement", "aligned_complement", "transcription_coding"}


def load_rl_rows():
    if not RL_CACHE.exists():
        sys.exit("Run finedna_demo.py first: it streams and caches the rl_mix sample.")
    return [json.loads(line) for line in RL_CACHE.open()]


def make_fresh(prompt, rnd):
    """Swap the prompt's DNA for a random sequence of the same length."""
    old = longest(prompt)
    new = "".join(rnd.choice("ACGT") for _ in old)
    return prompt.replace(old, new)


def has_final_answer(reply):
    """False when the reply ran out of tokens before an 'Answer:' line."""
    return any(ln.strip().lower().startswith("answer") for ln in reply.splitlines())


def extract_answer(reply, task_type):
    lines = [ln.strip() for ln in reply.strip().splitlines() if ln.strip()]
    ans = next((ln.split(":", 1)[1] for ln in reversed(lines) if ln.lower().startswith("answer")),
               lines[-1] if lines else "")
    ans = ans.strip().strip("*`. ")
    if task_type in SEQ_ANSWER:
        return longest(ans, min_len=1) or ans
    m = re.search(r"-?\d+(?:\.\d+)?", ans)
    if m and task_type not in ("translation_coding",):
        unit = "%" if "%" in ans else ""
        return m.group(0) + unit
    return ans


def build_questions(per_type, seed=0):
    rnd = random.Random(seed)
    by_type = collections.defaultdict(list)
    for r in load_rl_rows():
        if r["task_group"] == "basic_dna" and r["task_type"] in FRESH_OK:
            by_type[r["task_type"]].append(r)
    qs = []
    for k in range(per_type):  # round-robin so a partial run still covers every type
        for t in FRESH_OK:
            if k >= len(by_type[t]):
                continue
            r = by_type[t][k]
            p = r["prompt"][-1]["content"]
            qs.append({"qid": f"{t}-{k}-seen", "task_type": t, "kind": "seen",
                       "prompt": p, "gold": r["answer_text"]})
            fp = make_fresh(p, rnd)
            qs.append({"qid": f"{t}-{k}-fresh", "task_type": t, "kind": "fresh",
                       "prompt": fp, "gold": CHECKERS[t](fp)})
    return qs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--per-type", type=int, default=3, help="seen+fresh pairs per task type")
    ap.add_argument("--space", default=SPACE)
    ap.add_argument("--max-new-tokens", type=int, default=1024, help="the Space caps this at 1024")
    args = ap.parse_args()

    done = {}
    if OUT.exists():
        for line in OUT.open():
            d = json.loads(line)
            if d["correct"] or has_final_answer(d["reply"]):  # re-ask replies that were cut off
                done[d["qid"]] = d
    qs = [q for q in build_questions(args.per_type) if q["qid"] not in done]
    print(f"{len(done)} cached, {len(qs)} to ask on {args.space}", flush=True)

    if qs:
        from gradio_client import Client
        client = Client(args.space, token=os.environ.get("HF_TOKEN"), verbose=False)
        with OUT.open("a") as f:
            for q in qs:
                t0 = time.time()
                try:
                    reply = client.predict(q["prompt"], "", "", "", "extra", "", args.max_new_tokens, api_name=API)
                except Exception as e:  # quota exhausted, Space asleep, ...
                    print(f"  stopped at {q['qid']}: {type(e).__name__}: {str(e)[:200]}", flush=True)
                    break
                ans = extract_answer(reply, q["task_type"])
                rec = {**q, "reply": reply, "model_answer": ans, "correct": compare(ans, q["gold"]),
                       "truncated": not has_final_answer(reply),
                       "seconds": round(time.time() - t0, 1)}
                f.write(json.dumps(rec, ensure_ascii=False) + "\n"); f.flush()
                done[q["qid"]] = rec
                mark = "OK " if rec["correct"] else "CUT" if rec["truncated"] else "BAD"
                print(f"  {q['qid']:<28} {mark} gold={q['gold'][:30]!r} "
                      f"model={ans[:30]!r} ({rec['seconds']}s)", flush=True)

    recs = list(done.values())
    cut = [r["qid"] for r in recs if not r["correct"] and not has_final_answer(r["reply"])]
    table = collections.defaultdict(lambda: {"seen": [0, 0], "fresh": [0, 0]})
    for r in recs:
        if r["qid"] in cut:
            continue  # ran out of tokens: neither right nor wrong
        cell = table[r["task_type"]][r["kind"]]
        cell[0] += r["correct"]; cell[1] += 1
    print(f"\n{len(recs)} answered, {len(cut)} cut off before a final answer (excluded): {cut}")
    print(f"{'task_type':<22}{'seen':>8}{'fresh':>8}")
    tot = {"seen": [0, 0], "fresh": [0, 0]}
    for t in FRESH_OK:
        if t in table:
            s, fr = table[t]["seen"], table[t]["fresh"]
            print(f"{t:<22}{s[0]:>4}/{s[1]:<3}{fr[0]:>4}/{fr[1]:<3}")
            for k in tot:
                tot[k][0] += table[t][k][0]; tot[k][1] += table[t][k][1]
    for k, (c, n) in tot.items():
        print(f"{k:>6}: {c}/{n} correct ({100 * c / max(n, 1):.0f}%)")


if __name__ == "__main__":
    main()
    sys.stdout.flush()
    os._exit(0)
