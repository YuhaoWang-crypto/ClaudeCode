"""Quiz the live Carbon-X model on a Hugging Face Space with FineDNA questions.

    python finedna_demo/eval_carbonx_space.py --per-type 3

Carbon-X has no hosted inference endpoint, so this calls a public Gradio Space
that runs it on ZeroGPU. Two Spaces expose a usable API (see SPACES below); both
load HuggingFaceBio/Carbon-X with thinking off and sample at temperature 0.6,
the authors' setting. Set HF_TOKEN to use your own ZeroGPU quota, which is per
account and shared by every Space.

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

def _ask_mbcredits(client, prompt, max_new):
    # /ask_carbon(question, reference_dna, variant_dna, extra_dna, extra_name, context, max_new_tokens)
    return client.predict(prompt, "", "", "", "extra", "", max_new, api_name="/ask_carbon")


def _ask_mrfakename(client, prompt, max_new):
    # /chat(message, history, temperature, max_new_tokens) -> (history, textbox, timing)
    history, _, _ = client.predict(prompt, [], 0.6, max_new, api_name="/chat")
    content = history[-1]["content"]
    if isinstance(content, list):
        content = "".join(part.get("text", "") for part in content if isinstance(part, dict))
    return content


# name -> (caller, max_new_tokens cap). fffiloni/Carbon-X-Playground runs Carbon-X too,
# but exposes no API endpoints, so it can only be used in the browser.
SPACES = {
    "mrfakename/carbon-x": (_ask_mrfakename, 1536),
    "mbcredits/carbon-x-demo": (_ask_mbcredits, 1024),
}
DEFAULT_SPACE = "mrfakename/carbon-x"
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
    reply = reply.replace("`", "")  # mrfakename's Space wraps DNA-like spans in backticks
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
    ap.add_argument("--space", default=DEFAULT_SPACE, choices=list(SPACES))
    ap.add_argument("--max-new-tokens", type=int, help="default: the Space's cap")
    ap.add_argument("--skip-types", nargs="*", default=[],
                    help="task types not to ask, e.g. gc_content (its base-by-base count of a "
                         "~120 bp sequence needs ~1,900 tokens, more than either Space allows)")
    args = ap.parse_args()
    ask, cap = SPACES[args.space]
    max_new = min(args.max_new_tokens or cap, cap)

    done = {}
    if OUT.exists():
        for line in OUT.open():
            d = json.loads(line)
            # A reply that was cut off is asked again only with a bigger token budget.
            if d["correct"] or has_final_answer(d["reply"]) or d.get("max_new_tokens", 512) >= max_new:
                done[d["qid"]] = d
    qs = [q for q in build_questions(args.per_type)
          if q["qid"] not in done and q["task_type"] not in args.skip_types]
    print(f"{len(done)} cached, {len(qs)} to ask on {args.space}", flush=True)

    if qs:
        from gradio_client import Client
        client = Client(args.space, token=os.environ.get("HF_TOKEN"), verbose=False)
        with OUT.open("a") as f:
            for q in qs:
                t0 = time.time()
                try:
                    reply = ask(client, q["prompt"], max_new)
                except Exception as e:  # quota exhausted, Space asleep, ...
                    print(f"  stopped at {q['qid']}: {type(e).__name__}: {str(e)[:200]}", flush=True)
                    break
                ans = extract_answer(reply, q["task_type"])
                rec = {**q, "space": args.space, "max_new_tokens": max_new,
                       "reply": reply, "model_answer": ans, "correct": compare(ans, q["gold"]),
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
    by_space = collections.Counter(r.get("space", "mbcredits/carbon-x-demo") for r in recs)
    print(f"\n{len(recs)} answered ({dict(by_space)}), "
          f"{len(cut)} cut off before a final answer (excluded): {cut}")
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
