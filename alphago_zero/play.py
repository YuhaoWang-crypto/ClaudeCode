"""Play against the trained network in the terminal.

    python -m alphago_zero.play                         # you are black
    python -m alphago_zero.play --color white --sims 400
Moves: e.g.  C4, pass, undo, quit
"""
import argparse

import numpy as np

from .go import GoState, BLACK
from .mcts import Node, NetEvaluator, search
from .model import load


def ai_move(state, ev, sims, c_puct, rng):
    root = Node(state)
    search([root], ev, sims, c_puct, None, rng)
    # most visits, ties by mean value; never pass while passing back would lose
    n = state.n * state.n
    mine = state.score() * (1 if state.to_play == BLACK else -1)
    key = root.N + 1e-3 * root.W / np.maximum(root.N, 1)
    if mine < 0 and key[:n].max() > 0:
        key[n] = -np.inf
    a = int(np.argmax(key))
    return a, root


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ckpt", default="alphago_zero/runs/7x7/best.pt")
    ap.add_argument("--color", default="black", choices=["black", "white"])
    ap.add_argument("--sims", type=int, default=200)
    args = ap.parse_args()

    net, ck = load(args.ckpt)
    cfg = ck["cfg"]
    rng = np.random.default_rng()
    ev = NetEvaluator(net, cfg["history"], True, rng)
    human = BLACK if args.color == "black" else -BLACK
    state = GoState.new(cfg["board_size"], cfg["komi"])
    past = []
    print(f"AlphaGo Zero mini | {state.n}x{state.n}, komi {state.komi}, "
          f"generation {ck.get('generation')} | you are {args.color}")
    while not state.is_terminal():
        print("\n" + str(state))
        if state.to_play == human:
            txt = input("your move> ").strip()
            if txt in ("quit", "q"):
                return
            if txt == "undo" and len(past) >= 2:
                state = past[-2]
                past = past[:-2]
                continue
            try:
                mv = state.str_to_move(txt)
                if not state.legal_mask()[mv]:
                    raise ValueError
            except (ValueError, IndexError):
                print("illegal / unreadable move")
                continue
            past.append(state)
            state = state.play(mv)
        else:
            mv, root = ai_move(state, ev, args.sims, cfg["c_puct"], rng)
            v = root.value()
            print(f"AI plays {state.move_to_str(mv)}  "
                  f"(AI win-rate estimate {100 * (v + 1) / 2:.0f}%)")
            past.append(state)
            state = state.play(mv)
    print("\n" + str(state))
    s = state.score()
    print(f"final: {'Black' if s > 0 else 'White'} wins by {abs(s)} "
          f"(Tromp-Taylor, komi {state.komi})")


if __name__ == "__main__":
    main()
