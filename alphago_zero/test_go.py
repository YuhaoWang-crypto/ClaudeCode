"""Rule + symmetry tests.  Run:  python -m alphago_zero.test_go"""
import numpy as np

from .go import GoState, BLACK, WHITE, sym_policy, inv_sym_policy, sym_planes


def setup(n, moves):
    s = GoState.new(n, komi=0.5)
    for m in moves:
        s = s.play(s.str_to_move(m))
    return s


def test_capture():
    # black surrounds the white stone on B2 of a 5x5 board
    s = setup(5, ["B3", "B2", "A2", "E5", "C2", "E4", "B1"])
    assert s.board[s.str_to_move("B2")] == 0, "white stone should be captured"


def test_suicide_illegal():
    s = setup(5, ["B1", "E5", "A2", "E4"])          # black to move
    s = s.play(s.n * s.n)                           # black passes -> white to move
    s = s.play(s.str_to_move("D4"))                 # white elsewhere, black to move
    s = s.play(s.n * s.n)                           # black passes, white to move
    assert not s.legal_mask()[s.str_to_move("A1")], "A1 is suicide for white"


def test_ko():
    # classic ko shape on 5x5
    s = setup(5, ["B4", "C4", "A3", "D3", "B2", "C2", "E5", "B3", "C3"])
    # black C3 captured white B3; white may not retake at B3 immediately
    assert s.board[s.str_to_move("B3")] == 0
    assert s.to_play == WHITE
    assert not s.legal_mask()[s.str_to_move("B3")], "ko recapture must be illegal"
    s = s.play(s.str_to_move("E1"))                 # ko threat
    s = s.play(s.str_to_move("D1"))
    assert s.legal_mask()[s.str_to_move("B3")], "ko retake legal after exchange"


def test_scoring():
    s = setup(5, ["C1", "D1", "C2", "D2", "C3", "D3", "C4", "D4", "C5", "D5"])
    # black owns columns A-C (15), white D-E (10) -> black +5 - komi 0.5
    assert abs(s.score() - 4.5) < 1e-9, s.score()
    assert s.winner() == BLACK


def test_symmetry_roundtrip():
    n = 7
    pi = np.random.default_rng(0).random(n * n + 1)
    for k in range(8):
        assert np.allclose(inv_sym_policy(sym_policy(pi, n, k), n, k), pi)
        # board and policy transforms must agree
        board = np.arange(n * n).reshape(n, n)
        t = sym_planes(board, k).reshape(-1)
        assert np.allclose(sym_policy(pi, n, k)[: n * n], pi[t])


def test_terminal():
    s = GoState.new(5)
    s = s.play(25).play(25)
    assert s.is_terminal()


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok", name)
