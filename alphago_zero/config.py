"""Hyper-parameters.  Comments give the AlphaGo Zero paper value for comparison."""
from dataclasses import dataclass, asdict, field


@dataclass
class Config:
    # game
    board_size: int = 7            # paper: 19
    komi: float = 8.5              # paper: 7.5 (7x7 is a black win by 9 under area rules)
    history: int = 8               # paper: 8 past positions -> 17 input planes

    # network
    filters: int = 48              # paper: 256
    blocks: int = 4                # paper: 19 (20-block) / 39 (40-block)
    value_hidden: int = 48         # paper: 256

    # MCTS
    sims: int = 100                # paper: 1600 per move
    c_puct: float = 1.5            # paper: constant, value unpublished
    dir_alpha: float = 0.25        # paper: 0.03 for 19x19 (~10 / legal moves)
    dir_eps: float = 0.25          # paper: 0.25
    temp_moves: int = 8            # paper: tau=1 for first 30 moves, then tau->0

    # self-play
    games_per_iter: int = 128      # paper: ~25,000 between checkpoints
    workers: int = 4
    no_resign_frac: float = 0.1    # paper: 10% of games never resign
    resign_fp_target: float = 0.05 # paper: keep false positives < 5%
    resign_start_iter: int = 3     # gather calibration data before resigning

    # training
    window_games: int = 4000       # paper: most recent 500,000 games
    batch_size: int = 256          # paper: 2048
    sample_reuse: float = 4.0      # each new position is trained on ~4x (8x symmetries)
    lr: float = 0.02               # paper: 1e-2 -> 1e-3 -> 1e-4 (step annealing)
    lr_milestones: tuple = (0.65, 0.9)  # x0.1 at these fractions of the time budget
    momentum: float = 0.9          # paper: 0.9
    weight_decay: float = 1e-4     # paper: c = 1e-4

    # evaluation / gating
    eval_games: int = 40           # paper: 400
    eval_sims: int = 100           # paper: 1600
    gate_threshold: float = 0.55   # paper: 55%
    gating: bool = True            # AlphaZero (2018) later dropped gating
    random_eval_every: int = 2
    random_eval_games: int = 40

    # run
    iterations: int = 1000
    time_budget_min: float = 150.0
    eval_open_moves: int = 2       # arena: sample first 2 moves for game diversity
    out_dir: str = "alphago_zero/runs/7x7"
    seed: int = 0

    def to_dict(self):
        return asdict(self)
