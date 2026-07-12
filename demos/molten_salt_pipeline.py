"""
Molten-salt ML-potential pipeline on Modal.

  Stage 1  CP2K AIMD (CPU)      -> DFT trajectory (energies + forces + coords)
  Stage 2  convert (local)      -> DeepMD-kit training set (type.raw + set.000/*.npy)
  Stage 3  DeepMD training (GPU) -> trained potential (model + learning curve)

Demo system: molten NaCl. Uses a small supercell / few MD steps so the whole
loop finishes in minutes. Scaling to 郭硕's spec (50-100 atoms, 6 T points,
~10 ps equilibration) is just parameters + more steps -- see notes at the end.
"""
import io, json, os, re, subprocess, sys
import modal
# NOTE: numpy is imported lazily inside functions, NOT at module top level.
# Modal imports this whole module inside every container to hydrate functions,
# and the CP2K image's Python cannot import numpy (add_python quirk) -- but the
# CP2K function never needs numpy anyway.

app = modal.App("molten-salt-mlp-v2")

# ---------- images ----------
cp2k_image = modal.Image.from_registry("cp2k/cp2k:2024.1", add_python="3.11")

deepmd_image = (
    modal.Image.from_registry(
        "pytorch/pytorch:2.3.1-cuda12.1-cudnn8-runtime"
    )
    .pip_install("deepmd-kit", "numpy<2")
)

# ---------- unit conversions ----------
HARTREE_EV = 27.211386245988
BOHR_A = 0.529177210903
FORCE_AU_EVA = HARTREE_EV / BOHR_A      # Hartree/Bohr -> eV/Angstrom


# ======================================================================
# geometry (runs locally; sandbox has numpy)
# ======================================================================
def build_nacl(nx=2, ny=2, nz=2, a0=6.30):
    """Rocksalt NaCl supercell. a0 expanded toward molten density (~1.56 g/cc)."""
    import numpy as np
    basis = [  # (element, fractional pos in conventional cubic cell)
        ("Na", (0.0, 0.0, 0.0)), ("Na", (0.0, 0.5, 0.5)),
        ("Na", (0.5, 0.0, 0.5)), ("Na", (0.5, 0.5, 0.0)),
        ("Cl", (0.5, 0.5, 0.5)), ("Cl", (0.5, 0.0, 0.0)),
        ("Cl", (0.0, 0.5, 0.0)), ("Cl", (0.0, 0.0, 0.5)),
    ]
    syms, xyz = [], []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for el, (fx, fy, fz) in basis:
                    syms.append(el)
                    xyz.append([(i + fx) * a0, (j + fy) * a0, (k + fz) * a0])
    cell = np.array([nx * a0, ny * a0, nz * a0])
    return syms, np.array(xyz), cell


CP2K_TEMPLATE = """\
&GLOBAL
  PROJECT salt
  RUN_TYPE MD
  PRINT_LEVEL LOW
&END GLOBAL
&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    POTENTIAL_FILE_NAME GTH_POTENTIALS
    &MGRID
      CUTOFF 300
      REL_CUTOFF 50
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-10
    &END QS
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF 1.0E-6
      MAX_SCF 30
      &OT
        MINIMIZER DIIS
        PRECONDITIONER FULL_SINGLE_INVERSE
      &END OT
      &OUTER_SCF
        EPS_SCF 1.0E-6
        MAX_SCF 10
      &END OUTER_SCF
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC {Lx:.6f} {Ly:.6f} {Lz:.6f}
    &END CELL
    &COORD
{coords}
    &END COORD
    &KIND Na
      BASIS_SET DZVP-MOLOPT-SR-GTH
      POTENTIAL GTH-PBE-q9
    &END KIND
    &KIND Cl
      BASIS_SET DZVP-MOLOPT-SR-GTH
      POTENTIAL GTH-PBE-q7
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
&MOTION
  &MD
    ENSEMBLE NVT
    STEPS {steps}
    TIMESTEP 1.0
    TEMPERATURE {temp}
    &THERMOSTAT
      TYPE NOSE
      &NOSE
        TIMECON 50
      &END NOSE
    &END THERMOSTAT
  &END MD
  &PRINT
    &TRAJECTORY
      &EACH
        MD 1
      &END EACH
    &END TRAJECTORY
    &FORCES
      &EACH
        MD 1
      &END EACH
    &END FORCES
  &END PRINT
&END MOTION
"""


def make_cp2k_input(syms, xyz, cell, steps, temp):
    coords = "\n".join(
        f"      {el}  {x:.6f}  {y:.6f}  {z:.6f}"
        for el, (x, y, z) in zip(syms, xyz)
    )
    return CP2K_TEMPLATE.format(
        Lx=cell[0], Ly=cell[1], Lz=cell[2], coords=coords, steps=steps, temp=temp
    )


# ======================================================================
# STAGE 1 : CP2K AIMD  (CPU)
# ======================================================================
@app.function(image=cp2k_image, cpu=8.0, timeout=3600)
def run_cp2k(inp: str):
    import glob, shutil
    workdir = "/root/run"
    os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir)
    with open("salt.inp", "w") as f:
        f.write(inp)

    # locate CP2K data dir
    data_dir = os.environ.get("CP2K_DATA_DIR", "")
    if not (data_dir and os.path.exists(os.path.join(data_dir, "BASIS_MOLOPT"))):
        for c in ["/opt/cp2k/data", "/opt/cp2k-*/data"]:
            hits = glob.glob(os.path.join(c, "BASIS_MOLOPT"))
            if hits:
                data_dir = os.path.dirname(hits[0]); break
        else:
            hit = subprocess.run(["bash", "-lc",
                                  "find / -name BASIS_MOLOPT 2>/dev/null | head -1"],
                                 capture_output=True, text=True).stdout.strip()
            data_dir = os.path.dirname(hit)
    os.environ["CP2K_DATA_DIR"] = data_dir

    binary = shutil.which("cp2k.psmp") or shutil.which("cp2k.ssmp") or "cp2k.psmp"
    env = dict(os.environ, OMP_NUM_THREADS="8")
    p = subprocess.run([binary, "-i", "salt.inp", "-o", "cp2k.log"],
                       env=env, capture_output=True, text=True)

    def read(path):
        return open(path).read() if os.path.exists(path) else ""

    log = read("cp2k.log")
    return {
        "returncode": p.returncode,
        "data_dir": data_dir,
        "binary": binary,
        "log_tail": log[-4000:],
        "stderr_tail": p.stderr[-2000:],
        "pos_xyz": read("salt-pos-1.xyz"),
        "frc_xyz": read("salt-frc-1.xyz"),
    }


# ======================================================================
# STAGE 2 : parse CP2K -> arrays  (local)
# ======================================================================
def parse_xyz_frames(text):
    """Return list of (comment, np.array[N,3]) from a concatenated xyz."""
    import numpy as np
    lines = text.splitlines()
    frames, i = [], 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1; continue
        n = int(lines[i].split()[0])
        comment = lines[i + 1]
        arr = np.array([[float(v) for v in lines[i + 2 + k].split()[1:4]]
                        for k in range(n)])
        frames.append((comment, arr))
        i += 2 + n
    return frames


def cp2k_to_arrays(pos_xyz, frc_xyz, syms, cell):
    import numpy as np
    pos = parse_xyz_frames(pos_xyz)
    frc = parse_xyz_frames(frc_xyz)
    n = min(len(pos), len(frc))
    pos, frc = pos[:n], frc[:n]

    coord = np.array([p[1] for p in pos]).reshape(n, -1)          # Angstrom
    force = np.array([f[1] for f in frc]).reshape(n, -1) * FORCE_AU_EVA
    energy = np.array([float(re.search(r"E\s*=\s*(-?\d+\.\d+)", p[0]).group(1))
                       for p in pos]) * HARTREE_EV
    box = np.tile(np.diag(cell).reshape(1, 9), (n, 1))            # Angstrom

    type_map = ["Na", "Cl"]
    types = np.array([type_map.index(s) for s in syms])
    return dict(type_map=type_map, types=types,
                coord=coord, force=force, energy=energy, box=box)


# ======================================================================
# STAGE 3 : DeepMD training  (GPU)
# ======================================================================
DP_INPUT = {
    "model": {
        "type_map": ["Na", "Cl"],
        "descriptor": {"type": "se_e2_a", "sel": [80, 80],
                       "rcut_smth": 0.5, "rcut": 6.0,
                       "neuron": [25, 50, 100], "resnet_dt": False,
                       "axis_neuron": 12, "seed": 1},
        "fitting_net": {"neuron": [120, 120, 120], "resnet_dt": True, "seed": 1},
    },
    "learning_rate": {"type": "exp", "start_lr": 1e-3, "stop_lr": 3.5e-4,
                      "decay_steps": 500},
    "loss": {"type": "ener", "start_pref_e": 0.02, "limit_pref_e": 1,
             "start_pref_f": 1000, "limit_pref_f": 1},
    "training": {
        "training_data": {"systems": ["./data"], "batch_size": "auto"},
        "validation_data": {"systems": ["./data"], "batch_size": "auto"},
        "numb_steps": 2000, "seed": 10,
        "disp_file": "lcurve.out", "disp_freq": 100, "save_freq": 1000,
    },
}


@app.function(image=deepmd_image, gpu="T4", timeout=3600)
def train_deepmd(arrays: dict, dp_input: dict):
    import numpy as np
    work = "/root/train"; data = os.path.join(work, "data", "set.000")
    os.makedirs(data, exist_ok=True)
    os.chdir(work)

    d = arrays
    np.savetxt("data/type.raw", d["types"], fmt="%d")
    with open("data/type_map.raw", "w") as f:
        f.write("\n".join(d["type_map"]) + "\n")
    np.save(os.path.join(data, "coord.npy"),  d["coord"].astype(np.float64))
    np.save(os.path.join(data, "box.npy"),    d["box"].astype(np.float64))
    np.save(os.path.join(data, "energy.npy"), d["energy"].astype(np.float64))
    np.save(os.path.join(data, "force.npy"),  d["force"].astype(np.float64))

    with open("input.json", "w") as f:
        json.dump(dp_input, f, indent=2)

    gpu = subprocess.run(["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"],
                         capture_output=True, text=True).stdout.strip()

    # torch backend (base image ships torch+CUDA)
    train = subprocess.run(["dp", "--pt", "train", "input.json"],
                           capture_output=True, text=True)
    freeze = subprocess.run(["dp", "--pt", "freeze", "-o", "model.pth"],
                            capture_output=True, text=True)

    def read(p, b=False):
        if not os.path.exists(p):
            return b"" if b else ""
        return open(p, "rb").read() if b else open(p).read()

    return {
        "gpu": gpu,
        "train_rc": train.returncode,
        "freeze_rc": freeze.returncode,
        "train_tail": train.stdout[-3000:] + "\n--STDERR--\n" + train.stderr[-3000:],
        "lcurve": read("lcurve.out"),
        "model_bytes": read("model.pth", b=True),
        "nframes": int(d["energy"].shape[0]),
    }


# ======================================================================
# orchestration
# ======================================================================
@app.local_entrypoint()
def main(stage: str = "all", steps: int = 30, temp: float = 1200.0,
         nx: int = 2, ny: int = 2, nz: int = 2, dp_steps: int = 2000):
    import numpy as np
    outdir = "/home/user/ClaudeCode/demos/msalt_out"
    os.makedirs(outdir, exist_ok=True)

    syms, xyz, cell = build_nacl(nx, ny, nz)
    print(f"[geom] NaCl {len(syms)} atoms, cell {cell.round(2)} A")

    if stage in ("aimd", "all"):
        inp = make_cp2k_input(syms, xyz, cell, steps, temp)
        open(os.path.join(outdir, "salt.inp"), "w").write(inp)
        print(f"[cp2k] launching AIMD: {steps} steps @ {temp} K ...")
        res = run_cp2k.remote(inp)
        print(f"[cp2k] rc={res['returncode']} binary={res['binary']} data_dir={res['data_dir']}")
        if res["returncode"] != 0:
            print("[cp2k] LOG TAIL:\n", res["log_tail"][-1500:])
            print("[cp2k] STDERR:\n", res["stderr_tail"])
            return
        open(os.path.join(outdir, "salt-pos-1.xyz"), "w").write(res["pos_xyz"])
        open(os.path.join(outdir, "salt-frc-1.xyz"), "w").write(res["frc_xyz"])
        nfr = res["pos_xyz"].count("i =")
        print(f"[cp2k] done, ~{nfr} frames written")

    if stage in ("convert", "all"):
        pos = open(os.path.join(outdir, "salt-pos-1.xyz")).read()
        frc = open(os.path.join(outdir, "salt-frc-1.xyz")).read()
        arr = cp2k_to_arrays(pos, frc, syms, cell)
        np.savez(os.path.join(outdir, "arrays.npz"), **arr,
                 types=arr["types"], energy=arr["energy"])
        print(f"[conv] {arr['energy'].shape[0]} frames | "
              f"E range {arr['energy'].min():.2f}..{arr['energy'].max():.2f} eV | "
              f"|F|max {np.abs(arr['force']).max():.2f} eV/A")

    if stage in ("train", "all"):
        data = np.load(os.path.join(outdir, "arrays.npz"), allow_pickle=True)
        arrays = {k: data[k] for k in data.files}
        arrays["type_map"] = list(arrays["type_map"])
        dp_input = json.loads(json.dumps(DP_INPUT))
        dp_input["training"]["numb_steps"] = dp_steps
        print(f"[deepmd] training {dp_steps} steps on GPU ...")
        res = train_deepmd.remote(arrays, dp_input)
        print(f"[deepmd] gpu={res['gpu']} train_rc={res['train_rc']} "
              f"freeze_rc={res['freeze_rc']} nframes={res['nframes']}")
        if res["lcurve"]:
            open(os.path.join(outdir, "lcurve.out"), "w").write(res["lcurve"])
            print("[deepmd] learning curve (head/tail):")
            lines = res["lcurve"].strip().splitlines()
            for l in lines[:2] + ["..."] + lines[-3:]:
                print("   ", l)
        if res["model_bytes"]:
            open(os.path.join(outdir, "model.pth"), "wb").write(res["model_bytes"])
            print(f"[deepmd] saved model.pth ({len(res['model_bytes'])} bytes)")
        if res["train_rc"] != 0:
            print("[deepmd] TRAIN OUTPUT TAIL:\n", res["train_tail"][-2000:])
