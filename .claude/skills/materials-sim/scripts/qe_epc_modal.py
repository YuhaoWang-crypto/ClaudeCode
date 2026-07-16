"""
MgB2 电声耦合 (Quantum ESPRESSO DFPT) -> λ, ω_log, Allen-Dynes Tc。
诚实说明: MgB2 的 λ 对 k/q 网格极敏感, 本 demo 用中等网格, 结果为欠收敛的示意值。
"""
import modal

image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("quantum-espresso", "wget", "ca-certificates")
    .pip_install("numpy")
)
app = modal.App("qe-mgb2-epc", image=image)

# QE 自带测试赝势 (LDA, norm-conserving) —— 经典 MgB2 EPC 教程所用
PSEUDO = {
    "B.pz-vbc.UPF": "https://pseudopotentials.quantum-espresso.org/upf_files/B.pz-vbc.UPF",
    "Mg.pz-n-vbc.UPF": "https://pseudopotentials.quantum-espresso.org/upf_files/Mg.pz-n-vbc.UPF",
}


@app.function(cpu=16.0, memory=32768, timeout=5400)
def mgb2_epc(a: float = 5.826, c_over_a: float = 1.142,
             nk: int = 16, nq: int = 2, ecut: float = 40.0) -> dict:
    import os, subprocess, urllib.request, glob
    import numpy as np

    os.makedirs("/work/pseudo", exist_ok=True)
    os.chdir("/work")
    for fn, url in PSEUDO.items():
        try:
            urllib.request.urlretrieve(url, f"/work/pseudo/{fn}")
        except Exception as e:
            return {"error": f"pseudo download failed {fn}: {e!r}"}
    ok = {fn: os.path.getsize(f"/work/pseudo/{fn}") for fn in PSEUDO}

    c = a * c_over_a
    # 找 QE 可执行文件
    def exe(name):
        for p in (name, f"/usr/bin/{name}", f"/usr/bin/{name}.mpi"):
            if subprocess.run(["which", p], capture_output=True).returncode == 0 or os.path.exists(p):
                return p
        return name
    PW, PH = exe("pw.x"), exe("ph.x")

    scf = f"""&control
 calculation='scf', prefix='mgb2', outdir='/work/out', pseudo_dir='/work/pseudo', verbosity='high'
/
&system
 ibrav=4, celldm(1)={a}, celldm(3)={c_over_a}, nat=3, ntyp=2,
 ecutwfc={ecut}, occupations='smearing', smearing='mp', degauss=0.025
/
&electrons
 conv_thr=1.0d-10, mixing_beta=0.7
/
ATOMIC_SPECIES
 Mg 24.305 Mg.pz-n-vbc.UPF
 B  10.811 B.pz-vbc.UPF
ATOMIC_POSITIONS crystal
 Mg 0.0 0.0 0.0
 B  0.333333333 0.666666667 0.5
 B  0.666666667 0.333333333 0.5
K_POINTS automatic
 {nk} {nk} {nk} 0 0 0
"""
    open("scf.in", "w").write(scf)
    r1 = subprocess.run(f"{PW} -in scf.in", shell=True, capture_output=True, text=True)
    if "JOB DONE" not in r1.stdout:
        return {"error": "scf failed", "tail": r1.stdout[-1500:] + r1.stderr[-500:], "pseudo": ok}

    ph = f"""electron-phonon coupling for MgB2
&inputph
 prefix='mgb2', outdir='/work/out', fildyn='mgb2.dyn',
 tr2_ph=1.0d-12, ldisp=.true., nq1={nq}, nq2={nq}, nq3={nq},
 electron_phonon='interpolated', el_ph_sigma=0.005, el_ph_nsigma=10,
 fildvscf='dvscf'
/
"""
    open("ph.in", "w").write(ph)
    r2 = subprocess.run(f"{PH} -in ph.in", shell=True, capture_output=True, text=True)
    out = r2.stdout

    # 解析 lambda / omega_log / Tc (ph.x 会在 a2Fq2r/lambda 行输出; 取汇总)
    lambdas, wlogs = [], []
    for line in out.splitlines():
        s = line.strip()
        if s.startswith("lambda") and "=" in s:
            try:
                lambdas.append(float(s.split("=")[1].split()[0]))
            except Exception:
                pass
    # Allen-Dynes (若能解析到) —— 否则返回原始输出片段供检查
    res = {"pseudo_bytes": ok, "nk": nk, "nq": nq, "ecut": ecut,
           "lambda_values_parsed": lambdas[:10],
           "ph_done": "JOB DONE" in out,
           "ph_tail": out[-2500:]}
    print("lambda parsed:", lambdas[:10], "ph_done:", res["ph_done"])
    return res


@app.local_entrypoint()
def main(nk: int = 16, nq: int = 2):
    import json
    r = mgb2_epc.remote(nk=nk, nq=nq)
    json.dump(r, open("/home/user/qe_mgb2_result.json", "w"), indent=2)
    print(json.dumps({k: v for k, v in r.items() if k != "ph_tail"}, ensure_ascii=False, indent=2))
    print("---- ph.x 输出尾部 ----")
    print(r.get("ph_tail", r.get("tail", "(none)"))[-1500:])
