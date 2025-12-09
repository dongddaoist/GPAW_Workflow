#!/usr/bin/env python3
import argparse
import csv
import json
import sys
from pathlib import Path
import subprocess

from ase.io import read, write
from ase.units import Bohr
from ase.data import atomic_numbers
from ase.optimize import BFGS
from gpaw import GPAW, FermiDirac, mpi as gpaw_mpi
from gpaw import PoissonSolver

# MPI rank / communicator
rank = gpaw_mpi.rank
world = gpaw_mpi.world


# ---------- CSV / metadata helpers ----------

def load_row_for_mid(csv_path: Path, mid: str):
    """Return the row (dict) from csv_path matching material_id == mid."""
    with open(csv_path, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get("material_id", "").strip() == mid:
                return row
    return None


# ---------- Bader ACF parser ----------

def parse_acf(acf_path: Path):
    """
    Parse Henkelman-group ACF.dat file.

    Returns:
        rows: list of dicts with per-atom info
        total_e: total electrons in all basins
        vacuum_e: vacuum charge
    """
    rows = []
    total_e, vacuum_e = None, None
    with open(acf_path, "r") as f:
        lines = [l.rstrip("\n") for l in f]

    start = None
    for i, l in enumerate(lines):
        if l.strip().startswith("#") and "X" in l and "CHARGE" in l.upper():
            start = i + 1
            break
    if start is None:
        raise RuntimeError(f"Cannot find ACF header in {acf_path}")

    # Per-atom rows
    for l in lines[start:]:
        if not l.strip() or set(l.strip()) == set("-"):
            break
        parts = l.split()
        if len(parts) < 7:
            continue
        rows.append({
            "idx": int(parts[0]),
            "x": float(parts[1]),
            "y": float(parts[2]),
            "z": float(parts[3]),
            "electrons": float(parts[4]),
            "min_dist": float(parts[5]),
            "atomic_vol": float(parts[6]),
        })

    # Footer for total electrons and vacuum charge
    for l in lines[::-1]:
        up = l.upper()
        if "NUMBER OF ELECTRONS" in up and "VACUUM" not in up:
            try:
                total_e = float(l.split()[-1])
            except Exception:
                pass
        if "VACUUM CHARGE" in up:
            try:
                vacuum_e = float(l.split()[-1])
            except Exception:
                pass

    return rows, total_e, vacuum_e


# ---------- small utility to run external commands ----------

def run_cmd(cmd, cwd=None):
    p = subprocess.Popen(
        cmd,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    out, err = p.communicate()
    return p.returncode, out, err


# ---------- GPAW helpers ----------

def build_gpaw_calculator(mode="fd", h=0.18, ecut=600.0, kpts=None,
                          spinpol=True, txt_path=None):
    """
    Build a GPAW calculator with either finite-difference (fd) or PW mode.
    Uses Fermi-Dirac smearing and a Poisson solver with default settings.
    """
    if mode == "fd":
        calc = GPAW(
            mode="fd",
            h=h,
            xc="PBE",
            kpts=kpts,
            spinpol=spinpol,
            occupations=FermiDirac(0.1),
            poissonsolver=PoissonSolver(),  # no eps to avoid future warnings
            txt=str(txt_path) if txt_path else None,
            convergence={"energy": 1e-6, "density": 1e-5},
        )
    else:
        from gpaw import PW
        calc = GPAW(
            mode=PW(ecut),
            xc="PBE",
            kpts=kpts,
            spinpol=spinpol,
            occupations=FermiDirac(0.1),
            txt=str(txt_path) if txt_path else None,
            convergence={"energy": 1e-6, "density": 1e-5},
        )
    return calc


def guess_kmesh(atoms, target=0.25):
    """
    Guess a Gamma-centered Monkhorst-Pack kmesh from reciprocal cell lengths.

    target ~ 0.25 1/Å is a typical choice for ~0.2-0.3 Å^-1 spacing.
    """
    import numpy as np
    a, b, c = atoms.get_reciprocal_cell()
    lens = np.array([np.linalg.norm(v) for v in (a, b, c)])
    kmesh = np.maximum(1, (lens / target).astype(int))
    return tuple(int(x) for x in kmesh)


def relax_if_needed(atoms, calc, outdir: Path, do_relax: bool):
    """
    If do_relax == True: run a BFGS relaxation (fixed cell).
    Else: just run one SCF to converge the density at the current geometry.
    """
    atoms.calc = calc
    if not do_relax:
        E = atoms.get_potential_energy()
        return atoms, E

    traj_file = outdir / "relax.traj"
    log_file = outdir / "relax.log"
    dyn = BFGS(atoms, trajectory=str(traj_file), logfile=str(log_file))
    dyn.run(fmax=0.02)
    E = atoms.get_potential_energy()
    write(outdir / "POSCAR_relaxed", atoms, format="vasp")
    return atoms, E


# ---------- main workflow ----------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mid", required=True, help="material_id, e.g. mp-12345")
    ap.add_argument("--base-dir", default=".",
                    help="project root (contains poscar-input/)")
    ap.add_argument("--mode", choices=["fd", "pw"], default="fd",
                    help="GPAW mode: finite-difference (fd) or plane-wave (pw)")
    ap.add_argument("--h", type=float, default=0.18,
                    help="real-space grid spacing for fd mode")
    ap.add_argument("--ecut", type=float, default=600.0,
                    help="plane-wave energy cutoff (eV) for pw mode")
    ap.add_argument("--no_spin", action="store_true",
                    help="disable spin polarization")
    ap.add_argument("--gridref", type=int, default=4,
                    help="grid refinement for all-electron density")
    ap.add_argument("--bader_bin", default="bader",
                    help="path to Henkelman-group bader executable")
    args = ap.parse_args()

    base = Path(args.base_dir).resolve()
    poscar_input_dir = base / "poscar-input"
    csv_path = base / "poscar-input-list.csv"

    # Find metadata row for this material
    row = load_row_for_mid(csv_path, args.mid)
    if row is None:
        print(f"[ERR] {args.mid} not found in {csv_path}", file=sys.stderr)
        sys.exit(1)

    source_type = row.get("source_type", "").strip()
    if source_type not in ("optimized", "non_optimized"):
        print(f"[WARN] source_type for {args.mid} is '{source_type}', treating as non_optimized")
        source_type = "non_optimized"

    mp_dir = poscar_input_dir / args.mid
    if not mp_dir.is_dir():
        print(f"[ERR] directory {mp_dir} not found", file=sys.stderr)
        sys.exit(1)

    poscar_path = mp_dir / "POSCAR"
    if not poscar_path.exists():
        print(f"[ERR] POSCAR not found for {args.mid} at {poscar_path}", file=sys.stderr)
        sys.exit(1)

    outdir = mp_dir / "gpaw"
    outdir.mkdir(parents=True, exist_ok=True)

    # Mark job as running (only once)
    outdir = mp_dir / "gpaw"
    outdir.mkdir(parents=True, exist_ok=True)

    # Mark job as running (only once)
    if rank == 0:
        (mp_dir / "RUNNING").write_text("running\n")
        # Also update job_status in the CSV
        try:
            update_script = base / "update_job_status.py"
            subprocess.run(
                [
                "python3",
                str(update_script),
                "--csv",
                str(csv_path),
                "--mid",
                args.mid,
                "--status",
                "running",
                ],
                check=True,
            )
        except Exception as e:
            print(f"[WARN] Failed to set job_status=running for {args.mid}: {e}", file=sys.stderr)


    # 1. Read structure
    atoms = read(poscar_path)

    # 2. Build kmesh and calculator
    kmesh = guess_kmesh(atoms, target=0.25)
    spinpol = not args.no_spin
    calc = build_gpaw_calculator(
        mode=args.mode,
        h=args.h,
        ecut=args.ecut,
        kpts=kmesh,
        spinpol=spinpol,
        txt_path=outdir / "gpaw.log",
    )

    # 3. Relaxation if needed
    do_relax = (source_type == "non_optimized")
    atoms, E = relax_if_needed(atoms, calc, outdir, do_relax)

    # 4. Save final calc and structure
    calc.write(outdir / "calc.gpw", mode="all")
    write(outdir / "POSCAR_final", atoms, format="vasp")

    # 5. Densities: all-electron (for Bader), and optionally valence
    # All-electron density (Bohr^-3 -> Å^-3)
    rho_ae = calc.get_all_electron_density(gridrefinement=args.gridref)
    rho_ae_au = rho_ae * Bohr**3
    ae_cube = outdir / "density_ae.cube"
    write(ae_cube, atoms, data=rho_ae_au)

    # Valence/pseudo density (if available)
    try:
        rho_pseudo = calc.get_pseudo_density()
        rho_pseudo_au = rho_pseudo * Bohr**3
        val_cube = outdir / "density_valence.cube"
        write(val_cube, atoms, data=rho_pseudo_au)
    except Exception:
        val_cube = None

    # Ensure cube file is written before running Bader
    world.barrier()

    # 6. Run Bader on all-electron density (rank 0 only)
    if rank == 0:
        code, out, err = run_cmd(
            [args.bader_bin, "-p", "all_atom", "-p", "atom_index", str(ae_cube)],
            cwd=outdir,
        )
        with open(outdir / "bader_stdout.log", "w") as f:
            f.write(out)
            f.write("\n--- STDERR ---\n")
            f.write(err)

        if code != 0:
            print(f"[ERR] Bader failed for {args.mid}", file=sys.stderr)
            print("[ERR] Bader stdout/stderr:", file=sys.stderr)
            print(out, file=sys.stderr)
            print(err, file=sys.stderr)
            sys.exit(2)

    world.barrier()  # wait until Bader finishes

    # 7. Parse ACF.dat and write summary / CSV (rank 0 only)
    if rank == 0:
        acf = outdir / "ACF.dat"
        if not acf.exists():
            print(f"[ERR] ACF.dat not found after Bader for {args.mid}", file=sys.stderr)
            sys.exit(2)

        acf_rows, total_e, vacuum_e = parse_acf(acf)
        charges = []
        for r in acf_rows:
            idx0 = r["idx"] - 1
            Z = atomic_numbers[atoms[idx0].symbol]
            q = Z - r["electrons"]
            charges.append(q)

        # JSON summary
        summary = {
            "material_id": args.mid,
            "source_type": source_type,
            "energy_eV": float(E),
            "kmesh": kmesh,
            "spinpol": spinpol,
            "mode": args.mode,
            "h": args.h if args.mode == "fd" else None,
            "ecut_eV": args.ecut if args.mode == "pw" else None,
            "gridref": args.gridref,
            "total_electrons_in_partitions": total_e,
            "vacuum_electrons": vacuum_e,
            "atoms": [
                {
                    "index": r["idx"],
                    "symbol": atoms[r["idx"] - 1].symbol,
                    "Z": int(atomic_numbers[atoms[r["idx"] - 1].symbol]),
                    "x": r["x"],
                    "y": r["y"],
                    "z": r["z"],
                    "electrons": r["electrons"],
                    "bader_charge": charges[i],
                    "atomic_volume": r["atomic_vol"],
                    "min_dist": r["min_dist"],
                }
                for i, r in enumerate(acf_rows)
            ],
        }
        with open(outdir / "bader_summary.json", "w") as f:
            json.dump(summary, f, indent=2)

        # CSV with per-atom Bader charges
        with open(outdir / "bader_charges.csv", "w") as f:
            f.write("index,symbol,Z,x,y,z,electrons,bader_charge,volume,min_dist\n")
            for i, r in enumerate(acf_rows):
                a = atoms[r["idx"] - 1]
                f.write(
                    f"{r['idx']},{a.symbol},{atomic_numbers[a.symbol]},"
                    f"{r['x']},{r['y']},{r['z']},"
                    f"{r['electrons']},{charges[i]},"
                    f"{r['atomic_vol']},{r['min_dist']}\n"
                )

        # 8. Mark completion
        # 8. Mark completion (rank 0 only)
        if rank == 0:
            (mp_dir / "DONE").write_text("completed\n")
            try:
                (mp_dir / "RUNNING").unlink()
            except FileNotFoundError:
                pass

            # Update job_status in the CSV
            try:
                update_script = base / "update_job_status.py"
                subprocess.run(
                    [
                    "python3",
                    str(update_script),
                    "--csv",
                    str(csv_path),
                    "--mid",
                    args.mid,
                    "--status",
                    "completed",
                    ],
                    check=True,
                )
            except Exception as e:
                print(f"[WARN] Failed to set job_status=completed for {args.mid}: {e}", file=sys.stderr)

            print(f"[OK] {args.mid} finished.")


if __name__ == "__main__":
    main()

