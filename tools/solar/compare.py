"""
Solar system validation: simulation output vs JPL Horizons reference.

Reads tests/sim_output.bin and tests/jpl_reference.csv, reports per-body
max and RMS relative position error in two tables (relative %, absolute km),
and exports tests/comparison_results.json.

Usage:
    python tools/solar/compare.py
    python tools/solar/compare.py --bodies Mercury,Venus,Earth

Run tools/solar/fetch.py first to populate the reference CSV.
"""

import argparse, csv, json, struct
from collections import defaultdict
from pathlib import Path
import numpy as np


def load_sim_binary(path):
    """Read tests/sim_output.bin into per-body dicts of arrays."""
    data = defaultdict(lambda: {"t": [], "pos": [], "vel": []})
    with open(path, "rb") as f:
        num_bodies = struct.unpack("Q", f.read(8))[0]
        names = []
        for _ in range(num_bodies):
            raw = f.read(32)
            names.append(raw.split(b'\x00')[0].decode())

        doubles_per_frame = 2 + num_bodies * 6
        frame_size = doubles_per_frame * 8

        while True:
            raw = f.read(frame_size)
            if len(raw) < frame_size:
                break

            frame = np.frombuffer(raw, dtype=np.float64)
            time_s = frame[1]

            states = frame[2:].reshape(num_bodies, 6)
            for i, name in enumerate(names):
                data[name]["t"].append(time_s)
                data[name]["pos"].append(states[i, 0:3] / 1e3)  # m -> km
                data[name]["vel"].append(states[i, 3:6] / 1e3)

    return {n: {k: np.array(v) for k, v in d.items()} for n, d in data.items()}


def load_ref(path):
    """Read JPL reference CSV into per-body dicts of arrays."""
    data = defaultdict(lambda: {"jd": [], "pos": [], "vel": []})
    with open(path) as f:
        for row in csv.DictReader(f):
            n = row["name"]
            data[n]["jd"].append(float(row["jd"]))
            data[n]["pos"].append([float(row["x_km"]),  float(row["y_km"]),  float(row["z_km"])])
            data[n]["vel"].append([float(row["vx_kms"]), float(row["vy_kms"]), float(row["vz_kms"])])
    return {n: {k: np.array(v) for k, v in d.items()} for n, d in data.items()}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sim",    default="tests/sim_output.bin")
    p.add_argument("--ref",    default="tests/jpl_reference.csv")
    p.add_argument("--bodies", default=None,
                   help="comma-separated body names (default: all shared)")
    args = p.parse_args()

    ref = load_ref(args.ref)
    sim = load_sim_binary(args.sim)
    shared = sorted(set(sim) & set(ref))
    if args.bodies:
        shared = [b for b in args.bodies.split(",") if b in shared]
    if not shared:
        print(f"No shared bodies.\n  Sim: {sorted(sim)}\n  Ref: {sorted(ref)}")
        return

    print(f"Comparing {len(shared)} bodies\n")

    results = {}
    for name in shared:
        n = min(len(sim[name]["pos"]), len(ref[name]["pos"]))
        if n < 2: continue
        sp = sim[name]["pos"][:n]
        rp = ref[name]["pos"][:n]

        abs_err = np.linalg.norm(sp - rp, axis=1)
        r_ref = np.linalg.norm(rp, axis=1)
        rel = np.where(r_ref > 0, abs_err / r_ref, 0)

        results[name] = {
            "max_rel_pct":  rel[1:].max() * 100,
            "rms_rel_pct":  np.sqrt(np.mean(rel[1:]**2)) * 100,
            "max_abs_km":   abs_err[1:].max(),
            "rms_abs_km":   np.sqrt(np.mean(abs_err[1:]**2)),
            "mean_abs_km":  abs_err[1:].mean(),
        }

    print(f"{'Body':<14} {'Max Rel (%)':<16} {'RMS Rel (%)':<16}")
    print("=" * 48)
    for name in sorted(results, key=lambda n: results[n]["max_rel_pct"], reverse=True):
        r = results[name]
        print(f"{name:<14} {r['max_rel_pct']:<16.6f} {r['rms_rel_pct']:<16.6f}")

    rel_vals = [r["max_rel_pct"] for r in results.values()]
    rel_vals_no_sun = [r["max_rel_pct"] for n, r in results.items() if n != "Sun"]

    print("=" * 48)
    print(f"{'Worst:':<14} {max(rel_vals):.6f}%")
    print(f"{'Mean:':<14} {sum(rel_vals) / len(rel_vals):.6f}%")
    if len(rel_vals_no_sun) < len(rel_vals):
        print(f"{'Mean (ex Sun):':<14} {sum(rel_vals_no_sun) / len(rel_vals_no_sun):.6f}%")

    print(f"\n{'Body':<14} {'Max Abs (km)':<18} {'RMS Abs (km)':<18} {'Mean Abs (km)':<18}")
    print("=" * 70)
    for name in sorted(results, key=lambda n: results[n]["max_abs_km"], reverse=True):
        r = results[name]
        print(f"{name:<14} {r['max_abs_km']:<18.2f} {r['rms_abs_km']:<18.2f} {r['mean_abs_km']:<18.2f}")

    abs_vals = [r["max_abs_km"] for r in results.values()]
    abs_vals_no_sun = [r["max_abs_km"] for n, r in results.items() if n != "Sun"]

    print("=" * 70)
    print(f"{'Worst:':<14} {max(abs_vals):.2f} km")
    print(f"{'Mean:':<14} {sum(abs_vals) / len(abs_vals):.2f} km")
    if len(abs_vals_no_sun) < len(abs_vals):
        print(f"{'Mean (ex Sun):':<14} {sum(abs_vals_no_sun) / len(abs_vals_no_sun):.2f} km")

    summary = {
        "bodies": {n: r for n, r in results.items()},
        "summary": {
            "num_bodies": len(results),
            "mean_max_rel_pct": sum(rel_vals) / len(rel_vals),
            "mean_max_rel_pct_ex_sun": sum(rel_vals_no_sun) / len(rel_vals_no_sun) if rel_vals_no_sun else None,
            "worst_max_rel_pct": max(rel_vals),
            "worst_max_abs_km":  max(abs_vals),
            "mean_max_abs_km":   sum(abs_vals) / len(abs_vals),
        },
    }
    out_path = Path("tests/comparison_results.json")
    out_path.write_text(json.dumps(summary, indent=2, default=float))
    print(f"\n-> {out_path}")


if __name__ == "__main__":
    main()
