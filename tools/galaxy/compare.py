"""
MW-Andromeda separation validation.

Reads tests/sim_output.bin, computes the COM of MW (first N/2 particles) and
M31 (last N/2 particles) at every output frame, and plots their separation
vs time. First close passage predicted around t = 3.9 Gyr; final merger
around t = 5.9 Gyr per van der Marel et al. 2012.

Usage:
    python src/galaxy_compare.py
    python src/galaxy_compare.py --sim tests/sim_output.bin --out tests/galaxy_separation.png
"""

import struct, argparse
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

KPC = 3.086e19
GYR = 3.156e16


def load_positions(path):
    """Load positions only from sim_output.bin.

    Returns (positions, times, N) where positions is (n_frames, N, 3) in meters
    and times is in seconds. Format matches what main.cu writes via Output.cu.
    """
    with open(path, "rb") as f:
        N = struct.unpack("Q", f.read(8))[0]
        f.read(N * 32)  # skip 32-byte name records

        doubles_per_frame = 2 + N * 6
        frame_size = doubles_per_frame * 8

        frames = []
        times = []
        while True:
            raw = f.read(frame_size)
            if len(raw) < frame_size:
                break
            arr = np.frombuffer(raw, dtype=np.float64)
            times.append(arr[1])
            states = arr[2:].reshape(N, 6)
            frames.append(states[:, 0:3])

    return np.array(frames), np.array(times), N


def find_local_minima(y, min_separation=5):
    """Indices of local minima in y, with a minimum index separation between
    successive minima (avoids picking up noise-driven adjacent samples)."""
    out = []
    last = -min_separation
    for i in range(1, len(y) - 1):
        if y[i] < y[i - 1] and y[i] < y[i + 1] and i - last >= min_separation:
            out.append(i)
            last = i
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sim", default="tests/sim_output.bin")
    ap.add_argument("--out", default="tests/galaxy_separation.png")
    args = ap.parse_args()

    print(f"Loading {args.sim}...")
    positions, times, N = load_positions(args.sim)
    half = N // 2
    n_frames = len(positions)

    print(f"  {N} particles, {n_frames} frames")
    print(f"  Time span: {times[0] / GYR:.2f} -> {times[-1] / GYR:.2f} Gyr")

    # Equal-mass particles within each galaxy, so unweighted mean is the COM.
    mw_com  = positions[:, :half, :].mean(axis=1)
    m31_com = positions[:, half:, :].mean(axis=1)

    sep_m   = np.linalg.norm(mw_com - m31_com, axis=1)
    sep_kpc = sep_m / KPC
    t_gyr   = times / GYR

    minima = [m for m in find_local_minima(sep_kpc) if t_gyr[m] > 1.0]

    print(f"\nInitial separation: {sep_kpc[0]:.1f} kpc")
    labels = ["First passage", "Second passage", "Third passage"]
    if minima:
        for i, m in enumerate(minima[:3]):
            print(f"{labels[i]}: t = {t_gyr[m]:.2f} Gyr, min separation = {sep_kpc[m]:.1f} kpc")
    else:
        print("No local minima found (encounter has not occurred yet).")

    print(f"\nvan der Marel et al. 2012 predictions:")
    print(f"  First close passage:  ~3.9 Gyr")
    print(f"  Final merger:         ~5.9 Gyr")

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(t_gyr, sep_kpc, "b-", linewidth=2, label="simulation")
    ax.set_xlabel("Time (Gyr)")
    ax.set_ylabel("MW-M31 separation (kpc)")
    ax.set_title("Milky Way - Andromeda separation over time")
    ax.grid(True, alpha=0.3)
    ax.axvline(3.9, color="red",     linestyle="--", alpha=0.5, label="van der Marel 2012: first passage")
    ax.axvline(5.9, color="darkred", linestyle="--", alpha=0.5, label="van der Marel 2012: merger")

    for m in minima[:3]:
        ax.plot(t_gyr[m], sep_kpc[m], "ko", markersize=6)
        ax.annotate(f"{t_gyr[m]:.2f} Gyr\n{sep_kpc[m]:.0f} kpc",
                    xy=(t_gyr[m], sep_kpc[m]),
                    xytext=(10, 10), textcoords="offset points",
                    fontsize=8, alpha=0.8)

    ax.legend(loc="upper right")
    ax.set_ylim(bottom=0)

    out_path = Path(args.out)
    out_path.parent.mkdir(exist_ok=True, parents=True)
    plt.tight_layout()
    plt.savefig(out_path, dpi=100)
    plt.show()
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
