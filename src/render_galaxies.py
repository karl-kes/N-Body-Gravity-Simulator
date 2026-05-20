"""
MW-Andromeda Rerun visualizer.

Streams the simulation output as a 3D point cloud with a COM-separation
timeline alongside. Two colors: MW blue, M31 red. Assumes the first N/2
particles are MW and the last N/2 are M31 (matches generate_galaxies.py).

Usage:
    python src/render_galaxies.py                    # spawn viewer
    python src/render_galaxies.py --save sim.rrd     # save to file
    python src/render_galaxies.py --radius 0.5       # smaller markers
"""

import struct, argparse
import numpy as np
import rerun as rr
import rerun.blueprint as rrb

KPC = 3.086e19
GYR = 3.156e16
MYR = 3.156e13


def load(path):
    """Load positions only from sim_output.bin. Returns (pos_kpc, times_s, N)."""
    with open(path, "rb") as f:
        N = struct.unpack("Q", f.read(8))[0]
        f.read(N * 32)  # skip names

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

    return np.array(frames) / KPC, np.array(times), N


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sim", default="tests/sim_output.bin")
    ap.add_argument("--save", default=None,
                    help="save to .rrd file instead of spawning viewer")
    ap.add_argument("--radius", type=float, default=1.0,
                    help="point radius in kpc (default: 1.0)")
    args = ap.parse_args()

    print(f"Loading {args.sim}...")
    pos, times, N = load(args.sim)
    half = N // 2
    n_frames = len(pos)

    print(f"  {N} particles, {n_frames} frames")
    print(f"  Time span: {times[0]/GYR:.2f} -> {times[-1]/GYR:.2f} Gyr")

    # Per-frame COM separation (cheap O(N) diagnostic, drives the side panel).
    mw_com  = pos[:, :half, :].mean(axis=1)
    m31_com = pos[:, half:, :].mean(axis=1)
    separation_kpc = np.linalg.norm(mw_com - m31_com, axis=1)

    # Per-particle static color array.
    colors = np.empty((N, 4), dtype=np.uint8)
    colors[:half] = [74, 144, 217, 200]   # MW: blue
    colors[half:] = [230, 69, 69, 200]    # M31: red

    rr.init("MW-Andromeda Collision", spawn=(args.save is None))
    if args.save:
        rr.save(args.save)
        print(f"  Saving to {args.save}")

    rr.log("world", rr.ViewCoordinates.RIGHT_HAND_Z_UP, static=True)

    blueprint = rrb.Blueprint(
        rrb.Horizontal(
            rrb.Spatial3DView(origin="world", name="Collision"),
            rrb.TimeSeriesView(
                origin="diagnostics",
                name="MW-M31 separation (kpc)",
                contents=["diagnostics/separation_kpc"],
            ),
            column_shares=[3, 1],
        ),
        rrb.TimePanel(state="expanded"),
    )
    rr.send_blueprint(blueprint)

    print("Logging frames...")
    for f in range(n_frames):
        rr.set_time("timestep", sequence=f)
        rr.set_time("time_myr", sequence=int(times[f] / MYR))

        rr.log("world/galaxies", rr.Points3D(
            positions=pos[f].astype(np.float32),
            colors=colors,
            radii=args.radius,
        ))

        rr.log("diagnostics/separation_kpc", rr.Scalars(separation_kpc[f]))

        if (f + 1) % 50 == 0 or f == 0 or f == n_frames - 1:
            pct = 100.0 * (f + 1) / n_frames
            print(f"  Frame {f+1}/{n_frames} ({pct:.0f}%)")

    print("Done!")
    if args.save:
        print(f"Open with: rerun {args.save}")
    else:
        print("Scrub the timeline in the Rerun viewer to navigate.")


if __name__ == "__main__":
    main()
