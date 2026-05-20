import struct, re, argparse
import numpy as np
import rerun as rr
import rerun.blueprint as rrb
from pathlib import Path
from collections import OrderedDict

# Usage:
#   python render.py                              # -> opens Rerun viewer
#   python render.py --save sim.rrd               # -> saves to .rrd file
#   python render.py --sim tests/sim_output.bin   # -> custom binary path

AU = 1.496e11  # meters
G  = 6.6743e-11

# Body Metadata

PLANETS = {"Sun", "Mercury", "Venus", "Earth", "Mars",
           "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"}

COLORS = {
    "Sun":      [253, 184, 19],
    "Mercury":  [140, 120, 83],
    "Venus":    [230, 195, 92],
    "Earth":    [74, 144, 217],
    "Moon":     [170, 170, 170],
    "Mars":     [193, 68, 14],
    "Jupiter":  [200, 139, 58],
    "Saturn":   [232, 209, 145],
    "Uranus":   [115, 194, 198],
    "Neptune":  [59, 91, 165],
    "Pluto":    [210, 180, 140],
}
DEFAULT_COLOR = [120, 120, 120]

RADII = {
    "Sun": 0.20, "Jupiter": 0.12, "Saturn": 0.10, "Uranus": 0.07,
    "Neptune": 0.07, "Earth": 0.06, "Venus": 0.06, "Mars": 0.05,
    "Mercury": 0.04, "Pluto": 0.03,
}
DEFAULT_PLANET_RADIUS = 0.06

TRAIL_WIDTHS = {
    "Sun": 0.04, "Jupiter": 0.03, "Saturn": 0.03, "Uranus": 0.025,
    "Neptune": 0.025, "Earth": 0.02, "Venus": 0.02, "Mars": 0.02,
    "Mercury": 0.015, "Pluto": 0.015,
}
DEFAULT_TRAIL_WIDTH = 0.008

# Binary Loader

def load_binary(path):
    """Load simulation binary. Returns (names, pos_m, vel_ms).

    pos_m:  (n_frames, N, 3) in meters
    vel_ms: (n_frames, N, 3) in m/s
    """
    with open(path, "rb") as f:
        num_bodies = struct.unpack("Q", f.read(8))[0]
        names = []
        for _ in range(num_bodies):
            raw = f.read(32)
            names.append(raw.split(b'\x00')[0].decode())

        doubles_per_frame = 2 + num_bodies * 6
        frame_size = doubles_per_frame * 8

        raw = f.read()

    n_frames = len(raw) // frame_size
    all_data = np.frombuffer(raw[:n_frames * frame_size], dtype=np.float64)
    all_data = all_data.reshape(n_frames, doubles_per_frame)

    states = all_data[:, 2:].reshape(n_frames, num_bodies, 6)
    pos_m  = states[:, :, 0:3]
    vel_ms = states[:, :, 3:6]

    return names, pos_m, vel_ms


def load_masses(names):
    """Parse masses from Body.hpp, return array ordered by names."""
    body_hpp = Path("src/Body.hpp")
    masses = {}
    if body_hpp.exists():
        text = body_hpp.read_text()
        for m in re.finditer(r'\{\s*"(\w+)"\s*,\s*([0-9.eE+\-]+)', text):
            masses[m.group(1)] = float(m.group(2))
    return np.array([masses.get(n, 0.0) for n in names])


def read_config(path="src/Config.hpp"):
    try:
        text = Path(path).read_text()
        m = re.search(r'num_years\{\s*(\d+)\s*\}', text)
        return int(m.group(1)) if m else 100
    except FileNotFoundError:
        return 100


# Vectorized Diagnostics

def compute_all_diagnostics(pos, vel, mass, names):
    """Precompute all diagnostic timeseries (vectorized over frames).

    pos:  (F, N, 3) meters
    vel:  (F, N, 3) m/s
    mass: (N,) kg

    Returns dict of 1D arrays, each length F.
    """
    F, N, _ = pos.shape
    print("  Precomputing diagnostics...")

    # Kinetic energy: 0.5 * m * |v|^2 summed over bodies
    v_sq = np.sum(vel ** 2, axis=2)
    KE = 0.5 * np.sum(mass[None, :] * v_sq, axis=1)

    # Potential energy: -G * m_i * m_j / r_ij for all pairs
    PE = np.zeros(F)
    for i in range(N):
        for j in range(i + 1, N):
            rij = np.linalg.norm(pos[:, i, :] - pos[:, j, :], axis=1)
            rij = np.maximum(rij, 1e-10)
            PE -= G * mass[i] * mass[j] / rij

    total_E = KE + PE
    E0 = total_E[0]
    dE_rel = np.abs((total_E - E0) / E0) if E0 != 0 else np.zeros(F)

    # Angular momentum: L_i = m_i (r_i x v_i)
    cross = np.cross(pos, vel)
    L_vec = np.sum(mass[None, :, None] * cross, axis=1)
    L_mag = np.linalg.norm(L_vec, axis=1)
    L0 = L_mag[0]
    dL_rel = np.abs((L_mag - L0) / L0) if L0 != 0 else np.zeros(F)

    # Linear momentum: P_i = m_i * v_i
    P_vec = np.sum(mass[None, :, None] * vel, axis=1)
    P_mag = np.linalg.norm(P_vec, axis=1)
    P0 = P_mag[0]
    dP_rel = np.abs((P_mag - P0) / P0) if P0 != 0 else np.zeros(F)

    # Per-planet heliocentric distance (AU)
    sun_idx = names.index("Sun") if "Sun" in names else 0
    sun_pos = pos[:, sun_idx, :]
    distances = {}
    for i, name in enumerate(names):
        if name in PLANETS and name != "Sun":
            distances[name] = np.linalg.norm(pos[:, i, :] - sun_pos, axis=1) / AU

    print("  Diagnostics done.")
    return {
        "KE": KE, "PE": PE, "total_E": total_E, "dE_rel": dE_rel,
        "L_mag": L_mag, "dL_rel": dL_rel,
        "P_mag": P_mag, "dP_rel": dP_rel,
        "distances": distances,
    }


# Main

def main():
    parser = argparse.ArgumentParser(description="N-Body Rerun Visualizer")
    parser.add_argument("--sim", default="tests/sim_output.bin",
                        help="Path to simulation binary output")
    parser.add_argument("--save", default=None,
                        help="Save to .rrd file instead of spawning viewer")
    parser.add_argument("--trail-length", type=int, default=0,
                        help="Max trail points (0 = full trail)")
    args = parser.parse_args()

    # Load data
    print(f"Loading {args.sim}...")
    names, pos_m, vel_ms = load_binary(args.sim)
    mass = load_masses(names)
    num_years = read_config()
    n_frames = pos_m.shape[0]
    pos_au = pos_m / AU

    print(f"  {len(names)} bodies, {n_frames} frames, {num_years} years")
    if mass.sum() > 0:
        print(f"  Loaded masses for {np.count_nonzero(mass)} bodies from Body.hpp")

    # Precompute all diagnostics vectorized
    diag = compute_all_diagnostics(pos_m, vel_ms, mass, names)

    # Initialize Rerun
    rr.init("N-Body Solar System", spawn=(args.save is None))
    if args.save:
        rr.save(args.save)
        print(f"  Saving to {args.save}")

    # Static: coordinate system
    rr.log("world", rr.ViewCoordinates.RIGHT_HAND_Z_UP, static=True)

    # Default blueprint: all panels visible and laid out
    distance_planets = sorted(n for n in names if n in PLANETS and n != "Sun")

    blueprint = rrb.Blueprint(
        rrb.Horizontal(
            # Left: 3D view
            rrb.Spatial3DView(origin="world", name="Solar System"),
            # Right: stacked diagnostic plots
            rrb.Vertical(
                rrb.TimeSeriesView(
                    origin="energy",
                    name="Energy",
                    contents=["energy/total_energy", "energy/kinetic_energy", "energy/potential_energy"],
                ),
                rrb.TimeSeriesView(
                    origin="energy",
                    name="Energy Drift |ΔE/E₀|",
                    contents=["energy/dE_over_E0"],
                ),
                rrb.TimeSeriesView(
                    origin="momentum",
                    name="Momentum Drift",
                    contents=[
                        "momentum/angular_momentum_drift",
                        "momentum/linear_momentum_drift",
                    ],
                ),
                rrb.TimeSeriesView(
                    origin="momentum",
                    name="Momentum (absolute)",
                    contents=[
                        "momentum/angular_momentum",
                        "momentum/linear_momentum",
                    ],
                ),
                rrb.TimeSeriesView(
                    origin="distance",
                    name="Heliocentric Distance (AU)",
                    contents=[f"distance/{p}" for p in distance_planets],
                ),
            ),
            column_shares=[2, 1],
        ),
        rrb.TimePanel(fps=90.0),
    )
    rr.send_blueprint(blueprint)

    # Log each frame
    print("Logging frames...")
    for f_idx in range(n_frames):
        rr.set_time("timestep", sequence=f_idx)

        year = f_idx / max(n_frames - 1, 1) * num_years
        rr.set_time("year", sequence=int(year * 100))

        # 3D scene: planets only

        for i, name in enumerate(names):
            if name not in PLANETS:
                continue

            entity = f"world/{name}"
            x, y, z = pos_au[f_idx, i]

            color = COLORS.get(name, DEFAULT_COLOR)
            radius = RADII.get(name, DEFAULT_PLANET_RADIUS)

            rr.log(f"{entity}/position", rr.Points3D(
                positions=[[x, y, z]],
                colors=[color],
                radii=[radius],
                labels=[name],
            ))

            trail_start = 0
            if args.trail_length > 0:
                trail_start = max(0, f_idx - args.trail_length)

            if f_idx > trail_start:
                trail_pts = pos_au[trail_start:f_idx + 1, i, :].astype(np.float32)
                trail_color = color + [100]
                tw = TRAIL_WIDTHS.get(name, DEFAULT_TRAIL_WIDTH)

                rr.log(f"{entity}/trail", rr.LineStrips3D(
                    strips=[trail_pts],
                    colors=[trail_color],
                    radii=[tw],
                ))

        # Diagnostics: every frame

        # Energy
        rr.log("energy/total_energy",     rr.Scalars(diag["total_E"][f_idx]))
        rr.log("energy/kinetic_energy",   rr.Scalars(diag["KE"][f_idx]))
        rr.log("energy/potential_energy",  rr.Scalars(diag["PE"][f_idx]))
        rr.log("energy/dE_over_E0",       rr.Scalars(diag["dE_rel"][f_idx]))

        # Angular momentum
        rr.log("momentum/angular_momentum",       rr.Scalars(diag["L_mag"][f_idx]))
        rr.log("momentum/angular_momentum_drift",  rr.Scalars(diag["dL_rel"][f_idx]))

        # Linear momentum
        rr.log("momentum/linear_momentum",        rr.Scalars(diag["P_mag"][f_idx]))
        rr.log("momentum/linear_momentum_drift",   rr.Scalars(diag["dP_rel"][f_idx]))

        # Per-planet heliocentric distance
        for pname, d_arr in diag["distances"].items():
            rr.log(f"distance/{pname}", rr.Scalars(d_arr[f_idx]))

        # Progress
        if (f_idx + 1) % 100 == 0 or f_idx == 0 or f_idx == n_frames - 1:
            pct = 100.0 * (f_idx + 1) / n_frames
            print(f"  Frame {f_idx + 1}/{n_frames} ({pct:.0f}%)")

    print("Done!")
    if args.save:
        print(f"Open with: rerun {args.save}")
    else:
        print("Use the timeline slider in the Rerun viewer to scrub through timesteps.")


if __name__ == "__main__":
    main()