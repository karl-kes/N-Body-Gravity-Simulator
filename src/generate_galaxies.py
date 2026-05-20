import argparse, struct
import numpy as np
from pathlib import Path

# Usage:
#     python src/generate_galaxies.py --n 25000

# Binary IC format (matches main.cpp and jpl_compare.py):
#     Header:
#         uint64_t N
#         double   dt          (s)
#         double   total_time  (s)
#         double   output_dt   (s)
#         double   eps         (m)
#     Names (32 * N bytes, null-padded ASCII)
#     Particle data: 7 doubles per particle (mass, x, y, z, vx, vy, vz) in SI

# Layout: first N/2 particles are MW, last N/2 are M31. visualize_galaxies.py
# and galaxy_compare.py depend on this ordering.

# SI constants:
G_SI     = 6.6743e-11
M_SUN    = 1.989e30
KPC      = 3.086e19
KM_PER_S = 1e3

# Time constants for scenario parameters:
SEC_PER_YR = 365.25 * 86400.0
MYR_S      = 1.0e6 * SEC_PER_YR
GYR_S      = 1.0e9 * SEC_PER_YR

# Galaxy parameters (van der Marel et al. 2012, simplified to Plummer profile)
MW_MASS   = 1.5e12 * M_SUN
M31_MASS  = 1.5e12 * M_SUN
MW_SCALE  = 30.0 * KPC
M31_SCALE = 30.0 * KPC

# Initial conditions (van der Marel 2012):
SEPARATION     = 770.0 * KPC          # Current MW-M31 separation
RADIAL_VEL     = 109.0 * KM_PER_S     # Approach velocity along separation axis
TANGENTIAL_VEL = 17.0  * KM_PER_S     # Transverse component (M31 proper motion)

# Scenario parameters for the IC header:
GALAXY_DT          = 1.0 * MYR_S      # 1 Myr timestep
GALAXY_TOTAL_TIME  = 6.0 * GYR_S      # 6 Gyr integration
GALAXY_OUTPUT_DT   = 10.0 * MYR_S     # 10 Myr cadence
GALAXY_EPS         = 100.0 * 3.086e16  # 100 pc softening (kpc * 1e-3 in meters)


def sample_plummer_positions(n, a, rng):
    """Inverse-CDF sample n positions from a Plummer sphere of scale radius a.

    r = a / sqrt(u^(-2/3) - 1), with isotropic angles.
    Clamp u to avoid r = 0 (degenerate) and r >> a (outlier dominates the AABB).
    """
    u = rng.uniform(1e-6, 1.0 - 1e-3, n)
    r = a / np.sqrt(u**(-2.0 / 3.0) - 1.0)
    cos_theta = rng.uniform(-1.0, 1.0, n)
    sin_theta = np.sqrt(1.0 - cos_theta**2)
    phi = rng.uniform(0, 2 * np.pi, n)
    x = r * sin_theta * np.cos(phi)
    y = r * sin_theta * np.sin(phi)
    z = r * cos_theta
    return np.column_stack([x, y, z])


def sample_plummer_speeds(n, rng):
    """
    Rejection-sample n speeds q = v/v_esc from the Plummer DF:
        f(q) ~ q^2 (1 - q^2)^(7/2),  q in [0, 1]
    Maximum at q^2 = 2/9, value (2/9) * (7/9)^(7/2) ~= 0.0922.
    """
    g_max = (2.0 / 9.0) * (7.0 / 9.0)**3.5
    out = np.empty(n)
    filled = 0
    while filled < n:
        batch = max((n - filled) * 12, 1024)
        q = rng.uniform(0.0, 1.0, batch)
        u = rng.uniform(0.0, g_max, batch)
        accepted = q[u < q**2 * (1.0 - q**2)**3.5]
        take = min(len(accepted), n - filled)
        out[filled:filled + take] = accepted[:take]
        filled += take
    return out


def sample_plummer_velocities(positions, M, a, rng):
    """Self-consistent isotropic velocities given Plummer potential at each r.

    v_esc(r)^2 = 2GM/sqrt(r^2 + a^2). Speed from Plummer DF, direction isotropic.
    """
    n = len(positions)
    r = np.linalg.norm(positions, axis=1)
    v_esc = np.sqrt(2.0 * G_SI * M / np.sqrt(r**2 + a**2))
    v = sample_plummer_speeds(n, rng) * v_esc

    cos_theta = rng.uniform(-1.0, 1.0, n)
    sin_theta = np.sqrt(1.0 - cos_theta**2)
    phi = rng.uniform(0, 2 * np.pi, n)
    vx = v * sin_theta * np.cos(phi)
    vy = v * sin_theta * np.sin(phi)
    vz = v * cos_theta
    return np.column_stack([vx, vy, vz])


def make_galaxy(n, M, a, com_pos, com_vel, rng):
    """Self-consistent Plummer galaxy at given COM phase-space point."""
    pos_rel = sample_plummer_positions(n, a, rng)
    vel_rel = sample_plummer_velocities(pos_rel, M, a, rng)
    mass = np.full(n, M / n)
    return pos_rel + com_pos, vel_rel + com_vel, mass


def write_ic_file(path, mass, pos, vel, name_prefixes, dt, total_time, output_dt, eps):
    """Write the unified binary IC format."""
    N = len(mass)
    Path(path).parent.mkdir(exist_ok=True, parents=True)
    with open(path, "wb") as f:
        # Header
        f.write(struct.pack("Q", N))
        f.write(struct.pack("d", dt))
        f.write(struct.pack("d", total_time))
        f.write(struct.pack("d", output_dt))
        f.write(struct.pack("d", eps))

        # Names: prefix + zero-padded index, 32 bytes each.
        for i, prefix in enumerate(name_prefixes):
            name = f"{prefix}{i}"
            buf = name.encode("ascii", errors="replace")[:31]
            f.write(buf + b"\x00" * (32 - len(buf)))

        # Particle data
        data = np.column_stack([mass, pos, vel]).astype(np.float64)
        f.write(data.tobytes())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=25000,
                    help="Particles per galaxy (default: 25000, total 50000)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--output", default="tests/sim_ic.bin")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    n_per = args.n
    n_total = 2 * n_per

    print(f"Scenario: MW-Andromeda collision")
    print(f"  Particles: {n_total} ({n_per} per galaxy)")
    print(f"  dt: {GALAXY_DT / MYR_S:.1f} Myr")
    print(f"  total time: {GALAXY_TOTAL_TIME / GYR_S:.1f} Gyr")
    print(f"  output cadence: {GALAXY_OUTPUT_DT / MYR_S:.1f} Myr")
    print(f"  softening: {GALAXY_EPS / KPC * 1000:.0f} pc\n")

    # MW at -separation/2 moving toward +x; M31 mirrored.
    # Tangential split symmetrically so net angular momentum is consistent.
    mw_pos  = np.array([-SEPARATION / 2, 0.0, 0.0])
    mw_vel  = np.array([+RADIAL_VEL / 2, +TANGENTIAL_VEL / 2, 0.0])
    m31_pos = np.array([+SEPARATION / 2, 0.0, 0.0])
    m31_vel = np.array([-RADIAL_VEL / 2, -TANGENTIAL_VEL / 2, 0.0])

    mw_p,  mw_v,  mw_m  = make_galaxy(n_per, MW_MASS,  MW_SCALE,  mw_pos,  mw_vel,  rng)
    m31_p, m31_v, m31_m = make_galaxy(n_per, M31_MASS, M31_SCALE, m31_pos, m31_vel, rng)

    pos  = np.vstack([mw_p, m31_p])
    vel  = np.vstack([mw_v, m31_v])
    mass = np.concatenate([mw_m, m31_m])

    # Subtract total COM phase-space point so the binary stays centered.
    total_mass = mass.sum()
    com_p = (pos * mass[:, None]).sum(axis=0) / total_mass
    com_v = (vel * mass[:, None]).sum(axis=0) / total_mass
    pos -= com_p
    vel -= com_v

    # Name prefixes: first N/2 MW, last N/2 M31.
    name_prefixes = ["MW_"] * n_per + ["M31_"] * n_per

    write_ic_file(
        args.output, mass, pos, vel, name_prefixes,
        GALAXY_DT, GALAXY_TOTAL_TIME, GALAXY_OUTPUT_DT, GALAXY_EPS,
    )

    bytes_out = 40 + n_total * 32 + n_total * 7 * 8
    print(f"-> {args.output} ({n_total} particles, {bytes_out / 1e6:.1f} MB)")
    print(f"   Total mass:     {total_mass / M_SUN:.2e} M_sun")
    print(f"   Particle mass:  {mass[0] / M_SUN:.2e} M_sun")
    print(f"   Approach speed: {RADIAL_VEL / KM_PER_S:.0f} km/s")


if __name__ == "__main__":
    main()
