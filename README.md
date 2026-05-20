# N-Body Gravity Engine

A high-performance N-body gravitational simulator in C++, with two physical scenarios:

- **Solar system**: 35 bodies (Sun, 8 planets, Pluto, 25 natural satellites) over multi-century timescales, validated against NASA JPL Horizons DE441 ephemeris data.
- **Galaxy collision**: MW-Andromeda encounter, two Plummer galaxies (configurable N), validated against van der Marel et al. 2012 predicted timeline.

Both scenarios use a 4th-order Yoshida symplectic integrator and choose between a SIMD-vectorized direct O($N^2$) kernel or a Barnes-Hut O($N \log N$) octree kernel.

**[Technical Paper](docs/N_Body_Technical_Paper.pdf)** covers the symplectic integrator derivation, force model, validation methodology, and error analysis.

---

## Quick Start

```bash
# Solar system, full pipeline:
./scripts/run.sh --fetch --visualize

# Galaxy collision, Barnes-Hut, full pipeline:
./scripts/run.sh --galaxy --force bh --visualize

# Tests:
./scripts/test.sh
```

On Windows replace `./scripts/run.sh` with `.\scripts\run.ps1` and `--flags` with `-Flags` (`-Fetch`, `-Galaxy`, etc).

---

## Scripts

Two scripts handle everything. Both have Linux/macOS (`.sh`) and Windows (`.ps1`) versions with identical behavior and naming.

### `scripts/run.sh` / `scripts/run.ps1`

Full simulation pipeline: stage the initial conditions, build the code, run the simulation, validate against reference data (solar only), open viewers (optional).

| Flag (bash) | Flag (PowerShell) | Default | Description |
|---|---|---|---|
| `--fetch` | `-Fetch` | (off) | Fetch solar system IC from JPL Horizons before running |
| `--galaxy` | `-Galaxy` | (off) | Generate MW-Andromeda IC before running |
| `--n N` | `-N N` | `25000` | Particles per galaxy (galaxy only) |
| `--start DATE` | `-Start DATE` | `1950-01-01` | Solar fetch start date |
| `--no-moons` | `-NoMoons` | (off) | Solar: planets only, no moons (10 bodies instead of 35) |
| `--force MODE` | `-Force MODE` | `direct` | Force kernel: `direct` or `bh` |
| `--theta T` | `-Theta T` | `0.5` | Barnes-Hut opening angle (only used with `--force bh`) |
| `--visualize` | `-Visualize` | (off) | Open matplotlib viewer after the simulation completes |
| `--render` | `-Render` | (off) | Open Rerun viewer after the simulation completes |

**Scenario selection.** `--fetch` and `--galaxy` are mutually exclusive. If you pass neither, the script uses whatever IC is already in `tests/sim_ic.bin` (from a previous run).

**Common usage:**

```bash
# Solar, default everything (uses existing IC if present):
./scripts/run.sh

# Solar, refetch JPL data first:
./scripts/run.sh --fetch

# Solar, planets only, Barnes-Hut, both viewers:
./scripts/run.sh --fetch --no-moons --force bh --visualize --render

# Galaxy, 50k particles per galaxy, Barnes-Hut, matplotlib viewer:
./scripts/run.sh --galaxy --n 50000 --force bh --visualize

# Galaxy, Rerun viewer with separation timeline:
./scripts/run.sh --galaxy --force bh --render
```

### `scripts/test.sh` / `scripts/test.ps1`

Build and run the test suite and/or the scaling benchmark.

| Flag (bash) | Flag (PowerShell) | Default | Description |
|---|---|---|---|
| `--benchmark` | `-Benchmark` | (off) | Run the scaling benchmark after unit tests |
| `--bench-only` | `-BenchOnly` | (off) | Skip unit tests, run benchmark only |
| `--trials N` | `-Trials N` | `3` | Trials per N value in the benchmark |
| `--max-n N` | `-MaxN N` | `65536` | Maximum N for benchmark sweep |

**Common usage:**

```bash
./scripts/test.sh                                       # run all unit tests
./scripts/test.sh --benchmark                           # tests, then scaling benchmark
./scripts/test.sh --bench-only --trials 5 --max-n 16384 # benchmark only
```

---

## Prerequisites

- C++23 compiler (GCC / Clang / MSVC / MinGW)
- CMake 3.25+
- OpenMP
- Python 3.8+ with `numpy` and `matplotlib`
- Internet connection for `--fetch` (JPL Horizons API)
- Optional: `pip install rerun-sdk` for the Rerun viewer (`--render`)

---

## What the scripts call under the hood

The runner scripts orchestrate Python tools, the C++ binary, and viewers. Invoke any of these directly if you want finer control.

### Solar tools (`tools/solar/`)

```bash
# Fetch IC from JPL Horizons:
python tools/solar/fetch.py                        # planets only
python tools/solar/fetch.py --moons                # full 35 bodies
python tools/solar/fetch.py --start 1900-01-01     # custom start date

# Validate simulation against JPL reference:
python tools/solar/compare.py                              # all bodies
python tools/solar/compare.py --bodies Mercury,Venus,Earth # subset

# Viewers:
python tools/solar/visualize.py                    # matplotlib
python tools/solar/visualize.py --speed 4          # start at 4x playback speed
python tools/solar/render.py                       # Rerun
```

### Galaxy tools (`tools/galaxy/`)

```bash
# Generate MW-Andromeda IC:
python tools/galaxy/generate.py                    # 25000 per galaxy (50k total)
python tools/galaxy/generate.py --n 50000          # 50000 per galaxy (100k total)
python tools/galaxy/generate.py --seed 7           # deterministic with custom seed

# Validate against van der Marel 2012 predicted timeline:
python tools/galaxy/compare.py

# Viewers:
python tools/galaxy/visualize.py                   # matplotlib
python tools/galaxy/render.py                      # Rerun
python tools/galaxy/render.py --save sim.rrd       # save to .rrd file for sharing
```

### Simulation binary

```bash
./build/main                                  # default: direct, theta=0.5
./build/main --force direct
./build/main --force bh --theta 0.5
./build/main --force bh --theta 0.3           # tighter MAC, more accurate, slower
./build/main --ic some_other.bin              # use a different IC file
./build/main --out tests/run42.bin            # write output to a different path
```

Reads `tests/sim_ic.bin` by default, writes `tests/sim_output.bin` by default.

### Matplotlib viewer controls

| Key | Action |
|---|---|
| Space | Play / Pause |
| Right / Left | Step forward / backward |
| Up / Down | Speed up (2×) / slow down (0.5×) |
| R | Reset to frame 0 |
| Q | Quit |

---

## Key Results

**Solar system, 249-year simulation (1950–2199), 35 bodies, Δt = 360 s, Yoshida 4th-order:**

| Metric | Value |
|---|---|
| Mean max relative error (ex Sun) | 0.079% |
| Worst-case body | Callisto (0.411%) |
| Planet error range | 0.0007% – 0.417% |
| Energy conservation ΔE/E | 1.5 × 10⁻¹² |
| Angular momentum ΔL/L | 4.3 × 10⁻¹³ |

Residual errors are attributed to physics model differences (Newtonian gravity with 35 bodies vs. JPL DE441's ~300 bodies with post-Newtonian relativity, solar oblateness J₂, asteroid perturbations, and tidal effects), not numerical integration error. A timestep convergence study confirms that Δt ≤ 900 s lies on the floating-point precision floor.

**Galaxy collision, 6 Gyr simulation, 50k particles, Barnes-Hut θ = 0.5:**

First close passage predicted by van der Marel 2012 at $t \approx 3.9$ Gyr, final merger around $t \approx 5.9$ Gyr. Validation plot in `tests/galaxy_separation.png`.

---

## Project Structure

```
├── src/                            # C++ source (build inputs)
│   ├── main.cpp                    # Entry point (reads tests/sim_ic.bin)
│   ├── Config.hpp                  # Universal constants (G, time units)
│   ├── Force/                      # Direct O(N²) SIMD + Barnes-Hut O(N log N) octree
│   ├── Integrator/                 # Yoshida 4th-order + Velocity Verlet
│   ├── Particle/                   # SoA particle data (single contiguous allocation)
│   ├── Simulation/                 # Time-stepping loop, conservation diagnostics
│   └── Output/                     # Binary output format
│
├── tools/                          # Python tooling (auxiliary)
│   ├── solar/
│   │   ├── fetch.py                # JPL Horizons fetcher -> tests/sim_ic.bin
│   │   ├── compare.py              # Validate against JPL reference
│   │   ├── visualize.py            # Matplotlib viewer
│   │   └── render.py               # Rerun dashboard
│   └── galaxy/
│       ├── generate.py             # Plummer-sphere sampler -> tests/sim_ic.bin
│       ├── compare.py              # Separation-vs-time validation
│       ├── visualize.py            # Matplotlib viewer
│       └── render.py               # Rerun dashboard
│
├── scripts/
│   ├── run.sh / run.ps1            # Simulation pipeline runner
│   └── test.sh / test.ps1          # Test & benchmark runner
│
├── tests/
│   ├── unit_tests.cpp              # 23 unit tests
│   ├── benchmark.cpp               # Scaling benchmark
│   ├── sim_ic.bin                  # Initial conditions (generated, gitignored)
│   ├── sim_output.bin              # Simulation output (generated, gitignored)
│   ├── jpl_reference.csv           # JPL reference ephemeris (generated)
│   └── comparison_results.json     # Validation summary (generated)
│
├── docs/
│   └── N_Body_Technical_Paper.pdf
├── CMakeLists.txt
└── README.md
```

---

## Configuration

Scenario parameters (dt, total time, output cadence, softening) are determined by the IC file, not by `Config.hpp`. Each IC generator writes a header with the appropriate values for its scenario:

- Solar: dt = 900 s, 249 years, output every 487 hours, ε = 10⁻⁹ m
- Galaxy: dt = 1 Myr, 6 Gyr, output every 10 Myr, ε = 100 pc

To customize, edit the constants near the top of `tools/solar/fetch.py` (`SOLAR_*`) or `tools/galaxy/generate.py` (`GALAXY_*`), regenerate the IC, and rerun. `Config.hpp` only holds universal constants (G, time unit conversions, OpenMP threshold).

---

## Implementation Notes

**Yoshida 4th-order symplectic integrator.** Four position drifts interleaved with three force evaluations per timestep, achieving O(Δt⁴) local truncation error. The negative intermediate coefficient introduces a backward sub-step essential for error cancellation. Energy errors remain bounded and oscillatory over arbitrarily long integrations.

**Structure-of-Arrays memory layout.** All particle data occupies a single contiguous allocation with SIMD-aligned sub-arrays (AVX2: 32-byte, AVX-512: 64-byte). Each sub-array is padded to a SIMD-width boundary so that every array start is naturally aligned for vectorized loads/stores. `__restrict__`-qualified pointers enable SIMD auto-vectorization.

**Branchless force kernel.** Self-interaction is eliminated with a floating-point mask rather than a conditional branch, preserving SIMD vectorization. Newton's third law symmetry is intentionally not exploited; the doubled FLOP count is traded for regular memory access patterns and freedom from race conditions under OpenMP.

**Barnes-Hut octree.** A second `Force` implementation provides O(N log N) gravity for larger ensembles. The tree is rebuilt every timestep via top-down counting-sort partition into octants; storage capacity is sticky across calls so only the size resets. Each node stores center of mass, total mass, bounding-box geometry, eight child pointers, and an explicit `is_leaf` flag (disambiguating leaves from internal nodes with sparse octants).

**Binary IC format.** Both IC generators write a single binary file with a 40-byte header (N, dt, total_time, output_dt, eps) followed by 32-byte name records and 7-double-per-particle phase-space coordinates. The simulator reads this directly and configures itself from the header; no scenario-specific code paths exist in the C++ side.