# N-Body Gravity Engine

A high-performance N-body gravitational simulator in C++23, with two physical scenarios:

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

The native build targets Linux. On Windows, use the Bash scripts under WSL2.

---

## Scripts

The Bash scripts handle the supported Linux/WSL2 workflow. The older PowerShell/MinGW variants remain in the tree for reference, but they are not compatible with the required GCC 15/CUDA 13.3 toolchain.

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
| `--fp32` | — | (off) | Build with `fp_t` as `float` instead of `double` |

**Scenario selection.** `--fetch` and `--galaxy` are mutually exclusive. If you pass neither, the script uses whatever IC is already in `tests/sim_ic.bin` (from a previous run).

**Common usage:**

```bash
./scripts/run.sh                                              # use existing IC
./scripts/run.sh --fetch                                      # solar, refetch JPL
./scripts/run.sh --fetch --no-moons --force bh --visualize    # solar, BH, planets only
./scripts/run.sh --galaxy --n 50000 --force bh --visualize    # galaxy, 100k total
./scripts/run.sh --galaxy --force bh --render                 # galaxy, Rerun viewer
```

### `scripts/test.sh` / `scripts/test.ps1`

Build and run the test suite and/or the scaling benchmark.

| Flag (bash) | Flag (PowerShell) | Default | Description |
|---|---|---|---|
| `--benchmark` | `-Benchmark` | (off) | Run the scaling benchmark after unit tests |
| `--bench-only` | `-BenchOnly` | (off) | Skip unit tests, run benchmark only |
| `--trials N` | `-Trials N` | `3` | Trials per N value in the benchmark |
| `--max-n N` | `-MaxN N` | `65536` | Maximum N for benchmark sweep |
| `--fp32` | — | (off) | Build with `fp_t` as `float` instead of `double` |

**Common usage:**

```bash
./scripts/test.sh                                       # run all unit tests
./scripts/test.sh --benchmark                           # tests, then scaling benchmark
./scripts/test.sh --bench-only --trials 5 --max-n 16384 # benchmark only
```

---

## Prerequisites

- GCC 15.x (`g++-15`)
- NVIDIA CUDA Toolkit 13.3.x, with GCC 15 as its host compiler
- CMake 3.25+
- OpenMP
- Internet access during the first configure so CMake can fetch the pinned xpu dependency
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

### Solar system

249-year simulation (1950–2199), 35 bodies, Δt = 360 s, Yoshida 4th-order, direct kernel:

| Metric | Value |
|---|---|
| Mean max relative error (ex Sun) | 0.079% |
| Worst-case body | Callisto (0.411%) |
| Planet error range | 0.0007% – 0.417% |
| Energy conservation ΔE/E | 1.5 × 10⁻¹² |
| Angular momentum ΔL/L | 4.3 × 10⁻¹³ |

Residual errors are attributed to physics model differences (Newtonian gravity with 35 bodies vs. JPL DE441's ~300 bodies with post-Newtonian relativity, solar oblateness $J_2$, asteroid perturbations, and tidal effects), not numerical integration error. A timestep convergence study confirms that Δt ≤ 900 s lies on the floating-point precision floor.

### Galaxy collision

6 Gyr simulation, 50,000 particles ($N = 25{,}000$ per galaxy), Barnes-Hut at θ = 0.5, Plummer profile, $\epsilon = 100$ pc softening:

| Metric | Simulation | Reference (van der Marel 2012) |
|---|---|---|
| First close passage | 3.72 Gyr | ~3.9 Gyr |
| Initial separation | 769.8 kpc | 770 kpc |
| Energy drift ΔE/E | 1.44% | — |
| Angular momentum drift ΔL/L | 0.19% | — |
| Wall-clock runtime | 52.7 min | — |

First-passage time matches van der Marel's 2012 prediction within 5%, well inside the ~10-15% uncertainty band on the literature value (driven by observational uncertainty in M31's tangential velocity). Subsequent local minima in COM separation are post-merger oscillations of intermixed particle populations, not independent encounters; the galaxy cores have phase-mixed by t ≈ 5 Gyr. Validation plot saved to `tests/galaxy_separation.png` by `tools/galaxy/compare.py`.

### Scaling benchmark

12-thread OpenMP, Yoshida dt = 900 s, BH θ = 0.5. Single-threaded measurements are skipped above $N = 8192$ to keep wall-time tractable at large N:

| N | Steps | Direct OMP (ms) | BH OMP (ms) | BH speedup |
|---|---|---|---|---|
| 512 | 2020 | 1447.4 | 1316.1 | 1.10× |
| 1024 | 788 | 1930.5 | 1078.8 | 1.79× |
| 2048 | 206 | 1940.8 | 777.5 | 2.50× |
| 4096 | 54 | 1942.8 | 489.6 | 3.97× |
| 8192 | 14 | 2000.5 | 211.7 | 9.45× |
| 16384 | 3 | 1671.0 | 128.3 | 13.02× |
| 32768 | 3 | 6651.6 | 373.6 | 17.81× |
| 65536 | 3 | 27658.9 | 630.5 | 43.87× |
| 131072 | 3 | 106792.5 | 2030.5 | 52.59× |

Per-step direct kernel cost quadruples for every doubling of N, confirming textbook $O(N^2)$ scaling across more than two orders of magnitude in $N$. Barnes-Hut wins at every measured particle count, including the smallest tested ($N = 512$, 1.10× speedup); the crossover point sits below the practical range of the benchmark. At $N = 131{,}072$, Barnes-Hut delivers a **52× speedup** over the already-SIMD-vectorized OpenMP direct kernel — the practical justification for shipping the BH implementation in the codebase.

The galaxy collision demo ($N = 50{,}000$) sits firmly in the BH-dominant regime. Extrapolating from the table, direct mode at that particle count would have taken approximately 16-18 hours instead of the 52.7 minutes measured.

---

## Project Structure

```
├── src/                            # C++ source (build inputs)
│   ├── main.cu                     # Entry point (reads tests/sim_ic.bin)
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
│   ├── unit_tests/                 # 23 unit tests split by module
│   ├── benchmark.cu                # Scaling benchmark
│   ├── sim_ic.bin                  # Initial conditions (generated, gitignored)
│   ├── sim_output.bin              # Simulation output (generated, gitignored)
│   ├── jpl_reference.csv           # JPL reference ephemeris (generated)
│   ├── comparison_results.json     # Solar validation summary (generated)
│   └── galaxy_separation.png       # Galaxy validation plot (generated)
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
