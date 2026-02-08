# N-Body Physics Engine

A high-performance N-body gravitational simulator written in C++17, validated against NASA JPL Horizons ephemeris data. Simulates 35 solar system bodies (10 planets + 25 moons) over 100 years with 0.05% average positional deviation from NASA's DE441 planetary ephemeris and machine-precision energy conservation.

## Results

**100-year simulation, 35 bodies, dt = 3600s, Yoshida 4th-order symplectic integrator:**

| Metric | Value |
|---|---|
| Average relative error | 0.05% |
| Worst-case error (Ganymede) | 0.29% |
| Planet error range | 0.00003% – 0.10% |
| Energy conservation | 4.2e-8 % drift |
| Runtime (35 bodies, 100 yr) | ~7009 milliseconds |

**Error by body class:**

- **Gas giants** (Jupiter, Saturn, Uranus, Neptune): < 0.002% — essentially exact
- **Inner planets** (Venus, Earth, Mars): < 0.03%
- **Mercury**: 0.10% — highest planet error due to missing solar oblateness (J₂) and asteroid perturbations
- **Moons** (Io, Europa, Ganymede, etc.): 0.02% – 0.29% — phase drift from fast orbital periods relative to timestep

Residual errors are from physics model differences (Newtonian + 10 bodies vs. JPL DE441's ~300 bodies with full PPN relativity, solar oblateness, and tidal effects), not integration error.

## Features

- **Yoshida 4th-order symplectic integrator** with Velocity Verlet alternative
- **Automated JPL Horizons validation pipeline** — fetches reference ephemerides and compares against simulation output
- **35 solar system bodies** with JPL-sourced initial conditions (10 planets, 25 moons)
- **OpenMP parallelization** with SIMD vectorization
- **Interactive 3D visualization** via matplotlib
- **Single-file configuration** — change `dt`, `num_years`, or `output_hours` in `Config.hpp` and both the C++ simulator and Python validation tools adapt automatically

## Project Structure

```
├── main.cpp                    # Entry point
├── Config.hpp                  # Single-source configuration (dt, years, output interval)
├── Body.hpp                    # JPL Horizons initial conditions (auto-generated)
├── Output.hpp                  # CSV output handler
├── jpl_compare.py              # JPL fetch + validation (reads Config.hpp)
├── visualize.py                # Interactive 3D orbit viewer
├── Classes/
│   ├── Force/                  # Gravitational force computation
│   ├── Integrator/             # Yoshida 4th-order, Velocity Verlet
│   ├── Particle/               # SoA particle data structure
│   └── Simulation/             # Main simulation loop, energy diagnostics
└── validation/                 # Output directory for sim + reference data
    ├── sim_output.csv
    ├── jpl_reference.csv
    └── body_catalog.json
```

## Quick Start

### Prerequisites

- C++17 compiler (GCC / Clang / MSVC)
- Python 3.8+ with `numpy` and `matplotlib`
- Internet connection (for JPL Horizons fetch)

### Build & Run

```bash
# 1. Fetch JPL reference data (reads Config.hpp for step size and duration)
python jpl_compare.py fetch --moons

# 2. Compile
g++ -std=c++17 -O3 -march=native -fopenmp *.cpp Classes/Force/*.cpp Classes/Integrator/*.cpp Classes/Particle/*.cpp Classes/Simulation/*.cpp -o main.exe

# 3. Run simulation
./main.exe

# 4. Validate against JPL
python jpl_compare.py compare

# 5. Visualize
python visualize.py
```

## Configuration

All simulation parameters are controlled from `Config.hpp`:

```cpp
inline static constexpr double dt{ 3600.0 };              // Integration timestep (seconds)
inline static constexpr std::size_t num_years{ 100 };     // Simulation duration
inline static constexpr std::size_t output_hours{ 487 };  // Output interval (hours)
```

The Python scripts read `num_years` and `output_hours` directly from this file, so there's a single source of truth. A compile-time `static_assert` catches invalid `dt` / `output_hours` combinations.

`output_hours` must satisfy: `(output_hours × 3600) % dt == 0`, and JPL Horizons must accept it as an integer hour step.

## How It Works

### Integration

The simulator uses the **Yoshida 4th-order symplectic integrator**, which splits each timestep into a sequence of position and momentum sub-steps with specially chosen coefficients. This preserves the symplectic structure of Hamilton's equations, guaranteeing bounded energy error over arbitrarily long integrations — unlike Runge-Kutta methods which exhibit secular energy drift.

The 4th-order scheme achieves O(dt⁴) local truncation error with only 3 force evaluations per step, compared to 4 for classical RK4.

### Force Computation

Pairwise Newtonian gravity with O(N²) direct summation. The inner loop uses a branchless mask (`i == j ? 0 : 1`) to skip self-interaction without breaking SIMD vectorization. For N > 500 bodies, OpenMP parallelization activates automatically.

### Validation Pipeline

`jpl_compare.py` automates the full validation workflow:

1. **Fetch**: Queries NASA JPL Horizons API for state vectors of all 35 bodies at intervals matching the simulation output
2. **Compare**: Loads both datasets, computes position errors at each epoch, and reports max relative error per body

The comparison uses index-aligned matching — both the simulation and JPL reference output at exactly the same time intervals (`output_hours` hours), eliminating interpolation error.

### Visualization

`visualize.py` opens an interactive matplotlib 3D window with:

- Keyboard controls (Space = play/pause, arrows = step/speed, T = trails, R = reset)
- Colored trails with planet labels
- Automatic axis scaling from trajectory data

## Technical Notes

**Why 0.05% and not 0.00%?** The simulator solves Newtonian gravity for 35 bodies. JPL DE441 integrates ~300 bodies with full post-Newtonian relativity (PPN), solar oblateness (J₂), asteroid perturbations (Ceres, Pallas, Vesta), and tidal effects. The 0.05% residual is the gap between these two physics models, not integration error — as confirmed by the 8.8 × 10⁻¹¹% energy drift.

**Why not include relativistic corrections?** The EIH (Einstein-Infeld-Hoffmann) 1PN equations were implemented and tested. However, the velocity-dependent PN forces break the symplecticity of the Yoshida integrator, degrading energy conservation by 4 orders of magnitude (from 10⁻¹¹% to 10⁻⁷%). Over 100 years, this symplectic violation costs more accuracy than the 43 arcsec/century of Mercury precession it's intended to capture. A proper fix requires an implicit or splitting integrator for velocity-dependent potentials.

**Moon errors are phase drift, not integration failure.** Io orbits Jupiter every 1.77 days. With dt = 360s, that's ~425 steps per orbit — adequate for 4th-order accuracy per orbit, but phase errors accumulate over ~20,000 orbits in 100 years. Reducing dt to 60–120s or implementing adaptive timestepping would improve moon accuracy at the cost of runtime.