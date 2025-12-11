# N-Body Gravitational Simulation

A high-performance N-body gravitational simulation written in C++ with 3D visualization support using Python and VPython.

## Overview

This project simulates the gravitational interactions between multiple celestial bodies using Newtonian mechanics. It implements the Velocity Verlet integration method for accurate and energy-conserving orbital calculations, making it suitable for simulating planetary systems, star clusters, or any gravitationally-bound system.

## Features

- **Velocity Verlet Integration**: Second-order symplectic integrator that conserves energy over long simulations
- **Softening Parameter**: Prevents numerical instabilities during close encounters
- **Energy Drift Monitoring**: Tracks simulation accuracy by monitoring total energy conservation
- **CSV Import/Export**: Easy configuration of initial conditions and trajectory output
- **3D Visualization**: Real-time animated rendering with VPython
- **OpenMP Support**: Optional parallelization for large-scale simulations

## Project Structure

```
├── main.cpp                 # Entry point
├── Constants.hpp            # Physical constants (G, softening, etc.)
├── bodies.csv               # Input: initial positions, velocities, masses
├── trajectories.csv         # Output: simulation results
├── render.py                # 3D visualization script
└── Classes/
    ├── Vec_3D/
    │   ├── Vec_3D.hpp       # 3D vector class declaration
    │   └── Vec_3D.cpp       # 3D vector implementation
    ├── Body/
    │   ├── Body.hpp         # Celestial body class declaration
    │   └── Body.cpp         # Body physics implementation
    └── Simulation/
        ├── Simulation.hpp   # Simulation manager declaration
        └── Simulation.cpp   # Core simulation loop
```

## Getting Started

### Prerequisites

- C++ compiler with C++17 support (g++, clang++)
- Python 3.x (for visualization)
- VPython (`pip install vpython`)
- pandas (`pip install pandas`)
- numpy (`pip install numpy`)

### Compilation

```bash
g++ main.cpp Classes/Vec_3D/Vec_3D.cpp Classes/Body/Body.cpp Classes/Simulation/Simulation.cpp -o main.exe
```

For parallel execution with OpenMP (recommended for 1000+ bodies):

```bash
g++ main.cpp Classes/Vec_3D/Vec_3D.cpp Classes/Body/Body.cpp Classes/Simulation/Simulation.cpp -o main.exe -fopenmp
```

### Running the Simulation

1. **Configure initial conditions** in `bodies.csv`:

```csv
# x, y, z, vx, vy, vz, mass
0, 0, 0, 0, 0, 0, 1.989e30
1.496e11, 0, 0, 0, 29783, 0, 5.972e24
```

2. **Run the simulation**:

```bash
./main.exe
```

3. **Enter simulation parameters** when prompted:
   - Number of steps (e.g., 10000)

4. **Visualize the results**:

```bash
python render.py
```

## Input Format

The `bodies.csv` file defines initial conditions for each body:

| Column | Description | Units |
|--------|-------------|-------|
| x, y, z | Initial position | meters |
| vx, vy, vz | Initial velocity | m/s |
| mass | Body mass | kg |

### Example: Inner Solar System

```csv
# Sun
0, 0, 0, 0, 0, 0, 1.989e30
# Earth
1.496e11, 0, 0, 0, 29783, 0, 5.972e24
# Mars
2.279e11, 0, 0, 0, 24077, 0, 6.39e23
# Venus
1.082e11, 0, 0, 0, 35020, 0, 4.867e24
```

## Configuration

Physical constants can be adjusted in `Constants.hpp`:

| Constant | Default | Description |
|----------|---------|-------------|
| `G` | 6.67430e-11 | Gravitational constant (m³/kg/s²) |
| `EPSILON` | 1.0e3 | Softening parameter (m) |

**Timestep (`dt`)**: Default is 1000 seconds. Adjust in `Simulation` constructor for different time scales.

## Visualization Controls

The VPython renderer (`render.py`) provides interactive 3D visualization:

| Control | Action |
|---------|--------|
| **Spacebar** | Pause/Resume animation |
| **Mouse drag** | Rotate view |
| **Scroll** | Zoom in/out |
| **Right-click drag** | Pan view |

Bodies are colored automatically and sized logarithmically based on mass.

## Physics

### Gravitational Acceleration

The acceleration of body *i* due to all other bodies:

$$\vec{a}_i = \sum_{j \neq i} \frac{G \cdot m_j}{(|\vec{r}_{ij}|^2 + \epsilon^2)^{3/2}} \vec{r}_{ij}$$

Where ε is the softening parameter that prevents singularities.

### Velocity Verlet Integration

Position and velocity updates:

```
x(t+dt) = x(t) + v(t)·dt + 0.5·a(t)·dt²
v(t+dt) = v(t) + 0.5·(a(t) + a(t+dt))·dt
```

This method is time-reversible and provides excellent energy conservation.

## Performance Tips

| Body Count | Recommendation |
|------------|----------------|
| < 100 | Serial execution (no `-fopenmp`) |
| 100 - 1000 | Test both; serial often faster |
| > 1000 | Use OpenMP (`-fopenmp`) |

## Output

The simulation generates `trajectories.csv` with columns:

| Column | Description |
|--------|-------------|
| step | Simulation step number |
| body_id | Body index (0-based) |
| x, y, z | Position at this step (meters) |
