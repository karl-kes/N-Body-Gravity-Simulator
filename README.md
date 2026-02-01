# N-Body Gravitational Simulation

A high-performance N-body gravitational simulation written in C++17, featuring symplectic integrators for long-term orbital stability and energy conservation.

## Overview

This project simulates gravitational interactions between celestial bodies using Newtonian mechanics. It implements both Velocity Verlet and Yoshida 4th-order symplectic integrators, achieving energy drift below 10^-12% over century-scale simulations. The default configuration simulates the full solar system (Sun + 8 planets + Pluto) using JPL Horizons initial conditions.

## Features

- **Symplectic Integrators**: 
  - Velocity Verlet (2nd order)
  - Yoshida (4th order)
- **Structure-of-Arrays (SoA) Layout**: Cache-friendly memory layout for better performance
- **Energy Conservation Tracking**: Monitors total energy drift to validate simulation accuracy
- **Extensible Force Framework**: Polymorphic design supports adding new force types (Lennard-Jones, etc.)
- **OpenMP Support**: Optional parallelization for large-scale simulations (WIP)

## Project Structure

```
├── main.cpp                 # Entry point with solar system initial conditions
├── Constants.hpp            # Physical constants and simulation parameters
└── Classes/
    ├── Particle/
    │   ├── Particle.hpp
    │   └── Particle.cpp
    ├── Force/
    │   ├── Force.hpp        # Abstract force base class
    │   └── Force.cpp        # Force implementations
    ├── Integrator/
    │   ├── Integrator.hpp   # Integrator base class
    │   └── Integrator.cpp   # Velocity Verlet & Yoshida implementations
    └── Simulation/
        ├── Simulation.hpp   # Simulation manager
        └── Simulation.cpp   # Core loop and energy calculation
```

## Getting Started

### Prerequisites

- C++ compiler with C++17 support (g++, clang++)

### Compilation

```bash
g++ -std=c++17 -O3 -march=native *.cpp Classes/Force/*.cpp Classes/Integrator/*.cpp Classes/Particle/*.cpp Classes/Simulation/*.cpp -o main.exe
```

### Running

```bash
./main.exe
```

## Configuration

Simulation parameters in `Constants.hpp`:

| Constant | Default | Description |
|----------|---------|-------------|
| `G` | 6.6743e-11 | Gravitational constant (m³/kg/s²) |
| `EPS` | 1e-9 | Softening parameter (m) |
| `AU` | 1.496e11 | Astronomical unit (m) |
| `dt` | 3600 | Timestep (seconds) |
| `num_years` | 10 | Simulation duration |

## Initial Conditions

The default configuration uses approximate JPL Horizons state vectors for January 2025:

| Body | Mass (kg) | Position (km) | Velocity (km/s) |
|------|-----------|---------------|-----------------|
| Sun | 1.989e30 | Origin | Stationary |
| Mercury | 3.302e23 | (-1.478e7, -6.553e7, -3.898e6) | (36.2, -9.0, -4.0) |
| Venus | 4.869e24 | (-5.765e7, -9.361e7, 2.110e6) | (29.9, -18.5, -2.1) |
| Earth | 5.972e24 | (-2.627e7, 1.445e8, -1.049e4) | (-29.8, -5.4, 0.0) |
| Mars | 6.417e23 | (2.067e8, 4.500e7, -4.057e6) | (-3.9, 26.0, 0.65) |
| Jupiter | 1.898e27 | (5.765e8, 4.405e8, -1.493e7) | (-7.9, 10.7, 0.14) |
| Saturn | 5.683e26 | (1.357e9, -5.194e8, -4.480e7) | (2.9, 9.0, -0.26) |
| Uranus | 8.681e25 | (1.855e9, 2.233e9, -1.579e7) | (-5.2, 4.0, 0.08) |
| Neptune | 1.024e26 | (4.461e9, -2.705e8, -9.775e7) | (0.29, 5.5, -0.12) |
| Pluto | 1.303e22 | (2.595e9, -4.513e9, -2.816e8) | (4.8, 1.5, -1.6) |

## Physics

### Gravitational Acceleration

$$\vec{a}_i = \sum_{j \neq i} \frac{G \cdot m_j}{(|\vec{r}_{ij}|^2 + \epsilon^2)^{3/2}} \vec{r}_{ij}$$

### Yoshida 4th-Order Integrator

Coefficients derived from $w_1 = \frac{1}{2 - \sqrt[3]{2}}$, $w_0 = -\frac{\sqrt[3]{2}}{2 - \sqrt[3]{2}}$:

```
c₁ = c₄ = w₁/2 ≈ 0.6756
c₂ = c₃ = (w₀ + w₁)/2 ≈ -0.1756
d₁ = d₃ = w₁ ≈ 1.3512
d₂ = w₀ ≈ -1.7024
```

### Velocity Verlet Integration

```
x(t+dt) = x(t) + v(t)·dt + 0.5·a(t)·dt²
v(t+dt) = v(t) + 0.5·(a(t) + a(t+dt))·dt
```

## Performance

| Integrator | Timestep | 10-Year Energy Drift | Steps |
|------------|----------|----------------------|-------|
| Velocity Verlet | 1 hour | ~9e-5% | 87,660 |
| Yoshida | 1 hour | ~5e-12% | 87,660 |

## Sample Output

```
<--- Solar System Simulation --->
Bodies: 10
Duration: 10 years

Running simulation...
Progress: 100.0%
Max Energy Drift: 4.586857e-12%

Final distances from Sun:
Mercury   0.3051 AU
Venus     0.7528 AU
Earth     0.9829 AU
Mars      1.4655 AU
Jupiter   4.8432 AU
Saturn    8.9645 AU
Uranus    18.5853 AU
Neptune   29.8446 AU
Pluto     37.3802 AU
```

## Extending the Framework

### Adding a New Force

```cpp
class LennardJones : public Force {
    double epsilon_, sigma_;
public:
    LennardJones( double eps, double sig ) : epsilon_{ eps }, sigma_{ sig } {}
    
    void apply( Particles &particles ) const override {
        // Implementation
    }
};

// Usage:
sim.add_force( std::make_unique<LennardJones>( 1.0, 1.0 ) );
```

### Adding a New Integrator

```cpp
class RK4 : public Integrator {
public:
    RK4( double dt ) : Integrator{ dt } {}
    
    void integrate( Particles &p, std::vector<std::unique_ptr<Force>> const &f ) const override {
        // Implementation
    }
};
```

## References

- [JPL Horizons System](https://ssd.jpl.nasa.gov/horizons/) — Ephemeris data
- [JPL Planetary Fact Sheet](https://nssdc.gsfc.nasa.gov/planetary/factsheet/) — Physical parameters
- Yoshida, H. (1990). "Construction of higher order symplectic integrators". Physics Letters A.