#pragma once

#include "../Particle/Particle.hpp"
#include "../../Constants.hpp"

#define INLINE __forceinline

#include <vector>
#include <cmath>
#include <cstddef>
#include <omp.h>

class Force {
public:
    virtual ~Force() = default;
    virtual void apply( Particles &particles ) const = 0;
};

class Gravity : public Force {
public:
    Gravity();
    static void INLINE compute_forces(
        double const pxi, double const pyi, double const pzi, double const mi,
        double const vxi, double const vyi, double const vzi,
        double const pxj, double const pyj, double const pzj, double const mj,
        double const vxj, double const vyj, double const vzj,
        double &a_xi, double &a_yi, double &a_zi,
        double const G, double const eps_sq, double const c_sq, double const mask
    );

    void apply( Particles &particles ) const override;
};