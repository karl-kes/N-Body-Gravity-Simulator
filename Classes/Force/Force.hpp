#pragma once

#include "../Particle/Particle.hpp"
#include "../../Constants.hpp"

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
    inline void compute_forces(
        double pxi, double pyi, double pzi, double mi,
        double vxi, double vyi, double vzi,
        double pxj, double pyj, double pzj, double mj,
        double vxj, double vyj, double vzj,
        double& a_xi, double& a_yi, double& a_zi,
        double G, double eps_sq, double c_sq 
    ) const;
    void apply( Particles &particles ) const override;
};