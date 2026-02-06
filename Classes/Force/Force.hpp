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
    void apply( Particles &particles ) const override;
};