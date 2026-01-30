#pragma once

#include "Particle.hpp"
#include <vector>
#include <cmath>
#include <cstddef>

class Force {
public:
    virtual ~Force() = default;
    virtual void apply( Particles &particles ) const = 0;
};

class Gravity : public Force {
private:
    // Gravitational Constant:
    static constexpr double G_{ 6.6743e-11 };

public:
    Gravity();

    void apply( Particles &particles ) const override;
    [[nodiscard]] double G() const { return G_; }
};