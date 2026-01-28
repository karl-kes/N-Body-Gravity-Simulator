#pragma once

#include "Particle.hpp"
#include <vector>
#include <cmath>
#include <omp.h>

class Force_Law {
public:
    virtual ~Force_Law() = default;
    virtual void apply( std::vector<Particle> &particles ) const = 0;
};

class Gravity : public Force_Law {
private:
    double G_;

public:
    Gravity();

    void apply( std::vector<Particle> &particles ) const override;

    [[nodiscard]] double G() const { return G_; }
};