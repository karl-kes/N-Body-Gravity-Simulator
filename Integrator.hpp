#pragma once

#include "Particle.hpp"
#include "Force.hpp"
#include <vector>
#include <cstddef>
#include <memory>

class Integrator {
public:
    virtual ~Integrator() = default;
    virtual void integrate( std::vector<Particle> &particles, std::vector<std::unique_ptr<Force_Law>> const &forces ) const = 0;
};

class Velocity_Verlet : public Integrator {
private:
    double dt_;
public:
    Velocity_Verlet( double dt = 1.0 );

    void integrate( std::vector<Particle> &particles, std::vector<std::unique_ptr<Force_Law>> const &forces ) const override;
    [[nodiscard]] double dt() const { return dt_; }
};