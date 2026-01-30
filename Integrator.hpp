#pragma once

#include "Particle.hpp"
#include "Force.hpp"
#include <vector>
#include <cstddef>
#include <memory>

class Integrator {
protected:
    double dt_;
public:
    Integrator( double dt ) : dt_{ dt } {}
    virtual ~Integrator() = default;
    virtual void integrate( Particles &particles, std::vector<std::unique_ptr<Force>> const &forces ) const = 0;
    [[nodiscard]] double dt() const { return dt_; }
};

class Velocity_Verlet : public Integrator {
public:
    Velocity_Verlet( double dt = 1.0 );

    void integrate( Particles &particles, std::vector<std::unique_ptr<Force>> const &forces ) const override;
};

class Yoshida : public Integrator {
public:
    Yoshida( double dt = 1.0 );

    void integrate( Particles &particles, std::vector<std::unique_ptr<Force>> const &forces ) const override;
};