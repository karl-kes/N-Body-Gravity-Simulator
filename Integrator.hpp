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
private:
    double w_0_, w_1_;
    double c_1_, c_2_;
    double d_1_, d_2_;

public:
    Yoshida( double dt = 1.0 );
    void integrate( Particles &particles, std::vector<std::unique_ptr<Force>> const &forces ) const override;

    [[nodiscard]] double w_0() const { return w_0_; }
    [[nodiscard]] double w_1() const { return w_1_; }
    [[nodiscard]] double c_1() const { return c_1_; }
    [[nodiscard]] double c_2() const { return c_2_; }
    [[nodiscard]] double d_1() const { return d_1_; }
    [[nodiscard]] double d_2() const { return d_2_; }
};