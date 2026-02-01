#pragma once

#include "../Particle/Particle.hpp"
#include "../Force/Force.hpp"

#include <vector>
#include <cstddef>
#include <memory>
#include <string>
#include <omp.h>

class Integrator {
protected:
    double dt_;
    std::string name_;
public:
    Integrator( double dt, std::string const &name ) 
    : dt_{ dt }
    , name_{name }
    { }
    virtual ~Integrator() = default;
    virtual void integrate( Particles &particles, std::vector<std::unique_ptr<Force>> const &forces ) const = 0;
    [[nodiscard]] double dt() const { return dt_; }
    [[nodiscard]] const std::string &name() const { return name_; }
};

class Velocity_Verlet : public Integrator {
public:
    Velocity_Verlet( double dt = 1.0 );
    void integrate( Particles &particles, std::vector<std::unique_ptr<Force>> const &forces ) const override;
};

class Yoshida : public Integrator {
private:
    static constexpr double cbrt_2_{ 1.2599210498948732 };
    double w_0_, w_1_;
    double c_1_, c_2_, c_3_, c_4_;
    double d_1_, d_2_, d_3_;

public:
    Yoshida( double dt = 1.0 );
    void integrate( Particles &particles, std::vector<std::unique_ptr<Force>> const &forces ) const override;

    [[nodiscard]] double cbrt_2() const { return cbrt_2_; }

    [[nodiscard]] double w_0() const { return w_0_; }
    [[nodiscard]] double w_1() const { return w_1_; }

    [[nodiscard]] double c_1() const { return c_1_; }
    [[nodiscard]] double c_2() const { return c_2_; }
    [[nodiscard]] double c_3() const { return c_3_; }
    [[nodiscard]] double c_4() const { return c_4_; }

    [[nodiscard]] double d_1() const { return d_1_; }
    [[nodiscard]] double d_2() const { return d_2_; }
    [[nodiscard]] double d_3() const { return d_3_; }
};