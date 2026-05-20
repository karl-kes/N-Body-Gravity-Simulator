#pragma once

#include "../Particle/Particle.hpp"
#include "../Force/Force.hpp"
#include "../Integrator/Integrator.hpp"
#include "../Config.hpp"
#include "../Output/Output.hpp"

#include <iostream>
#include <string>
#include <memory>
#include <vector>
#include <cstddef>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <stdexcept>

#if defined(__GNUC__) || defined(__clang__)
    #define RESTRICT __restrict__
#elif defined(_MSC_VER)
    #define RESTRICT __restrict
#else
    #define RESTRICT
#endif

class Simulation {
private:
    Particles particles_;
    std::vector<std::unique_ptr<Force>> forces_;
    std::unique_ptr<Integrator> integrator_;
    std::size_t num_bodies_;
    std::size_t num_steps_;
    std::size_t output_interval_;
    double dt_;
    double eps_;
    std::vector<std::string> body_names_;
    std::string output_path_;

    double total_energy() const;
    double total_ang_momentum() const;
    double total_lin_momentum() const;
    
    void total_ang_momentum_vec( double &Lx, double &Ly, double &Lz ) const;
    void total_lin_momentum_vec( double &Px, double &Py, double &Pz ) const;

    void print_progress( std::size_t const current, std::size_t const total ) const {
        double const percent{ 100.0 * current / total };

        std::cout << "\rProgress: " << std::fixed << std::setprecision( 0 )
                  << percent << "%" << std::flush;
    }

public:
    explicit Simulation(
        std::size_t const num_particles,
        std::size_t const steps,
        std::size_t const output_interval,
        double const dt,
        double const eps,
        std::vector<std::string> names,
        std::string output_path
    );

    Particles &particles() { return particles_; }
    [[nodiscard]] Particles const &particles() const { return particles_; }

    [[nodiscard]] std::size_t steps() const { return num_steps_; }
    [[nodiscard]] std::size_t output_interval() const { return output_interval_; }
    [[nodiscard]] std::size_t num_bodies() const { return num_bodies_; }
    [[nodiscard]] double dt() const { return dt_; }
    [[nodiscard]] double eps() const { return eps_; }

    [[nodiscard]] std::vector<std::unique_ptr<Force>> &forces() { return forces_; }
    [[nodiscard]] std::unique_ptr<Integrator> &integrator() { return integrator_; }

    void run();
    void add_force( std::unique_ptr<Force> force );
    void set_integrator( std::unique_ptr<Integrator> sim_integrator );
    void initial_output();
};