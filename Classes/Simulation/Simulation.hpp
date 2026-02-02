#pragma once

#include "../Particle/Particle.hpp"
#include "../Force/Force.hpp"
#include "../Integrator/Integrator.hpp"
#include "../../Constants.hpp"
#include "../../Body.hpp"

#include <memory>
#include <vector>
#include <cstddef>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <omp.h>

class Simulation {
private:
    Particles particles_;
    std::vector<std::unique_ptr<Force>> forces_;
    std::unique_ptr<Integrator> integrator_;
    std::size_t num_bodies_;
    std::size_t num_steps_;
    std::size_t output_interval_;

    double total_energy() const;

    void print_progress( std::size_t current, std::size_t total ) const {
        double percent{ 100.0 * current / total };

        std::cout << "\rProgress: " << std::fixed << std::setprecision(0) 
                  << percent << "%" << std::flush;
    }

public:
    // Constructor, Destructor, Move, Copy:
    explicit Simulation( std::size_t num_particles, std::size_t steps, std::size_t output_interval );
    ~Simulation() = default;
    Simulation( const Simulation& ) = delete;
    Simulation& operator=( const Simulation& ) = delete;
    Simulation( Simulation&& ) = delete;
    Simulation& operator=( Simulation&& ) = delete;

    // Getters & Setters:
    // Particles:
    Particles &particles() { return particles_; }
    [[nodiscard]] const Particles &particles() const { return particles_; }

    [[nodiscard]] std::size_t steps() const { return num_steps_; }
    [[nodiscard]] std::size_t output_interval() const { return output_interval_; }
    [[nodiscard]] std::size_t num_bodies() const { return num_bodies_; }

    std::vector<std::unique_ptr<Force>> &forces() { return forces_; }
    std::unique_ptr<Integrator> &integrator() { return integrator_; }

    // Simulation Functions:
    void run();
    void add_force( std::unique_ptr<Force> force );
    void set_integrator( std::unique_ptr<Integrator> sim_integrator );
    void final_output( Body const *bodies ) const;
    void initial_output();
};