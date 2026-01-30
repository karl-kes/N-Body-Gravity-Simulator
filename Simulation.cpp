#include "Simulation.hpp"

Simulation::Simulation( std::size_t num_particles, std::size_t steps, std::size_t output_interval )
: particles_{ num_particles }
, forces_{}
, integrator_{ nullptr }
, num_steps_{ steps }
, output_interval_{ output_interval }
{ }

void Simulation::run() {
    double inital_energy{ total_energy() };
    double max_energy{ inital_energy };

    for ( std::size_t curr_step{}; curr_step < steps(); ++curr_step ) {
        integrator()->integrate( particles(), forces() );
        max_energy = std::max( total_energy(), max_energy );
        print_progress( curr_step, steps() );
    }
}

void Simulation::add_force( std::unique_ptr<Force> force ) {
    forces().emplace_back( force );
}

void Simulation::set_integrator( std::unique_ptr<Integrator> sim_integrator ) {
    integrator() = std::move( sim_integrator );
}

double Simulation::total_energy() const {
    std::size_t const N{ particles().num_particles() };
    double energy{};

    for ( std::size_t i{}; i < N; ++i ) {
        auto &a{ particles() };
        double vel_sq{ a.vel_x()[i]*a.vel_x()[i] + a.vel_y()[i]*a.vel_y()[i] + a.vel_z()[i]*a.vel_z()[i]  };

        // Kinetic Energy:
        // K = 1/2 * m * v^2
        energy += 0.5 * a.mass()[i] * vel_sq;

        for ( std::size_t j{i+1}; j < N; ++j ) {
            auto &b{ particles() };

            double dist_x{ b.pos_x()[j] - a.pos_x()[i] };
            double dist_y{ b.pos_y()[j] - a.pos_y()[i] };
            double dist_z{ b.pos_z()[j] - a.pos_z()[i] };

            double const inv_R{ 1.0 / std::sqrt( dist_x*dist_x + dist_y*dist_y + dist_z*dist_z ) };

            // Potential Energy:
            // U = -G * m_1 * m_2 / R
            energy -= G() * a.mass()[i] * b.mass()[j] * inv_R;
        }
    }
    return energy;
}