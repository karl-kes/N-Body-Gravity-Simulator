#include "Simulation.hpp"

Simulation::Simulation( std::size_t num_particles, std::size_t steps, std::size_t output_interval )
: particles_{ num_particles }
, forces_{}
, integrator_{ nullptr }
, num_bodies_{ num_particles }
, num_steps_{ steps }
, output_interval_{ output_interval }
{ }

void Simulation::run() {
    double inital_energy{ total_energy() };
    double max_energy{ inital_energy };

    auto start_time{ std::chrono::high_resolution_clock::now() };

    for ( std::size_t curr_step{}; curr_step < steps(); ++curr_step ) {
        integrator()->integrate( particles(), forces() );
        max_energy = std::max( total_energy(), max_energy );

        if ( curr_step % output_interval() == 0 ) { print_progress( curr_step, steps() ); }
    }
    std::cout << "\rProgress: 100%" << std::flush;

    auto end_time{ std::chrono::high_resolution_clock::now() };
    auto duration{ std::chrono::duration_cast<std::chrono::milliseconds>( end_time - start_time ) };

    double drift{ std::abs( 100.0 * ( max_energy - inital_energy ) / inital_energy ) };
    std::cout << "\nMax Energy Drift: " << std::scientific << std::setprecision( 6 ) << drift << "%" << std::endl;
    std::cout << "Duration of Simulation: " << duration.count() << " ms" << std::endl;
}

void Simulation::add_force( std::unique_ptr<Force> force ) {
    forces().emplace_back( std::move( force ) );
}

void Simulation::set_integrator( std::unique_ptr<Integrator> sim_integrator ) {
    integrator() = std::move( sim_integrator );
}

double Simulation::total_energy() const {
    std::size_t const N{ particles().num_particles() };
    double energy{};

    #pragma omp parallel for reduction( +:energy ) schedule( static )
    for ( std::size_t i = 0; i < N; ++i ) {
        double const vel_sq{ particles().vel_x()[i]*particles().vel_x()[i] +
                             particles().vel_y()[i]*particles().vel_y()[i] +
                             particles().vel_z()[i]*particles().vel_z()[i]  };

        // Kinetic Energy:
        // K = 1/2 * m * v^2
        energy += 0.5 * particles().mass()[i] * vel_sq;

        for ( std::size_t j = i + 1; j < N; ++j ) {
            double const dist_x{ particles().pos_x()[j] - particles().pos_x()[i] };
            double const dist_y{ particles().pos_y()[j] - particles().pos_y()[i] };
            double const dist_z{ particles().pos_z()[j] - particles().pos_z()[i] };

            double const dist_sq{ dist_x*dist_x + dist_y*dist_y + dist_z*dist_z + constant::EPS*constant::EPS };

            double const inv_R{ 1.0 / std::sqrt( dist_sq ) };

            // Potential Energy:
            // U = -G * m_1 * m_2 / R
            energy -= constant::G * particles().mass()[i] * particles().mass()[j] * inv_R;
        }
    }
    
    return energy;
}

void Simulation::final_output( Body const *bodies ) const {
    std::cout << "\nFinal distances from Sun:" << std::endl;
    for ( std::size_t i{ 1 }; i < num_bodies(); ++i ) {
        double const dist_x{ particles().pos_x(i) - particles().pos_x(0) };
        double const dist_y{ particles().pos_y(i) - particles().pos_y(0) };
        double const dist_z{ particles().pos_z(i) - particles().pos_z(0) };

        double R{ std::sqrt( dist_x*dist_x + dist_y*dist_y + dist_z*dist_z ) };

        std::cout << std::left << std::setw( 10 ) << bodies[i].name
                  << std::fixed << std::setprecision( 4 )
                  << R / constant::AU << " AU" << std::endl;
    }
}

void Simulation::initial_output() {
    std::cout << "<--- Solar System Simulation --->" << std::endl;
    std::cout << "Bodies: " << num_bodies() << std::endl;
    std::cout << "Integrator: " << integrator()->name() << std::endl;
    std::cout << "Duration: " << constant::num_years << " years" << std::endl;
    std::cout << std::endl;
}