#include "Simulation.hpp"

#define RESTRICT __restrict

Simulation::Simulation( std::size_t const num_particles, std::size_t const steps, std::size_t const output_interval )
: particles_{ num_particles }
, forces_{}
, integrator_{ nullptr }
, num_bodies_{ num_particles }
, num_steps_{ steps }
, output_interval_{ output_interval }
{ }

void Simulation::run() {
    double const initial_energy{ total_energy() };
    double max_energy{ initial_energy };
    double min_energy{ initial_energy };

    auto const start_time{ std::chrono::high_resolution_clock::now() };

    for ( std::size_t curr_step{}; curr_step < steps(); ++curr_step ) {
        integrator()->integrate( particles(), forces() );

        if ( curr_step % ( constant::steps_per_year / 12 ) == 0 ) {
            double const E{ total_energy() };
            max_energy = std::max( E, max_energy );
            min_energy = std::min( E, min_energy );
        }
        if ( curr_step % output_interval() == 0 ) { print_progress( curr_step, steps() ); }
    }
    std::cout << "\rProgress: 100%" << std::flush;

    auto const end_time{ std::chrono::high_resolution_clock::now() };
    auto const duration{ std::chrono::duration_cast<std::chrono::milliseconds>( end_time - start_time ) };

    double const drift{ std::abs( 100.0 * ( max_energy - min_energy ) / initial_energy ) };
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

    double* RESTRICT px{ particles().pos_x().get() };
    double* RESTRICT py{ particles().pos_y().get() };
    double* RESTRICT pz{ particles().pos_z().get() };

    double* RESTRICT vx{ particles().vel_x().get() };
    double* RESTRICT vy{ particles().vel_y().get() };
    double* RESTRICT vz{ particles().vel_z().get() };

    double const* RESTRICT mass{ particles().mass().get() };

    constexpr double eps_sq{ constant::EPS*constant::EPS };
    constexpr double G{ constant::G };
    constexpr double OMP_THRESHOLD{ constant::OMP_THRESHOLD };

    double kinetic_energy{};
    #pragma omp parallel for reduction( +:kinetic_energy ) schedule( static )
    for ( std::size_t i = 0; i < N; ++i ) {
        double const vel_sq{ vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]  };

        kinetic_energy += 0.5 * mass[i] * vel_sq;
    }

    double potential_energy{};
    #pragma omp parallel for reduction( +:potential_energy ) schedule( guided )
    for ( std::size_t i = 0; i < N; ++i ) {
        double const pxi{ px[i] }, pyi{ py[i] }, pzi{ pz[i] };
        double const mi{ mass[i] };
        double U_i{};

        #pragma omp simd reduction( +:U_i )
        for ( std::size_t j = i + 1; j < N; ++j ) {
            double const dx{ px[j] - pxi };
            double const dy{ py[j] - pyi };
            double const dz{ pz[j] - pzi };

            double const dist_sq{ dx*dx + dy*dy + dz*dz + eps_sq };
            double const inv_R{ 1.0 / std::sqrt( dist_sq ) };

            U_i += mass[j] * inv_R;
        }

        potential_energy -= G * mi * U_i;
    }
    
    return kinetic_energy + potential_energy;
}

void Simulation::final_output( Body const *bodies ) const {
    std::cout << "\nFinal distances from Sun:" << std::endl;
    for ( std::size_t i{ 1 }; i < num_bodies(); ++i ) {
        double const dist_x{ particles().pos_x(i) - particles().pos_x(0) };
        double const dist_y{ particles().pos_y(i) - particles().pos_y(0) };
        double const dist_z{ particles().pos_z(i) - particles().pos_z(0) };

        double const R{ std::sqrt( dist_x*dist_x + dist_y*dist_y + dist_z*dist_z ) };

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
    std::cout << "Post-Newtonian: " << ( constant::ENABLE_PN ? "Enabled" : "Disabled" ) << std::endl;
    std::cout << std::endl;
}

void Simulation::output_positions( Body const *bodies, std::size_t const curr_time ) const {
    std::size_t const size{ sizeof(bodies) / sizeof(bodies[0]) };


}