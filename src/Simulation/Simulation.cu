#include "Simulation.hpp"

Simulation::Simulation(
    std::size_t const num_particles,
    std::size_t const steps,
    std::size_t const output_interval,
    double const dt,
    double const eps,
    std::vector<std::string> names,
    std::string output_path )
: particles_{ num_particles }
, forces_{}
, integrator_{ nullptr }
, num_bodies_{ num_particles }
, num_steps_{ steps }
, output_interval_{ output_interval }
, dt_{ dt }
, eps_{ eps }
, body_names_{ std::move( names ) }
, output_path_{ std::move( output_path ) }
{ }

void Simulation::run() {
    if ( !integrator_ ) {
        throw std::runtime_error( "No integrator set. Call set_integrator() before run()." );
    }
    if ( forces_.empty() ) {
        throw std::runtime_error( "No forces added. Call add_force() before run()." );
    }

    // Compute initial accelerations so that integrators which read acceleration
    // on their first step (e.g. Velocity Verlet) start from a valid state.
    // Yoshida does not need this (it computes forces internally), but it is
    // harmless and makes the API contract explicit.
    for ( std::size_t i{}; i < num_bodies(); ++i ) {
        particles().acc_x()[i] = 0.0;
        particles().acc_y()[i] = 0.0;
        particles().acc_z()[i] = 0.0;
    }
    for ( auto const &force : forces() ) {
        force->apply( particles() );
    }

    double const initial_energy{ total_energy() };
    double min_energy{ initial_energy };
    double max_energy{ initial_energy };

    double L0_x{}, L0_y{}, L0_z{};
    total_ang_momentum_vec( L0_x, L0_y, L0_z );
    double const initial_ang_momentum{ std::sqrt( L0_x*L0_x + L0_y*L0_y + L0_z*L0_z ) };
    double max_ang_momentum_vec_drift{};

    double P0_x{}, P0_y{}, P0_z{};
    total_lin_momentum_vec( P0_x, P0_y, P0_z );
    double const initial_lin_momentum{ std::sqrt( P0_x*P0_x + P0_y*P0_y + P0_z*P0_z ) };
    double max_lin_momentum_vec_drift{};

    initial_output();

    Binary_Output bin{ output_path_, body_names_, num_bodies() };
    bin.write( particles(), 0, 0.0 );

    auto const start_time{ std::chrono::high_resolution_clock::now() };

    for ( std::size_t curr_step{1}; curr_step <= steps(); ++curr_step ) {
        integrator()->integrate( particles(), forces() );

        // Sample conservation quantities at 10x the output cadence.
        if ( curr_step % ( 10*output_interval() ) == 0 ) {
            double const E{ total_energy() };
            max_energy = std::max( E, max_energy );
            min_energy = std::min( E, min_energy );

            double Lx{}, Ly{}, Lz{};
            total_ang_momentum_vec( Lx, Ly, Lz );
            double const dLx{ Lx - L0_x }, dLy{ Ly - L0_y }, dLz{ Lz - L0_z };
            double const L_vec_drift{ std::sqrt( dLx*dLx + dLy*dLy + dLz*dLz ) };
            max_ang_momentum_vec_drift = std::max( max_ang_momentum_vec_drift, L_vec_drift );

            double Px{}, Py{}, Pz{};
            total_lin_momentum_vec( Px, Py, Pz );
            double const dPx{ Px - P0_x }, dPy{ Py - P0_y }, dPz{ Pz - P0_z };
            double const P_vec_drift{ std::sqrt( dPx*dPx + dPy*dPy + dPz*dPz ) };
            max_lin_momentum_vec_drift = std::max( max_lin_momentum_vec_drift, P_vec_drift );

            print_progress( curr_step, steps() );
        }

        if ( curr_step % output_interval() == 0 ) {
            bin.write( particles(), curr_step, curr_step * dt_ );
        }
    }
    std::cout << "\rProgress: 100%" << std::flush;

    auto const end_time{ std::chrono::high_resolution_clock::now() };
    auto const duration{ std::chrono::duration_cast<std::chrono::milliseconds>( end_time - start_time ) };

    double const energy_drift{ std::abs( 100.0 * ( max_energy - min_energy ) / initial_energy ) };

    double const ang_momentum_drift{ initial_ang_momentum > 0.0
        ? 100.0 * max_ang_momentum_vec_drift / initial_ang_momentum : 0.0 };
    double const lin_momentum_drift{ initial_lin_momentum > 0.0
        ? 100.0 * max_lin_momentum_vec_drift / initial_lin_momentum : 0.0 };

    std::cout << "\nMax Energy Drift: " << std::scientific << std::setprecision( 6 ) << energy_drift << "%" << std::endl;
    std::cout << "Max Angular Momentum Drift: " << std::scientific << std::setprecision( 6 ) << ang_momentum_drift << "%"
              << "  (max ||L(t)-L(0)|| = " << max_ang_momentum_vec_drift << ")" << std::endl;
    std::cout << "Max Linear Momentum Drift: " << std::scientific << std::setprecision( 6 ) << lin_momentum_drift << "%"
              << "  (max ||P(t)-P(0)|| = " << max_lin_momentum_vec_drift << ")" << std::endl;
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

    double const* RESTRICT px{ particles().pos_x() };
    double const* RESTRICT py{ particles().pos_y() };
    double const* RESTRICT pz{ particles().pos_z() };
    double const* RESTRICT vx{ particles().vel_x() };
    double const* RESTRICT vy{ particles().vel_y() };
    double const* RESTRICT vz{ particles().vel_z() };
    double const* RESTRICT mass{ particles().mass() };

    double const eps_sq{ eps_ * eps_ };
    constexpr double G{ config::G };
    constexpr std::size_t OMP_THRESHOLD{ config::OMP_THRESHOLD };

    double kinetic_energy{};
    auto kinetic_kernel = [vx, vy, vz, mass]( std::size_t i ) -> double {
        double const v_sq{ vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i] };
        return 0.5 * mass[i] * v_sq;
    };

    double potential_energy{};
    auto potential_kernel = [px, py, pz, mass, N, eps_sq]( std::size_t i ) -> double {
        double const pxi{ px[i] }, pyi{ py[i] }, pzi{ pz[i] };
        double const mi{ mass[i] };
        double row_pot{ 0.0 };

        #pragma omp simd reduction( +:row_pot )
        for ( std::size_t j = i + 1; j < N; ++j ) {
            double const dx{ px[j] - pxi };
            double const dy{ py[j] - pyi };
            double const dz{ pz[j] - pzi };

            double const R_sq{ dx*dx + dy*dy + dz*dz + eps_sq };
            double const inv_R{ 1.0 / std::sqrt( R_sq ) };

            row_pot -= mass[j] * inv_R;
        }

        return G * mi * row_pot;
    };

    if ( N >= OMP_THRESHOLD ) {
        #pragma omp parallel for simd reduction( +:kinetic_energy ) schedule( static )
        for ( std::size_t i = 0; i < N; ++i ) {
            kinetic_energy += kinetic_kernel( i );
        }
    } else {
        #pragma omp simd reduction( +:kinetic_energy )
        for ( std::size_t i = 0; i < N; ++i ) {
            kinetic_energy += kinetic_kernel( i );
        }
    }

    if ( N >= OMP_THRESHOLD ) {
        #pragma omp parallel for reduction( +:potential_energy ) schedule( guided )
        for ( std::size_t i = 0; i < N; ++i ) {
            potential_energy += potential_kernel( i );
        }
    } else {
        for ( std::size_t i = 0; i < N; ++i ) {
            potential_energy += potential_kernel( i );
        }
    }

    return kinetic_energy + potential_energy;
}

void Simulation::total_ang_momentum_vec( double &Lx, double &Ly, double &Lz ) const {
    std::size_t const N{ particles().num_particles() };
    double const* RESTRICT m{ particles().mass() };

    double const* RESTRICT px{ particles().pos_x() };
    double const* RESTRICT py{ particles().pos_y() };
    double const* RESTRICT pz{ particles().pos_z() };

    double const* RESTRICT vx{ particles().vel_x() };
    double const* RESTRICT vy{ particles().vel_y() };
    double const* RESTRICT vz{ particles().vel_z() };

    Lx = 0.0; Ly = 0.0; Lz = 0.0;
    #pragma omp simd reduction( +:Lx, Ly, Lz )
    for ( std::size_t i = 0; i < N; ++i ) {
        Lx += m[i] * ( py[i]*vz[i] - pz[i]*vy[i] );
        Ly += m[i] * ( pz[i]*vx[i] - px[i]*vz[i] );
        Lz += m[i] * ( px[i]*vy[i] - py[i]*vx[i] );
    }
}

void Simulation::total_lin_momentum_vec( double &Px, double &Py, double &Pz ) const {
    std::size_t const N{ particles().num_particles() };

    double const* RESTRICT m{ particles().mass() };

    double const* RESTRICT vx{ particles().vel_x() };
    double const* RESTRICT vy{ particles().vel_y() };
    double const* RESTRICT vz{ particles().vel_z() };

    Px = 0.0; Py = 0.0; Pz = 0.0;
    #pragma omp simd reduction( +:Px, Py, Pz )
    for ( std::size_t i = 0; i < N; ++i ) {
        Px += m[i] * vx[i];
        Py += m[i] * vy[i];
        Pz += m[i] * vz[i];
    }
}

double Simulation::total_ang_momentum() const {
    double Lx{}, Ly{}, Lz{};
    total_ang_momentum_vec( Lx, Ly, Lz );
    return std::sqrt( Lx*Lx + Ly*Ly + Lz*Lz );
}

double Simulation::total_lin_momentum() const {
    double Px{}, Py{}, Pz{};
    total_lin_momentum_vec( Px, Py, Pz );
    return std::sqrt( Px*Px + Py*Py + Pz*Pz );
}

void Simulation::initial_output() {
    std::cout << "\n<--- N-Body Simulation --->" << std::endl;
    std::cout << "Bodies: " << num_bodies() << std::endl;
    std::cout << "Integrator: " << integrator()->name() << std::endl;
    std::cout << "Dt: " << dt_ << " seconds" << std::endl;
    std::cout << "Steps: " << num_steps_ << std::endl;
    std::cout << "Eps: " << eps_ << " m" << std::endl;
    std::cout << "Parallelization: " << ( ( config::OMP_THRESHOLD <= num_bodies() ) ? "Enabled" : "Disabled" ) << std::endl;
    std::cout << std::endl;
}