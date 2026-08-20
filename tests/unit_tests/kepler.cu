// Two-body Kepler and three-body integration: conservation laws over real orbits.

#include "test_harness.hpp"
#include "helpers.hpp"

#include "Force/Force.hpp"
#include "Integrator/Integrator.hpp"
#include "Particle/Particle.hpp"
#include "Config.hpp"

#include <cmath>
#include <memory>
#include <numbers>
#include <vector>


TEST( kepler_circular_orbit_energy_conservation ) {
    double const M{ 1.989e30 };
    double const m{ 5.972e24 };
    double const r{ 1.496e11 };

    Particles p{ 2 };
    setup_two_body( p, M, m, r );

    double const T_orb{ 2.0 * std::numbers::pi * std::sqrt( r*r*r / ( config::G * ( M + m ) ) ) };
    double const dt{ 360.0 };
    std::size_t const steps{ static_cast<std::size_t>( 10.0 * T_orb / dt ) };
    double const E0{ total_energy( p ) };

    std::vector<std::unique_ptr<Force>> forces;
    forces.push_back( std::make_unique<Gravity>() );
    Yoshida integ{ dt };
    step_n( p, integ, forces, steps );

    double const dE_rel{ std::abs( ( total_energy( p ) - E0 ) / E0 ) };
    ASSERT_LT( dE_rel, 1e-10 );
    ++g_pass;
}

TEST( kepler_circular_orbit_returns_to_start ) {
    double const M{ 1.989e30 };
    double const m{ 5.972e24 };
    double const r{ 1.496e11 };

    Particles p{ 2 };
    setup_two_body( p, M, m, r );
    double const x1_init{ p.pos_x()[1] };
    double const y1_init{ p.pos_y()[1] };

    double const T_orb{ 2.0 * std::numbers::pi * std::sqrt( r*r*r / ( config::G * ( M + m ) ) ) };
    double const dt{ 360.0 };
    std::size_t const steps{ static_cast<std::size_t>( T_orb / dt ) };

    std::vector<std::unique_ptr<Force>> forces;
    forces.push_back( std::make_unique<Gravity>() );
    Yoshida integ{ dt };
    step_n( p, integ, forces, steps );

    double const dx{ p.pos_x()[1] - x1_init };
    double const dy{ p.pos_y()[1] - y1_init };
    double const pos_err{ std::sqrt( dx*dx + dy*dy ) };
    ASSERT_LT( pos_err / r, 1e-4 );
    ++g_pass;
}

TEST( kepler_angular_momentum_conservation ) {
    double const M{ 1.989e30 };
    double const m{ 5.972e24 };
    double const r{ 1.496e11 };

    Particles p{ 2 };
    setup_two_body( p, M, m, r );
    double const L0{ angular_momentum( p ) };

    double const T_orb{ 2.0 * std::numbers::pi * std::sqrt( r*r*r / ( config::G * ( M + m ) ) ) };
    double const dt{ 360.0 };
    std::size_t const steps{ static_cast<std::size_t>( 5.0 * T_orb / dt ) };

    std::vector<std::unique_ptr<Force>> forces;
    forces.push_back( std::make_unique<Gravity>() );
    Yoshida integ{ dt };
    step_n( p, integ, forces, steps );

    double const dL_rel{ std::abs( ( angular_momentum( p ) - L0 ) / L0 ) };
    ASSERT_LT( dL_rel, 1e-12 );
    ++g_pass;
}

TEST( kepler_linear_momentum_conservation ) {
    double const M{ 1.989e30 };
    double const m{ 5.972e24 };
    double const r{ 1.496e11 };

    Particles p{ 2 };
    setup_two_body( p, M, m, r );

    double const T_orb{ 2.0 * std::numbers::pi * std::sqrt( r*r*r / ( config::G * ( M + m ) ) ) };
    double const dt{ 360.0 };
    std::size_t const steps{ static_cast<std::size_t>( 5.0 * T_orb / dt ) };

    std::vector<std::unique_ptr<Force>> forces;
    forces.push_back( std::make_unique<Gravity>() );
    Yoshida integ{ dt };
    step_n( p, integ, forces, steps );

    double const P_final{ linear_momentum( p ) };
    double const body_momentum_scale{ m * std::sqrt( config::G * ( M + m ) / r ) };
    ASSERT_LT( P_final / body_momentum_scale, 1e-10 );
    ++g_pass;
}

TEST( three_body_energy_conservation ) {
    Particles p{ 3 };

    p.pos_x()[0] = 0.0; p.pos_y()[0] = 0.0; p.pos_z()[0] = 0.0;
    p.vel_x()[0] = 0.0; p.vel_y()[0] = 0.0; p.vel_z()[0] = 0.0;
    p.mass()[0]  = 1.989e30;

    double const rJ{ 5.2 * config::AU };
    double const vJ{ std::sqrt( config::G * p.mass()[0] / rJ ) };
    p.pos_x()[1] = rJ;  p.pos_y()[1] = 0.0; p.pos_z()[1] = 0.0;
    p.vel_x()[1] = 0.0; p.vel_y()[1] = vJ;  p.vel_z()[1] = 0.0;
    p.mass()[1]  = 1.898e27;

    double const rS{ 9.5 * config::AU };
    double const vS{ std::sqrt( config::G * p.mass()[0] / rS ) };
    p.pos_x()[2] = 0.0;  p.pos_y()[2] = rS;  p.pos_z()[2] = 0.0;
    p.vel_x()[2] = -vS;  p.vel_y()[2] = 0.0; p.vel_z()[2] = 0.0;
    p.mass()[2]  = 5.683e26;

    for ( std::size_t i{}; i < 3; ++i ) {
        p.acc_x()[i] = 0.0; p.acc_y()[i] = 0.0; p.acc_z()[i] = 0.0;
    }

    double const E0{ total_energy( p ) };
    double const T_jup{ 2.0 * std::numbers::pi * std::sqrt( rJ*rJ*rJ / ( config::G * p.mass()[0] ) ) };
    double const dt{ 900.0 };
    std::size_t const steps{ static_cast<std::size_t>( T_jup / dt ) };

    std::vector<std::unique_ptr<Force>> forces;
    forces.push_back( std::make_unique<Gravity>() );
    Yoshida integ{ dt };
    step_n( p, integ, forces, steps );

    double const dE_rel{ std::abs( ( total_energy( p ) - E0 ) / E0 ) };
    ASSERT_LT( dE_rel, 1e-10 );
    ++g_pass;
}
