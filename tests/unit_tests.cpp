// tests/unit_tests.cpp
//
// Self-contained unit tests for the N-body engine.
// No external test framework required - builds standalone.
//
// Build:  cmake --build build --target tests
// Run:    ./build/tests

#include "../src/Particle/Particle.hpp"
#include "../src/Force/Force.hpp"
#include "../src/Force/BarnesHut.hpp"
#include "../src/Integrator/Integrator.hpp"
#include "../src/Config.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <memory>
#include <functional>
#include <numbers>
#include <random>

// Minimal test harness

static int g_pass{ 0 };
static int g_fail{ 0 };

#define TEST( name ) static void name(); \
    static bool name##_reg = ( tests().push_back( { #name, name } ), true ); \
    static void name()

#define ASSERT_NEAR( actual, expected, tol ) \
    do { \
        double const _a{ (actual) }; \
        double const _e{ (expected) }; \
        double const _t{ (tol) }; \
        if ( std::abs( _a - _e ) > _t ) { \
            std::cerr << "    FAIL: " << #actual << " = " << std::scientific \
                      << std::setprecision( 10 ) << _a << ", expected " << _e \
                      << " +/- " << _t << "\n"; \
            ++g_fail; return; \
        } \
    } while ( 0 )

#define ASSERT_LT( actual, bound ) \
    do { \
        double const _a{ (actual) }; \
        double const _b{ (bound) }; \
        if ( _a >= _b ) { \
            std::cerr << "    FAIL: " << #actual << " = " << std::scientific \
                      << std::setprecision( 10 ) << _a << " >= " << _b << "\n"; \
            ++g_fail; return; \
        } \
    } while ( 0 )

#define ASSERT_TRUE( cond ) \
    do { \
        if ( !(cond) ) { \
            std::cerr << "    FAIL: " << #cond << " is false\n"; \
            ++g_fail; return; \
        } \
    } while ( 0 )

struct TestEntry {
    std::string name;
    std::function<void()> fn;
};

static std::vector<TestEntry> &tests() {
    static std::vector<TestEntry> t;
    return t;
}

// Helpers

static double kinetic_energy( Particles const &p ) {
    double ke{};
    for ( std::size_t i{}; i < p.num_particles(); ++i ) {
        double const v_sq{ p.vel_x()[i]*p.vel_x()[i]
                         + p.vel_y()[i]*p.vel_y()[i]
                         + p.vel_z()[i]*p.vel_z()[i] };
        ke += 0.5 * p.mass()[i] * v_sq;
    }
    return ke;
}

static double potential_energy( Particles const &p ) {
    double pe{};
    std::size_t const N{ p.num_particles() };
    for ( std::size_t i{}; i < N; ++i ) {
        for ( std::size_t j{ i + 1 }; j < N; ++j ) {
            double const dx{ p.pos_x()[j] - p.pos_x()[i] };
            double const dy{ p.pos_y()[j] - p.pos_y()[i] };
            double const dz{ p.pos_z()[j] - p.pos_z()[i] };
            double const r{ std::sqrt( dx*dx + dy*dy + dz*dz ) };
            pe -= config::G * p.mass()[i] * p.mass()[j] / r;
        }
    }
    return pe;
}

static double total_energy( Particles const &p ) {
    return kinetic_energy( p ) + potential_energy( p );
}

static double linear_momentum( Particles const &p ) {
    double px{}, py{}, pz{};
    for ( std::size_t i{}; i < p.num_particles(); ++i ) {
        px += p.mass()[i] * p.vel_x()[i];
        py += p.mass()[i] * p.vel_y()[i];
        pz += p.mass()[i] * p.vel_z()[i];
    }
    return std::sqrt( px*px + py*py + pz*pz );
}

static double angular_momentum( Particles const &p ) {
    double lx{}, ly{}, lz{};
    for ( std::size_t i{}; i < p.num_particles(); ++i ) {
        double const m{ p.mass()[i] };
        lx += m * ( p.pos_y()[i]*p.vel_z()[i] - p.pos_z()[i]*p.vel_y()[i] );
        ly += m * ( p.pos_z()[i]*p.vel_x()[i] - p.pos_x()[i]*p.vel_z()[i] );
        lz += m * ( p.pos_x()[i]*p.vel_y()[i] - p.pos_y()[i]*p.vel_x()[i] );
    }
    return std::sqrt( lx*lx + ly*ly + lz*lz );
}

// Set up a circular two-body Kepler orbit (Sun + Earth-like).
// Body 0 and body 1 placed at the barycenter with circular velocity.
static void setup_two_body( Particles &p, double const M, double const m, double const r ) {
    double const mu{ m / ( M + m ) };
    double const r0{ -mu * r };
    double const r1{ ( 1.0 - mu ) * r };

    double const v_orb{ std::sqrt( config::G * ( M + m ) / r ) };
    double const v0{ -mu * v_orb };
    double const v1{ ( 1.0 - mu ) * v_orb };

    p.pos_x()[0] = r0;  p.pos_y()[0] = 0.0; p.pos_z()[0] = 0.0;
    p.vel_x()[0] = 0.0; p.vel_y()[0] = v0;  p.vel_z()[0] = 0.0;
    p.mass()[0]  = M;

    p.pos_x()[1] = r1;  p.pos_y()[1] = 0.0; p.pos_z()[1] = 0.0;
    p.vel_x()[1] = 0.0; p.vel_y()[1] = v1;  p.vel_z()[1] = 0.0;
    p.mass()[1]  = m;

    for ( std::size_t i{}; i < 2; ++i ) {
        p.acc_x()[i] = 0.0;
        p.acc_y()[i] = 0.0;
        p.acc_z()[i] = 0.0;
    }
}

static void step_n( Particles &p, Integrator &integ,
                    std::vector<std::unique_ptr<Force>> &forces, std::size_t const n ) {
    for ( std::size_t s{}; s < n; ++s ) {
        integ.integrate( p, forces );
    }
}

// Sample N particles from a Plummer sphere with scale radius a and total mass M.
// r = a / sqrt(u^(-2/3) - 1), inverse CDF of the enclosed-mass profile.
// Isotropic angular sampling. Velocities and accelerations zeroed.
static void populate_plummer( Particles &p, std::size_t const N,
                              double const a, double const M_tot,
                              unsigned const seed ) {
    std::mt19937_64 rng{ seed };
    std::uniform_real_distribution<double> u_dist{ 0.0, 1.0 };
    double const m_per{ M_tot / static_cast<double>( N ) };

    for ( std::size_t i{}; i < N; ++i ) {
        double const u{ std::clamp( u_dist( rng ), 1e-6, 1.0 - 1e-6 ) };
        double const r{ a / std::sqrt( std::pow( u, -2.0 / 3.0 ) - 1.0 ) };

        double const cos_theta{ 2.0 * u_dist( rng ) - 1.0 };
        double const sin_theta{ std::sqrt( 1.0 - cos_theta * cos_theta ) };
        double const phi{ 2.0 * std::numbers::pi * u_dist( rng ) };

        p.pos_x()[i] = r * sin_theta * std::cos( phi );
        p.pos_y()[i] = r * sin_theta * std::sin( phi );
        p.pos_z()[i] = r * cos_theta;
        p.vel_x()[i] = 0.0; p.vel_y()[i] = 0.0; p.vel_z()[i] = 0.0;
        p.acc_x()[i] = 0.0; p.acc_y()[i] = 0.0; p.acc_z()[i] = 0.0;
        p.mass()[i]  = m_per;
    }
}

// Uniform random positions in [-L, L]^3 with equal mass per particle.
static void populate_uniform_cube( Particles &p, std::size_t const N,
                                   double const L, double const M_tot,
                                   unsigned const seed ) {
    std::mt19937_64 rng{ seed };
    std::uniform_real_distribution<double> dist{ -L, L };
    double const m_per{ M_tot / static_cast<double>( N ) };

    for ( std::size_t i{}; i < N; ++i ) {
        p.pos_x()[i] = dist( rng );
        p.pos_y()[i] = dist( rng );
        p.pos_z()[i] = dist( rng );
        p.vel_x()[i] = 0.0; p.vel_y()[i] = 0.0; p.vel_z()[i] = 0.0;
        p.acc_x()[i] = 0.0; p.acc_y()[i] = 0.0; p.acc_z()[i] = 0.0;
        p.mass()[i]  = m_per;
    }
}

struct ForceErrorStats {
    double mean_rel;
    double max_rel;
};

// Apply direct and BH forces to identical configurations, return per-particle
// relative force error |a_BH - a_direct| / |a_direct|, aggregated as mean and max.
static ForceErrorStats compare_forces( Particles &p_direct, Particles &p_bh,
                                       double const theta ) {
    Gravity{}.apply( p_direct );
    Gravity_BarnesHut{ theta }.apply( p_bh );

    std::size_t const N{ p_direct.num_particles() };
    double sum_rel{};
    double max_rel{};
    std::size_t counted{};

    for ( std::size_t i{}; i < N; ++i ) {
        double const ax{ p_direct.acc_x()[i] };
        double const ay{ p_direct.acc_y()[i] };
        double const az{ p_direct.acc_z()[i] };
        double const a_mag{ std::sqrt( ax*ax + ay*ay + az*az ) };
        if ( a_mag == 0.0 ) continue;

        double const dx{ p_bh.acc_x()[i] - ax };
        double const dy{ p_bh.acc_y()[i] - ay };
        double const dz{ p_bh.acc_z()[i] - az };
        double const err{ std::sqrt( dx*dx + dy*dy + dz*dz ) };

        double const rel{ err / a_mag };
        sum_rel += rel;
        max_rel = std::max( max_rel, rel );
        ++counted;
    }

    return { sum_rel / static_cast<double>( counted ), max_rel };
}


// 1. Yoshida coefficients

TEST( yoshida_coefficients_sum_to_unity ) {
    Yoshida const y{ 1.0 };
    double const c_sum{ y.c_1() + y.c_2() + y.c_3() + y.c_4() };
    double const d_sum{ y.d_1() + y.d_2() + y.d_3() };
    ASSERT_NEAR( c_sum, 1.0, 1e-14 );
    ASSERT_NEAR( d_sum, 1.0, 1e-14 );
    ++g_pass;
}

TEST( yoshida_coefficients_symmetry ) {
    Yoshida const y{ 1.0 };
    ASSERT_NEAR( y.c_1(), y.c_4(), 1e-15 );
    ASSERT_NEAR( y.c_2(), y.c_3(), 1e-15 );
    ASSERT_NEAR( y.d_1(), y.d_3(), 1e-15 );
    ++g_pass;
}

TEST( yoshida_d2_is_negative ) {
    Yoshida const y{ 1.0 };
    ASSERT_TRUE( y.d_2() < 0.0 );
    ++g_pass;
}

TEST( yoshida_coefficients_dt_independent ) {
    // Raw weights are dt-independent; scaling happens inside integrate().
    Yoshida const y1{ 1.0 };
    Yoshida const y2{ 7.0 };
    ASSERT_NEAR( y2.c_1() / y2.c_4(), 1.0, 1e-14 );
    ASSERT_NEAR( y1.c_1(), y2.c_1(), 1e-14 );
    ASSERT_NEAR( y1.d_2(), y2.d_2(), 1e-14 );
    ++g_pass;
}


// 2. Force kernel

TEST( gravity_two_body_inverse_square ) {
    Particles p{ 2 };
    double const r{ 1e10 };
    p.pos_x()[0] = 0.0;   p.pos_y()[0] = 0.0; p.pos_z()[0] = 0.0;
    p.pos_x()[1] = r;     p.pos_y()[1] = 0.0; p.pos_z()[1] = 0.0;
    p.mass()[0] = 1e30;   p.mass()[1] = 1e24;
    for ( std::size_t i{}; i < 2; ++i ) {
        p.acc_x()[i] = 0.0; p.acc_y()[i] = 0.0; p.acc_z()[i] = 0.0;
    }

    Gravity g;
    g.apply( p );

    double const expected_a0{ config::G * p.mass()[1] / ( r * r ) };
    ASSERT_NEAR( p.acc_x()[0], expected_a0, expected_a0 * 1e-9 );
    ASSERT_NEAR( p.acc_y()[0], 0.0, 1e-30 );

    double const expected_a1{ -config::G * p.mass()[0] / ( r * r ) };
    ASSERT_NEAR( p.acc_x()[1], expected_a1, std::abs( expected_a1 ) * 1e-9 );
    ++g_pass;
}

TEST( gravity_self_interaction_zero ) {
    Particles p{ 1 };
    p.pos_x()[0] = 1e11; p.pos_y()[0] = 2e11; p.pos_z()[0] = 3e11;
    p.mass()[0] = 1.989e30;
    p.acc_x()[0] = 0.0; p.acc_y()[0] = 0.0; p.acc_z()[0] = 0.0;

    Gravity g;
    g.apply( p );

    ASSERT_NEAR( p.acc_x()[0], 0.0, 1e-30 );
    ASSERT_NEAR( p.acc_y()[0], 0.0, 1e-30 );
    ASSERT_NEAR( p.acc_z()[0], 0.0, 1e-30 );
    ++g_pass;
}

TEST( gravity_newton_third_law ) {
    Particles p{ 2 };
    p.pos_x()[0] = 0.0;     p.pos_y()[0] = 0.0;   p.pos_z()[0] = 0.0;
    p.pos_x()[1] = 1.5e11;  p.pos_y()[1] = 1e11;  p.pos_z()[1] = 0.0;
    p.mass()[0] = 1.989e30; p.mass()[1] = 5.972e24;
    for ( std::size_t i{}; i < 2; ++i ) {
        p.acc_x()[i] = 0.0; p.acc_y()[i] = 0.0; p.acc_z()[i] = 0.0;
    }

    Gravity g;
    g.apply( p );

    double const fx_total{ p.mass()[0]*p.acc_x()[0] + p.mass()[1]*p.acc_x()[1] };
    double const fy_total{ p.mass()[0]*p.acc_y()[0] + p.mass()[1]*p.acc_y()[1] };
    double const fz_total{ p.mass()[0]*p.acc_z()[0] + p.mass()[1]*p.acc_z()[1] };

    double const f_scale{ std::abs( p.mass()[0] * p.acc_x()[0] ) };
    ASSERT_LT( std::abs( fx_total ) / f_scale, 1e-12 );
    ASSERT_LT( std::abs( fy_total ) / f_scale, 1e-12 );
    ASSERT_NEAR( fz_total, 0.0, 1e-30 );
    ++g_pass;
}

TEST( gravity_accumulates_not_overwrites ) {
    Particles p{ 2 };
    p.pos_x()[0] = 0.0;   p.pos_y()[0] = 0.0; p.pos_z()[0] = 0.0;
    p.pos_x()[1] = 1e11;  p.pos_y()[1] = 0.0; p.pos_z()[1] = 0.0;
    p.mass()[0] = 1e30;   p.mass()[1] = 1e30;

    double const init_ax{ 42.0 };
    p.acc_x()[0] = init_ax; p.acc_y()[0] = 0.0; p.acc_z()[0] = 0.0;
    p.acc_x()[1] = 0.0;    p.acc_y()[1] = 0.0; p.acc_z()[1] = 0.0;

    Gravity g;
    g.apply( p );

    ASSERT_TRUE( std::abs( p.acc_x()[0] - init_ax ) > 1e-20 );
    ASSERT_TRUE( p.acc_x()[0] > init_ax );
    ++g_pass;
}


// 3. Two-body Kepler orbit

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


// 4. Yoshida vs Verlet convergence order

TEST( yoshida_fourth_order_convergence ) {
    // Compare position error at two timesteps covering the same physical duration.
    // dt ratio = 2, so position error ratio should be near 2^4 = 16.
    double const M{ 1.989e30 };
    double const m{ 5.972e24 };
    double const r{ 1.496e11 };

    double const dt_coarse{ 14400.0 };
    double const dt_fine{ 7200.0 };
    double const dt_ref{ 900.0 };
    std::size_t const n_coarse{ 2000 };
    double const T_phys{ dt_coarse * n_coarse };
    std::size_t const n_fine{ static_cast<std::size_t>( T_phys / dt_fine ) };
    std::size_t const n_ref{ static_cast<std::size_t>( T_phys / dt_ref ) };

    auto run = [&]( double const dt, std::size_t const steps ) -> std::pair<double, double> {
        Particles p{ 2 };
        setup_two_body( p, M, m, r );
        std::vector<std::unique_ptr<Force>> forces;
        forces.push_back( std::make_unique<Gravity>() );
        Yoshida integ{ dt };
        step_n( p, integ, forces, steps );
        return { p.pos_x()[1], p.pos_y()[1] };
    };

    auto [rx, ry] = run( dt_ref, n_ref );
    auto [cx, cy] = run( dt_coarse, n_coarse );
    auto [fx, fy] = run( dt_fine, n_fine );

    double const err_coarse{ std::sqrt( (cx-rx)*(cx-rx) + (cy-ry)*(cy-ry) ) };
    double const err_fine{ std::sqrt( (fx-rx)*(fx-rx) + (fy-ry)*(fy-ry) ) };

    double const ratio{ err_coarse / err_fine };
    ASSERT_TRUE( ratio > 6.0 );
    ASSERT_TRUE( ratio < 40.0 );
    ++g_pass;
}

TEST( verlet_second_order_convergence ) {
    // Same approach for Velocity Verlet: dt ratio = 2, expect ratio near 2^2 = 4.
    double const M{ 1.989e30 };
    double const m{ 5.972e24 };
    double const r{ 1.496e11 };

    double const dt_coarse{ 14400.0 };
    double const dt_fine{ 7200.0 };
    double const dt_ref{ 900.0 };
    std::size_t const n_coarse{ 2000 };
    double const T_phys{ dt_coarse * n_coarse };
    std::size_t const n_fine{ static_cast<std::size_t>( T_phys / dt_fine ) };
    std::size_t const n_ref{ static_cast<std::size_t>( T_phys / dt_ref ) };

    auto run = [&]( double const dt, std::size_t const steps ) -> std::pair<double, double> {
        Particles p{ 2 };
        setup_two_body( p, M, m, r );
        std::vector<std::unique_ptr<Force>> forces;
        forces.push_back( std::make_unique<Gravity>() );
        Velocity_Verlet integ{ dt };
        step_n( p, integ, forces, steps );
        return { p.pos_x()[1], p.pos_y()[1] };
    };

    auto [rx, ry] = run( dt_ref, n_ref );
    auto [cx, cy] = run( dt_coarse, n_coarse );
    auto [fx, fy] = run( dt_fine, n_fine );

    double const err_coarse{ std::sqrt( (cx-rx)*(cx-rx) + (cy-ry)*(cy-ry) ) };
    double const err_fine{ std::sqrt( (fx-rx)*(fx-rx) + (fy-ry)*(fy-ry) ) };

    double const ratio{ err_coarse / err_fine };
    ASSERT_TRUE( ratio > 2.0 );
    ASSERT_TRUE( ratio < 8.0 );
    ++g_pass;
}


// 5. Three-body conservation

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


// 6. Particle SoA layout

TEST( particle_soa_contiguous_memory ) {
    // Verify that SoA sub-arrays are evenly strided in memory.
    // With aligned allocation, the stride may be larger than N due to
    // SIMD-width padding (e.g. N=10 rounds up to 12 on AVX2).
    Particles p{ 10 };
    double* px{ p.pos_x() };
    double* py{ p.pos_y() };

    // pos_y should start exactly one stride after pos_x.
    std::ptrdiff_t const stride{ py - px };
    ASSERT_TRUE( stride >= 10 );

    // All subsequent arrays should be at the same stride interval.
    double* pz{ p.pos_z() };
    ASSERT_TRUE( pz - py == stride );

    // mass is sub-array index 12, so it should be 12 strides from pos_x.
    double* mass{ p.mass() };
    ASSERT_TRUE( mass == px + 12 * stride );
    ++g_pass;
}


// 7. Barnes-Hut

TEST( bh_self_interaction_zero ) {
    Particles p{ 1 };
    p.pos_x()[0] = 1e11; p.pos_y()[0] = 2e11; p.pos_z()[0] = 3e11;
    p.mass()[0] = 1.989e30;
    p.acc_x()[0] = 0.0; p.acc_y()[0] = 0.0; p.acc_z()[0] = 0.0;

    Gravity_BarnesHut bh{ 0.5 };
    bh.apply( p );

    ASSERT_NEAR( p.acc_x()[0], 0.0, 1e-30 );
    ASSERT_NEAR( p.acc_y()[0], 0.0, 1e-30 );
    ASSERT_NEAR( p.acc_z()[0], 0.0, 1e-30 );
    ++g_pass;
}

TEST( bh_two_body_matches_direct ) {
    // With N=2, BH never has internal nodes to summarize, so it must reduce
    // to the direct kernel up to summation-order rounding.
    auto setup = []( Particles &p ) {
        p.pos_x()[0] = 0.0;     p.pos_y()[0] = 0.0;   p.pos_z()[0] = 0.0;
        p.pos_x()[1] = 1.5e11;  p.pos_y()[1] = 1e11;  p.pos_z()[1] = 0.0;
        p.mass()[0] = 1.989e30; p.mass()[1] = 5.972e24;
        for ( std::size_t i{}; i < 2; ++i ) {
            p.acc_x()[i] = 0.0; p.acc_y()[i] = 0.0; p.acc_z()[i] = 0.0;
        }
    };

    Particles p_dir{ 2 }; setup( p_dir );
    Particles p_bh{ 2 };  setup( p_bh );

    Gravity{}.apply( p_dir );
    Gravity_BarnesHut{ 0.5 }.apply( p_bh );

    for ( std::size_t i{}; i < 2; ++i ) {
        double const scale{ std::sqrt( p_dir.acc_x()[i]*p_dir.acc_x()[i]
                                     + p_dir.acc_y()[i]*p_dir.acc_y()[i]
                                     + p_dir.acc_z()[i]*p_dir.acc_z()[i] ) };
        double const dx{ p_bh.acc_x()[i] - p_dir.acc_x()[i] };
        double const dy{ p_bh.acc_y()[i] - p_dir.acc_y()[i] };
        double const dz{ p_bh.acc_z()[i] - p_dir.acc_z()[i] };
        double const err{ std::sqrt( dx*dx + dy*dy + dz*dz ) };
        ASSERT_LT( err / scale, 1e-10 );
    }
    ++g_pass;
}

TEST( bh_small_cluster_matches_direct_at_low_theta ) {
    // theta close to zero forces the MAC to fail almost everywhere, so BH
    // opens nearly all internal nodes and converges to direct summation.
    // Tolerance is generous (1e-9) because summation order differs.
    constexpr std::size_t N{ 20 };
    std::mt19937_64 rng{ 12345 };
    std::uniform_real_distribution<double> pos_dist{ -1e11, 1e11 };
    std::uniform_real_distribution<double> mass_dist{ 1e23, 1e26 };

    auto populate = [&]( Particles &p ) {
        rng.seed( 12345 );
        for ( std::size_t i{}; i < N; ++i ) {
            p.pos_x()[i] = pos_dist( rng );
            p.pos_y()[i] = pos_dist( rng );
            p.pos_z()[i] = pos_dist( rng );
            p.mass()[i] = mass_dist( rng );
            p.acc_x()[i] = 0.0; p.acc_y()[i] = 0.0; p.acc_z()[i] = 0.0;
        }
    };

    Particles p_dir{ N }; populate( p_dir );
    Particles p_bh{ N };  populate( p_bh );

    Gravity{}.apply( p_dir );
    Gravity_BarnesHut{ 0.05, 1 }.apply( p_bh );  // tight MAC, single-particle leaves

    for ( std::size_t i{}; i < N; ++i ) {
        double const scale{ std::sqrt( p_dir.acc_x()[i]*p_dir.acc_x()[i]
                                     + p_dir.acc_y()[i]*p_dir.acc_y()[i]
                                     + p_dir.acc_z()[i]*p_dir.acc_z()[i] ) };
        double const dx{ p_bh.acc_x()[i] - p_dir.acc_x()[i] };
        double const dy{ p_bh.acc_y()[i] - p_dir.acc_y()[i] };
        double const dz{ p_bh.acc_z()[i] - p_dir.acc_z()[i] };
        double const err{ std::sqrt( dx*dx + dy*dy + dz*dz ) };
        ASSERT_LT( err / scale, 1e-9 );
    }
    ++g_pass;
}

TEST( bh_kepler_two_body_returns_to_start ) {
    // Same orbit setup as the direct version but driven by the BH force.
    // Tolerance is looser (1e-3) because BH approximates, though at N=2 the
    // tree has nothing to summarize. The looser bound guards against future
    // implementation changes (different schedules, etc.).
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
    forces.push_back( std::make_unique<Gravity_BarnesHut>( 0.5 ) );
    Yoshida integ{ dt };
    step_n( p, integ, forces, steps );

    double const dx{ p.pos_x()[1] - x1_init };
    double const dy{ p.pos_y()[1] - y1_init };
    double const pos_err{ std::sqrt( dx*dx + dy*dy ) };
    ASSERT_LT( pos_err / r, 1e-3 );
    ++g_pass;
}

TEST( bh_energy_drift_bounded ) {
    // BH breaks the symplectic guarantee: the multipole approximation is not
    // a conservative force at the discrete level. Energy is no longer
    // conserved at the 1e-12 floor of the direct kernel, but it should
    // remain bounded over short integrations. 1e-3 over 200 steps is the
    // regression bound; observed drift is typically much tighter.
    constexpr std::size_t N{ 50 };
    std::mt19937_64 rng{ 777 };
    std::uniform_real_distribution<double> pos_dist{ -5e11, 5e11 };
    std::uniform_real_distribution<double> vel_dist{ -1e3, 1e3 };
    std::uniform_real_distribution<double> mass_dist{ 1e23, 1e26 };

    Particles p{ N };
    for ( std::size_t i{}; i < N; ++i ) {
        p.pos_x()[i] = pos_dist( rng );
        p.pos_y()[i] = pos_dist( rng );
        p.pos_z()[i] = pos_dist( rng );
        p.vel_x()[i] = vel_dist( rng );
        p.vel_y()[i] = vel_dist( rng );
        p.vel_z()[i] = vel_dist( rng );
        p.mass()[i] = mass_dist( rng );
        p.acc_x()[i] = 0.0; p.acc_y()[i] = 0.0; p.acc_z()[i] = 0.0;
    }

    double const E0{ total_energy( p ) };

    std::vector<std::unique_ptr<Force>> forces;
    forces.push_back( std::make_unique<Gravity_BarnesHut>( 0.5 ) );
    Yoshida integ{ 900.0 };
    step_n( p, integ, forces, 200 );

    double const dE_rel{ std::abs( ( total_energy( p ) - E0 ) / E0 ) };
    ASSERT_LT( dE_rel, 1e-3 );
    ++g_pass;
}

TEST( bh_force_accuracy_plummer ) {
    // Plummer sphere: centrally concentrated, stresses adaptive tree refinement.
    // Bounds are loose initial regression values. After first run on production
    // hardware, observe the actual values and tighten to ~2-3x observed.
    // Reference: Hernquist 1987 reports RMS relative force error of ~0.5% at
    // theta=0.5 for monopole BH; bounds allow headroom for N=2000 finite-size
    // effects and summation-order rounding.
    constexpr std::size_t N{ 2000 };
    double const a{ 1e11 };
    double const M_tot{ 1e30 };

    Particles p_direct{ N };
    Particles p_bh{ N };
    populate_plummer( p_direct, N, a, M_tot, 42 );
    populate_plummer( p_bh,     N, a, M_tot, 42 );

    ForceErrorStats const stats{ compare_forces( p_direct, p_bh, 0.5 ) };

    ASSERT_LT( stats.mean_rel, 1e-2 );
    ASSERT_LT( stats.max_rel,  1e-1 );
    ++g_pass;
}

TEST( bh_force_accuracy_uniform_cube ) {
    // Uniform cube: balanced tree, isolates bulk-traversal correctness from
    // adaptive-refinement behavior. Divergence between this and the Plummer
    // test localizes bugs to one path or the other.
    constexpr std::size_t N{ 2000 };
    double const L{ 1e11 };
    double const M_tot{ 1e30 };

    Particles p_direct{ N };
    Particles p_bh{ N };
    populate_uniform_cube( p_direct, N, L, M_tot, 42 );
    populate_uniform_cube( p_bh,     N, L, M_tot, 42 );

    ForceErrorStats const stats{ compare_forces( p_direct, p_bh, 0.5 ) };

    // Cube errors should be smaller than Plummer (no deep clusters).
    ASSERT_LT( stats.mean_rel, 5e-3 );
    ASSERT_LT( stats.max_rel,  5e-2 );
    ++g_pass;
}

TEST( bh_force_error_decreases_with_theta ) {
    // Verify BH error decreases monotonically as theta tightens. This is a
    // weaker claim than the theta^2 leading-order scaling but it's the one
    // that actually holds for the mean-relative-error metric. The per-cell
    // theta^2 scaling is real but doesn't survive aggregation across particles
    // (cancellations from symmetric distributions produce steeper empirical
    // scaling, closer to theta^3 or theta^4 depending on the theta range).
    constexpr std::size_t N{ 2000 };
    double const L{ 1e11 };
    double const M_tot{ 1e30 };

    auto run = [&]( double const theta ) -> ForceErrorStats {
        Particles p_d{ N }, p_b{ N };
        populate_uniform_cube( p_d, N, L, M_tot, 42 );
        populate_uniform_cube( p_b, N, L, M_tot, 42 );
        return compare_forces( p_d, p_b, theta );
    };

    ForceErrorStats const s_10{  run( 1.0 ) };
    ForceErrorStats const s_05{  run( 0.5 ) };
    ForceErrorStats const s_025{ run( 0.25 ) };

    std::cout << "\n    [diag] mean rel err: theta=1.0 " << std::scientific << std::setprecision( 3 )
              << s_10.mean_rel << ", theta=0.5 " << s_05.mean_rel
              << ", theta=0.25 " << s_025.mean_rel << "\n   ";

    // Monotonic improvement (the actually-true scaling property).
    ASSERT_TRUE( s_10.mean_rel > s_05.mean_rel );
    ASSERT_TRUE( s_05.mean_rel > s_025.mean_rel );

    // Significant total dynamic range across the theta sweep. Empirically
    // ~165x on uniform cube N=2000; assert > 20x to catch "BH ignores theta"
    // failure modes while leaving slack for hardware and seed variation.
    ASSERT_TRUE( s_10.mean_rel / s_025.mean_rel > 20.0 );
    ++g_pass;
}

// Main

int main() {
    std::cout << "\n=== N-Body Unit Tests ===\n\n";

    for ( auto const &t : tests() ) {
        std::cout << "  " << t.name << "... " << std::flush;
        int const fail_before{ g_fail };
        t.fn();
        if ( g_fail == fail_before ) {
            std::cout << "PASS\n";
        }
    }

    std::cout << "\n" << g_pass << " passed, " << g_fail << " failed.\n\n";
    return g_fail > 0 ? 1 : 0;
}