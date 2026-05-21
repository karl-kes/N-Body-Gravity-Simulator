// Barnes-Hut octree: correctness, convergence to direct at low theta,
// energy bounds, and per-cell error monotonicity with theta.

#include "test_harness.hpp"
#include "helpers.hpp"

#include "Force/Force.hpp"
#include "Force/BarnesHut.hpp"
#include "Integrator/Integrator.hpp"
#include "Particle/Particle.hpp"
#include "Config.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numbers>
#include <random>
#include <vector>


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
    // Hernquist 1987 reports RMS relative force error of ~0.5% at theta=0.5 for
    // monopole BH. Bounds allow headroom for finite-N (N=2000) and summation
    // order variation.
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

    // Monotonic improvement.
    ASSERT_TRUE( s_10.mean_rel > s_05.mean_rel );
    ASSERT_TRUE( s_05.mean_rel > s_025.mean_rel );

    // Significant total dynamic range across the theta sweep.
    // Empirically ~150x on uniform cube N=2000; assert > 20x to catch
    // "BH ignores theta" failure modes while leaving slack for hardware
    // and seed variation.
    ASSERT_TRUE( s_10.mean_rel / s_025.mean_rel > 20.0 );
    ++g_pass;
}
