// Direct gravity kernel correctness: pair interactions, masks, accumulation.

#include "test_harness.hpp"
#include "helpers.hpp"

#include "Force/Force.hpp"
#include "Particle/Particle.hpp"
#include "Config.hpp"

#include <cmath>


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
