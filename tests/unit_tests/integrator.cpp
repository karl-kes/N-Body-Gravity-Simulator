// Yoshida and Velocity Verlet: coefficient invariants and empirical convergence order.

#include "test_harness.hpp"
#include "helpers.hpp"

#include "Force/Force.hpp"
#include "Integrator/Integrator.hpp"
#include "Config.hpp"

#include <cmath>
#include <numbers>


// Yoshida coefficient invariants

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


// Empirical convergence order via Richardson extrapolation

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
