#include "Integrator.hpp"

#include <omp.h>

Velocity_Verlet::Velocity_Verlet( double dt )
: Integrator{ dt, "Velocity Verlet" }
{ }

void Velocity_Verlet::integrate( Particles &particles, std::vector<std::unique_ptr<Force>> const &forces ) const {
    std::size_t const N{ particles.num_particles() };

    double* RESTRICT px{ particles.pos_x() };
    double* RESTRICT py{ particles.pos_y() };
    double* RESTRICT pz{ particles.pos_z() };

    double* RESTRICT vx{ particles.vel_x() };
    double* RESTRICT vy{ particles.vel_y() };
    double* RESTRICT vz{ particles.vel_z() };

    double* RESTRICT ax{ particles.acc_x() };
    double* RESTRICT ay{ particles.acc_y() };
    double* RESTRICT az{ particles.acc_z() };

    double* RESTRICT o_ax{ particles.old_acc_x() };
    double* RESTRICT o_ay{ particles.old_acc_y() };
    double* RESTRICT o_az{ particles.old_acc_z() };

    double const dt_local{ dt() };

    #pragma omp parallel for schedule( static ) if ( N >= config::OMP_THRESHOLD )
    for ( std::size_t i = 0; i < N; ++i ) {
        px[i] += dt_local * ( vx[i] + 0.5 * ax[i] * dt_local );
        py[i] += dt_local * ( vy[i] + 0.5 * ay[i] * dt_local );
        pz[i] += dt_local * ( vz[i] + 0.5 * az[i] * dt_local );

        o_ax[i] = ax[i];
        o_ay[i] = ay[i];
        o_az[i] = az[i];

        ax[i] = 0.0;
        ay[i] = 0.0;
        az[i] = 0.0;
    }

    for ( auto const &force : forces ) {
        force->apply( particles );
    }

    #pragma omp parallel for schedule( static ) if ( N >= config::OMP_THRESHOLD )
    for ( std::size_t i = 0; i < N; ++i ) {
        vx[i] += 0.5 * ( o_ax[i] + ax[i] ) * dt_local;
        vy[i] += 0.5 * ( o_ay[i] + ay[i] ) * dt_local;
        vz[i] += 0.5 * ( o_az[i] + az[i] ) * dt_local;
    }
}

Yoshida::Yoshida( double const dt )
: Integrator{ dt, "Yoshida" }
, c_1_{ w_1() / 2.0 }
, c_2_{ ( w_0() + w_1() ) / 2.0 }
, c_3_{ ( w_0() + w_1() ) / 2.0 }
, c_4_{ w_1() / 2.0 }
, d_1_{ w_1() }
, d_2_{ w_0() }
, d_3_{ w_1() }
{ }

void Yoshida::integrate( Particles &particles, std::vector<std::unique_ptr<Force>> const &forces ) const {
    std::size_t const N{ particles.num_particles() };

    double* RESTRICT px{ particles.pos_x() };
    double* RESTRICT py{ particles.pos_y() };
    double* RESTRICT pz{ particles.pos_z() };

    double* RESTRICT vx{ particles.vel_x() };
    double* RESTRICT vy{ particles.vel_y() };
    double* RESTRICT vz{ particles.vel_z() };

    double* RESTRICT ax{ particles.acc_x() };
    double* RESTRICT ay{ particles.acc_y() };
    double* RESTRICT az{ particles.acc_z() };

    auto calculate_pos = [this, px, py, pz, vx, vy, vz, N]( double const c ) {
        double const c_dt{ c * dt() };

        #pragma omp simd
        for ( std::size_t i = 0; i < N; ++i ) {
            px[i] += c_dt * vx[i];
            py[i] += c_dt * vy[i];
            pz[i] += c_dt * vz[i];
        }
    };

    auto apply_force = [&forces, &particles, ax, ay, az, N]() {
        #pragma omp simd
        for ( std::size_t i = 0; i < N; ++i ) {
            ax[i] = 0.0;
            ay[i] = 0.0;
            az[i] = 0.0;
        }

        for ( auto const &force : forces ) {
            force->apply( particles );
        };
    };

    auto calculate_vel = [this, vx, vy, vz, ax, ay, az, N]( double const d ) {
        double const d_dt{ d * dt() };

        #pragma omp simd
        for ( std::size_t i = 0; i < N; ++i ) {
            vx[i] += d_dt * ax[i];
            vy[i] += d_dt * ay[i];
            vz[i] += d_dt * az[i];
        }
    };

    // d2 (= w0) is negative: the backward sub-step cancels lower-order error terms.
    calculate_pos( c_1() );
    apply_force();
    calculate_vel( d_1() );

    calculate_pos( c_2() );
    apply_force();
    calculate_vel( d_2() );

    calculate_pos( c_3() );
    apply_force();
    calculate_vel( d_3() );

    calculate_pos( c_4() );
}