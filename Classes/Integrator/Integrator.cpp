#include "Integrator.hpp"

#define RESTRICT __restrict

Velocity_Verlet::Velocity_Verlet( double dt )
: Integrator{ dt, "Velocity Verlet" }
{ }

void Velocity_Verlet::integrate( Particles &particles, std::vector<std::unique_ptr<Force>> const &forces ) const {
    std::size_t const N{ particles.num_particles() };

    #pragma omp parallel for schedule( static )
    for ( std::size_t i = 0; i < N; ++i ) {
        particles.pos_x()[i] += dt() * ( particles.vel_x()[i] + 0.5 * particles.acc_x()[i] * dt() );
        particles.pos_y()[i] += dt() * ( particles.vel_y()[i] + 0.5 * particles.acc_y()[i] * dt() );
        particles.pos_z()[i] += dt() * ( particles.vel_z()[i] + 0.5 * particles.acc_z()[i] * dt() );

        particles.old_acc_x()[i] = particles.acc_x()[i];
        particles.old_acc_y()[i] = particles.acc_y()[i];
        particles.old_acc_z()[i] = particles.acc_z()[i];

        particles.acc_x()[i] = 0.0;
        particles.acc_y()[i] = 0.0;
        particles.acc_z()[i] = 0.0;
    }

    for ( auto const &force : forces ) {
        force->apply( particles );
    }

    #pragma omp parallel for schedule( static )
    for ( std::size_t i = 0; i < N; ++i ) {
        particles.vel_x()[i] += 0.5 * ( particles.old_acc_x()[i] + particles.acc_x()[i] ) * dt();
        particles.vel_y()[i] += 0.5 * ( particles.old_acc_y()[i] + particles.acc_y()[i] ) * dt();
        particles.vel_z()[i] += 0.5 * ( particles.old_acc_z()[i] + particles.acc_z()[i] ) * dt();
    }
}

Yoshida::Yoshida( double const dt )
: Integrator{ dt, "Yoshida" }
, w_0_{ -cbrt_2() / ( 2.0 - cbrt_2() ) }
, w_1_{ 1.0 / ( 2.0 - cbrt_2() ) }
, c_1_{ w_1() / 2.0 }
, c_2_{ ( w_0() + w_1() ) / 2.0 }
, c_3_{ c_2() }
, c_4_{ c_1() }
, d_1_{ w_1() }
, d_2_{ w_0() }
, d_3_{ w_1() }
{ }

void Yoshida::integrate( Particles &particles, std::vector<std::unique_ptr<Force>> const &forces ) const {
    std::size_t const N{ particles.num_particles() };

    double* RESTRICT px{ particles.pos_x().get() };
    double* RESTRICT py{ particles.pos_y().get() };
    double* RESTRICT pz{ particles.pos_z().get() };

    double* RESTRICT vx{ particles.vel_x().get() };
    double* RESTRICT vy{ particles.vel_y().get() };
    double* RESTRICT vz{ particles.vel_z().get() };

    double* RESTRICT ax{ particles.acc_x().get() };
    double* RESTRICT ay{ particles.acc_y().get() };
    double* RESTRICT az{ particles.acc_z().get() };

    auto calculate_pos = [&]( double const c ) {
        double const c_dt{ c * dt() };

        #pragma omp parallel for schedule( static )
        for ( std::size_t i = 0; i < N; ++i ) {
            px[i] += c_dt * vx[i];
            py[i] += c_dt * vy[i];
            pz[i] += c_dt * vz[i];
        }
    };

    auto apply_force = [&](){
        #pragma omp parallel for schedule( static )
        for ( std::size_t i = 0; i < N; ++i ) {
            ax[i] = 0.0;
            ay[i] = 0.0;
            az[i] = 0.0;
        }

        for ( auto const &force : forces ) {
            force->apply( particles );
        };
    };

    auto calculate_vel = [&]( double const d ) {
        double const d_dt{ d * dt() };
        #pragma omp parallel for schedule( static )
        for ( std::size_t i = 0; i < N; ++i ) {
            vx[i] += d_dt * ax[i];
            vy[i] += d_dt * ay[i];
            vz[i] += d_dt * az[i];
        }
    };

    // Step 1:
    calculate_pos( c_1() );
    apply_force();
    calculate_vel( d_1() );

    // Step 2:
    calculate_pos( c_2() );
    apply_force();
    calculate_vel( d_2() );

    // Step 3:
    calculate_pos( c_3() );
    apply_force();
    calculate_vel( d_3() );

    // Step 4:
    calculate_pos( c_4() );
}