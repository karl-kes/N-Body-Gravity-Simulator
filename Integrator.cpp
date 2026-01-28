#include "Integrator.hpp"

Velocity_Verlet::Velocity_Verlet( double dt )
: dt_{ dt }
{ }

void Velocity_Verlet::integrate( std::vector<Particle> &particles, std::vector<std::unique_ptr<Force_Law>> const &forces ) const {
    std::size_t const N{ particles.size() };

    std::vector<double> old_acc_x{}, old_acc_y{}, old_acc_z{};
    old_acc_x.resize( N );
    old_acc_y.resize( N );
    old_acc_z.resize( N );

    for ( std::size_t i = 0; i < N; ++i ) {
        auto &a{ particles[i] };

        a.pos_x() += dt() * ( a.vel_x() + 0.5 * a.acc_x() * dt() );
        a.pos_y() += dt() * ( a.vel_y() + 0.5 * a.acc_y() * dt() );
        a.pos_z() += dt() * ( a.vel_z() + 0.5 * a.acc_z() * dt() );

        old_acc_x[i] = a.acc_x();
        old_acc_y[i] = a.acc_y();
        old_acc_z[i] = a.acc_z();
    }

    for ( auto &particle : particles ) {
        particle.clear_acc();
    }

    for ( auto const &force : forces ) {
        force->apply( particles );
    }

    for ( std::size_t i = 0; i < N; ++i ) {
        auto &a{ particles[i] };

        a.vel_x() += 0.5 * ( old_acc_x[i] + a.acc_x() ) * dt();
        a.vel_y() += 0.5 * ( old_acc_y[i] + a.acc_y() ) * dt();
        a.vel_z() += 0.5 * ( old_acc_z[i] + a.acc_z() ) * dt();
    }
}