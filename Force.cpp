#include "Force.hpp"

Gravity::Gravity()
: G_{ 6.6743e-11 }
{ }

void Gravity::apply( std::vector<Particle> &particles ) const {
    std::size_t const N{ particles.size() };
    double const EPS{ 1e-9 };

    for ( std::size_t i = 0; i < N; ++i ) {
        auto &a{ particles[i] };
        for ( std::size_t j = i + 1; j < N; ++j ) {
            auto &b{ particles[j] };

            double dx{ b.pos_x() - a.pos_x() };
            double dy{ b.pos_y() - a.pos_y() };
            double dz{ b.pos_z() - a.pos_z() };

            double R_sq{ dx*dx + dy*dy + dz*dz + EPS };
            double inv_R_cb{ 1.0 / ( std::sqrt( R_sq ) * R_sq ) };

            a.acc_x() += G() * b.mass() * inv_R_cb * dx;
            a.acc_y() += G() * b.mass() * inv_R_cb * dy;
            a.acc_z() += G() * b.mass() * inv_R_cb * dz;

            b.acc_x() -= G() * a.mass() * inv_R_cb * dx;
            b.acc_y() -= G() * a.mass() * inv_R_cb * dy;
            b.acc_z() -= G() * a.mass() * inv_R_cb * dz;
        }
    }
}