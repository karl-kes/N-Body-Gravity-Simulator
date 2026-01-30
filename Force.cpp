#include "Force.hpp"

Gravity::Gravity()
{ }

void Gravity::apply( Particles &particles ) const {
    std::size_t const N{ particles.num_particles() };
    double const EPS{ 1e-9 };

    for ( std::size_t i{}; i < N; ++i ) {
        auto &a{ particles };
        for ( std::size_t j{i + 1}; j < N; ++j ) {
            auto &b{ particles };

            double dx{ b.pos_x()[j] - a.pos_x()[i] };
            double dy{ b.pos_y()[j] - a.pos_y()[i] };
            double dz{ b.pos_z()[j] - a.pos_z()[i] };

            double R_sq{ dx*dx + dy*dy + dz*dz + EPS*EPS };
            double inv_R_cb{ 1.0 / ( std::sqrt( R_sq ) * R_sq ) };

            a.acc_x()[i] += G() * b.mass()[j] * inv_R_cb * dx;
            a.acc_y()[i] += G() * b.mass()[j] * inv_R_cb * dy;
            a.acc_z()[i] += G() * b.mass()[j] * inv_R_cb * dz;

            b.acc_x()[j] -= G() * a.mass()[i] * inv_R_cb * dx;
            b.acc_y()[j] -= G() * a.mass()[i] * inv_R_cb * dy;
            b.acc_z()[j] -= G() * a.mass()[i] * inv_R_cb * dz;
        }
    }
}