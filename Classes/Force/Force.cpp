#include "Force.hpp"

Gravity::Gravity()
{ }

void Gravity::apply( Particles &particles ) const {
    std::size_t const N{ particles.num_particles() };

    #pragma omp parallel for schedule( static )
    for ( std::size_t i = 0; i < N; ++i ) {
        double a_x{}, a_y{}, a_z{};

        for ( std::size_t j = 0; j < N; ++j ) {
            double const dx{ particles.pos_x()[j] - particles.pos_x()[i] };
            double const dy{ particles.pos_y()[j] - particles.pos_y()[i] };
            double const dz{ particles.pos_z()[j] - particles.pos_z()[i] };

            double const R_sq{ dx*dx + dy*dy + dz*dz + constant::EPS*constant::EPS };
            double const inv_R_cb{ 1.0 / ( std::sqrt( R_sq ) * R_sq ) };

            double const factor{ constant::G * particles.mass()[j] * inv_R_cb };

            a_x += factor * dx;
            a_y += factor * dy;
            a_z += factor * dz;
        }

        particles.acc_x()[i] = a_x;
        particles.acc_y()[i] = a_y;
        particles.acc_z()[i] = a_z;
    }
}