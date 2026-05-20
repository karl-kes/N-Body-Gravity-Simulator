#include "Force.hpp"

#include <omp.h>

Gravity::Gravity( double const eps )
: eps_{ eps }
{ }

void Gravity::apply( Particles &particles ) const {
    std::size_t const N{ particles.num_particles() };

    double const* RESTRICT px{ particles.pos_x() };
    double const* RESTRICT py{ particles.pos_y() };
    double const* RESTRICT pz{ particles.pos_z() };

    double* RESTRICT ax{ particles.acc_x() };
    double* RESTRICT ay{ particles.acc_y() };
    double* RESTRICT az{ particles.acc_z() };

    double const* RESTRICT mass{ particles.mass() };

    double const eps_sq{ eps_ * eps_ };
    constexpr double G{ config::G };

    auto apply_kernel = [px, py, pz, ax, ay, az, mass, N, eps_sq]( std::size_t i ) {
        double const pxi{ px[i] }, pyi{ py[i] }, pzi{ pz[i] };

        double a_xi{}, a_yi{}, a_zi{};

        #pragma omp simd reduction( +:a_xi, a_yi, a_zi )
        for ( std::size_t j = 0; j < N; ++j ) {
            double const mask{ ( i == j ) ? 0.0 : 1.0 };

            accumulate_pairwise (
                pxi, pyi, pzi,
                px[j], py[j], pz[j], mass[j],
                a_xi, a_yi, a_zi,
                G, eps_sq, mask
            );
        }

        ax[i] += a_xi;
        ay[i] += a_yi;
        az[i] += a_zi;
    };

    if ( N >= config::OMP_THRESHOLD ) {
        #pragma omp parallel for schedule( static )
        for ( std::size_t i = 0; i < N; ++i ) {
            apply_kernel( i );
        }
    } else {
        for ( std::size_t i = 0; i < N; ++i ) {
            apply_kernel( i );
        }
    }
}