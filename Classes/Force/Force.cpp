#include "Force.hpp"

#define RESTRICT __restrict

Gravity::Gravity()
{ }

void Gravity::apply( Particles &particles ) const {
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

    double const* RESTRICT mass{ particles.mass().get() };

    constexpr double eps_sq{ constant::EPS * constant::EPS };
    constexpr double G{ constant::G };
    constexpr double c_sq{ constant::C_SQ };

    constexpr std::size_t sun_idx{0};
    double const GM_sun{ G * mass[sun_idx] };
    double const sun_x{ px[sun_idx] }, sun_y{ py[sun_idx] }, sun_z{ pz[sun_idx] };

    #pragma omp parallel for schedule( static )
    for ( std::size_t i = 0; i < N; ++i ) {
        double const pxi{ px[i] }, pyi{ py[i] }, pzi{ pz[i] };
        double const vxi{ vx[i] }, vyi{ vy[i] }, vzi{ vz[i] };
        double a_x{}, a_y{}, a_z{};

        for ( std::size_t j = 0; j < N ; ++j ) {
            double const dx{ px[j] - pxi };
            double const dy{ py[j] - pyi };
            double const dz{ pz[j] - pzi };

            double const R_sq{ dx*dx + dy*dy + dz*dz + eps_sq };
            double const R_inv{ 1.0 / std::sqrt( R_sq ) };
            double const R_inv_cb{ R_inv*R_inv*R_inv };

            double const factor{ G * mass[j] * R_inv_cb };

            a_x += factor * dx;
            a_y += factor * dy;
            a_z += factor * dz;
        }

        if ( constant::ENABLE_PN && i != sun_idx ) {
            double const dx{ pxi - sun_x };
            double const dy{ pyi - sun_y };
            double const dz{ pzi - sun_z };

            double const R_sq{ dx*dx + dy*dy + dz*dz + eps_sq };
            double const R{ std::sqrt( R_sq ) };
            double const R_inv{ 1.0 / R };

            double const v_sq{ vxi*vxi + vyi*vyi + vzi*vzi };
            double const r_dot_v{ dx*vxi + dy*vyi + dz*vzi };

            double const coef_term{ GM_sun * R_inv * R_inv / c_sq };
            double const rad_term{ 4.0 * GM_sun * R_inv - v_sq };
            double const vel_term{ 4.0 * r_dot_v * R_inv };

            a_x += coef_term * ( rad_term * dx * R_inv + vel_term * vxi );
            a_y += coef_term * ( rad_term * dy * R_inv + vel_term * vyi );
            a_z += coef_term * ( rad_term * dz * R_inv + vel_term * vzi );   
        }

        ax[i] = a_x;
        ay[i] = a_y;
        az[i] = a_z;
    }
}