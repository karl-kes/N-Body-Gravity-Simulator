#include "Force.hpp"

#define RESTRICT __restrict

Gravity::Gravity()
{ }

void INLINE Gravity::compute_forces (
        double const pxi, double const pyi, double const pzi, double const mi,
        double const vxi, double const vyi, double const vzi,
        double const pxj, double const pyj, double const pzj, double const mj,
        double const vxj, double const vyj, double const vzj,
        double &a_xi, double &a_yi, double &a_zi,
        double const G, double const eps_sq, double const c_sq, double const mask
    ) {

    double const dx{ pxj - pxi };
    double const dy{ pyj - pyi };
    double const dz{ pzj - pzi };

    double const R_sq{ dx*dx + dy*dy + dz*dz + eps_sq };
    double const R_inv{ 1.0 / std::sqrt( R_sq ) };
    double const R_inv_cb{ R_inv * R_inv * R_inv };

    double const G_mj_R_inv_cb{ G * mj * R_inv_cb };

    double f_x{ G_mj_R_inv_cb * dx };
    double f_y{ G_mj_R_inv_cb * dy };
    double f_z{ G_mj_R_inv_cb * dz };

    if ( constant::ENABLE_PN ) {
        double const v_sq_i{ vxi*vxi + vyi*vyi + vzi*vzi };
        double const r_dot_v_i{ dx*vxi + dy*vyi + dz*vzi };
        double const coef_i{ G * mj * R_inv * R_inv * R_inv / c_sq };
        double const rad_term_i{ 4.0 * G * mj * R_inv - v_sq_i };
        double const vel_term_i{ 4.0 * r_dot_v_i * R_inv };

        f_x += coef_i * ( rad_term_i * dx + vel_term_i * vxi );
        f_y += coef_i * ( rad_term_i * dy + vel_term_i * vyi );
        f_z += coef_i * ( rad_term_i * dz + vel_term_i * vzi );
    }

    a_xi += mask * f_x;
    a_yi += mask * f_y;
    a_zi += mask * f_z;
}

void Gravity::apply( Particles &particles ) const {
    std::size_t const N{ particles.num_particles() };

    double const* RESTRICT px{ particles.pos_x().get() };
    double const* RESTRICT py{ particles.pos_y().get() };
    double const* RESTRICT pz{ particles.pos_z().get() };

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
    constexpr double OMP_THRESHOLD{ constant::OMP_THRESHOLD };

    auto apply_kernel = [=]( std::size_t i ) {
        double const pxi{ px[i] }, pyi{ py[i] }, pzi{ pz[i] };
        double const vxi{ vx[i] }, vyi{ vy[i] }, vzi{ vz[i] };
        double const mi{ mass[i] };
        
        double a_xi{}, a_yi{}, a_zi{};

        #pragma omp simd reduction( +:a_xi, a_yi, a_zi )
        for ( std::size_t j = 0; j < N; ++j ) {
            double const mask{ ( i == j ) ? 0.0 : 1.0 };

            compute_forces (
                pxi, pyi, pzi, mi, vxi, vyi, vzi,
                px[j], py[j], pz[j], mass[j],
                vx[j], vy[j], vz[j],
                a_xi, a_yi, a_zi,
                G, eps_sq, c_sq, mask
            );
        }

        ax[i] += a_xi;
        ay[i] += a_yi;
        az[i] += a_zi;
    };

    if ( N >= OMP_THRESHOLD ) {
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