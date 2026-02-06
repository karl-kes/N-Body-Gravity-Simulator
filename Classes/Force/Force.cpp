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
    constexpr double OMP_THRESHOLD{ constant::OMP_THRESHOLD };

    #pragma omp parallel for if( N >= OMP_THRESHOLD ) schedule( static )
    for ( std::size_t i = 0; i < N; ++i ) {
        double const pxi{ px[i] }, pyi{ py[i] }, pzi{ pz[i] };
        double const vxi{ vx[i] }, vyi{ vy[i] }, vzi{ vz[i] };
        double const mi{ mass[i] };
        
        double a_xi{}, a_yi{}, a_zi{};

        for ( std::size_t j = i + 1; j < N; ++j ) {
            double const dx{ px[j] - pxi };
            double const dy{ py[j] - pyi };
            double const dz{ pz[j] - pzi };

            double const R_sq{ dx*dx + dy*dy + dz*dz + eps_sq };
            double const R_inv{ 1.0 / std::sqrt( R_sq ) };
            double const R_inv_cb{ R_inv * R_inv * R_inv };

            double const mj{ mass[j] };

            double const G_mj_R_inv_cb{ G * mj * R_inv_cb };
            double const G_mi_R_inv_cb{ G * mi * R_inv_cb };

            double f_xi{ G_mj_R_inv_cb * dx };
            double f_yi{ G_mj_R_inv_cb * dy };
            double f_zi{ G_mj_R_inv_cb * dz };

            double f_xj{ -G_mi_R_inv_cb * dx };
            double f_yj{ -G_mi_R_inv_cb * dy };
            double f_zj{ -G_mi_R_inv_cb * dz };

            if ( constant::ENABLE_PN ) {
                double const vxj{ vx[j] }, vyj{ vy[j] }, vzj{ vz[j] };

                double const v_sq_i{ vxi*vxi + vyi*vyi + vzi*vzi };
                double const r_dot_v_i{ dx*vxi + dy*vyi + dz*vzi };
                double const coef_i{ G * mj * R_inv * R_inv * R_inv / c_sq };
                double const rad_term_i{ 4.0 * G * mj * R_inv - v_sq_i };
                double const vel_term_i{ 4.0 * r_dot_v_i * R_inv };

                f_xi += coef_i * ( rad_term_i * dx + vel_term_i * vxi );
                f_yi += coef_i * ( rad_term_i * dy + vel_term_i * vyi );
                f_zi += coef_i * ( rad_term_i * dz + vel_term_i * vzi );

                double const v_sq_j{ vxj*vxj + vyj*vyj + vzj*vzj };
                double const r_dot_v_j{ -dx*vxj - dy*vyj - dz*vzj };
                double const coef_j{ G * mi * R_inv * R_inv * R_inv / c_sq };
                double const rad_term_j{ 4.0 * G * mi * R_inv - v_sq_j };
                double const vel_term_j{ 4.0 * r_dot_v_j * R_inv };

                f_xj -= coef_j * ( rad_term_j * dx + vel_term_j * vxj );
                f_yj -= coef_j * ( rad_term_j * dy + vel_term_j * vyj );
                f_zj -= coef_j * ( rad_term_j * dz + vel_term_j * vzj );
            }

            a_xi += f_xi;
            a_yi += f_yi;
            a_zi += f_zi;

            #pragma omp atomic
            ax[j] += f_xj;
            #pragma omp atomic
            ay[j] += f_yj;
            #pragma omp atomic
            az[j] += f_zj;
        }

        #pragma omp atomic
        ax[i] += a_xi;
        #pragma omp atomic
        ay[i] += a_yi;
        #pragma omp atomic
        az[i] += a_zi;
    }
}