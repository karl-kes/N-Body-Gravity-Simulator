#pragma once

#include "../Particle/Particle.hpp"
#include "../Config.hpp"

#include <vector>
#include <cmath>
#include <cstddef>

#if defined(__GNUC__) || defined(__clang__)
    #define RESTRICT __restrict__
#elif defined(_MSC_VER)
    #define RESTRICT __restrict
#else
    #define RESTRICT
#endif

class Force {
public:
    virtual ~Force() = default;
    virtual void apply( Particles &particles ) const = 0;
};

class Gravity : public Force {
public:
    explicit Gravity( double const eps = 1e-9 );

    static inline void accumulate_pairwise (
        double const pxi, double const pyi, double const pzi,
        double const pxj, double const pyj, double const pzj, double const mj,
        double &a_xi, double &a_yi, double &a_zi,
        double const G, double const eps_sq, double const mask ) {

        double const dx{ pxj - pxi };
        double const dy{ pyj - pyi };
        double const dz{ pzj - pzi };

        double const R_sq{ dx*dx + dy*dy + dz*dz + eps_sq };
        double const R_inv{ 1.0 / std::sqrt( R_sq ) };
        double const R_inv_cb{ R_inv * R_inv * R_inv };

        double const G_mj_R_inv_cb{ G * mj * R_inv_cb };

        double const f_x{ G_mj_R_inv_cb * dx };
        double const f_y{ G_mj_R_inv_cb * dy };
        double const f_z{ G_mj_R_inv_cb * dz };

        a_xi += mask * f_x;
        a_yi += mask * f_y;
        a_zi += mask * f_z;
    }

    void apply( Particles &particles ) const override;

    [[nodiscard]] double eps() const { return eps_; }

private:
    double eps_;
};