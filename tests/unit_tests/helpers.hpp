#pragma once

// Shared test utilities: conservation diagnostics, IC samplers, force comparison.
// All functions are inline so multiple test TUs can include this header without
// linker conflicts.

#include "Particle/Particle.hpp"
#include "Force/Force.hpp"
#include "Force/BarnesHut.hpp"
#include "Integrator/Integrator.hpp"
#include "Config.hpp"

#include <cmath>
#include <vector>
#include <memory>
#include <random>
#include <algorithm>
#include <numbers>
#include <cstddef>

inline double kinetic_energy( Particles const &p ) {
    double ke{};
    for ( std::size_t i{}; i < p.num_particles(); ++i ) {
        double const v_sq{ p.vel_x()[i]*p.vel_x()[i]
                         + p.vel_y()[i]*p.vel_y()[i]
                         + p.vel_z()[i]*p.vel_z()[i] };
        ke += 0.5 * p.mass()[i] * v_sq;
    }
    return ke;
}

inline double potential_energy( Particles const &p ) {
    double pe{};
    std::size_t const N{ p.num_particles() };
    for ( std::size_t i{}; i < N; ++i ) {
        for ( std::size_t j{ i + 1 }; j < N; ++j ) {
            double const dx{ p.pos_x()[j] - p.pos_x()[i] };
            double const dy{ p.pos_y()[j] - p.pos_y()[i] };
            double const dz{ p.pos_z()[j] - p.pos_z()[i] };
            double const r{ std::sqrt( dx*dx + dy*dy + dz*dz ) };
            pe -= config::G * p.mass()[i] * p.mass()[j] / r;
        }
    }
    return pe;
}

inline double total_energy( Particles const &p ) {
    return kinetic_energy( p ) + potential_energy( p );
}

inline double linear_momentum( Particles const &p ) {
    double px{}, py{}, pz{};
    for ( std::size_t i{}; i < p.num_particles(); ++i ) {
        px += p.mass()[i] * p.vel_x()[i];
        py += p.mass()[i] * p.vel_y()[i];
        pz += p.mass()[i] * p.vel_z()[i];
    }
    return std::sqrt( px*px + py*py + pz*pz );
}

inline double angular_momentum( Particles const &p ) {
    double lx{}, ly{}, lz{};
    for ( std::size_t i{}; i < p.num_particles(); ++i ) {
        double const m{ p.mass()[i] };
        lx += m * ( p.pos_y()[i]*p.vel_z()[i] - p.pos_z()[i]*p.vel_y()[i] );
        ly += m * ( p.pos_z()[i]*p.vel_x()[i] - p.pos_x()[i]*p.vel_z()[i] );
        lz += m * ( p.pos_x()[i]*p.vel_y()[i] - p.pos_y()[i]*p.vel_x()[i] );
    }
    return std::sqrt( lx*lx + ly*ly + lz*lz );
}

// Circular two-body Kepler orbit, placed at the barycenter with circular velocity.
inline void setup_two_body( Particles &p, double const M, double const m, double const r ) {
    double const mu{ m / ( M + m ) };
    double const r0{ -mu * r };
    double const r1{ ( 1.0 - mu ) * r };

    double const v_orb{ std::sqrt( config::G * ( M + m ) / r ) };
    double const v0{ -mu * v_orb };
    double const v1{ ( 1.0 - mu ) * v_orb };

    p.pos_x()[0] = r0;  p.pos_y()[0] = 0.0; p.pos_z()[0] = 0.0;
    p.vel_x()[0] = 0.0; p.vel_y()[0] = v0;  p.vel_z()[0] = 0.0;
    p.mass()[0]  = M;

    p.pos_x()[1] = r1;  p.pos_y()[1] = 0.0; p.pos_z()[1] = 0.0;
    p.vel_x()[1] = 0.0; p.vel_y()[1] = v1;  p.vel_z()[1] = 0.0;
    p.mass()[1]  = m;

    for ( std::size_t i{}; i < 2; ++i ) {
        p.acc_x()[i] = 0.0;
        p.acc_y()[i] = 0.0;
        p.acc_z()[i] = 0.0;
    }
}

inline void step_n( Particles &p, Integrator &integ,
                    std::vector<std::unique_ptr<Force>> &forces, std::size_t const n ) {
    for ( std::size_t s{}; s < n; ++s ) {
        integ.integrate( p, forces );
    }
}


// Plummer sphere position sampler: r = a / sqrt(u^(-2/3) - 1), isotropic angles.
inline void populate_plummer( Particles &p, std::size_t const N,
                              double const a, double const M_tot,
                              unsigned const seed ) {
    std::mt19937_64 rng{ seed };
    std::uniform_real_distribution<double> u_dist{ 0.0, 1.0 };
    double const m_per{ M_tot / static_cast<double>( N ) };

    for ( std::size_t i{}; i < N; ++i ) {
        double const u{ std::clamp( u_dist( rng ), 1e-6, 1.0 - 1e-6 ) };
        double const r{ a / std::sqrt( std::pow( u, -2.0 / 3.0 ) - 1.0 ) };

        double const cos_theta{ 2.0 * u_dist( rng ) - 1.0 };
        double const sin_theta{ std::sqrt( 1.0 - cos_theta * cos_theta ) };
        double const phi{ 2.0 * std::numbers::pi * u_dist( rng ) };

        p.pos_x()[i] = r * sin_theta * std::cos( phi );
        p.pos_y()[i] = r * sin_theta * std::sin( phi );
        p.pos_z()[i] = r * cos_theta;
        p.vel_x()[i] = 0.0; p.vel_y()[i] = 0.0; p.vel_z()[i] = 0.0;
        p.acc_x()[i] = 0.0; p.acc_y()[i] = 0.0; p.acc_z()[i] = 0.0;
        p.mass()[i]  = m_per;
    }
}

// Uniform random positions in [-L, L]^3 with equal mass per particle.
inline void populate_uniform_cube( Particles &p, std::size_t const N,
                                   double const L, double const M_tot,
                                   unsigned const seed ) {
    std::mt19937_64 rng{ seed };
    std::uniform_real_distribution<double> dist{ -L, L };
    double const m_per{ M_tot / static_cast<double>( N ) };

    for ( std::size_t i{}; i < N; ++i ) {
        p.pos_x()[i] = dist( rng );
        p.pos_y()[i] = dist( rng );
        p.pos_z()[i] = dist( rng );
        p.vel_x()[i] = 0.0; p.vel_y()[i] = 0.0; p.vel_z()[i] = 0.0;
        p.acc_x()[i] = 0.0; p.acc_y()[i] = 0.0; p.acc_z()[i] = 0.0;
        p.mass()[i]  = m_per;
    }
}

struct ForceErrorStats {
    double mean_rel;
    double max_rel;
};

// Apply direct and BH forces to identical configurations, return per-particle
// relative force error |a_BH - a_direct| / |a_direct|, aggregated as mean and max.
inline ForceErrorStats compare_forces( Particles &p_direct, Particles &p_bh,
                                       double const theta ) {
    Gravity{}.apply( p_direct );
    Gravity_BarnesHut{ theta }.apply( p_bh );

    std::size_t const N{ p_direct.num_particles() };
    double sum_rel{};
    double max_rel{};
    std::size_t counted{};

    for ( std::size_t i{}; i < N; ++i ) {
        double const ax{ p_direct.acc_x()[i] };
        double const ay{ p_direct.acc_y()[i] };
        double const az{ p_direct.acc_z()[i] };
        double const a_mag{ std::sqrt( ax*ax + ay*ay + az*az ) };
        if ( a_mag == 0.0 ) continue;

        double const dx{ p_bh.acc_x()[i] - ax };
        double const dy{ p_bh.acc_y()[i] - ay };
        double const dz{ p_bh.acc_z()[i] - az };
        double const err{ std::sqrt( dx*dx + dy*dy + dz*dz ) };

        double const rel{ err / a_mag };
        sum_rel += rel;
        max_rel = std::max( max_rel, rel );
        ++counted;
    }

    return { sum_rel / static_cast<double>( counted ), max_rel };
}
