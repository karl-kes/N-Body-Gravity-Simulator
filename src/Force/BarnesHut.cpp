#include "BarnesHut.hpp"

#include <algorithm>
#include <cmath>

#include <omp.h>

Gravity_BarnesHut::Gravity_BarnesHut( double const theta,
                                      std::size_t const leaf_bucket,
                                      double const eps )
: theta_{ theta }
, leaf_bucket_{ leaf_bucket }
, eps_{ eps }
{ }


void Gravity_BarnesHut::build_tree( Particles const &particles ) const {
    std::size_t const N{ particles.num_particles() };

    double const* RESTRICT px{ particles.pos_x() };
    double const* RESTRICT py{ particles.pos_y() };
    double const* RESTRICT pz{ particles.pos_z() };
    double const* RESTRICT mass{ particles.mass() };

    nodes_.clear();
    indices_.resize( N );
    for ( std::size_t i{}; i < N; ++i ) indices_[i] = static_cast<int>( i );
    scratch_.resize( N );

    if ( N == 0 ) return;

    // Cubic root box centered on the AABB midpoint, sized to enclose every
    // body. A cube keeps octant subdivisions geometrically clean.
    double min_x{ px[0] }, max_x{ px[0] };
    double min_y{ py[0] }, max_y{ py[0] };
    double min_z{ pz[0] }, max_z{ pz[0] };
    for ( std::size_t i{ 1 }; i < N; ++i ) {
        min_x = std::min( min_x, px[i] ); max_x = std::max( max_x, px[i] );
        min_y = std::min( min_y, py[i] ); max_y = std::max( max_y, py[i] );
        min_z = std::min( min_z, pz[i] ); max_z = std::max( max_z, pz[i] );
    }
    double const cx{ 0.5 * ( min_x + max_x ) };
    double const cy{ 0.5 * ( min_y + max_y ) };
    double const cz{ 0.5 * ( min_z + max_z ) };
    double half{ 0.5 * std::max( { max_x - min_x, max_y - min_y, max_z - min_z } ) };
    if ( half <= 0.0 ) half = 1.0;  // All particles coincident; degenerate but well-defined.
    half *= 1.0 + 1e-12;            // Pad so positions on the boundary fall inside.

    nodes_.emplace_back();
    build_recursive( 0, 0, static_cast<int>( N ),
                     cx, cy, cz, half, 0,
                     px, py, pz, mass );
}


void Gravity_BarnesHut::build_recursive(
    int const node_id,
    int const begin, int const end,
    double const bcx, double const bcy, double const bcz,
    double const half, int const depth,
    double const* RESTRICT px, double const* RESTRICT py, double const* RESTRICT pz,
    double const* RESTRICT mass ) const
{
    int const count{ end - begin };

    {
        BHNode &n{ nodes_[node_id] };
        n.box_cx = bcx; n.box_cy = bcy; n.box_cz = bcz;
        n.half_width = half;
        for ( int k{}; k < 8; ++k ) n.children[k] = -1;
    }

    if ( count <= static_cast<int>( leaf_bucket_ ) || depth >= MAX_DEPTH ) {
        double m_sum{};
        double cx_sum{}, cy_sum{}, cz_sum{};
        for ( int k{ begin }; k < end; ++k ) {
            int const i{ indices_[k] };
            double const m{ mass[i] };
            m_sum  += m;
            cx_sum += m * px[i];
            cy_sum += m * py[i];
            cz_sum += m * pz[i];
        }
        BHNode &n{ nodes_[node_id] };
        n.total_mass = m_sum;
        if ( m_sum > 0.0 ) {
            n.com_x = cx_sum / m_sum;
            n.com_y = cy_sum / m_sum;
            n.com_z = cz_sum / m_sum;
        } else {
            n.com_x = bcx; n.com_y = bcy; n.com_z = bcz;
        }
        n.children[0] = -2;
        n.children[1] = begin;
        n.children[2] = count;
        return;
    }

    // Counting-sort indices_[begin..end) into 8 octants.
    // Octant index = (x>=cx) | ((y>=cy)<<1) | ((z>=cz)<<2).
    int counts[8]{};
    for ( int k{ begin }; k < end; ++k ) {
        int const i{ indices_[k] };
        int const oct{ ( px[i] >= bcx ? 1 : 0 )
                     | ( py[i] >= bcy ? 2 : 0 )
                     | ( pz[i] >= bcz ? 4 : 0 ) };
        ++counts[oct];
    }

    int starts[8]{};
    int offsets[8]{};
    starts[0] = begin;
    for ( int k{ 1 }; k < 8; ++k ) starts[k] = starts[k-1] + counts[k-1];
    for ( int k{}; k < 8; ++k ) offsets[k] = starts[k];

    for ( int k{ begin }; k < end; ++k ) {
        int const i{ indices_[k] };
        int const oct{ ( px[i] >= bcx ? 1 : 0 )
                     | ( py[i] >= bcy ? 2 : 0 )
                     | ( pz[i] >= bcz ? 4 : 0 ) };
        scratch_[ offsets[oct]++ ] = i;
    }
    for ( int k{ begin }; k < end; ++k ) indices_[k] = scratch_[k];

    double const child_half{ 0.5 * half };
    for ( int k{}; k < 8; ++k ) {
        if ( counts[k] == 0 ) continue;

        double const ccx{ bcx + ( ( k & 1 ) ? child_half : -child_half ) };
        double const ccy{ bcy + ( ( k & 2 ) ? child_half : -child_half ) };
        double const ccz{ bcz + ( ( k & 4 ) ? child_half : -child_half ) };

        // Reserve child node id BEFORE the recursive call. Re-access parent
        // through nodes_[node_id] AFTER the call, since emplace_back may
        // have reallocated the underlying buffer.
        int const child_id{ static_cast<int>( nodes_.size() ) };
        nodes_.emplace_back();
        nodes_[node_id].children[k] = child_id;

        build_recursive( child_id,
                         starts[k], starts[k] + counts[k],
                         ccx, ccy, ccz, child_half, depth + 1,
                         px, py, pz, mass );
    }

    double m_sum{};
    double cx_sum{}, cy_sum{}, cz_sum{};
    for ( int k{}; k < 8; ++k ) {
        int const child_id{ nodes_[node_id].children[k] };
        if ( child_id < 0 ) continue;
        BHNode const &c{ nodes_[child_id] };
        m_sum  += c.total_mass;
        cx_sum += c.total_mass * c.com_x;
        cy_sum += c.total_mass * c.com_y;
        cz_sum += c.total_mass * c.com_z;
    }
    BHNode &n{ nodes_[node_id] };
    n.total_mass = m_sum;
    if ( m_sum > 0.0 ) {
        n.com_x = cx_sum / m_sum;
        n.com_y = cy_sum / m_sum;
        n.com_z = cz_sum / m_sum;
    } else {
        n.com_x = bcx; n.com_y = bcy; n.com_z = bcz;
    }
}


void Gravity_BarnesHut::traverse_for_particle(
    std::size_t const i,
    double const* RESTRICT px, double const* RESTRICT py, double const* RESTRICT pz,
    double const* RESTRICT mass,
    double &a_xi, double &a_yi, double &a_zi ) const
{
    double const eps_sq{ eps_ * eps_ };
    constexpr double G{ config::G };

    double const theta_sq{ theta_ * theta_ };

    double const pxi{ px[i] }, pyi{ py[i] }, pzi{ pz[i] };

    // A balanced octree at N=131k has depth ~6, so 256 stack frames covers
    // MAX_DEPTH=32 with a wide margin (worst case 7 siblings per level).
    int stack[256];
    int sp{ 0 };
    stack[sp++] = 0;

    while ( sp > 0 ) {
        int const node_id{ stack[--sp] };
        BHNode const &n{ nodes_[node_id] };

        if ( n.children[0] == -2 ) {
            // Leaf bucket: sum each particle, masking the self-term.
            int const first{ n.children[1] };
            int const cnt{ n.children[2] };
            for ( int b{}; b < cnt; ++b ) {
                int const j{ indices_[ first + b ] };
                double const mask{ ( static_cast<std::size_t>( j ) == i ) ? 0.0 : 1.0 };
                Gravity::accumulate_pairwise(
                    pxi, pyi, pzi,
                    px[j], py[j], pz[j], mass[j],
                    a_xi, a_yi, a_zi,
                    G, eps_sq, mask
                );
            }
            continue;
        }

        // MAC: open if (2 * half_width)^2 >= theta^2 * d^2, else approximate.
        double const dx{ n.com_x - pxi };
        double const dy{ n.com_y - pyi };
        double const dz{ n.com_z - pzi };
        double const d_sq{ dx*dx + dy*dy + dz*dz };
        double const s{ 2.0 * n.half_width };

        if ( s * s < theta_sq * d_sq ) {
            Gravity::accumulate_pairwise(
                pxi, pyi, pzi,
                n.com_x, n.com_y, n.com_z, n.total_mass,
                a_xi, a_yi, a_zi,
                G, eps_sq, 1.0
            );
        } else {
            for ( int k{}; k < 8; ++k ) {
                int const child_id{ n.children[k] };
                if ( child_id >= 0 ) stack[sp++] = child_id;
            }
        }
    }
}


void Gravity_BarnesHut::apply( Particles &particles ) const {
    std::size_t const N{ particles.num_particles() };
    if ( N == 0 ) return;

    build_tree( particles );

    double const* RESTRICT px{ particles.pos_x() };
    double const* RESTRICT py{ particles.pos_y() };
    double const* RESTRICT pz{ particles.pos_z() };
    double const* RESTRICT mass{ particles.mass() };

    double* RESTRICT ax{ particles.acc_x() };
    double* RESTRICT ay{ particles.acc_y() };
    double* RESTRICT az{ particles.acc_z() };

    if ( N >= config::OMP_THRESHOLD ) {
        #pragma omp parallel for schedule( dynamic, 32 )
        for ( std::size_t i = 0; i < N; ++i ) {
            double a_xi{}, a_yi{}, a_zi{};
            traverse_for_particle( i, px, py, pz, mass, a_xi, a_yi, a_zi );
            ax[i] += a_xi;
            ay[i] += a_yi;
            az[i] += a_zi;
        }
    } else {
        for ( std::size_t i = 0; i < N; ++i ) {
            double a_xi{}, a_yi{}, a_zi{};
            traverse_for_particle( i, px, py, pz, mass, a_xi, a_yi, a_zi );
            ax[i] += a_xi;
            ay[i] += a_yi;
            az[i] += a_zi;
        }
    }
}