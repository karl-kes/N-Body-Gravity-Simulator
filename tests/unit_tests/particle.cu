// SoA memory layout: verify aligned, evenly-strided sub-arrays.

#include "test_harness.hpp"
#include "Particle/Particle.hpp"

#include <cstddef>


TEST( particle_soa_contiguous_memory ) {
    // Verify that SoA sub-arrays are evenly strided in memory.
    // With aligned allocation, the stride may be larger than N due to
    // SIMD-width padding (e.g. N=10 rounds up to 12 on AVX2).
    Particles p{ 10 };
    double* px{ p.pos_x() };
    double* py{ p.pos_y() };

    // pos_y should start exactly one stride after pos_x.
    std::ptrdiff_t const stride{ py - px };
    ASSERT_TRUE( stride >= 10 );

    // All subsequent arrays should be at the same stride interval.
    double* pz{ p.pos_z() };
    ASSERT_TRUE( pz - py == stride );

    // mass is sub-array index 12, so it should be 12 strides from pos_x.
    double* mass{ p.mass() };
    ASSERT_TRUE( mass == px + 12 * stride );
    ++g_pass;
}
