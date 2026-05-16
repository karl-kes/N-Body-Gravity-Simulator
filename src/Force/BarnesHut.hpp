#pragma once

#include "Force.hpp"

#include <vector>
#include <cstddef>

// One node of the Barnes-Hut octree.
// Internal nodes: children[k] is a node index for octant k, or -1 if empty.
// Leaf nodes:     children[0] = -2 (sentinel), children[1] = first index into
//                 indices_, children[2] = bucket count.
// Disambiguation: a node is a leaf iff children[0] == -2.
// (-2 is distinct from the "empty child" sentinel -1; using -1 for both
//  would cause internal nodes with empty octant 0 to be misread as leaves.)
struct BHNode {
    double com_x, com_y, com_z;
    double total_mass;
    double box_cx, box_cy, box_cz;
    double half_width;
    int children[8];
};

class Gravity_BarnesHut : public Force {
public:
    explicit Gravity_BarnesHut( double const theta = 0.5,
                                std::size_t const leaf_bucket = 8 );

    void apply( Particles &particles ) const override;

    [[nodiscard]] double theta() const { return theta_; }
    [[nodiscard]] std::size_t leaf_bucket() const { return leaf_bucket_; }

private:
    double theta_;
    std::size_t leaf_bucket_;

    // Tree state. Rebuilt every apply(); capacity is sticky to avoid
    // re-allocation across integration steps.
    mutable std::vector<BHNode> nodes_;
    mutable std::vector<int> indices_;
    mutable std::vector<int> scratch_;

    // Hard recursion cap. Degenerate clustered input collapses into a single
    // bucket leaf at this depth; the only correctness fallback in the build.
    static constexpr int MAX_DEPTH{ 32 };

    void build_tree( Particles const &particles ) const;

    // Always re-access nodes via nodes_[id]; never hold a reference across
    // recursive calls, since vector growth invalidates references.
    void build_recursive( int const node_id,
                          int const begin, int const end,
                          double const bcx, double const bcy, double const bcz,
                          double const half, int const depth,
                          double const* RESTRICT px, double const* RESTRICT py, double const* RESTRICT pz,
                          double const* RESTRICT mass ) const;

    void traverse_for_particle( std::size_t const i,
                                double const* RESTRICT px, double const* RESTRICT py, double const* RESTRICT pz,
                                double const* RESTRICT mass,
                                double &a_xi, double &a_yi, double &a_zi ) const;
};
