#pragma once

#include "../../Constants.hpp"
#include "../Vec_3D/Vec_3D.hpp"
#include <vector>

class Body {
private:
    Vec_3D pos_;
    Vec_3D vel_;
    Vec_3D acc_;
    Vec_3D old_acc_;
    double mass_;

public:
    Body ( Vec_3D new_pos = 0.0, Vec_3D new_vel = 0.0, double new_mass = 0.0 );
    
    // Calculates new acceleration based on forces from other bodies.
    void calculate_new_acc( std::vector<Body> const &other_bodies, std::size_t const &self_idx );

    // Updates body.
    void update( double const &dt );

    // Body getters:
    const Vec_3D &get_pos() const;
    const Vec_3D &get_vel() const;
    const Vec_3D &get_acc() const;
    const Vec_3D &get_old_acc() const;
    double get_mass() const;

    // Body setters:
    void set_pos( Vec_3D const &new_pos );
    void set_vel( Vec_3D const &new_vel );
    void set_acc( Vec_3D const &new_acc );
    void set_old_acc( Vec_3D const &new_old_acc );
};