#include "Body.hpp"
#include "Constants.hpp"

Body::Body ( Vec_3D new_pos,
             Vec_3D new_vel,
             double new_mass ):
pos_( new_pos ),
vel_( new_vel ), 
acc_{ 0, 0, 0 }, 
old_acc_{ 0, 0, 0 }, 
mass_( new_mass ) {

}

// Calculates new acceleration based on forces from other bodies.
void Body::calculate_new_acc( std::vector<Body> const &other_bodies, std::size_t const &self_idx ) {
    set_old_acc( acc_ );
    Vec_3D total_acc{};

    for ( std::size_t idx = 0; idx < other_bodies.size(); ++idx ) {
        if ( idx == self_idx ) continue;
        
        Vec_3D R{ other_bodies[idx].get_pos() - get_pos() };
        double dist_squared{ R.norm_squared() + EPSILON * EPSILON };
        if ( dist_squared > 1e24 ) continue;
        double dist{ std::sqrt( dist_squared ) };

        // Acceleration from law of gravitation: a = R_vector * (GM / r^3) 
        total_acc += R * ( ( G * other_bodies[idx].get_mass() ) / ( dist * dist * dist ) );
    }
    set_acc( total_acc );
}

// Updates body.
void Body::update( double const &dt ) {
    set_pos( get_pos() + ( get_vel() * dt + get_acc() * ( 0.5 * dt * dt ) ) );
    set_vel( get_vel() + ( ( get_acc() + get_old_acc() ) * ( 0.5 * dt ) ) );
}

// Body getters:
const Vec_3D &Body::get_pos() const { 
    return pos_;
}
const Vec_3D &Body::get_vel() const {
    return vel_;
}
const Vec_3D &Body::get_acc() const {
    return acc_;
}
const Vec_3D &Body::get_old_acc() const {
    return old_acc_;
}
double Body::get_mass() const {
    return mass_;
}

// Body setters:
void Body::set_pos( Vec_3D new_pos ) {
    pos_ = new_pos;
}
void Body::set_vel( Vec_3D new_vel ) {
    vel_ = new_vel;
}
void Body::set_acc( Vec_3D new_acc ) {
    acc_ = new_acc;
}
void Body::set_old_acc( Vec_3D new_old_acc ) {
    old_acc_ = new_old_acc;
}