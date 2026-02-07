#include "Particle.hpp"

Particles::Particles( std::size_t const num_particles )
: num_particles_{ num_particles } {
    auto allocate = [size = num_particles]() { return std::make_unique<double[]>( size ); };

    mass_ = allocate();

    pos_x_ = allocate();
    pos_y_ = allocate();
    pos_z_ = allocate();

    vel_x_ = allocate();
    vel_y_ = allocate();
    vel_z_ = allocate();

    acc_x_ = allocate();
    acc_y_ = allocate();
    acc_z_ = allocate();

    old_acc_x_ = allocate();
    old_acc_y_ = allocate();
    old_acc_z_ = allocate();

    potential_energy_ = allocate();
}