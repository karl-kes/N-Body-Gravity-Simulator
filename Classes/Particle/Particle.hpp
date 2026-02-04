#pragma once

#include <memory>
#include <iostream>

class Particles {
private:
    // Number of Particles:
    std::size_t num_particles_;

    // Mass:
    std::unique_ptr<double[]> mass_; // kg

    // Position, Velocity, and Acceleration of particle:
    std::unique_ptr<double[]> pos_x_, pos_y_, pos_z_; // m
    std::unique_ptr<double[]> vel_x_, vel_y_, vel_z_; // m/s
    std::unique_ptr<double[]> acc_x_, acc_y_, acc_z_; // m/s^2
    std::unique_ptr<double[]> old_acc_x_, old_acc_y_, old_acc_z_; // m/s^2

public:
    // Constructor, Destructor, Copy, Move:
    explicit Particles( std::size_t const num_particles );
    
    // Getters & Setters:
    // Num Particles:
    [[nodiscard]] std::size_t num_particles() const { return num_particles_; }

    // Mass
    [[nodiscard]] double mass( std::size_t const idx ) const { return mass_[idx]; }

    // Getters:
    // Positions:
    [[nodiscard]] double pos_x( std::size_t const idx ) const { return pos_x_[idx]; }
    [[nodiscard]] double pos_y( std::size_t const idx ) const { return pos_y_[idx]; }
    [[nodiscard]] double pos_z( std::size_t const idx ) const { return pos_z_[idx]; }

    // Velocities:
    [[nodiscard]] double vel_x( std::size_t const idx ) const { return vel_x_[idx]; }
    [[nodiscard]] double vel_y( std::size_t const idx ) const { return vel_y_[idx]; }
    [[nodiscard]] double vel_z( std::size_t const idx ) const { return vel_z_[idx]; }

    // Accelerations:
    [[nodiscard]] double acc_x( std::size_t const idx ) const { return acc_x_[idx]; }
    [[nodiscard]] double acc_y( std::size_t const idx ) const { return acc_y_[idx]; }
    [[nodiscard]] double acc_z( std::size_t const idx ) const { return acc_z_[idx]; }

    [[nodiscard]] double old_acc_x( std::size_t const idx ) const { return old_acc_x_[idx]; }
    [[nodiscard]] double old_acc_y( std::size_t const idx ) const { return old_acc_y_[idx]; }
    [[nodiscard]] double old_acc_z( std::size_t const idx ) const { return old_acc_z_[idx]; }

    // Pointer:
    [[nodiscard]] const std::unique_ptr<double[]> &pos_x() const { return pos_x_; }
    [[nodiscard]] const std::unique_ptr<double[]> &pos_y() const { return pos_y_; }
    [[nodiscard]] const std::unique_ptr<double[]> &pos_z() const { return pos_z_; }

    [[nodiscard]] const std::unique_ptr<double[]> &vel_x() const { return vel_x_; }
    [[nodiscard]] const std::unique_ptr<double[]> &vel_y() const { return vel_y_; }
    [[nodiscard]] const std::unique_ptr<double[]> &vel_z() const { return vel_z_; }

    [[nodiscard]] const std::unique_ptr<double[]> &acc_x() const { return acc_x_; }
    [[nodiscard]] const std::unique_ptr<double[]> &acc_y() const { return acc_y_; }
    [[nodiscard]] const std::unique_ptr<double[]> &acc_z() const { return acc_z_; }

    [[nodiscard]] const std::unique_ptr<double[]> &old_acc_x() const { return old_acc_x_; }
    [[nodiscard]] const std::unique_ptr<double[]> &old_acc_y() const { return old_acc_y_; }
    [[nodiscard]] const std::unique_ptr<double[]> &old_acc_z() const { return old_acc_z_; }

    [[nodiscard]] const std::unique_ptr<double[]> &mass() const { return mass_; }

    // Setters:
    // Mass:
    std::unique_ptr<double[]> &mass() { return mass_; }

    // Positions:
    std::unique_ptr<double[]> &pos_x() { return pos_x_; }
    std::unique_ptr<double[]> &pos_y() { return pos_y_; }
    std::unique_ptr<double[]> &pos_z() { return pos_z_; }

    // Velocities:
    std::unique_ptr<double[]> &vel_x() { return vel_x_; }
    std::unique_ptr<double[]> &vel_y() { return vel_y_; }
    std::unique_ptr<double[]> &vel_z() { return vel_z_; }

    // Accelerations:
    std::unique_ptr<double[]> &acc_x() { return acc_x_; }
    std::unique_ptr<double[]> &acc_y() { return acc_y_; }
    std::unique_ptr<double[]> &acc_z() { return acc_z_; }

    std::unique_ptr<double[]> &old_acc_x() { return old_acc_x_; }
    std::unique_ptr<double[]> &old_acc_y() { return old_acc_y_; }
    std::unique_ptr<double[]> &old_acc_z() { return old_acc_z_; }
};