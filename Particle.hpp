#pragma once

class Particle {
private:
    // Mass:
    double mass_; // kg

    // Position, Velocity, and Acceleration of particle:
    double pos_[3]{}; // m
    double vel_[3]{}; // m/s
    double acc_[3]{}; // m/s^2

public:
    // Constructor:
    explicit Particle(
        double mass,
        double init_pos_x = 0.0, double init_pos_y = 0.0, double init_pos_z = 0.0,
        double init_vel_x = 0.0, double init_vel_y = 0.0, double init_vel_z = 0.0,
        double init_acc_x = 0.0, double init_acc_y = 0.0, double init_acc_z = 0.0 );

    // Getters & Setters:
    // Mass
    [[nodiscard]] double mass() const { return mass_; }

    // Getters:
    // Positions:
    [[nodiscard]] double pos_x() const { return pos_[0]; }
    [[nodiscard]] double pos_y() const { return pos_[1]; }
    [[nodiscard]] double pos_z() const { return pos_[2]; }

    // Velocities:
    [[nodiscard]] double vel_x() const { return vel_[0]; }
    [[nodiscard]] double vel_y() const { return vel_[1]; }
    [[nodiscard]] double vel_z() const { return vel_[2]; }

    // Accelerations:
    [[nodiscard]] double acc_x() const { return acc_[0]; }
    [[nodiscard]] double acc_y() const { return acc_[1]; }
    [[nodiscard]] double acc_z() const { return acc_[2]; }

    // Pointers:
    [[nodiscard]] const double *pos() const { return pos_; }
    [[nodiscard]] const double *vel() const { return vel_; }
    [[nodiscard]] const double *acc() const { return acc_; }

    // Setters:
    // Positions:
    double &pos_x() { return pos_[0]; }
    double &pos_y() { return pos_[1]; }
    double &pos_z() { return pos_[2]; }

    // Velocities:
    double &vel_x() { return vel_[0]; }
    double &vel_y() { return vel_[1]; }
    double &vel_z() { return vel_[2]; }

    // Accelerations:
    double &acc_x() { return acc_[0]; }
    double &acc_y() { return acc_[1]; }
    double &acc_z() { return acc_[2]; }

    // Pointers:
    double *pos() noexcept { return pos_; }
    double *vel() noexcept { return vel_; }
    double *acc() noexcept { return acc_; }

    // Helpers:
    void clear_acc() { acc_x() = 0; acc_y() = 0; acc_z() = 0; };
};