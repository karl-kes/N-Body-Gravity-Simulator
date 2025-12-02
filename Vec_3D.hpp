#pragma once

#include <cmath>

class Vec_3D {
private:
    double x_;
    double y_;
    double z_;

public:
    Vec_3D ( double new_x = 0.0,
             double new_y = 0.0, 
             double new_z = 0.0 );

    // Getters:
    const double &get_x() const;
    const double &get_y() const;
    const double &get_z() const;

    // Setters:
    void set_x( double new_x );
    void set_y( double new_y );
    void set_z( double new_z );

    // Helper functions:
    double norm() const;
    double norm_squared() const;

    // Operator overloads:
    Vec_3D operator+( Vec_3D const &other_vec ) const;
    Vec_3D operator-( Vec_3D const &other_vec ) const;
    Vec_3D operator*( double const &constant ) const;
    Vec_3D &operator+=( Vec_3D const &other_vec );
};

