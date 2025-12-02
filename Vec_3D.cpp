#include "Vec_3D.hpp"

// Constructor
Vec_3D::Vec_3D ( double new_x, 
                 double new_y, 
                 double new_z ):
x_{ new_x },
y_{ new_y },
z_{ new_z } {

}

// Getters:
const double &Vec_3D::get_x() const {
    return x_;
}
const double &Vec_3D::get_y() const {
    return y_;
}
const double &Vec_3D::get_z() const {
    return z_;
}

// Setters:
void Vec_3D::set_x( double new_x ) {
    x_ = new_x;
}
void Vec_3D::set_y( double new_y ) {
    y_ = new_y;
}
void Vec_3D::set_z( double new_z ) {
    z_ = new_z;
}

// Helper functions:
double Vec_3D::norm() const {
    return std::sqrt( get_x()*get_x() + get_y()*get_y() + get_z()*get_z() );
}
double Vec_3D::norm_squared() const {
    return get_x()*get_x() + get_y()*get_y() + get_z()*get_z();
}

// Operator overloads:
Vec_3D Vec_3D::operator+( Vec_3D const &other_vec ) const {
    return { get_x() + other_vec.get_x(), get_y() + other_vec.get_y(), get_z() + other_vec.get_z() };
}
Vec_3D Vec_3D::operator-( Vec_3D const &other_vec ) const {
    return { get_x() - other_vec.get_x(), get_y() - other_vec.get_y(), get_z() - other_vec.get_z() };
}
Vec_3D Vec_3D::operator*( double const &constant ) const {
    return { get_x()*constant, get_y()*constant, get_z()*constant };
}
Vec_3D &Vec_3D::operator+=( Vec_3D const &other_vec ) {
    set_x( get_x() + other_vec.get_x() );
    set_y( get_y() + other_vec.get_y() );
    set_z( get_z() + other_vec.get_z() );
    return *this;
}