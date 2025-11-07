#include <iostream>
#include <cassert>
#include <cmath>
#include <iomanip>

double G{ 6.67e-11 };

class Vec_3D {
public:
    double x_, y_, z_;
    double magnitude() {
        return std::sqrt( x_*x_ + y_*y_ + z_*z_ );
    }
};

Vec_3D operator+( Vec_3D const &vec, Vec_3D const &other ) {
    return{ vec.x_ + other.x_, vec.y_ + other.y_, vec.z_ + other.z_ };
}

Vec_3D operator-( Vec_3D const &vec, Vec_3D const &other ) {
    return{ vec.x_ - other.x_, vec.y_ - other.y_, vec.z_ - other.z_ };
}

Vec_3D operator*( Vec_3D const &vec, double const &constant ) {
    return{ vec.x_ * constant, vec.y_ * constant, vec.z_ * constant };
}

Vec_3D &operator+=( Vec_3D &vec, Vec_3D const &other ) {
    vec.x_ += other.x_;
    vec.y_ += other.y_;
    vec.z_ += other.z_;
    return vec;
}

class Body {
private:
    Vec_3D pos_;
    Vec_3D vel_;
    Vec_3D acc_;
    double mass_;

public:
    Body ( Vec_3D pos, Vec_3D vel, double mass ) 
         : pos_( pos ), vel_( vel ), acc_{ 0, 0, 0 }, mass_( mass ) {}
    
    Vec_3D calculate_force( Body const &other ) {
        Vec_3D R{ other.pos_ - pos_ };
        double dist{ R.magnitude() };

        if ( dist < 1.0e-10 ) return{ 0, 0, 0 };
        double force_mag{ G*mass_*other.mass_ / ( dist * dist ) };

        return R * ( force_mag / dist );
    }

    void reset_acceleration() { acc_ = { 0, 0, 0 }; }
    void apply_force( Vec_3D const &force ) {
        acc_ += ( force * ( 1.0 / mass_ ) );
    }

    void update( double dt ) {
        vel_ += acc_ * dt;
        pos_ += vel_ * dt;
    }

    Vec_3D get_pos() const { return pos_; }
    Vec_3D get_vel() const { return vel_; }
    double get_mass() const { return mass_; }
};

int main() {
    static constexpr double MASS{ 1.0e20 };
    static constexpr double POS{ 1.0e4 };
    static constexpr double dt{ 0.1 };
    int steps{ 0 };
    int current_step{ 0 };

    Body B_1{ { -1.0*POS, POS, POS }, { 0, 0, 0 }, MASS };
    Body B_2{ { POS, -1.0*POS, POS }, { 0, 0, 0 }, MASS };
    Body B_3{ { POS, POS, -1.0*POS }, { 0, 0, 0 }, MASS };

    std::cout << "<--- 3-Body Simulation --->" << std::endl;

    std::cout << "\nNumber of simulation steps: ";
    std::cin >> steps;

    std::cout << "\nStarting 3-Body Simulation..." << std::endl;
    std::cout << std::fixed << std::setprecision(2);

    for( steps; steps > 0; --steps ) {
        Vec_3D force_1_2 = B_1.calculate_force( B_2 );
        Vec_3D force_1_3 = B_1.calculate_force( B_3 );
        Vec_3D force_2_3 = B_2.calculate_force( B_3 );
        
        B_1.reset_acceleration();
        B_2.reset_acceleration();
        B_3.reset_acceleration();
        
        B_1.apply_force( force_1_2 );
        B_1.apply_force( force_1_3 );
        
        B_2.apply_force( force_1_2 * -1.0 );
        B_2.apply_force( force_2_3 );
        
        B_3.apply_force( force_1_3 * -1.0 );
        B_3.apply_force( force_2_3 * -1.0 );
        
        B_1.update( dt );
        B_2.update( dt );
        B_3.update( dt );

        std::cout << "\n<--- Step: " << ( current_step + 1 ) << " --->" << std::endl;
        
        Vec_3D pos1 = B_1.get_pos();
        Vec_3D vel1 = B_1.get_vel();
        std::cout << "Body 1: Pos(" << pos1.x_/1000.0 << ", " 
                                    << pos1.y_/1000.0 << ", " 
                                    << pos1.z_/1000.0 << ") km, ";
        std::cout << "Vel: " << vel1.magnitude() << " m/s" << std::endl;
        
        Vec_3D pos2 = B_2.get_pos();
        Vec_3D vel2 = B_2.get_vel();
        std::cout << "Body 2: Pos(" << pos2.x_/1000.0 << ", " 
                                    << pos2.y_/1000.0 << ", " 
                                    << pos2.z_/1000.0 << ") km, ";
        std::cout << "Vel: " << vel2.magnitude() << " m/s" << std::endl;
        
        Vec_3D pos3 = B_3.get_pos();
        Vec_3D vel3 = B_3.get_vel();
        std::cout << "Body 3: Pos(" << pos3.x_/1000.0 << ", " 
                                    << pos3.y_/1000.0 << ", " 
                                    << pos3.z_/1000.0 << ") km, ";
        std::cout << "Vel: " << vel3.magnitude() << " m/s" << std::endl;
        
        std::cout << "Distances: B1-B2: " << ( B_1.get_pos() - B_2.get_pos() ).magnitude() / 1000.0 << " km, ";
        std::cout << "B1-B3: " << ( B_1.get_pos() - B_3.get_pos() ).magnitude() / 1000.0 << " km, ";
        std::cout << "B2-B3: " << ( B_2.get_pos() - B_3.get_pos() ).magnitude() / 1000.0 << " km" << std::endl;

        ++current_step;
    }

    std::cout << "\n<--- End of Simulation --->" << std::endl;

    return 0;
}