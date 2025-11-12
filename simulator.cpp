#include <iostream>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

static constexpr double G{ 6.67e-11 };
static constexpr double CONVERT_TO_KMS{ 1.0 / 1000.0 };
static constexpr double CONVERT_TO_KM{ 1.0 / 1000.0 };

struct Vec_3D {
    double x_, y_, z_;
    double norm() {
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
    
    void calculate_new_acc( std::vector<Body> const &other_bodies ) {
        Vec_3D total_force{};

        for ( std::size_t idx{ 0 }; idx < other_bodies.size(); ++idx ) {
            Vec_3D R{ other_bodies[idx].get_pos() - pos_ };
            double dist{ R.norm() };

            if ( dist < 1.0e-10 ) continue;
            double force_mag{ ( G * mass_ * other_bodies[idx].get_mass() ) / ( dist * dist ) };

            total_force += R * ( force_mag / dist );
        }
        acc_ = total_force * ( 1.0 / mass_ );
    }

    void update( double dt ) {
        vel_ += acc_ * dt;
        pos_ += vel_ * dt;
    }

    Vec_3D get_pos() const { return pos_; }
    Vec_3D get_vel() const { return vel_; }
    Vec_3D get_acc() const { return acc_; }
    double get_mass() const { return mass_; }
};

// Calculates total energy of all bodies.
double calculate_total_energy( std::vector<Body> const &bodies ) {
    double K{ 0.0 };
    double V{ 0.0 };

    for (size_t i = 0; i < bodies.size(); ++i) {
        K += 0.5 * bodies[i].get_mass() * bodies[i].get_vel().norm() * bodies[i].get_vel().norm();
    }

    for (size_t i = 0; i < bodies.size(); ++i) {
        for (size_t j = i + 1; j < bodies.size(); ++j) {
            Vec_3D R = bodies[i].get_pos() - bodies[j].get_pos();
            double dist = R.norm();

            if (dist < 1.0e-10) continue;

            V -= G * bodies[i].get_mass() * bodies[j].get_mass() / dist;
        }
    }
    
    return K + V;
}
// Load bodies.
void load_csv_bodies( std::string file_name, std::vector<Body> &bodies ) {
    std::ifstream file( file_name );
    if ( !file.is_open() ) {
        std::cout << "Error opening " << file_name << ".csv." << std::endl;
        return;
    }

    std::string line{};
    while ( std::getline( file, line ) ) {
        if ( line.empty() || line[0] == '#' ) continue;

        std::stringstream string_stream( line );
        std::string segment{};
        std::vector<double> values{};

        while ( std::getline( string_stream, segment, ',' ) ) {
            values.push_back( std::stod( segment ) );
        }

        if ( values.size() == 7 ) {
            bodies.emplace_back( Vec_3D{ values[0], values[1], values[2] },   // Positions
                                   Vec_3D{ values[3], values[4], values[5] }, // Velocities
                                   values[6]                                  // Mass
                                );
        }
    }
    file.close();
}

// Load the dt and number of steps.
void load_simulation_settings( std::string file_name, double &dt, int &steps ) {
    std::ifstream file( file_name );
    if ( !file.is_open() ) {
        std::cout << "Error opening " << file_name << ".csv." << std::endl;
        return;
    }

    std::string line{};
    while ( std::getline( file, line ) ) {
        if ( line.empty() || line[0] == '#' ) continue;

        std::stringstream string_stream( line );
        std::string segment{};
        std::vector<double> values{};

        while ( std::getline( string_stream, segment, ',' ) ) {
            values.push_back( std::stod( segment ) );
        }

        if ( values.size() == 2 ) {
            steps = values[0];
            dt = values[1];
        }
    }
    file.close();
}

int main() {
    double dt{};
    int steps{};
    load_simulation_settings( "config.csv", dt, steps );

    std::vector<Body> bodies{};
    load_csv_bodies( "bodies.csv", bodies );

    std::cout << "<--- N-Body Simulation --->" << std::endl;
    std::cout << "\nStarting N-Body Simulation..." << std::endl;
    std::cout << std::fixed << std::setprecision( 4 );

    double initial_energy{ calculate_total_energy( bodies ) };
    double max_energy_drift{};

    for( int current_step{}; current_step < steps; ++current_step ) {
        // Calculates new acceleration.
        for ( std::size_t idx{ 0 }; idx < bodies.size(); ++idx ) {
            bodies[idx].calculate_new_acc( bodies );
        }

        // Updates position and velocity for all bodies.
        for ( std::size_t idx{ 0 }; idx < bodies.size(); ++idx ) {
            bodies[idx].update( dt );
        }

        std::cout << "\n<--- Step: " << ( current_step + 1 ) << " --->" << std::endl;

        // Displays the current position and velocity for all bodies.
        for ( std::size_t idx{ 0 }; idx < bodies.size(); ++idx ) {
            Vec_3D curr_body_pos{ bodies[idx].get_pos() };
            Vec_3D curr_body_vel{ bodies[idx].get_vel() };

            std::cout << "Body " << idx + 1 << ": Pos(" << curr_body_pos.x_ * CONVERT_TO_KM << ", " 
                                                    << curr_body_pos.y_ * CONVERT_TO_KM << ", " 
                                                    << curr_body_pos.z_ * CONVERT_TO_KM << ") km, ";

            std::cout << "Vel(" << curr_body_vel.x_ * CONVERT_TO_KMS << ", " 
                                << curr_body_vel.y_ * CONVERT_TO_KMS << ", " 
                                << curr_body_vel.z_ * CONVERT_TO_KMS << ") km/s, ";
            std::cout << "Speed: " << curr_body_vel.norm() * CONVERT_TO_KMS << " km/s" << std::endl;
        }
        std::cout << std::endl;

        // Displays distance form all bodies to the other.
        for ( std::size_t idx{ 0 }; idx < bodies.size(); ++idx ) {
            for ( std::size_t second_idx{ idx + 1 }; second_idx < bodies.size(); ++second_idx ) {
                Vec_3D R{ bodies[idx].get_pos() - bodies[second_idx].get_pos() };
                std::cout << "Distance: B" << idx + 1 << "-B"
                                           << second_idx + 1 << ": "
                                           << R.norm() * CONVERT_TO_KM << " km" << std::endl;
            }
        }

        double current_energy{ calculate_total_energy( bodies ) };
        double energy_drift_percent{ 100.0 * std::abs( current_energy - initial_energy ) / std::abs( initial_energy ) };
        max_energy_drift = std::max( max_energy_drift, energy_drift_percent );
    }

    std::cout << "Max Energy Drift: " << max_energy_drift << "%." << std::endl;
    std::cout << "\n<--- End of Simulation --->" << std::endl;

    return 0;
}