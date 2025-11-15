#include <iostream>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <chrono>
#include <omp.h>

static constexpr double G{ 6.67430e-11 };
static constexpr double EPSILON{ 1.0e2 };
static constexpr double CONVERT_TO_KMS{ 1.0 / 1000.0 };
static constexpr double CONVERT_TO_KM{ 1.0 / 1000.0 };
static constexpr double CONVERT_TO_SEC{ 1.0e-9 };

struct Vec_3D {
    double x_, y_, z_;
    double norm() const {
        return std::sqrt( x_*x_ + y_*y_ + z_*z_ );
    }
    double norm_squared() const {
        return x_*x_ + y_*y_ + z_*z_;
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
    Vec_3D old_acc_;
    double mass_;

public:
    Body ( Vec_3D pos, Vec_3D vel, double mass ) 
         : pos_( pos ), vel_( vel ), acc_{ 0, 0, 0 }, old_acc_{ 0, 0, 0 }, mass_( mass ) {}
    
    // Calculates new acceleration based on forces from other bodies.
    void calculate_new_acc( std::vector<Body> const &other_bodies, std::size_t const &self_idx ) {
        old_acc_ = acc_;
        Vec_3D total_force{};

        for ( std::size_t idx = 0; idx < other_bodies.size(); ++idx ) {
            if ( idx == self_idx ) continue;
            
            Vec_3D R{ other_bodies[idx].get_pos() - pos_ };
            double dist_squared{ R.norm_squared() + EPSILON * EPSILON };
            double force_mag{ ( G * mass_ * other_bodies[idx].get_mass() ) / ( dist_squared ) };

            total_force += R * ( force_mag / std::sqrt( dist_squared ) );
        }
        acc_ = total_force * ( 1.0 / mass_ );
    }

    // Updates body.
    void update( double const &dt ) {
        pos_ += vel_ * dt + acc_ * ( 0.5 * dt * dt );
        vel_ += ( acc_ + old_acc_ ) * ( 0.5 * dt );
    }

    // Body getters.
    Vec_3D get_pos() const { return pos_; }
    Vec_3D get_vel() const { return vel_; }
    Vec_3D get_acc() const { return acc_; }
    double get_mass() const { return mass_; }
};

// Calculates total energy of all bodies.
double calculate_total_energy( std::vector<Body> const &bodies ) {
    double K{ 0.0 };
    double V{ 0.0 };

    for ( std::size_t i = 0; i < bodies.size(); ++i ) {
        for ( size_t j = i + 1; j < bodies.size(); ++j ) {
            Vec_3D R{ bodies[i].get_pos() - bodies[j].get_pos() };
            double dist{ R.norm() + EPSILON };

            V -= G * bodies[i].get_mass() * bodies[j].get_mass() / dist;
        }
        K += 0.5 * bodies[i].get_mass() * bodies[i].get_vel().norm_squared();
    }
    
    return K + V;
}

// Load bodies.
void load_csv_bodies( std::string const &file_name, std::vector<Body> &bodies ) {
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
            bodies.emplace_back( Vec_3D{ values[0], values[1], values[2] }, // Positions
                                 Vec_3D{ values[3], values[4], values[5] }, // Velocities
                                 values[6]                                  // Mass
                                );
        }
    }
    file.close();
}

int main() {
    double dt{ 1000 };
    int steps{};
    int num_outputs{};

    std::vector<Body> bodies{};
    load_csv_bodies( "bodies.csv", bodies );

    std::cout << "<--- N-Body Simulation --->" << std::endl;

    std::cout << "Enter number of steps: ";
    std::cin >> steps;
    std::cout << "Enter number of outputs: ";
    std::cin >> num_outputs;

    std::cout << "\nStarting N-Body Simulation..." << std::endl;
    std::ofstream out_file( "trajectories.csv" );
    out_file << "step,body_id,x,y,z\n";

    // Keep values in scientific notation to 3 sig figs.
    std::cout << std::scientific << std::setprecision( 3 );
    double initial_energy{ calculate_total_energy( bodies ) };
    double max_energy_drift{};

    auto start_time{ std::chrono::high_resolution_clock::now() };

    #pragma omp parallel 
    {
        for( int current_step{}; current_step < steps; ++current_step ) {
            // Calculates new acceleration.
            #pragma omp for
            for ( std::size_t idx = 0; idx < bodies.size(); ++idx ) {
                bodies[idx].calculate_new_acc( bodies, idx );
            }

            // Updates position and velocity for all bodies.
            #pragma omp for
            for ( std::size_t idx = 0; idx < bodies.size(); ++idx ) {
                bodies[idx].update( dt );
            }

            // Outputs information in specific number of outputs.
            int output_interval{};
            output_interval = steps / num_outputs;

            #pragma omp single
            {
                if (current_step % output_interval == 0 || current_step == steps - 1) {
                    // Calculates maximum energy drift in system.
                    double current_energy{ calculate_total_energy( bodies ) };
                    double energy_drift_percent{ 100.0 * std::abs( current_energy - initial_energy ) / std::abs( initial_energy ) };
                    max_energy_drift = std::max( max_energy_drift, energy_drift_percent );

                    // Outputs the current position for all bodies.
                    for ( std::size_t idx = 0; idx < bodies.size(); ++idx ) {
                        Vec_3D curr_body_pos{ bodies[idx].get_pos() };
                        Vec_3D curr_body_vel{ bodies[idx].get_vel() };

                        out_file << current_step << "," 
                                 << idx << ","
                                 << curr_body_pos.x_ << ","
                                 << curr_body_pos.y_ << ","
                                 << curr_body_pos.z_ << "\n";
                    }
                }

                if ( current_step % (steps / 100) == 0 ) {
                    std::cout << "Progress: " << (current_step * 100 / steps) << "%\r" << std::flush;
                }
            }
        }
    }
    auto time_elapsed{ ( std::chrono::high_resolution_clock::now() - start_time ).count() * CONVERT_TO_SEC };

    std::cout << std::endl << std::fixed << std::scientific << std::setprecision( 4 );
    std::cout << "\nMax Energy Drift: " << max_energy_drift << "%." << std::endl;
    std::cout << "Time elapsed: " << time_elapsed << " seconds." << std::endl;
    std::cout << "\n<--- End of Simulation --->" << std::endl;

    out_file.close();

    return 0;
}